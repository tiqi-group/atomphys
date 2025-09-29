#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 09/2023
# Author: Wojciech Adamczyk <wadamczyk@phys.ethz.ch>

import pint
import numpy as np
from atomphys.calc.rabi_frequency import Rabi_Frequency
from atomphys.electric_field import ElectricField
from atomphys.state import State
from atomphys.utils.utils import spherical_basis_vector
from atomphys.transition import TransitionType
from atomphys.calc.matrix_element import electric_dipole_matrix_element
from atomphys.calc.matrix_element import dipole_matrix_element_basis


def ac_stark_shift(
    state: State,
    mJ: float,
    El_field: ElectricField,
    _ureg: pint.UnitRegistry,
    wavelengths: np.ndarray | None = None,
    #B: pint.Quantity | float = 0.0,
):
    #B = set_default_units(B, "tesla")

    if wavelengths is None:
        omega_field = El_field.angular_frequency
        delta_E = 0 * _ureg("k*mK")
    else:
        omega_field = 2 * np.pi * _ureg("c") / wavelengths
        # Ensure array-shaped accumulator even if no transitions contribute
        delta_E = np.zeros_like(omega_field.magnitude, dtype=float) * _ureg("k*mK")

    for transition in state.transitions_from:
        #zeeman_shifts = transition.sublevels_zeeman_shift(B)
        state_i = state
        state_f = transition.state_f

        mJ_i = mJ
        for mJ_f in state_f.sublevels:
            try:
                tr_omega = transition.angular_frequency  # + zeeman_shifts[(mJ_i, mJ_f)]
                Omega = complex(
                    (
                        Rabi_Frequency(
                            E_field=El_field,
                            transition=transition,
                            mJ_i=mJ_i,
                            mJ_f=mJ_f,
                            _ureg=_ureg,
                        )
                    )
                    .to("MHz")
                    .magnitude
                ) * _ureg("MHz")
                stark_contribution = (
                    _ureg("hbar")
                    / 4
                    * (
                        (Omega * np.conj(Omega)) / (-tr_omega - omega_field)
                        + (Omega * np.conj(Omega)) / (-tr_omega + omega_field)
                    )
                )
                # Ensure we only add the real part to avoid dtype casting issues
                delta_E += stark_contribution.real
            except Exception as e:
                print(f"Exception in transitions_from loop: {e}")
                pass

    for transition in state.transitions_to:
        #zeeman_shifts = transition.sublevels_zeeman_shift(B)
        state_i = transition.state_i
        state_f = state

        mJ_f = mJ
        for mJ_i in state_i.sublevels:
            try:
                tr_omega = transition.angular_frequency  # + zeeman_shifts[(mJ_i, mJ_f)]
                Omega = complex(
                    (
                        Rabi_Frequency(
                            E_field=El_field,
                            transition=transition,
                            mJ_i=mJ_i,
                            mJ_f=mJ_f,
                            _ureg=_ureg,
                        )
                    )
                    .to("MHz")
                    .magnitude
                ) * _ureg("MHz")
                stark_contribution = (
                    _ureg("hbar")
                    / 4
                    * (
                        (Omega * np.conj(Omega)) / (-tr_omega - omega_field)
                        + (Omega * np.conj(Omega)) / (-tr_omega + omega_field)
                    )
                )
                # Ensure we only add the real part to avoid dtype casting issues
                delta_E += stark_contribution.real
            except Exception as e:
                print(f"Exception in transitions_from loop: {e}")
                pass

    if wavelengths is None:
        return delta_E.to("k*mK").magnitude.real * _ureg("k*mK")
    else:
        return (delta_E.to("k*mK").magnitude.real) * _ureg("k*mK")

# This is work in progress, but so far it reproduces Sec. III.A. from
# https://journals.aps.org/pra/pdf/10.1103/PhysRevA.97.063419
# and simplifies https://arxiv.org/pdf/2505.22466 by not considering ladder schemes.
def off_resonant_scattering_rate(
    state_I: State,
    state_F: State,
    mJ_i: float,
    mJ_f: float,
    El_field: ElectricField,
    _ureg: pint.UnitRegistry,
    wavelengths: np.ndarray | None = None,
):
    """
    Off-resonant scattering rate for E1 transitions only.

    Implements the expression from the docs (see ac_stark_shifts.md) specialized to electric-dipole (E1) transitions,
    summing over intermediate states k and spherical components q in {-1,0,1}.
    """
    if wavelengths is None:
        omega_l = El_field.angular_frequency
        Gamma = 0 * _ureg("Hz")
    else:
        omega_l = 2 * np.pi * _ureg("c") / wavelengths
        # Ensure array-shaped accumulator even if no transitions contribute
        Gamma = np.zeros_like(omega_l.magnitude, dtype=float) * _ureg("Hz")

    omega_if = (state_I.energy - state_F.energy).to("THz", "sp") * _ureg("_2pi")
    omega_prime = omega_l + omega_if
    prefactor = (omega_l ** 2) * (omega_prime ** 3) / _ureg("12*pi*vacuum_permittivity*c^3*hbar")
    J_i = state_I.quantum_numbers["J"]
    J_f = state_F.quantum_numbers["J"]
    
    Gamma = 0 * _ureg('Hz')
    
    for q in (-1, 0, 1):
        D_q = 0 * _ureg('e*a0*s')
        for transition in state_I.transitions_from:
            if transition.type != TransitionType.E1:
                continue
            state_i = state_I
            state_k = transition.state_f
            J_k = state_k.quantum_numbers["J"]
            omega_ki = (state_k.energy - state_i.energy).to("THz", "sp") * _ureg("_2pi")
            omega_fk = (state_F.energy - state_k.energy).to("THz", "sp") * _ureg("_2pi")

            for mJ_k in state_k.sublevels:
                try:
                    Omega = complex(
                        (
                            Rabi_Frequency(
                                E_field=El_field,
                                transition=transition,
                                mJ_i=mJ_i,
                                mJ_f=mJ_k,
                                _ureg=_ureg,
                            )
                        )
                        .to("Hz")
                        .magnitude
                    ) * _ureg("Hz")
                    
                    fdqk = dipole_matrix_element_basis(transition.A, transition.k, J_k, J_f, mJ_k, mJ_f, q, _ureg)*_ureg('e')
                    kdqi = dipole_matrix_element_basis(transition.A, transition.k, J_i, J_k, mJ_i, mJ_k, q, _ureg)*_ureg('e')
                    
                    D_q += fdqk * Omega / (omega_ki - omega_l) / omega_ki + kdqi * Omega / (omega_ki + omega_prime) / omega_fk
                except Exception as e:
                    print(f"Exception in transitions_from loop: {e}")
                    pass
        
        for transition in state_I.transitions_to:
            if transition.type != TransitionType.E1:
                continue
            state_i = state_I
            state_k = transition.state_i
            omega_ki = (state_k.energy - state_i.energy).to("THz", "sp") * _ureg("_2pi")
            omega_fk = (state_F.energy - state_k.energy).to("THz", "sp") * _ureg("_2pi")
            J_k = state_k.quantum_numbers["J"]

            for mJ_k in state_k.sublevels:
                try:
                    Omega = complex(
                        (
                            Rabi_Frequency(
                                E_field=El_field,
                                transition=transition,
                                mJ_i=mJ_k,
                                mJ_f=mJ_i,
                                _ureg=_ureg,
                            )
                        )
                        .to("Hz")
                        .magnitude
                    ) * _ureg("Hz")
                    
                    fdqk = dipole_matrix_element_basis(transition.A, transition.k, J_k, J_f, mJ_k, mJ_f, q, _ureg)*_ureg('e')
                    kdqi = dipole_matrix_element_basis(transition.A, transition.k, J_i, J_k, mJ_i, mJ_k, q, _ureg)*_ureg('e')
                    
                    D_q += fdqk * Omega / (omega_ki - omega_l) / omega_ki + kdqi * Omega / (omega_ki + omega_prime) / omega_fk
                except Exception as e:
                    print(f"Exception in transitions_to loop: {e}")
                    pass
        
        Gamma += prefactor * D_q * D_q.conj()
        
    return Gamma.to("Hz")
    
    
                


    # Build polarization unit vector epsilon_l (cartesian Jones) from the field
    E_vec = El_field.field()  # [V/m]
    E0 = np.linalg.norm(E_vec)
    if E0 == 0:
        return 0 * _ureg("1/s")
    # epsilon (unit vector) not explicitly used further because we compute d·epsilon via Rabi_Frequency

    # Helper to compute <a| d_q | b> in units of e a0 by projecting cartesian vector onto spherical basis
    def d_q(a: State, mJ_a: float, b: State, mJ_b: float, q: int):
        # Only include transitions that are E1-allowed between a and b
        # Find matching Transition object if it exists
        candidates = []
        for tr in a.transitions_from:
            if tr.state_f == b and tr.type == TransitionType.E1:
                candidates.append(tr)
        for tr in a.transitions_to:
            if tr.state_i == b and tr.type == TransitionType.E1:
                candidates.append(tr)
        if len(candidates) == 0:
            return 0.0 * _ureg("e*a0")
        # Use the first matching E1 transition (A, k define magnitude; sign via w3j inside electric_dipole_matrix_element)
        tr = candidates[0]
        J_i = tr.state_i.quantum_numbers["J"]
        J_f = tr.state_f.quantum_numbers["J"]
        # Determine direction of matrix element based on ordering of states in the Transition
        if tr.state_i == a and tr.state_f == b:
            val = electric_dipole_matrix_element(tr.A, tr.k, J_i, J_f, mJ_a, mJ_b, _ureg)
        else:
            # Element <a|d|b> = conjugate of <b|d|a>; for real Wigner-3j-based elements this is equal
            val = electric_dipole_matrix_element(tr.A, tr.k, J_f, J_i, mJ_b, mJ_a, _ureg)
        # Project onto spherical component q via cartesian basis vector c_q
        c_q = spherical_basis_vector(q)
        # electric_dipole_matrix_element returns cartesian vector valued quantity expressed via spherical expansion inside
        # Here, treat val as cartesian vector magnitude along c_q direction by dot product
        # If val is scalar in units e a0, embed along c_q direction
        # If val is vector-like (np.ndarray), project; otherwise return scalar
        proj = (np.dot(c_q, val) if hasattr(val, "__len__") else val)
        return proj.to("e*a0")

    # Build the Rabi coupling <k| d·epsilon | i> and <f| d·epsilon | k>
    def d_dot_eps(a: State, mJ_a: float, b: State, mJ_b: float):
        # d·epsilon = sum_mu epsilon^mu d_mu (cartesian). Using electric_dipole_matrix_element which already returns vector via w3j and spherical basis
        # Reuse Rabi_Frequency definition: Omega_ab = (E0 e / ħ) epsilon · d_ab
        # => epsilon · d_ab = (ħ/eE0) Omega_ab
        # We can compute Omega via existing helper restricted to E1 automatically by transition.type
        # Find relevant E1 transition between a and b
        tr_match = None
        for tr in a.transitions_from:
            if tr.state_f == b and tr.type == TransitionType.E1:
                tr_match = tr
                break
        if tr_match is None:
            for tr in a.transitions_to:
                if tr.state_i == b and tr.type == TransitionType.E1:
                    tr_match = tr
                    break
        if tr_match is None:
            return 0 * _ureg("e*a0")
        Omega = Rabi_Frequency(El_field, tr_match, mJ_a, mJ_b, _ureg).to("1/s")
        return (Omega * _ureg("hbar/e") / E0).to("e*a0")

    # Sum over intermediate states k that are E1-coupled to i and f
    amplitude_sum = 0 * _ureg("e^2*a0^2*s^2")
    for tr_k in state_i.transitions_from:
        if tr_k.type != TransitionType.E1:
            continue
        k = tr_k.state_f
        omega_ki = tr_k.angular_frequency
        for mJ_k in k.sublevels:
            # term 1: <f| d_q | k> * <k| d·epsilon | i> / ((omega_ki - omega_l) * omega_ki)
            num1 = 0 + 0j
            for q in (-1, 0, 1):
                d_fk_q = d_q(state_f, mJ_f, k, mJ_k, q)
                num1 += d_fk_q
            dkei = d_dot_eps(k, mJ_k, state_i, mJ_i)
            term1 = (num1 * dkei) / (omega_ki - omega_l) / omega_ki

            # term 2: <f| d·epsilon | k> * <k| d_q | i> / ((omega_fk) * (omega_ki + omega_l + omega_if))
            # omega_fk is transition frequency between f and k
            # Find transition f<->k
            tr_fk = None
            for tr in state_f.transitions_from:
                if tr.state_f == k and tr.type == TransitionType.E1:
                    tr_fk = tr
                    break
            if tr_fk is None:
                for tr in state_f.transitions_to:
                    if tr.state_i == k and tr.type == TransitionType.E1:
                        tr_fk = tr
                        break
            if tr_fk is None:
                omega_fk = 0 * _ureg("1/s")
            else:
                omega_fk = tr_fk.angular_frequency
            dfke = d_dot_eps(state_f, mJ_f, k, mJ_k)
            num2 = 0 + 0j
            for q in (-1, 0, 1):
                d_ki_q = d_q(k, mJ_k, state_i, mJ_i, q)
                num2 += d_ki_q
            denom2 = (omega_fk if omega_fk.m != 0 else 1 * _ureg("1/s")) * (omega_ki + omega_l + omega_if)
            term2 = (dfke * num2) / denom2

            amplitude_sum += (term1 + term2)

    amplitude_sq = (amplitude_sum * np.conj(amplitude_sum))
    Gamma = (prefactor * amplitude_sq).to("1/s")
    return Gamma
    


def polarizability(
    state: State,
    mJ: float,
    El_field: ElectricField,
    _ureg: pint.UnitRegistry,
    wavelengths: np.ndarray | None = None,
):
    AC_stark = ac_stark_shift(
        state=state, mJ=mJ, El_field=El_field, _ureg=_ureg, wavelengths=wavelengths
    )
    E0 = El_field.field_amplitude
    return (-4 * AC_stark / E0**2).to("e^2*a_0^2/E_h")

