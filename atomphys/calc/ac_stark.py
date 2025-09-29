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
from atomphys.atom import Atom
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
    atom: Atom,
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
        for transition_ik in state_I.transitions_from:
            state_i = state_I
            state_k = transition_ik.state_f
            transition_fk = atom.transition_between(state_F, state_k) if state_F.energy < state_k.energy else atom.transition_between(state_k, state_F)
            omega_ki = (state_k.energy - state_i.energy).to("THz", "sp") * _ureg("_2pi")
            omega_fk = (state_F.energy - state_k.energy).to("THz", "sp") * _ureg("_2pi")
            J_k = state_k.quantum_numbers["J"]

            if transition_fk==None:
                continue
            if transition_ik.type != TransitionType.E1 or transition_ik.state_f != state_k:
                continue

            for mJ_k in state_k.sublevels:
                try:
                    Omega_ik = complex(
                        (
                            Rabi_Frequency(
                                E_field=El_field,
                                transition=transition_ik,
                                mJ_i=mJ_i,
                                mJ_f=mJ_k,
                                _ureg=_ureg,
                            )
                        )
                        .to("Hz")
                        .magnitude
                    ) * _ureg("Hz")
                    
                    Omega_fk = complex(
                        (
                            Rabi_Frequency(
                                E_field=El_field,
                                transition=transition_fk,
                                mJ_i=mJ_k,
                                mJ_f=mJ_f,
                                _ureg=_ureg,
                            )
                        )
                        .to("Hz")
                        .magnitude
                    ) * _ureg("Hz")
                    
                    fdqk = dipole_matrix_element_basis(transition_fk.A, transition_fk.k, J_k, J_f, mJ_k, mJ_f, q, _ureg)*_ureg('e')
                    kdqi = dipole_matrix_element_basis(transition_ik.A, transition_ik.k, J_i, J_k, mJ_i, mJ_k, q, _ureg)*_ureg('e')
                    D_q += fdqk * Omega_ik / (omega_ki - omega_l) / omega_ki + kdqi * Omega_fk / (omega_ki + omega_prime) / omega_fk
                except Exception as e:
                    print(f"Exception in transitions_to loop: {e}")
                    pass
        
        
        for transition_ik in state_I.transitions_to:
            state_i = state_I
            state_k = transition_ik.state_i
            transition_fk = atom.transition_between(state_F, state_k) if state_F.energy < state_k.energy else atom.transition_between(state_k, state_F)
            omega_ki = (state_k.energy - state_i.energy).to("THz", "sp") * _ureg("_2pi")
            omega_fk = (state_F.energy - state_k.energy).to("THz", "sp") * _ureg("_2pi")
            J_k = state_k.quantum_numbers["J"]
            if transition_fk==None:
                continue
            if transition_ik.type != TransitionType.E1 or transition_ik.state_f != state_k:
                continue
            

            for mJ_k in state_k.sublevels:
                try:
                    Omega_ik = complex(
                        (
                            Rabi_Frequency(
                                E_field=El_field,
                                transition=transition_ik,
                                mJ_i=mJ_i,
                                mJ_f=mJ_k,
                                _ureg=_ureg,
                            )
                        )
                        .to("Hz")
                        .magnitude
                    ) * _ureg("Hz")
                    
                    Omega_fk = complex(
                        (
                            Rabi_Frequency(
                                E_field=El_field,
                                transition=transition_fk,
                                mJ_i=mJ_k,
                                mJ_f=mJ_f,
                                _ureg=_ureg,
                            )
                        )
                        .to("Hz")
                        .magnitude
                    ) * _ureg("Hz")
                    
                    fdqk = dipole_matrix_element_basis(transition_fk.A, transition_fk.k, J_k, J_f, mJ_k, mJ_f, q, _ureg)*_ureg('e')
                    kdqi = dipole_matrix_element_basis(transition_ik.A, transition_ik.k, J_i, J_k, mJ_i, mJ_k, q, _ureg)*_ureg('e')
                    D_q += fdqk * Omega_ik / (omega_ki - omega_l) / omega_ki + kdqi * Omega_fk / (omega_ki + omega_prime) / omega_fk
                except Exception as e:
                    print(f"Exception in transitions_to loop: {e}")
                    pass
        
        Gamma += prefactor * D_q * D_q.conj()    
    return Gamma.to("Hz")
    


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

