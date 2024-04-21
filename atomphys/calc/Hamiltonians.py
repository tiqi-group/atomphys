#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 09/2023
# Authors:
#   Wojciech Adamczyk <wadamczyk@phys.ethz.ch>
#   Carmelo Mordini <cmordini@phys.ethz.ch>

from .rabi_frequency import Rabi_Frequency
from .zeeman import zeeman_shift
import qutip
from ..atom import Atom
from ..state import State
import pint
from ..electric_field import ElectricField
from .util import find_rotating_frame
from ..transition import Transition
from .lindblad_operators import sqrt_lindblad_operator

import numpy as np
from itertools import combinations


def kets(states: list[State]):
    """
    For given list of states, returns a dictionary that maps all posible (state, mJ) pairs to a ket in the basis of all states (sorted by energy
    and then by mJ)

    Args:
        states (list[State]): List of states to be included in the hamiltonian

    Returns:
        dict[(State, int), qutip.Qobj]: Dictionary that maps all posible (state, mJ) pairs to a ket in the basis of all states (sorted by energy
        and then by mJ)
    
    """

    states.sort(key=lambda state: state.energy)
    states_m = []
    for state in states:
        for m in state.sublevels:
            states_m.append((state, m))

    total_number_of_states = len(states_m)
    kets_dict = {states_m[i]: qutip.basis(total_number_of_states, i) for i in range(total_number_of_states)}
    return kets_dict

def zero_hamiltonian(states):
    n_states = sum(len(s.sublevels) for s in states)
    return qutip.qzero(n_states)

def H0(states: list[State], _ureg: pint.UnitRegistry):
    """"
    Returns the atomic hamiltonian with no B field.

    H0 = sum_i (E_i |i><i|)

    Args:
        states (list[State]): List of states to be included in the hamiltonian
        _ureg (pint.UnitRegistry): Unit registry
    Returns:
        qutip.Qobj: Hamiltonian
    """
    H = 0
    kets_dict = kets(states)
    for (state, m), ket in kets_dict.items():
        H += (state.angular_frequency.to("MHz")).m * ket * ket.dag()
    return H

def H_zeeman(states: list[State], B_field: pint.Quantity):
    """
    Returns extra zeeman shift due to the B field. This is not the full hamiltonian, but rather a correction to the H0 hamiltonian.

    H_zeeman = sum_i (gJ_i mJ_i mu_b B |i><i|)
    """

    H = 0
    kets_dict = kets(states)
    for (state, m), ket in kets_dict.items():
        H += (zeeman_shift(state.g, m, B_field, state._ureg).to("MHz")).magnitude * ket * ket.dag()
    return H




def H_int(
    atom: Atom,
    states: list[State],
    fields: dict[ElectricField, list[Transition]],
    _ureg: pint.UnitRegistry,
):
    """Returns the interaction hamiltonian in rotating frame of reference for a given atom and electric fields, with states decided by the user
    The idea is that if the electric fields do not form a closed roop, one can enter rotating frame of reference (interaction picture), such
    that fields either rotate at 0 frequency or the frequency much higher than the rabi frequency of any transition. This leaves the hamiltonian
    time independent, and so is relatively cheap and computationally easy to solve it using qutip, and hence the interesting case for it.

    I am aware that for some cases, like raman transitions my algorithm woulc claim that the rotating frame of reference is not possible, but it is.
    I dont have time right now to implement it now.

    Args:
        atom (Atom): Atom object
        states (list[State]): List of states to be included in the hamiltonian
        fields (dict[ElectricField, list[Transition]]): Dictionary of electric fields and transitions that are driven by them
        _ureg (pint.UnitRegistry): Unit registry

    Returns:
        qutip.Qobj: Hamiltonian in rotating frame of reference plus correction on the base hamiltonian
    """
    H_total = []
    kets_dict = kets(states)

    for field in fields:
        H = 0
        H_c = 0
        for (state_i, m_i), ket_i in kets_dict.items():
            for (state_f, m_f), ket_j in kets_dict.items():

                transition = atom.transition_between(state_i, state_f)
                if transition is None:
                    continue
                Omega_ij = abs(Rabi_Frequency(field, transition, mJ_i=m_i, mJ_f=m_f).to("MHz").m)

                H_td_coeff = f'exp(1j * {field.angular_frequency.to("MHz").m} * t)'

                H_td_coeff_c = f'exp(-1j * {field.angular_frequency.to("MHz").m} * t)'
                
                H += -1 / 2 * Omega_ij * ket_i * ket_j.dag()
                H_c += -1 / 2 * np.conj(Omega_ij) * ket_i * ket_j.dag()

                H += -1 / 2 * np.conj(Omega_ij) * ket_j * ket_i.dag()
                H_c += -1 / 2 * Omega_ij * ket_j * ket_i.dag()

        H_total.append([H, H_td_coeff, field.angular_frequency])
        H_total.append([H_c, H_td_coeff_c, -field.angular_frequency])
    return H_total


def merge_redundancies(H_total):
    """
        Given the list of hamiltonians with their associated time dependence, merges the redundant terms.
        It merges them by comparing the frequency of the time dependence. If the frequencies are the same, it merges the terms into one.

        Merges terms in the hamiltonian that have the same time dependence. 
        The idea is that



        Merges redundant terms in the hamiltonian. This is done by adding the hamiltonians that have the same time dependence.
        This is done to reduce the number of terms in the hamiltonian, which is important for the performance of the solver.

        Args:
            H_total (list): List of hamiltonians to be merged

        Returns:
            list: List of merged hamiltonians
    """
    H_total_merged = []
    for H, H_td_coeff, omega in H_total:
        found = False
        i=0
        for H_merged, H_td_coeff_merged, omega_merged in H_total_merged:
            if omega == omega_merged:
                H_total_merged[i] += H
                found = True
                break
            i+=1
        if not found:
            H_total_merged.append([H, H_td_coeff, omega])
    return H_total_merged


def kets(atom: Atom, states: list[State]):
    states.sort(key=lambda state: state.energy)
    states_m = []
    for state in states:
        for m in state.sublevels:
            states_m.append((state, m))

    total_number_of_states = len(states_m)
    kets_dict = {states_m[i]: qutip.basis(total_number_of_states, i) for i in range(total_number_of_states)}
    return kets_dict


def zero_hamiltonian(states):
    n_states = sum(len(s.sublevels) for s in states)
    return qutip.qzero(n_states)


def H0(atom: Atom, states: list[State], _ureg: pint.UnitRegistry):
    """Returns the atomic hamiltonian with no B field"""
    H = 0
    kets_dict = kets(atom, states)
    for (state, m), ket in kets_dict.items():
        H += (state.angular_frequency.to("MHz")).magnitude * ket * ket.dag()
    return H


def H_zeeman(atom: Atom, states: list[State], B_field: pint.Quantity):
    """Returns the atomic hamiltonian with no B field"""
    H = 0
    kets_dict = kets(atom, states)
    for (state, m), ket in kets_dict.items():
        H += (zeeman_shift(state.g, m, B_field, state._ureg).to("MHz")).magnitude * ket * ket.dag()
    return H


def H_int(
    atom: Atom,
    states: list[State],
    fields: dict[ElectricField, list[Transition]],
    _ureg: pint.UnitRegistry,
    rwa_shifts=True,
):
    """Returns the interaction hamiltonian in rotating frame of reference for a given atom and electric fields, with states decided by the user
    The idea is that if the electric fields do not form a closed roop, one can enter rotating frame of reference (interaction picture), such
    that fields either rotate at 0 frequency or the frequency much higher than the rabi frequency of any transition. This leaves the hamiltonian
    time independent, and so is relatively cheap and computationally easy to solve it using qutip, and hence the interesting case for it.

    I am aware that for some cases, like raman transitions my algorithm woulc claim that the rotating frame of reference is not possible, but it is.
    I dont have time right now to implement it now.

    Args:
        atom (Atom): Atom object
        states (list[State]): List of states to be included in the hamiltonian
        fields (dict[ElectricField, list[Transition]]): Dictionary of electric fields and transitions that are driven by them
        _ureg (pint.UnitRegistry): Unit registry

    Returns:
        qutip.Qobj: Hamiltonian in rotating frame of reference plus correction on the base hamiltonian
    """
    H = 0
    kets_dict = kets(atom, states)
    for field, transitions in fields.items():
        for transition in transitions:
            state_i = transition.state_i
            state_f = transition.state_f
            for mJ_i in state_i.sublevels:
                for mJ_f in state_f.sublevels:
                    ket_i = kets_dict[(state_i, mJ_i)]
                    ket_f = kets_dict[(state_f, mJ_f)]

                    Omega_ij = Rabi_Frequency(field, transition, mJ_i=mJ_i, mJ_f=mJ_f).to("MHz")
                    h = 1 / 2 * Omega_ij.magnitude * (ket_i * ket_f.dag())
                    H += h + h.dag()
        # transition_graph = H.full()
    if rwa_shifts:
        # TODO: perhaps this could be another dedicated function?
        RWA_shifts, nodes = find_rotating_frame(fields)
        for i, state in enumerate(nodes):
            for m in state.sublevels:
                ket = kets_dict[(state, m)]
                H += complex(RWA_shifts[i]) * ket * ket.dag()
    return H


def collapse_operators(atom: Atom, states: list[State], _ureg: pint.UnitRegistry, print_added_operators=False):
    # Technically I am not doing it fully right, as I dont take under account if the photons are exactly the same
    """Returns the atomic hamiltonian with no B field"""

    list_all_operators = []
    kets_dict = kets(atom, states)

    states_mJ = list(kets_dict.keys())

    for (state_i, mJ_i), (state_f, mJ_f) in combinations(states_mJ, 2):
        tr = atom.transition_between(state_i, state_f)
        if tr is not None:
            ket_i = kets_dict[(state_i, mJ_i)]
            ket_f = kets_dict[(state_f, mJ_f)]
            lo = sqrt_lindblad_operator(tr, mJ_i, mJ_f, _ureg).to("MHz**0.5")
            c_ij = complex(lo.magnitude) * (ket_i * ket_f.dag())
            # c_ij = np.sqrt(tr.Gamma.to("MHz").m) * (ket_i * ket_f.dag())
            list_all_operators.append(c_ij)
            if print_added_operators:
                print(
                    f"Added {(lo**2).to('MHz'):.3g~P} c operator from {(state_f.term, mJ_f)} to {(state_i.term, mJ_i)}"
                )

    return list_all_operators


def adiabatic_elimination(relevant_states, irrelevant_states, kets, H):
    """adiabatic_elimination

    refs:
    https://www.quantuminfo.physik.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaaajiobd pag 22
    https://link.springer.com/article/10.1140/epjp/i2014-14012-8

    Args:
        relevant_states (_type_): _description_
        irrelevant_states (_type_): _description_
        kets (_type_): _description_
        H (_type_): _description_

    Returns:
        _type_: _description_
    """
    _states = list(kets.keys())
    r_ix = [_states.index(s) for s in relevant_states]
    ir_ix = [_states.index(s) for s in irrelevant_states]
    nr = len(r_ix)

    w = qutip.Qobj(H[np.ix_(r_ix, r_ix)])
    omega = qutip.Qobj(H[np.ix_(r_ix, ir_ix)])
    delta = qutip.Qobj(H[np.ix_(ir_ix, ir_ix)])

    H_eff = w - omega * delta.inv() * omega.dag()
    rkets = {s: qutip.basis(nr, j) for j, s in enumerate(relevant_states)}
    return rkets, H_eff
