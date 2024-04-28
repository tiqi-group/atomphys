#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 09/2023
# Authors:
#   Wojciech Adamczyk <wadamczyk@phys.ethz.ch>
#   Carmelo Mordini <cmordini@phys.ethz.ch>

from itertools import combinations
import qutip
import pint
from atomphys.atom import Atom
from atomphys.state import State
from atomphys.electric_field import ElectricField
from atomphys.transition import Transition
from atomphys.calc.rabi_frequency import Rabi_Frequency
from atomphys.calc.zeeman import zeeman_shift
from atomphys.calc.lindblad_operators import sqrt_lindblad_operator


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
    kets_dict = {
        states_m[i]: qutip.basis(total_number_of_states, i)
        for i in range(total_number_of_states)
    }
    return kets_dict


def zero_hamiltonian(states):
    n_states = sum(len(s.sublevels) for s in states)
    return qutip.qzero(n_states)


def H0(states: list[State], _ureg: pint.UnitRegistry):
    """ "
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


def H_detuning(states, dets_vector, _ureg: pint.UnitRegistry):
    total = 0
    for state in states:
        J = state.quantum_numbers["J"]
        for i in range(int(2 * J) + 1):
            total += 1

    counter = 0

    ket = [qutip.basis(total, i) for i in range(total)]

    H = 0

    for s, state in enumerate(states):
        J = state.quantum_numbers["J"]
        for i in range(int(2 * J) + 1):
            for si in range(s + 1):
                H += (
                    complex((dets_vector[si] * _ureg("_2pi")).to("MHz").magnitude)
                    * ket[counter]
                    * ket[counter].dag()
                )
            counter += 1
    return H


def H_zeeman(states: list[State], B_field: pint.Quantity):
    """
    Returns extra zeeman shift due to the B field. This is not the full hamiltonian, but rather a correction to the H0 hamiltonian.

    H_zeeman = sum_i (gJ_i mJ_i mu_b B |i><i|)
    """

    H = 0
    kets_dict = kets(states)
    for (state, m), ket in kets_dict.items():
        H += (
            (zeeman_shift(state.g, m, B_field, state._ureg).to("MHz")).magnitude
            * ket
            * ket.dag()
        )
    return H


def H_int_offdiagonal(
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
    H = 0
    kets_dict = kets(states)
    for field, transitions in fields.items():
        for transition in transitions:
            state_i = transition.state_i
            state_f = transition.state_f
            for mJ_i in state_i.sublevels:
                for mJ_f in state_f.sublevels:
                    ket_i = kets_dict[(state_i, mJ_i)]
                    ket_f = kets_dict[(state_f, mJ_f)]

                    Omega_ij = Rabi_Frequency(
                        field, transition, mJ_i=mJ_i, mJ_f=mJ_f
                    ).to("MHz")
                    h = 1 / 2 * Omega_ij.magnitude * (ket_i * ket_f.dag())
                    H += h + h.dag()
    return H


def collapse_operators(
    atom: Atom,
    states: list[State],
    _ureg: pint.UnitRegistry,
    print_added_operators=False,
):
    # Technically I am not doing it fully right, as I dont take under account if the photons are exactly the same
    """Returns the atomic hamiltonian with no B field"""

    list_all_operators = []
    kets_dict = kets(states)

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
