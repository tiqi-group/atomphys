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

from itertools import combinations


def kets(atom: Atom, states: list[State]):
    states.sort(key=lambda state: state.energy)
    states_with_mJ = []
    for state in states:
        for m in state.sublevels:
            states_with_mJ.append((m, state))

    total_number_of_states = len(states_with_mJ)
    kets_dict = {states_with_mJ[i]: qutip.basis(total_number_of_states, i) for i in range(total_number_of_states)}
    return kets_dict


def H0(atom: Atom, states: list[State], _ureg: pint.UnitRegistry):
    """Returns the atomic hamiltonian with no B field"""
    H = 0
    kets_dict = kets(atom, states)
    for state_mJ, ket in kets_dict.items():
        H += (state_mJ[1].angular_frequency.to("MHz")).magnitude * ket * ket.dag()
    return H


def H_zeeman(atom: Atom, states: list[State], B_field: pint.Quantity):
    """Returns the atomic hamiltonian with no B field"""
    H = 0
    kets_dict = kets(atom, states)
    for state_mJ, ket in kets_dict.items():
        m, s = state_mJ
        H += (zeeman_shift(s.g, m, B_field, s._ureg).to("MHz")).magnitude * ket * ket.dag()

    return H


def H_int(atom: Atom, states: list[State], fields: dict[ElectricField, list[Transition]], _ureg: pint.UnitRegistry):
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
                    ket_i = kets_dict[(mJ_i, state_i)]
                    ket_f = kets_dict[(mJ_f, state_f)]

                    Omega_ij = Rabi_Frequency(field, transition, mJ_i=mJ_i, mJ_f=mJ_f).to("MHz")
                    h = 1 / 2 * complex(abs(Omega_ij).magnitude) * (ket_i * ket_f.dag())
                    H += h + h.dag()
        # transition_graph = H.full()

    RWA_shifts, nodes = find_rotating_frame(fields)
    for i, state in enumerate(nodes):
        for mJ in state.sublevels:
            ket = kets_dict[(mJ, state)]
            H += complex(RWA_shifts[i]) * ket * ket.dag()

    return H


def collapse_operators(atom: Atom, states: list[State], _ureg: pint.UnitRegistry, print_added_operators=False):
    # Technically I am not doing it fully right, as I dont take under account if the photons are exactly the same
    """Returns the atomic hamiltonian with no B field"""

    list_all_operators = []
    kets_dict = kets(atom, states)

    states_mJ = list(kets_dict.keys())

    for (mJ_i, state_i), (mJ_f, state_f) in combinations(states_mJ, 2):
        tr = atom.transition_between(state_i, state_f)
        if tr is not None:
            ket_i = kets_dict[(mJ_i, state_i)]
            ket_f = kets_dict[(mJ_f, state_f)]
            lo = sqrt_lindblad_operator(tr, mJ_i, mJ_f, _ureg)
            c_ij = complex(lo.magnitude) * (ket_i * ket_f.dag())
            # c_ij = np.sqrt(tr.Gamma.to("MHz").m) * (ket_i * ket_f.dag())
            list_all_operators.append(c_ij)
            if print_added_operators:
                print(
                    f"Added {complex(lo.magnitude)**2} MHz c operator from {(mJ_f, state_f.term)} to {(mJ_i, state_i.term)}"
                )

    return list_all_operators
