from .rabi_frequency import Rabi_Frequency
import qutip
from ..atom import Atom
from ..state import State
import pint
from ..electric_field import ElectricField
from .util import find_rotating_frame
from ..transition import Transition
from .lindblad_operators import sqrt_lindblad_operator



def H0(atom: Atom, states: list[State], _ureg: pint.UnitRegistry):
    """Returns the atomic hamiltonian with no B field"""
    H = 0
    

    states.sort(key=lambda state: state.energy)
    states_with_mJ = []
    for state in states:
        J = state.quantum_numbers['J']
        for i in range(int(2*J+1)):
            mJ = -J+i
            states_with_mJ.append((mJ, state))
    
    total_number_of_states = len(states_with_mJ)
    ket = [qutip.basis(total_number_of_states, i) for i in range(total_number_of_states)]

    for i, state in enumerate(states_with_mJ):
        H += (state[1].angular_frequency.to('MHz')).magnitude * ket[i] * ket[i].dag()

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

    
    states.sort(key=lambda state: state.energy)
    states_with_mJ = []
    dict_states_with_mJ = {}
    acc=0
    for state in states:
        J = state.quantum_numbers['J']
        for i in range(int(2*J+1)):
            mJ = -J+i
            states_with_mJ.append((mJ, state))
            dict_states_with_mJ[(mJ, state)] = acc
            acc+=1
    
    total_number_of_states = len(states_with_mJ)
    ket = [qutip.basis(total_number_of_states, i) for i in range(total_number_of_states)]

    H = 0
    for c, field in enumerate(fields.keys()):
        for transition in fields[field]:
            state_i = transition.state_i
            state_f = transition.state_f
            J_i = state_i.quantum_numbers['J']
            J_f = state_f.quantum_numbers['J']
            for i in range(int(2*J_i+1)):
                mJ_i = -J_i+i
                for f in range(int(2*J_f+1)):
                    mJ_f = -J_f+f
                    index_i = dict_states_with_mJ[(mJ_i, state_i)]
                    index_f = dict_states_with_mJ[(mJ_f, state_f)]

                    Omega_ij = Rabi_Frequency(field, transition, mJ_i=mJ_i, mJ_f=mJ_f).to('MHz')
                    H += -1/2*complex(Omega_ij.magnitude)*(ket[index_i]*ket[index_f].dag())
                    H += -1/2*complex(Omega_ij.magnitude)*(ket[index_f]*ket[index_i].dag())
        transition_graph = H.full()
    
    RWA_shifts, nodes = find_rotating_frame(fields)
    for i, state in enumerate(nodes):
        J = state.quantum_numbers['J']
        for j in range(int(2*J+1)):
            mJ = -J+j
            index = dict_states_with_mJ[(mJ, state)]
            H += complex(RWA_shifts[i]) * ket[index] * ket[index].dag()
        
    return H


def collapse_operators(atom: Atom, states: list[State], _ureg: pint.UnitRegistry, print_added_operators=False):
    #Technically I am not doing it fully right, as I dont take under account if the photons are exactly the same
    
    """Returns the atomic hamiltonian with no B field"""
    H = 0
    
    list_all_operators = []
    states.sort(key=lambda state: state.energy)
    states_with_mJ = []
    dict_states_with_mJ = {}
    acc=0
    for state in states:
        J = state.quantum_numbers['J']
        for i in range(int(2*J+1)):
            mJ = -J+i
            states_with_mJ.append((mJ, state))
            dict_states_with_mJ[(mJ, state)] = acc
            acc+=1
    
    total_number_of_states = len(states_with_mJ)
    ket = [qutip.basis(total_number_of_states, i) for i in range(total_number_of_states)]
    
    for i in range(len(states_with_mJ)):
        mJ_f, state_f = states_with_mJ[-(i+1)]
        for j in range(len(states_with_mJ)):
            if j>i:
                mJ_i, state_i = states_with_mJ[-(j+1)]
                tr = atom.transition_between(state_i, state_f)
                if tr != None:
                    index_i = dict_states_with_mJ[(mJ_i, state_i)]
                    index_f = dict_states_with_mJ[(mJ_f, state_f)]
                    lo = sqrt_lindblad_operator(tr, mJ_i, mJ_f, _ureg)
                    c_ij = complex(lo.magnitude) * (ket[index_i]*ket[index_f].dag())
                    list_all_operators.append(c_ij)
                    if print_added_operators:
                        print(f'Added {complex(lo.magnitude)**2} MHz c operator from {(mJ_f, state_f.term)} to {(mJ_i, state_i.term)}')
                    
    return list_all_operators