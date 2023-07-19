from .rabi_frequency import Rabi_Frequency
import qutip
from ..atom import Atom
from ..state import State
import pint
from ..electric_field import ElectricField
from .util import find_rotating_frame
from ..transition import Transition



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
    truncated_atom = atom.remove_all_but_states(states, copy=True)
    ket = [qutip.basis(total_number_of_states, i) for i in range(total_number_of_states)]

    for i, state in enumerate(states_with_mJ):
        H += (state[1].angular_frequency.to('MHz')).magnitude * ket[i] * ket[i].dag()

    return H

def H_int_graph(atom: Atom, states: list[State], fields: dict[ElectricField, list[Transition]], _ureg: pint.UnitRegistry):
    
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
    





def sqrt_Γ_ij(state_i, state_f, transition):
    ΔL=abs(state_f.L-state_i.L)
    prefactor = np.sqrt((2*state_f.J + 1))
    summation = 0
    for q in range(-ΔL, ΔL+1): summation += wigner_3j(state_i.J, ΔL, state_f.J, -state_i.m_j, q, state_f.m_j)
    return summation * prefactor * np.sqrt(transition.Γ)


def Linblat_operators(kets, states, lasers):
    #Technically I am not doing it fully right, as I dont take under account if the photons are exactly the same
    
    L_op = []
    for c, laser in enumerate(lasers.keys()):
        states_i, states_f = lasers[laser]
        
        for i in states_f:
            for j in states_i:
                transitions = [i for i in filter(lambda x: (x.state_f.name == states[j].name and x.state_f.configuration==states[j].configuration), states[i].transitions)]
                if transitions != []:
                    transition = transitions[0]
                    if j>i:
                        _sqrt_Γ_ij = (sqrt_Γ_ij(states[i], states[j], transition).to('MHz**0.5')).magnitude
                        c_ij = float(abs(_sqrt_Γ_ij)) * (ket[i]*ket[j].dag())
                    
                        L_op.append(c_ij)
    return L_op

def H_interaction(kets, states, lasers, v):
    H = 0
    for c, laser in enumerate(lasers.keys()):
        states_i, states_f = lasers[laser]
        
        for i in states_f:
            for j in states_i:
                transitions = [i for i in filter(lambda x: (x.state_f.name == states[j].name and x.state_f.configuration==states[j].configuration), states[i].transitions)]
                if transitions != []:
                    transition = transitions[0]
                    Ω_ij = (laser.Ω_ij(states[i], states[j], transition).to('MHz'))
                    #H += (-u('hbar')/2*Ω_ij)*(ket[i]*ket[j].dag() + ket[i].dag()*ket[j])
                    #H += (-u('hbar')/2*Ω_ij)*(ket[i]*ket[j].dag())
                    #H += (-u('hbar')/2*Ω_ij)*(ket[j]*ket[i].dag())
                    H += -1/2*float(abs(Ω_ij.magnitude))*(ket[i]*ket[j].dag())
                    H += -1/2*float(abs(Ω_ij.magnitude))*(ket[j]*ket[i].dag())
        


        min_val = min(states_i)
        for i in range(len(kets)):
            if i>=min_val:
                dH = -((laser.ω-np.dot(laser.k,v)).to('MHz')).magnitude
                H += dH*(ket[i]*ket[i].dag())
            
    H_tot = H
    
    return H_tot
     