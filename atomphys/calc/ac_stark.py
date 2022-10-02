from math import cos
from math import pi as π

import pint

from atomphys.calc.wigner import ishalfint, isint
from sympy.physics.wigner import wigner_6j as w6j
from sympy.physics.quantum.cg import CG
from sympy import S
import numpy as np

def scalar(state, omega, mj):
    """Calculate the scalar polarizability of a state |Ψ⟩

    """
    ω = omega
    ħ = state._ureg.ħ
    j = state.J
    α = 0
    k=0

    for transition in state.up:
        #print(transition.matrix_element)
        d = transition.matrix_element
        other_state = transition.state_f
        jd = other_state.quantum_numbers['J']
        ω0 = transition.ω
        w6 = w6j(j, k, j, 1, jd, 1, prec=10)
        α = α + np.sqrt(2*k+1)*(-1)**(1+j+jd)*w6*d**2 * ((-1)**(k)/(ω0*ħ+ω*ħ) + 1/(ω0*ħ-ω*ħ))
        
    for transition in state.down:
        #print(transition.matrix_element)
        d = transition.matrix_element
        other_state = transition.state_f
        jd = other_state.quantum_numbers['J']
        ω0 = -transition.ω
        w6 = w6j(j, k, j, 1, jd, 1, prec=10)
        α = α + np.sqrt(2*k+1)*(-1)**(1+j+jd)*w6*d**2 * ((-1)**(k)/(ω0*ħ+ω*ħ) + 1/(ω0*ħ-ω*ħ))
        
    return CG(j, mj, k, 0, j, mj).doit()/np.sqrt(2*j+1), α.to_base_units()


def vector(state, omega, mj):
    """Calculate the scalar polarizability of a state |Ψ⟩

    """
    ω = omega
    ħ = state._ureg.ħ
    j = state.J
    α = 0
    k=1

    for transition in state.up:
        #print(transition.matrix_element)
        d = transition.matrix_element
        other_state = transition.state_f
        jd = other_state.quantum_numbers['J']
        ω0 = transition.ω
        w6 = w6j(j, k, j, 1, jd, 1, prec=10)
        α = α + np.sqrt(2*k+1)*(-1)**(1+j+jd)*w6*d**2 * ((-1)**(k)/(ω0*ħ+ω*ħ) + 1/(ω0*ħ-ω*ħ))
        
    for transition in state.down:
        #print(transition.matrix_element)
        d = transition.matrix_element
        other_state = transition.state_f
        jd = other_state.quantum_numbers['J']
        ω0 = -transition.ω
        w6 = w6j(j, k, j, 1, jd, 1, prec=10)
        α = α + np.sqrt(2*k+1)*(-1)**(1+j+jd)*w6*d**2 * ((-1)**(k)/(ω0*ħ+ω*ħ) + 1/(ω0*ħ-ω*ħ))
        
    return CG(j, mj, k, 0, j, mj).doit()/np.sqrt(2*j+1), α.to_base_units()




def tensor(state, omega, mj):
    """Calculate the scalar polarizability of a state |Ψ⟩
    """
    ω = omega
    ħ = state._ureg.ħ
    j = state.J
    α = 0
    k=0

    for transition in state.up:
        #print(transition.matrix_element)
        d = transition.matrix_element
        other_state = transition.state_f
        jd = other_state.quantum_numbers['J']
        ω0 = transition.ω
        w6 = w6j(j, k, j, 1, jd, 1, prec=10)
        α = α + np.sqrt(2*k+1)*(-1)**(1+j+jd)*w6*d**2 * ((-1)**(k)/(ω0*ħ+ω*ħ) + 1/(ω0*ħ-ω*ħ))
        
    for transition in state.down:
        #print(transition.matrix_element)
        d = transition.matrix_element
        other_state = transition.state_f
        jd = other_state.quantum_numbers['J']
        ω0 = -transition.ω
        w6 = w6j(j, k, j, 1, jd, 1, prec=10)
        α = α + np.sqrt(2*k+1)*(-1)**(1+j+jd)*w6*d**2 * ((-1)**(k)/(ω0*ħ+ω*ħ) + 1/(ω0*ħ-ω*ħ))
        
    return CG(j, mj, k, 0, j, mj).doit()/np.sqrt(2*j+1), α.to_base_units()




def total_ACshift(
    state,
    mJ,
    omega,
    eps,
    e_z,
    I,
    u, 
):
    """Calculate the polarizability of a state for a given field polarization

    Calculates the dynamical polarizability α(ω) by taking the sum of the
    scalar, vector, and tensor poarts according to

    \\[
        \\alpha(\\omega) = \\alpha^0(\\omega) +
        A \\cos(\\theta_k) \\frac{m_J}{J} \\alpha^1(\\omega) +
        \\frac{1}{2}\\left(3 \\cos^2(\\theta_p) - 1\\right)\\frac{3m_J^2 - J(J+1)}{J(2J-1)}\\alpha^2(\\omega)
    \\]

    Arguments:
        state (State): state |Ψ⟩ to calculate the polarizability
        mJ: Zeeman sublevel to calculate the vector and tensor
            polarizability. If `None` only calculates the scalar
            component. Must be an integer if J is an integer, or
            a half-integer if J is a half integer.
        omega: angular frequency of the field to calculate
            the dynamical polarizability. Defaults to zero,
            in which case it calculates the static polarizability.
            Must have units of `Hz`.
        A: degree of circular polarization, with ±1 being circular
            polarization and 0 being linear polarization
        theta_k: angle between wave vector and quantization axis
        theta_p: angle between polarization vector and quantization axis

    Returns:
        Polarizability. Has units of `C m / (V/m)`, or dipole moment / electric field.

    Raises:
        ZeroDivisionError
    """
    
    J = state.J
    if not (isint(mJ) or (ishalfint(mJ) and ishalfint(J) and not isint(J))):
        raise ValueError(
            "mJ must be an integer, or a half-integer if J is a half-integer"
        )
    
    prefix0, α0 = scalar(state, omega, mJ)
    prefix1, α1 = vector(state, omega, mJ)
    prefix2, α2 = tensor(state, omega, mJ)
    
    aeps = np.sqrt(abs(eps[0])**2 + abs(eps[1])**2 + abs(eps[2])**2)
    eps_n = (eps[0]/aeps, eps[1]/aeps, eps[2]/aeps)
    
    ae_z = np.sqrt(abs(e_z[0])**2 + abs(e_z[1])**2 + abs(e_z[2])**2)
    e_z_n = (e_z[0]/ae_z, e_z[1]/ae_z, e_z[2]/ae_z)
    
    c0 = 1/np.sqrt(3)
    c1 = 1j/np.sqrt(2) * np.dot(np.cross(eps, np.conjugate(eps)), e_z)
    c2 = 1/np.sqrt(6) * (3*abs(np.dot(eps, e_z))**2-1)
                            

    E_square = 2*I/(u.ε_0*u.c)    
    
    return E_square*(c0*prefix0*α0 + c1*prefix1*α1 + c2*prefix2*α2).to_base_units()
