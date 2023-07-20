#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>
#
# Reference for this calculation is James 1998: Quantum dynamics of cold trappedions, with application to quantum computation
# https://arxiv.org/abs/quant-ph/9702053
# 
# Whilst reading it, be careful about the definitions of j' and j. It is not consistent throughout the paper. 
# I took j' to be the lower level states, as it agrees with other sources and makes sense intuitively 
# -> it increases spontaneous decay rate the more possibilities of decays to the lower states there are.


import pint
from math import pi
from sympy.physics.wigner import wigner_6j as w6j
from sympy.physics.wigner import wigner_3j as w3j
from sympy.physics.quantum.cg import CG
from .util import spherical_basis_second_rank_tensor, spherical_basis_vector


def reduced_dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, _ureg: pint.UnitRegistry):
    """
    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        _ureg: Unit registry

    Returns:
        Reduced dipole matrix element [a0]
    """

    return (_ureg('3/(4*c*alpha)')*(2*J_i+1)/k**3 * A)**(1/2)

def reduced_electric_dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, _ureg: pint.UnitRegistry):
    """
    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        _ureg: Unit registry

    Returns:
        Reduced electric dipole matrix element [e a0]
    """
    
    return _ureg('e') * reduced_dipole_matrix_element(A, k, J_i, _ureg)

def reduced_quadrupole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, _ureg: pint.UnitRegistry):
    """
    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        _ureg: Unit registry

    Returns:
        Reduced dipole matrix element [a0**2]
    """
 
    return (_ureg('15/(c*alpha)')*(2*J_i+1)/k**5 * A)**(1/2)

def reduced_electric_quadrupole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry):
    """
    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        _ureg: Unit registry

    Returns:
        Reduced dipole matrix element [e a0**2]
    """

    return _ureg('e') * reduced_quadrupole_matrix_element(A, k, J_f, _ureg)

def dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    """
    Dipole Matrix Element is reduced dipole matrix element decorated with clebsh gordan coefficients (wigner 3j sybol), 
    hence expressing the coupling between different magnetic sublevels.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry

    Returns:
        Dipole matrix element [a0]
    """
    
    dipole_matrix_element = 0
    d = reduced_dipole_matrix_element(A, k, J_i, _ureg)
    for q in range(-1, 2):
        w3j_coeff = w3j(J_f, 1, J_i, -mJ_f, q, mJ_i)
        dipole_matrix_element += w3j_coeff * spherical_basis_vector(q) * d
    return dipole_matrix_element.to('a0').magnitude*_ureg('a0')

def quadrupole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    """
    Quadrupole Matrix Element is reduced quadrupole matrix element decorated with clebsh gordan coefficients (wigner 3j sybol), 
    hence expressing the coupling between different magnetic sublevels.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry

    Returns:
        Quadrupole matrix element [a0**2]
    """
    
    quadrupole_matrix_element = 0
    red_q = reduced_quadrupole_matrix_element(A, k, J_i, _ureg)
    for q in range(-2, 3):
        w3j_coeff = float(w3j(J_f, 2, J_i, -mJ_f, q, mJ_i))
        sbt = spherical_basis_second_rank_tensor(q)
        quadrupole_matrix_element += w3j_coeff * sbt * red_q
    return abs(quadrupole_matrix_element.to('a0**2').magnitude)*_ureg('a0**2')

def electric_dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    """
    See dipole matrix element for details how is it different from reduced matrix elements.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry
        
    Returns:
        Electric dipole matrix element [e a0]
    """
    
    return (_ureg('e') * dipole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg))

def electric_quadrupole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    """
    See quadrupole matrix element for details how is it different from reduced matrix elements.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry
        
    Returns:
        Electric quadrupole matrix element [e a0]
    """

    return (_ureg('e') * quadrupole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)).to('e*a0**2')
        
