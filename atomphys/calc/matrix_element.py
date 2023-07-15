#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>
#
# # Reference for the calculation is 
#         - Thesis of Christoph Fisher, page 34.
#         - Quantum dynamics of cold trappedions, with application to quantum computation by Daniel F. V. James
#             - eq (5.9.) 


import pint
from math import pi
from sympy.physics.wigner import wigner_6j as w6j
from sympy.physics.wigner import wigner_3j as w3j
from sympy.physics.quantum.cg import CG
from .util import spherical_basis_second_rank_tensor, spherical_basis_vector


def reduced_dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry):

    return (_ureg('3/(4*c*alpha)')*(2*J_f+1)/k**3 * A)**(1/2)

def reduced_electric_dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry):
    
    return _ureg('e') * reduced_dipole_matrix_element(A, k, J_f, _ureg)

def reduced_quadrupole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry):
 
    return (_ureg('15/(c*alpha)')*(2*J_f+1)/k**5 * A)**(1/2)

def reduced_electric_quadrupole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry):

    return _ureg('e') * reduced_quadrupole_matrix_element(A, k, J_f, _ureg)

def dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    
    dipole_matrix_element = 0
    d = reduced_dipole_matrix_element(A, k, J_f, _ureg)
    for q in range(-1, 2):
        w3j_coeff = w3j(J_i, 1, J_f, -mJ_i, q, mJ_f)
        dipole_matrix_element += w3j_coeff * spherical_basis_vector(q) * d
    return dipole_matrix_element

def quadrupole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    
    quadrupole_matrix_element = 0
    reduced_quadrupole_matrix_element = reduced_quadrupole_matrix_element(A, k, J_f, _ureg)
    for q in range(-2, 3):
        w3j_coeff = w3j(J_i, 2, J_f, -mJ_i, q, mJ_f)
        quadrupole_matrix_element += w3j_coeff * spherical_basis_second_rank_tensor(q) * reduced_quadrupole_matrix_element
    return quadrupole_matrix_element

def electric_dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    
    return (_ureg('e') * dipole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg))

def electric_quadrupole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):

    return (_ureg('e') * quadrupole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)).to('e*a0**2')
        
