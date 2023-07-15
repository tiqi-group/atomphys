#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 13/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch> & Wojciech Adamczyk <wadamczyk@phys.ethz.ch> 
#
# TODO: find reference and document
#
"""
 \begin{aligned}
\Omega_{i j} & =\left\langle i\left|\frac{e}{\hbar} \dot{\vec{r}} \vec{A}_{0, l} e^{-i \vec{k}_l \vec{r}}\right| j\right\rangle \\
& =-i \omega_l\left\langle i\left|\vec{r} \vec{A}_{0, l}\right| j\right\rangle-\omega_l \frac{1}{2}\left\langle i\left|\left(\vec{r} \vec{k}_l\right)\left(\vec{r} \vec{A}_{0, l}\right)\right| j\right\rangle+\cdots
\end{aligned}
"""

from sympy.physics.wigner import wigner_6j as w6j
from sympy.physics.wigner import wigner_3j as w3j
import pint

def dipole_matrix_element(
        A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry | None = None):
    
    dipole_matrix_element = 0
    reduced_dipole_matrix_element = reduced_electric_dipole_matrix_element(A, k, J_f, _ureg)
    for q in range(-1, 2):
        w3j_coeff = w3j(J_i, 1, J_f, -mJ_i, q, mJ_f)
        dipole_matrix_element += w3j_coeff * spherical_basis_vector(q) * reduced_dipole_matrix_element
    return dipole_matrix_element


def dipole_Rabi_Frequency(
        dipole_matrix_element: pint.Quantity, electric_field: ElectricField, _ureg: pint.UnitRegistry | None = None):
    
    _ureg = pint.get_application_registry() if _ureg is None else _ureg

    return electric_field.amplitude * dipole_matrix_element

def dipole_Rabi_Frequency(
        transition: Transition, electric_field: ElectricField, _ureg: pint.UnitRegistry | None = None):
    
    _ureg = pint.get_application_registry() if _ureg is None else _ureg

    return electric_field.amplitude * dipole_matrix_element(


