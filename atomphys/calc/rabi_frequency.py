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
from .matrix_element import dipole_matrix_element, quadrupole_matrix_element
import numpy as np

def dipole_Rabi_Frequency(
        E: pint.Quantity, A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry | None = None):

    return np.dot(E, dipole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg))*_ureg('e/hbar')

def quadrupole_Rabi_Frequency(E_gradient: pint.Quantity, A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry | None = None):
    qme = quadrupole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)
    Eg = E_gradient

    return 1/2*np.sum(Eg*qme)*_ureg('e/hbar')

        





