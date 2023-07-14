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


def dipole_bare_Rabi_Frequency(transition, I):


