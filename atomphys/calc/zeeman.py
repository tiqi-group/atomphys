#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>
"""
Lande factors and Zeeman shift for fine and hyperfine states

References:
[1] Hans A. Bethe and Edwin E. Salpeter, Quantum Mechanics of One- and Two-Electron Atoms (SpringerVerlag, Berlin, 1957).
"""

from pint import Quantity

g_spin = 2.00231930436256


def g_lande_fine_LS(L, S, J, gL=1):
    """
    Calculates the Landé g-factor for fine structure states with LS coupling. [1]

    Args:
        L (float): Total orbital angular momentum.
        S (float): Total electron spin.
        J (float): Total angular momentum.
        gL (float, optional): Landé g-factor for the orbital angular momentum. Default is 1.

    Returns:
        float: Landé g-factor for the total angular momentum J.

    """
    if J == 0:
        return 0
    gJ = gL * (J * (J + 1) - S * (S + 1) + L * (L + 1)) / (2 * J * (J + 1)) + \
        g_spin * (J * (J + 1) + S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1))
    return gJ


def g_lande_hyperfine(J, I, F, gJ, gI=0):  # noqa: E741
    """
    Calculates the Landé g-factor for hyperfine structure. [1]

    Args:
        J (float): Total angular momentum.
        I (float): Nuclear spin.
        F (float): Total angular momentum (including electron and nuclear spins).
        gJ (float): Landé g-factor for the total angular momentum J.
        gI (float, optional): Landé g-factor for the nuclear spin. Default is 0.

    Returns:
        float: Landé g-factor for the total angular momentum F.

    """
    if F == 0:
        return 0
    gF = gJ * (F * (F + 1) - I * (I + 1) + J * (J + 1)) / (2 * F * (F + 1))
    if gI != 0:
        gF += gI * (F * (F + 1) + I * (I + 1) - J * (J + 1)) / (2 * F * (F + 1))
    return gF


def zeeman_shift(g, m, B: Quantity):
    mu = B._REGISTRY.bohr_magneton
    return mu * g * m * B.to('tesla', 'Gaussian')
