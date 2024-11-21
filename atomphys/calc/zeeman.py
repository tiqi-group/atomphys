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

import numpy as np
import pint

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
    gJ = gL * (J * (J + 1) - S * (S + 1) + L * (L + 1)) / (2 * J * (J + 1)) + g_spin * (
        J * (J + 1) + S * (S + 1) - L * (L + 1)
    ) / (2 * J * (J + 1))
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


def zeeman_shift(g, m, B: pint.Quantity, _ureg: pint.UnitRegistry | None = None):
    _ureg = pint.get_application_registry() if _ureg is None else _ureg
    mu = _ureg.bohr_magneton
    return mu * g * m * B.to("tesla", "Gaussian")


def quadratic_clock_shift(gJ, dE_hf, B, _ureg=None):
    _ureg = pint.get_application_registry() if _ureg is None else _ureg
    mu = _ureg.bohr_magneton
    return gJ**2 * mu**2 * B.to("tesla", "Gaussian") ** 2 / 2 / dE_hf


def field_sensitivity(g, m, _ureg=None):
    _ureg = pint.get_application_registry() if _ureg is None else _ureg
    mu = (1 * _ureg.bohr_magneton).to("J/gauss", "Gaussian")
    return (mu * g * m / _ureg.planck_constant).to("MHz/gauss")


def breit_rabi(B, Ahf, I, mJ, mI, gJ, gI=0, _ureg=None):
    """F = 1, F = 2 Na B-field shift from Breit rabi formula.
    freq = breit_rabi(B, F, mF)

    Args:
        B : B field (Gauss)
        F : Hyperfine level F. Must be 1 or 2 (raises AssertionError instead)
        mF: Zeeman sublevel. Must be abs(mF) <= F (raises AssertionError)

    Returns:
        freq: frequency shift (MHz) wrt 3S_{1/2} fine structure energy.
              Note that mF = 0 levels have nonzero energy in this way.
    """
    # B (Gauss) -> freq (MHz)
    assert np.abs(mI) <= I
    assert mJ in (-0.5, 0.5)
    _ureg = pint.get_application_registry() if _ureg is None else _ureg

    Ehf = Ahf * (I + 1 / 2)
    mu = _ureg.bohr_magneton
    B = B.to("tesla", "Gaussian")
    sign = np.sign(mJ)
    if np.abs(mI) == I and np.sign(mI) == sign:  # stretched states
        dE = Ehf * I / (2 * I + 1) + sign * 0.5 * (gJ + 2 * I * gI) * mu * B
    else:
        mF = mI + mJ
        x = (gJ - gI) * mu * B / Ehf
        dE = -Ehf / 2 / (2 * I + 1) + sign * Ehf / 2 * np.sqrt(
            1 + 4 * mF * x / (2 * I + 1) + x**2
        )
        if gI != 0:
            dE += gI * mu * mF * B
    return dE


def breit_rabi_dEdB(B, Ahf, I, mJ, mI, gJ, gI=0, _ureg=None):
    """F = 1, F = 2 Na B-field shift from Breit rabi formula.
    freq = breit_rabi(B, F, mF)

    Args:
        B : B field (Gauss)
        F : Hyperfine level F. Must be 1 or 2 (raises AssertionError instead)
        mF: Zeeman sublevel. Must be abs(mF) <= F (raises AssertionError)

    Returns:
        freq: frequency shift (MHz) wrt 3S_{1/2} fine structure energy.
              Note that mF = 0 levels have nonzero energy in this way.
    """
    # B (Gauss) -> freq (MHz)
    assert np.abs(mI) <= I
    assert mJ in (-0.5, 0.5)
    _ureg = pint.get_application_registry() if _ureg is None else _ureg

    Ehf = Ahf * (I + 1 / 2)
    mu = _ureg.bohr_magneton
    B = B.to("tesla", "Gaussian")
    sign = np.sign(mJ)
    if np.abs(mI) == I and np.sign(mI) == sign:  # stretched states
        dEdB = sign * 0.5 * (gJ + 2 * I * gI) * mu * np.ones_like(B)
    else:
        mF = mI + mJ
        x = (gJ - gI) * mu * B / Ehf
        dEdB = (
            sign
            / 4
            * (2 * x + 4 * mF / (2 * I + 1))
            / np.sqrt(1 + 4 * mF * x / (2 * I + 1) + x**2)
            * (gJ - gI)
            * mu
        )
        if gI != 0:
            dEdB += gI * mu * mF
    return dEdB
