#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>
"""
Hyperfine energy shift

References:
[1] Alan Corney, Atomic and Laser Spectroscopy (Oxford, 1977)
"""

from pint import Quantity


def hyperfine_shift(J: float, I: float, F: float, Ahf: Quantity, Bhf: Quantity = 0) -> Quantity:
    """
    Calculates the hyperfine energy shift for an atomic or molecular fine structure state [1].

    Args:
        J (float): Fine structure angular momentum.
        I (float): Nuclear spin.
        F (float): Total angular momentum.
        Ahf (Quantity): Hyperfine magnetic dipole constant.
        Bhf (Quantity, optional): Hyperfine electric quadrupole constant. Default is 0.

    Returns:
        Quantity: Energy shift of the hyperfine state.

    """
    K = F * (F + 1) - I * (I + 1) - J * (J + 1)
    Ehf = 1 / 2 * Ahf * K
    if Bhf != 0:
        den = 4 * I * (2 * I - 1) * J * (2 * J - 1)
        if den != 0:
            Ehf += Bhf * (3 / 2 * K * (K + 1) - 2 * I * (I + 1) * J * (J + 1)) / den
    return Ehf
