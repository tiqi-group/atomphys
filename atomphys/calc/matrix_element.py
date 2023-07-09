#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>
#
# TODO: find reference and document


import pint
from math import pi


def reduced_dipole_matrix_element(
        A: pint.Quantity, wavelength: pint.Quantity, _ureg: pint.UnitRegistry | None = None):
    _ureg = pint.get_application_registry() if _ureg is None else _ureg
    C = 3 * pi * _ureg.epsilon_0 * _ureg.hbar
    return (C * (wavelength / 2 / pi)**3 * A)**(1 / 2)


def reduced_quadrupole_matrix_element(
        A: pint.Quantity, wavelength: pint.Quantity, _ureg: pint.UnitRegistry | None = None):
    _ureg = pint.get_application_registry() if _ureg is None else _ureg
    beta = 1  # TODO: find the correct constant
    C = 4 * beta * pi * _ureg.epsilon_0 * _ureg.hbar
    return (C * (wavelength / 2 / pi)**5 * A)**(1 / 2)
