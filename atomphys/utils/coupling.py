#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>
#
# https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy-2

from enum import Enum
from .quantum_numbers import QuantumNumbers


class Coupling(Enum):
    LS = "LS"  # Russell-Saunders coupling
    jj = "jj"
    KS = "KS"  # pair coupling


def get_coupling(quantum_numbers: QuantumNumbers) -> Coupling:
    if quantum_numbers.L is not None and quantum_numbers.S is not None:
        return Coupling.LS
    if quantum_numbers.J1 is not None and quantum_numbers.J2 is not None:
        return Coupling.jj
    if quantum_numbers.K is not None and quantum_numbers.S2 is not None:
        return Coupling.KS
