#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>
#
# https://en.wikipedia.org/wiki/Selection_rule#Summary_table

from enum import Enum
from atomphys.quantum_numbers import QuantumNumbers


class TransitionType(Enum):
    E1 = "E1"    # electric dipole
    M1 = "M1"    # magnetic dipole
    E2 = "E2"    # electric quadrupole
    M2 = "M2"    # magnetic quadrupole
    # E3 = "E3"    # electric octupole
    NONE = None  # Not found


def get_transition_type_LS(qn1: QuantumNumbers, qn2: QuantumNumbers) -> TransitionType:
    deltaL = qn1.L - qn2.L
    parity = qn1.parity * qn2.parity
    if abs(deltaL) > 2:
        return TransitionType.NONE
    elif abs(deltaL) == 2:
        return TransitionType.E2
    elif abs(deltaL) == 1:
        return TransitionType.E1
    elif parity == 1:
        return TransitionType.M1
    else:
        return TransitionType.M2


# def is_electric_dipole(J1: float, J2: float, parity1: int, parity2: int):
#     parity = parity1 * parity2
#     deltaJ = J2 - J1
#     forbidden = J1 == 0 and J2 == 0
#     return abs(deltaJ) <= 1 and not forbidden and parity == -1


# def is_magnetic_dipole(J1: float, J2: float, parity1: int, parity2: int):
#     parity = parity1 * parity2
#     deltaJ = J2 - J1
#     forbidden = J1 == 0 and J2 == 0
#     return abs(deltaJ) <= 1 and not forbidden and parity == 1


# def is_electric_quadrupole(J1: float, J2: float, parity1: int, parity2: int):
#     parity = parity1 * parity2
#     deltaJ = J2 - J1
#     forbidden = (J1 == 0 and J2 == 0) or \
#         (J1 == 0 and J2 == 1) or \
#         (J1 == 0.5 and J2 == 0.5)
#     return abs(deltaJ) <= 2 and not forbidden and parity == 1
