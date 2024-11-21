#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>


def couple_angular_momenta(L1: float, L2: float) -> list[float]:
    return [abs(L2 - L1) + m for m in range(int(2 * min(L2, L1)) + 1)]


def magnetic_sublevels(L: float) -> list[float]:
    return [-L + m for m in range(int(2 * L) + 1)]
