#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Reference for this calculation is James 1998: Quantum dynamics of cold trappedions, with application to quantum computation
# https://arxiv.org/abs/quant-ph/9702053
# Also very usefull whilst coding this up was
# -> PhD Thesis of Lindenfelser, Frieder from 2017 (ETH Zürich) - Trapped Ion Quantum Informaiton Group
# -> Master Thesis of Beck, Gillenhall from 2020 (ETH Zürich) - Trapped Ion Quantum Informaiton Group

import pint
import numpy as np
from atomphys.transition import Transition
from atomphys.calc.linewidths import transition_specific_linewidth


def sqrt_lindblad_operator(
    transition: Transition, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry
):
    """
    Args:
        transition: Transition object
        mJ_i: mJ of the lower state
        mJ_f: mJ of the upper state
        _ureg: Unit registry

    Returns:
        sqrt of the Linewidth between two zeeman states [MHz**0.5]
    """

    lo = transition_specific_linewidth(transition, mJ_i, mJ_f, _ureg)

    return np.sqrt(complex((lo.to("MHz")).magnitude)) * _ureg("MHz**0.5")
