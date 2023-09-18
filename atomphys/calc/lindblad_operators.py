#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 13/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch> & Wojciech Adamczyk <wadamczyk@phys.ethz.ch> 
#
# Reference for this calculation is James 1998: Quantum dynamics of cold trappedions, with application to quantum computation
# https://arxiv.org/abs/quant-ph/9702053
# Also very usefull whilst coding this up was 
# -> PhD Thesis of Lindenfelser, Frieder from 2017 (ETH Zürich) - Trapped Ion Quantum Informaiton Group
# -> Master Thesis of Beck, Gillenhall from 2020 (ETH Zürich) - Trapped Ion Quantum Informaiton Group

import pint
from ..transition import Transition, TransitionType
from sympy.physics.wigner import wigner_3j as w3j
import numpy as np
from .matrix_element import reduced_dipole_matrix_element, reduced_quadrupole_matrix_element

def sqrt_lindblad_operator(transition: Transition, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    """
    Args:
        transition: Transition object
        mJ_i: mJ of the lower state
        mJ_f: mJ of the upper state
        _ureg: Unit registry
    
    Returns:
        sqrt of the Linewidth between two zeeman states [MHz**0.5]
    """

    state_i = transition.state_i
    state_f = transition.state_f
    J_f = state_f.quantum_numbers['J']
    J_i = state_i.quantum_numbers['J']
    k = transition.k


    if transition.type == TransitionType.E1:
        rd_sq = reduced_dipole_matrix_element(A = transition.A, k=k, J_i=J_i, _ureg=_ureg)**2
        pre = _ureg('4*c*alpha/3')*k**3
        lo = 0*_ureg('MHz')
        for q in [-1,0,1]: lo += pre* rd_sq * w3j(J_f, 1, J_i, -mJ_f, q, mJ_i)**2
    elif transition.type == TransitionType.E2:
        rd_sq = reduced_quadrupole_matrix_element(transition.A, k, J_f, _ureg)**2
        pre = _ureg('c*alpha/15')*k**5
        lo = 0*_ureg('MHz')
        for q in [-2, -1,0,1, 2]: lo += pre* rd_sq * w3j(J_f, 2, J_i, -mJ_f, q, mJ_i)**2
    else:
        raise NotImplementedError(f"Transition type {transition.type} not implemented")

    return np.sqrt(complex((lo.to('MHz')).magnitude))*_ureg('MHz**0.5')