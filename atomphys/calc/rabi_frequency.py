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
from .matrix_element import (
    dipole_matrix_element, quadrupole_matrix_element,
    reduced_dipole_matrix_element, w3j
)
import numpy as np
from ..electric_field import ElectricField
from ..transition import Transition, TransitionType


def dipole_Rabi_Frequency(E_field: pint.Quantity, A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    """
    Function will only work for a dipole transitions

    Args:
        E_field: Electric field strength [V/m]
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry

    Returns:
        Rabi frequency for a dipole transition [2pi*MHz]
    """
    d = dipole_matrix_element(A=A, k=k, J_i=J_i, J_f=J_f, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
    return (np.dot(E_field, d) * _ureg('e/hbar')).to('MHz')


def electric_field_from_pi_dipole_rabi_frequency(
        rabi_frequency: pint.Quantity,
        A: pint.Quantity, k: pint.Quantity,
        J_i: float, _ureg: pint.UnitRegistry):
    """electric_field_from_pi_dipole_rabi_frequency

    Return the electric field giving a requested rabi frequency on a pi transition
    (J_i, 0) -> (J_i + 1, 0)
    Don't forget to include 2*pi in the rabi_frequency!

    Args:
        rabi_frequency: Rabi frequency [2pi * MHz]
        A: Einstein coefficient [1/s]
        k: Transition wavenumber [1/m]
        J_i: Angular momentum of the lower state
        _ureg: Unit registry. Defaults to None.

    Returns:
        _type_: _description_
    """
    J_f = J_i + 1
    d = reduced_dipole_matrix_element(A=A, k=k, J_f=J_f, _ureg=_ureg)
    w = float(w3j(J_f, 1, J_i, 0, 0, 0))
    return (rabi_frequency / abs(d * w) * _ureg('hbar / e')).to('V/m')


def quadrupole_Rabi_Frequency(E_gradient: pint.Quantity, A: pint.Quantity, k: pint.Quantity, J_i: float, J_f: float, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry):
    """
    Funciton will only work for a quadrupole transitions

    Args:
        E_gradient: Electric field gradient [V/m**2]
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry

    Returns:
        Rabi frequency for a quadrupole transition [2pi*Hz]
    """
    qme = quadrupole_matrix_element(A=A, k=k, J_i=J_i, J_f=J_f, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
    return (1 / 2 * np.sum(E_gradient * qme) * _ureg('e/hbar')).to('MHz')


def Rabi_Frequency(E_field: ElectricField, transition: Transition, mJ_i: float, mJ_f: float, _ureg: pint.UnitRegistry | None = None):
    """
    Args:
        E_field: ElectricField object -> can be for instance a laser beam
        transition: Transition object
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry

    Returns:
        Rabi frequency for a transition [2pi*MHz]
    """
    _ureg = E_field._ureg if _ureg is None else _ureg
    A = transition.A
    k = transition.k
    J_i = transition.state_i.quantum_numbers['J']
    J_f = transition.state_f.quantum_numbers['J']
    if transition.type == TransitionType.E1:
        return dipole_Rabi_Frequency(E_field.field(), A, k, J_i, J_f, mJ_i, mJ_f, _ureg)
    elif transition.type == TransitionType.E2:
        return quadrupole_Rabi_Frequency(E_field.gradient(), A, k, J_i, J_f, mJ_i, mJ_f, _ureg)
    else:
        raise NotImplementedError(f"Transition type {transition.type} not implemented")
