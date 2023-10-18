#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Authors: Carmelo Mordini <cmordini@phys.ethz.ch> & Wojciech Adamczyk <wadamczyk@phys.ethz.ch>


import numpy as np
import pint

from atomphys.calc.rabi_frequency import Rabi_Frequency
from atomphys.electric_field import ElectricField
from atomphys.state import State
from atomphys.util import set_default_units


def AC_stark_shift(
    state: State,
    mJ: float,
    El_field: ElectricField,
    B: pint.Quantity = 0,
    _ureg: pint.UnitRegistry | None = None,
):

    ΔE = 0
    B = set_default_units(B, 'tesla')

    omega_field = El_field.ω

    for transition in state.transitions_from:
        zeeman_shifts = transition.sublevels_zeeman_shift(B)
        state_i = state
        state_f = transition.state_f

        mJ_i = mJ
        for mJ_f in state_f.sublevels:

            Ω = Rabi_Frequency(E_field=El_field, transition=transition, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
            if Ω != 0:
                tr_omega = transition.omega + zeeman_shifts[(mJ_i, mJ_f)]
                ΔE += _ureg('hbar') / 4 * ((Ω * np.conj(Ω)) / (-tr_omega - omega_field) +
                                           (Ω * np.conj(Ω)) / (-tr_omega + omega_field))

    for transition in state.transitions_to:
        zeeman_shifts = transition.sublevels_zeeman_shift(B)
        state_i = transition.state_i
        state_f = state

        mJ_f = mJ
        for mJ_i in state_i.sublevels:

            Ω = Rabi_Frequency(E_field=El_field, transition=transition, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
            if Ω != 0:
                tr_omega = transition.omega + zeeman_shifts[(mJ_i, mJ_f)]
                ΔE += _ureg('hbar') / 4 * ((Ω * np.conj(Ω)) / (tr_omega - omega_field) +
                                           (Ω * np.conj(Ω)) / (tr_omega + omega_field))
    return ΔE  # .to('k*mK').magnitude
