#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 09/2023
# Author: Wojciech Adamczyk <wadamczyk@phys.ethz.ch>

import pint
import numpy as np
from atomphys.calc.rabi_frequency import Rabi_Frequency
from atomphys.electric_field import ElectricField
from atomphys.state import State


def AC_stark_shift(
    state: State,
    mJ: float, 
    El_field: ElectricField,
    _ureg: pint.UnitRegistry,
):

    delta_E = 0 * _ureg('k*mK')

    omega_field = El_field.angular_frequency
    

    for transition in state.transitions_from:
        state_i = state
        state_f = transition.state_f

        for mJ_f in range(int(state_f.quantum_numbers['J']*2+1)):
            mJ_i = mJ

            Omega = Rabi_Frequency(E_field=ElectricField, transition=transition, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
            delta_E += _ureg('hbar')/4*((Omega*np.conj(Omega))/(-transition.angular_frequency - omega_field)+(Omega*np.conj(Omega))/(-transition.angular_frequency + omega_field))

    for transition in state.transitions_to:
        state_i = transition.state_i
        state_f = state
        
        for mJ_i in range(int(state_i.quantum_numbers['J']*2+1)):
            mJ_f = mJ

            Omega = Rabi_Frequency(E_field=ElectricField, transition=transition, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
            delta_E += _ureg('hbar')/4*((Omega*np.conj(Omega))/(transition.ω - omega_field) + (Omega*np.conj(Omega))/(transition.angular_frequency + omega_field))
    return delta_E.to('k*mK')

def polarizability(
    state: State,
    mJ: float, 
    El_field: ElectricField,
    _ureg: pint.UnitRegistry,
):
    AC_stark = AC_stark_shift(state=state, mJ=mJ, El_field=El_field, _ureg=_ureg)
    E0 = El_field.amplitude
    return (-4*AC_stark/E0**2).to('e^2*a_0^2/E_h')