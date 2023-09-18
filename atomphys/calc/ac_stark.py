from math import cos

import pint

from atomphys.calc.wigner import ishalfint, isint
from sympy.physics.wigner import wigner_6j as w6j
from sympy.physics.quantum.cg import CG
import numpy as np
from atomphys.calc.rabi_frequency import Rabi_Frequency
from atomphys.electric_field import ElectricField
from atomphys.transition import Transition, TransitionType
from atomphys.state import State


def AC_stark_shift(
    state: State,
    mJ: float,
    El_field: ElectricField,
    _ureg: pint.UnitRegistry | None = None,
):

    ΔE = 0

    omega_field = El_field.ω

    for transition in state.transitions_from:
        state_i = state
        state_f = transition.state_f

        for mJ_f in state_f.sublevels:
            mJ_i = mJ

            Ω = Rabi_Frequency(E_field=El_field, transition=transition, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
            ΔE += _ureg('hbar') / 4 * ((Ω * np.conj(Ω)) / (-transition.ω - omega_field) +
                                       (Ω * np.conj(Ω)) / (-transition.ω + omega_field))

    for transition in state.transitions_to:
        state_i = transition.state_i
        state_f = state

        for mJ_i in state_i.sublevels:
            mJ_f = mJ

            Ω = Rabi_Frequency(E_field=El_field, transition=transition, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
            ΔE += _ureg('hbar') / 4 * ((Ω * np.conj(Ω)) / (transition.ω - omega_field) +
                                       (Ω * np.conj(Ω)) / (transition.ω + omega_field))
    return ΔE.to('k*mK').magnitude
