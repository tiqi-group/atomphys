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
from atomphys.util import set_default_units


def ac_stark_shift(
    state: State,
    mJ: float,
    El_field: ElectricField,
    wavelengths: np.ndarray | None = None,
    B: pint.Quantity = 0,
    _ureg: pint.UnitRegistry | None = None,
):
    delta_E = 0 * _ureg("k*mK")

    B = set_default_units(B, 'tesla')

    if wavelengths is None:
        omega_field = El_field.angular_frequency
    else:
        omega_field = np.pi * 2 * _ureg("c") / wavelengths

    for transition in state.transitions_from:
        zeeman_shifts = transition.sublevels_zeeman_shift(B)
        state_i = state
        state_f = transition.state_f

        mJ_i = mJ
        for mJ_f in state_f.sublevels:

            try:
                tr_omega = transition.angular_frequency + zeeman_shifts[(mJ_i, mJ_f)]

                Omega = complex(
                    (
                        Rabi_Frequency(
                            E_field=El_field,
                            transition=transition,
                            mJ_i=mJ_i,
                            mJ_f=mJ_f,
                            _ureg=_ureg,
                        )
                    )
                    .to("MHz")
                    .magnitude
                ) * _ureg("MHz")
                delta_E += (
                    _ureg("hbar")
                    / 4
                    * (
                        (Omega * np.conj(Omega))
                        / (-tr_omega - omega_field)
                        + (Omega * np.conj(Omega))
                        / (-tr_omega + omega_field)
                    )
                )
            except Exception as e:
                pass
                # Continue with the next operation or iteration here, if needed.

    for transition in state.transitions_to:
        zeeman_shifts = transition.sublevels_zeeman_shift(B)
        state_i = transition.state_i
        state_f = state

        mJ_f = mJ
        for mJ_i in state_i.sublevels:

            try:
                tr_omega = transition.angular_frequency + zeeman_shifts[(mJ_i, mJ_f)]

                Omega = complex(
                    (
                        Rabi_Frequency(
                            E_field=El_field,
                            transition=transition,
                            mJ_i=mJ_i,
                            mJ_f=mJ_f,
                            _ureg=_ureg,
                        )
                    )
                    .to("MHz")
                    .magnitude
                ) * _ureg("MHz")
                delta_E += (
                    _ureg("hbar")
                    / 4
                    * (
                        (Omega * np.conj(Omega))
                        / (tr_omega - omega_field)
                        + (Omega * np.conj(Omega))
                        / (tr_omega + omega_field)
                    )
                )
            except Exception as e:
                pass

    if wavelengths is None:
        return delta_E.to("k*mK").magnitude.real* _ureg("k*mK")
    else:
        return np.array(
            [
                ((delta_single_E.to("k*mK")).magnitude).real
                for delta_single_E in delta_E
            ]
        ) * _ureg("k*mK")


def polarizability(
    state: State,
    mJ: float,
    El_field: ElectricField,
    _ureg: pint.UnitRegistry,
    wavelengths: np.ndarray | None = None,
):
    AC_stark = ac_stark_shift(
        state=state, mJ=mJ, El_field=El_field, _ureg=_ureg, wavelengths=wavelengths
    )
    E0 = El_field.field_amplitude
    return (-4 * AC_stark / E0**2).to("e^2*a_0^2/E_h")
