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


def ac_stark_shift(
    state: State,
    mJ: float,
    El_field: ElectricField,
    _ureg: pint.UnitRegistry,
    wavelengths: np.ndarray | None = None,
):
    delta_E = 0 * _ureg("k*mK")

    if wavelengths is None:
        omega_field = El_field.angular_frequency
        print(omega_field)
    else:
        omega_field = np.pi * 2 * _ureg("c") / wavelengths

    for transition in state.transitions_from:
        state_i = state
        state_f = transition.state_f

        for i in range(int(state_f.quantum_numbers["J"] * 2 + 1)):
            mJ_i = mJ
            mJ_f = -state_f.quantum_numbers["J"] + i

            try:
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
                        / (-transition.angular_frequency - omega_field)
                        + (Omega * np.conj(Omega))
                        / (-transition.angular_frequency + omega_field)
                    )
                )
            except Exception as e:
                print(
                    f"Error encountered: {e}. Transition type might not be implemented."
                )
                # Continue with the next operation or iteration here, if needed.

    for transition in state.transitions_to:
        state_i = transition.state_i
        state_f = state

        for i in range(int(state_i.quantum_numbers["J"] * 2 + 1)):
            mJ_f = mJ
            mJ_i = -state_i.quantum_numbers["J"] + i

            try:
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
                        / (transition.angular_frequency - omega_field)
                        + (Omega * np.conj(Omega))
                        / (transition.angular_frequency + omega_field)
                    )
                )
            except Exception as e:
                print(
                    f"Error encountered: {e}. Transition type might not be implemented."
                )
                # Continue with the next operation or iteration here, if needed.

    if wavelengths is None:
        return delta_E
    else:
        return np.array(
            [
                complex((delta_single_E.to("k*mK")).magnitude)
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
