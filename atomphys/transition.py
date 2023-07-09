#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 06/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import pint
from math import pi

from .state import State
from .util import default_units
from .calc.coupling import Coupling
from .calc.selection_rules import get_transition_type_LS, TransitionType
from .calc.matrix_element import reduced_dipole_matrix_element, reduced_quadrupole_matrix_element


class Transition:
    state_i: State
    state_f: State
    _A: pint.Quantity
    _ureg: pint.UnitRegistry

    def __init__(self, state_i: State, state_f: State, A: pint.Quantity,
                 _ureg: pint.UnitRegistry | None = None):
        self._ureg = _ureg if _ureg is not None else pint.get_application_registry()
        # sort states
        self.state_i, self.state_f = sorted([state_i, state_f])
        self.A = A

    def __repr__(self) -> str:
        return f"Transition({self.state_i.name} --> {self.state_f.name} {self.wavelength} ({self.type.value}))"

    @property
    def A(self) -> pint.Quantity:
        """Einstein coefficient"""
        return self._A

    @A.setter
    @default_units('rad / s')
    def A(self, value: pint.Quantity):
        self._A = value

    @property
    def energy(self) -> pint.Quantity:
        return self.state_f.energy - self.state_i.energy

    @property
    def frequency(self) -> pint.Quantity:
        return self.energy.to('THz', 'sp')

    @property
    def wavelength(self) -> pint.Quantity:
        try:
            return self.energy.to('nm', 'sp')
        except ZeroDivisionError:
            return self._ureg.Quantity("inf nm")

    @property
    def Gamma(self) -> pint.Quantity:
        return self.A.to('MHz') / 2 / pi

    @property
    def type(self) -> TransitionType:
        # TODO: get the correct transition type from quantum numbers
        if self.state_i.coupling == Coupling.LS and self.state_f.coupling == Coupling.LS:
            qn_i = self.state_i.quantum_numbers
            qn_f = self.state_f.quantum_numbers
            return get_transition_type_LS(qn_i, qn_f)
        else:
            # print(
            #     f"Transition type calculation is implemented only between states with LS coupling, but transition has {self.state_i.coupling} --> {self.state_f.coupling}")
            return TransitionType.NONE

    @property
    def matrix_element(self):
        if self.type == TransitionType.E1:
            return reduced_dipole_matrix_element(self.A, self.wavelength, self._ureg)
        elif self.type == TransitionType.E2:
            return reduced_quadrupole_matrix_element(self.A, self.wavelength, self._ureg)
        else:
            raise NotImplementedError(
                f"Matrix element calculation is implemented only for type E1/E2, but transition has {self.type}")
