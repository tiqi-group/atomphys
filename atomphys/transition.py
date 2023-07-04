#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 06/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

from dataclasses import dataclass
from math import inf

from .state import State
import pint


@dataclass(frozen=True, repr=False)
class Transition:
    # TODO: not sure dataclasses are good, this is ust a prototype
    state_i: State
    state_f: State
    _ureg: pint.UnitRegistry = None

    def __post_init__(self):
        if self._ureg is None:
            super().__setattr__('_ureg', pint.get_application_registry())
        # sort states
        si, sf = sorted([self.state_i, self.state_f])
        super().__setattr__('state_i', si)
        super().__setattr__('state_f', sf)

    def __repr__(self) -> str:
        return f"Transition({self.state_i.name} --> {self.state_f.name} {self.wavelength})"

    @property
    def A(self):
        # TODO this is just a mock value, get it from input
        return self._ureg.Quantity('inf s^-1')

    @property
    def wavelength(self):
        try:
            return (self.state_f.energy - self.state_i.energy).to('nm', 'sp')
        except ZeroDivisionError:
            return self._ureg.Quantity(inf, "nm")
