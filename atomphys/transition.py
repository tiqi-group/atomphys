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
from .calc.matrix_element import reduced_electric_dipole_matrix_element, reduced_electric_quadrupole_matrix_element


def make_alias(attr_name: str, get_unit: str = None):
    @property
    def prop(self):
        if get_unit is None:
            return getattr(self, attr_name)
        else:
            return getattr(self, attr_name).to(get_unit)
    return prop

class Transition:
    """
    This class represents a transition between two states in atomic physics.
    
    Properties:
        state_i (State): The initial state of the transition.
        state_f (State): The final state of the transition.
        A (pint.Quantity): The Einstein coefficient of the transition.
        energy (pint.Quantity): The energy of the transition.
        frequency (pint.Quantity): The frequency of the transition.
        wavelength (pint.Quantity): The wavelength of the transition.
        Gamma (pint.Quantity): The decay rate of the transition.
        type (TransitionType): The type of the transition.
    """
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
        return self._A.to('1/s')

    @A.setter
    @default_units('1 / s')
    def A(self, value: pint.Quantity):
        self._A = value

    @property
    def energy(self) -> pint.Quantity:
        return (self.state_f.energy - self.state_i.energy).to('Ry')
    
    @property
    def wavelength(self) -> pint.Quantity:
        try:
            return self.energy.to('nm', 'sp')
        except ZeroDivisionError:
            return self._ureg.Quantity("inf nm")
        
    @property
    def k(self) -> pint.Quantity:
        return (2*pi/self.wavelength).to('1/m')
        
    @property
    def frequency(self) -> pint.Quantity:
        return self.energy.to('THz', 'sp')

    @property
    def angular_frequency(self) -> pint.Quantity:
        return self.energy.to('1/s', 'sp')*self._ureg('2*pi')
      
    @property
    def reduced_matrix_element(self):
        if self.type == TransitionType.E1:
            return (reduced_electric_dipole_matrix_element(self.A, self.wavelength, self.state_f.quantum_numbers.J, self._ureg)).to('e * a0')
        elif self.type == TransitionType.E2:
            return reduced_electric_quadrupole_matrix_element(self.A, self.wavelength, self._ureg)
        else:
            raise NotImplementedError(
                f"Matrix element calculation is implemented only for type E1/E2, but transition has {self.type}")

    @property
    def saturation_intensity(self):
        h = self._ureg.planck_constant
        c = self._ureg.c
        return self._ureg('pi*planck_constant*c/3')* self.Γ / (self.λ ** 3)
    
    @property
    def cross_section(self):
        return self._ureg('hbar/2') * self.ω * self.Γ / (self.Isat)
    
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
    
    Einstein_coefficient = make_alias(attr_name='_A', get_unit='1/s')
    Γ = make_alias(attr_name='_A', get_unit='2*pi*MHz')
    Gamma = make_alias(attr_name='_A', get_unit='2*pi*MHz')
    ν = make_alias(attr_name='frequency', get_unit='THz')
    nu = make_alias(attr_name='frequency', get_unit='THz')
    ω = make_alias(attr_name='angular_frequency', get_unit='2*pi*1/s')
    omega = make_alias(attr_name='angular_frequency', get_unit='2*pi*1/s') 
    λ = make_alias('wavelength', 'nm')
    d = make_alias('reduced_matrix_element')
    I_sat = make_alias(attr_name='saturation_intensity', get_unit='W/cm^2')
    Isat = make_alias(attr_name='saturation_intensity', get_unit='W/cm^2')
    I_s = make_alias(attr_name='saturation_intensity', get_unit='W/cm^2')
    σ0 = make_alias(attr_name='cross_section', get_unit='cm^2')

    