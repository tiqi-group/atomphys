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


def make_property(attr_name: str, get_unit: str = None, conversion: str = None, extra_factor: pint.Quantity = 1):
    @property
    def prop(self):
        if get_unit is None:
            return getattr(self, attr_name)
        if conversion is None:
            return (getattr(self, attr_name) * extra_factor).to(get_unit)
        else:
            return (getattr(self, attr_name) * extra_factor).to(get_unit, conversion)
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

    A = make_property(attr_name='_A', get_unit='1/s')
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
    def reduced_matrix_element(self):
        if self.type == TransitionType.E1:
            return (reduced_dipole_matrix_element(self.A, self.wavelength, self._ureg)).to('e * a0')
        elif self.type == TransitionType.E2:
            return reduced_quadrupole_matrix_element(self.A, self.wavelength, self._ureg)
        else:
            raise NotImplementedError(
                f"Matrix element calculation is implemented only for type E1/E2, but transition has {self.type}")


    Einstein_coefficient = make_property(attr_name='_A', get_unit='1/s')
    Γ = make_property(attr_name='_A', get_unit='2*pi*MHz')
    Gamma = make_property(attr_name='_A', get_unit='2*pi*MHz')
    
    frequency = make_property(attr_name='energy', get_unit='THz', conversion='sp')
    ν = make_property(attr_name='energy', get_unit='THz', conversion='sp')
    nu = make_property(attr_name='energy', get_unit='THz', conversion='sp')
    
    ω = make_property(attr_name='energy', get_unit='2*pi*THz', conversion='sp', extra_factor=pi*2)
    omega = make_property(attr_name='energy', get_unit='2*pi*THz', conversion='sp', extra_factor=pi*2)
    angular_frequency = make_property(attr_name='energy', get_unit='2*pi*THz', conversion='sp', extra_factor=pi*2)
    
    λ = make_property('wavelength', 'nm')
    
    d = make_property('reduced_matrix_element')


    @property
    def saturation_intensity(self):
        h = self._ureg.planck_constant
        c = self._ureg.c
        return π * h * c * self.Γ / (3 * self.λ ** 3)
    

    I_sat = def make_property(attr_name='saturation_intensity', get_unit='W/cm^2')
    Isat = def make_property(attr_name='saturation_intensity', get_unit='W/cm^2')
    I_s = def make_property(attr_name='saturation_intensity', get_unit='W/cm^2')
        
    @property
    def cross_section(self):
        ħ = self._ureg.ħ
        return ħ * self.ω * self.Γ / (2 * self.Isat)
    
    σ0 = def make_property(attr_name='cross_section', get_unit='cm^2')

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

    

    def _reduced_dipole_matrix_element(self):
        """
        Calculates the reduced dipole matrix element of the transition.
        
        Reference for the calculation is Thesis of Christoph Fisher, page 34.

        Returns:
            pint.Quantity: The reduced dipole matrix element of the transition.
        """
        C = self._ureg('3 * pi * epsilon_0 * hba * c^3')
        return (C/self.omega**3 * self.A * (2*self.state_f.quantum_numbers.J+1))**(1/2)


    def _reduced_quadrupole_matrix_element(self):
        beta = 1  # TODO: find the correct constant
        C = 4 * beta * pi * _ureg.epsilon_0 * _ureg.hbar
        return (C * (wavelength / 2 / pi)**5 * A)**(1 / 2)