#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 06/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import pint
from math import pi

from .state import State
from .util import default_units, make_alias
from .calc.coupling import Coupling
from .calc.selection_rules import get_transition_type_LS, TransitionType
from .calc.matrix_element import reduced_electric_dipole_matrix_element, reduced_electric_quadrupole_matrix_element, electric_dipole_matrix_element, electric_quadrupole_matrix_element, dipole_matrix_element, quadrupole_matrix_element
from .calc.zeeman import field_sensitivity


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
        return f"Transition({self.state_i.name} --> {self.state_f.name} {self.wavelength:~.1fP} ({self.type.value}))"

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
        return (2 * pi / self.wavelength).to('1/m')

    @property
    def frequency(self) -> pint.Quantity:
        return self.energy.to('THz', 'sp')

    @property
    def angular_frequency(self) -> pint.Quantity:
        return self.energy.to('1/s', 'sp') * self._ureg('_2pi')

    @property
    def reduced_electric_matrix_element(self):
        if self.type == TransitionType.E1:
            J_f = self.state_f.quantum_numbers['J']
            return (reduced_electric_dipole_matrix_element(self.A, self.k, J_f, self._ureg)).to('e a0')
        elif self.type == TransitionType.E2:
            J_f = self.state_f.quantum_numbers['J']
            return reduced_electric_quadrupole_matrix_element(self.A, self.k, J_f, self._ureg).to('e a0**2')
        else:
            raise NotImplementedError(
                f"Matrix element calculation is implemented only for type E1/E2, but transition has {self.type}")

    def matrix_element(self, mJ_i: int, mJ_f: int):
        if self.type == TransitionType.E1:
            J_i = self.state_i.quantum_numbers['J']
            J_f = self.state_f.quantum_numbers['J']
            return dipole_matrix_element(self.A, self.k, J_i, J_f, mJ_i, mJ_f, self._ureg).to('a0')
        elif self.type == TransitionType.E2:
            J_i = self.state_i.quantum_numbers['J']
            J_f = self.state_f.quantum_numbers['J']
            return quadrupole_matrix_element(self.A, self.k, J_i, J_f, mJ_i, mJ_f, self._ureg).to('a0**2')

    def electric_matrix_element(self, mJ_i: int, mJ_f: int):
        if self.type == TransitionType.E1:
            J_i = self.state_i.quantum_numbers['J']
            J_f = self.state_f.quantum_numbers['J']
            return electric_dipole_matrix_element(self.A, self.k, J_i, J_f, mJ_i, mJ_f, self._ureg).to('e a0')
        elif self.type == TransitionType.E2:
            J_i = self.state_i.quantum_numbers['J']
            J_f = self.state_f.quantum_numbers['J']
            return electric_quadrupole_matrix_element(self.A, self.k, J_i, J_f, mJ_i, mJ_f, self._ureg).to('e a0**2')

    @property
    def delta_m(self):
        if self.type == TransitionType.E1:
            return 1
        elif self.type == TransitionType.E2:
            return 2
        else:
            return 0

    @property
    def sublevels(self):
        mis = self.state_i.sublevels
        mfs = self.state_f.sublevels
        ss = []
        for mi in mis:
            for mf in mfs:
                if abs(mf - mi) <= self.delta_m:
                    ss.append((mi, mf))
        return ss

    @property
    def sublevels_field_sentitivity(self):
        return {(mi, mf): field_sensitivity(self.state_f.g, mf, self._ureg) - field_sensitivity(self.state_i.g, mi, self._ureg)
                for (mi, mf) in self.sublevels}

    @property
    def saturation_intensity(self):
        return (self._ureg('pi*planck_constant*c/3') * self.A / (self.wavelength ** 3)).to('mW/cm^2')

    @property
    def cross_section(self):
        return (self._ureg('hbar/2') * self.ω * self.A / (self.Isat)).to('cm^2')

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

    def to_json(self):
        return {
            "A": f"{self.A.to('1/s').m} s^-1",
            "state_i": {
                "term": self.state_i.term,
                "energy": f"{self.state_i.energy.to('Ry'):f~P}"
            },
            "state_f": {
                "term": self.state_f.term,
                "energy": f"{self.state_f.energy.to('Ry'):f~P}"
            }
        }

    Einstein_coefficient = make_alias(attr_name='_A', get_unit='1/s')
    Γ = make_alias(attr_name='_A', get_unit='_2pi*MHz')
    Gamma = make_alias(attr_name='_A', get_unit='_2pi*MHz')
    ν = make_alias(attr_name='frequency', get_unit='THz')
    nu = make_alias(attr_name='frequency', get_unit='THz')
    ω = make_alias(attr_name='angular_frequency', get_unit='_2pi*1/s')
    omega = make_alias(attr_name='angular_frequency', get_unit='_2pi*1/s')
    λ = make_alias('wavelength', 'nm')
    d = make_alias('reduced_matrix_element')
    I_sat = make_alias(attr_name='saturation_intensity', get_unit='mW/cm^2')
    Isat = make_alias(attr_name='saturation_intensity', get_unit='mW/cm^2')
    I_s = make_alias(attr_name='saturation_intensity', get_unit='mW/cm^2')
    σ0 = make_alias(attr_name='cross_section', get_unit='cm^2')
