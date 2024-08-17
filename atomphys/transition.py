#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pint
from .state import State
from .utils.utils import default_units
from .utils.coupling import Coupling
from .utils.selection_rules import get_transition_type_LS, TransitionType
from .calc.matrix_element import (
    reduced_electric_dipole_matrix_element,
    reduced_electric_quadrupole_matrix_element,
    electric_dipole_matrix_element,
    electric_quadrupole_matrix_element,
    dipole_matrix_element,
    quadrupole_matrix_element,
)
from .calc.zeeman import field_sensitivity, zeeman_shift


class Transition:
    state_i: State
    state_f: State
    _A: pint.Quantity
    _ureg: pint.UnitRegistry

    def __init__(self, state_i: State, state_f: State, A: pint.Quantity, _ureg: pint.UnitRegistry | None = None):
        self._ureg = _ureg if _ureg is not None else pint.get_application_registry()
        self.state_i, self.state_f = sorted([state_i, state_f])
        self.A = A

    def __repr__(self) -> str:
        return f"Transition({self.state_i.name} --> {self.state_f.name} {self.wavelength:~.1fP} ({self.type.value}))"

    @property
    def A(self) -> pint.Quantity:
        return self._A.to("1/s")

    @A.setter
    @default_units("1 / s")
    def A(self, value: pint.Quantity):
        self._A = value
        
    @property
    def Gamma(self) -> pint.Quantity:
        return self.A.to('_2pi*MHz')
    
    @Gamma.setter
    @default_units("_2pi*MHz")
    def Gamma(self, value: pint.Quantity):
        self.A = value.to('1/s')

    @property
    def energy(self) -> pint.Quantity:
        return (self.state_f.energy - self.state_i.energy).to("Ry")

    @property
    def wavelength(self) -> pint.Quantity:
        try:
            return self.energy.to("nm", "sp")
        except ZeroDivisionError:
            return self._ureg.Quantity("inf nm")

    @property
    def k(self) -> pint.Quantity:
        return (self._ureg('_2pi') * 1 / self.wavelength).to("1/m")

    @property
    def frequency(self) -> pint.Quantity:
        return self.energy.to("THz", "sp")

    @property
    def angular_frequency(self) -> pint.Quantity:
        return self.energy.to('1/s', 'sp') * self._ureg('_2pi')

    @property
    def reduced_electric_matrix_element(self):
        """
        Reduced electric matrix element in units of e a0

        Returns:
            pint.Quantity: reduced electric matrix element
        """
        if self.type == TransitionType.E1:
            J_f = self.state_f.quantum_numbers["J"]
            return (reduced_electric_dipole_matrix_element(self.A, self.k, J_f, self._ureg)).to("e a0")
        elif self.type == TransitionType.E2:
            J_f = self.state_f.quantum_numbers["J"]
            return reduced_electric_quadrupole_matrix_element(self.A, self.k, J_f, self._ureg).to("e a0**2")
        else:
            raise NotImplementedError(
                f"Matrix element calculation is implemented only for type E1/E2, but transition has {self.type}"
            )

    def matrix_element(self, mJ_i: int, mJ_f: int):
        if self.type == TransitionType.E1:
            J_i = self.state_i.quantum_numbers["J"]
            J_f = self.state_f.quantum_numbers["J"]
            return dipole_matrix_element(self.A, self.k, J_i, J_f, mJ_i, mJ_f, self._ureg).to("a0")
        elif self.type == TransitionType.E2:
            J_i = self.state_i.quantum_numbers["J"]
            J_f = self.state_f.quantum_numbers["J"]
            return quadrupole_matrix_element(self.A, self.k, J_i, J_f, mJ_i, mJ_f, self._ureg).to("a0**2")
        else:
            raise NotImplementedError(
                f"Matrix element calculation is implemented only for type E1/E2, but transition has {self.type}"
            )
    
    def electric_matrix_element(self, mJ_i: int, mJ_f: int):
        if self.type == TransitionType.E1:
            J_i = self.state_i.quantum_numbers["J"]
            J_f = self.state_f.quantum_numbers["J"]
            return electric_dipole_matrix_element(self.A, self.k, J_i, J_f, mJ_i, mJ_f, self._ureg).to("e a0")
        elif self.type == TransitionType.E2:
            J_i = self.state_i.quantum_numbers["J"]
            J_f = self.state_f.quantum_numbers["J"]
            return electric_quadrupole_matrix_element(self.A, self.k, J_i, J_f, mJ_i, mJ_f, self._ureg).to("e a0**2")
        else:
            raise NotImplementedError(
                f"Matrix element calculation is implemented only for type E1/E2, but transition has {self.type}"
            )

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
        return {
            (mi, mf): (field_sensitivity(self.state_f.g, mf, self._ureg) -
                       field_sensitivity(self.state_i.g, mi, self._ureg))
            for (mi, mf) in self.sublevels
        }

    def sublevels_zeeman_shift(self, B: pint.Quantity) -> dict[float: pint.Quantity]:
        return {
            (mi, mf): (zeeman_shift(self.state_f.g, mf, B, self._ureg) -
                       zeeman_shift(self.state_i.g, mi, B, self._ureg)).to('MHz')
            for (mi, mf) in self.sublevels
        }

    @property
    def saturation_intensity(self):
        return (self._ureg('pi*planck_constant*c/3') * self.A / (self.wavelength ** 3)).to('mW/cm^2')
    
    @property
    def cross_section(self):
        return (self._ureg('hbar/2') * self.angular_frequency * self.A / (self.saturation_intensity)).to('cm^2')
    
    @property
    def type(self) -> TransitionType:
        # TODO: get the correct transition type from quantum numbers
        if self.state_i.coupling == Coupling.LS and self.state_f.coupling == Coupling.LS:
            qn_i = self.state_i.quantum_numbers
            qn_f = self.state_f.quantum_numbers
            return get_transition_type_LS(qn_i, qn_f)
        else:
            return TransitionType.NONE

    def to_json(self):
        return {
            "A": f"{self.A.to('1/s').m} s^-1",
            "state_i": {"term": self.state_i.term, "energy": f"{self.state_i.energy.to('Ry'):f~P}"},
            "state_f": {"term": self.state_f.term, "energy": f"{self.state_f.energy.to('Ry'):f~P}"},
        }
