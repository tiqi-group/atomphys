#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 06/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>


import pint
from math import pi
import warnings
from fractions import Fraction

from .quantum_numbers import QuantumNumbers
from .util import default_units

from .calc.hyperfine import hyperfine_shift
from .calc.zeeman import g_lande_fine_LS, g_lande_hyperfine
from .calc.angular_momentum import couple_angular_momenta, magnetic_sublevels
from .calc.coupling import get_coupling, Coupling


class State:

    configuration: str
    quantum_numbers: QuantumNumbers
    _energy: pint.Quantity
    _ureg: pint.UnitRegistry

    def __init__(self, configuration: str, term: str, energy: pint.Quantity,
                 atom=None,
                 _ureg: pint.UnitRegistry | None = None):

        self._atom = atom
        if atom is not None:
            self._ureg = atom._ureg
        elif _ureg is not None:
            self._ureg = _ureg
        else:
            self._ureg = pint.get_application_registry()

        self.configuration = configuration
        self.quantum_numbers = QuantumNumbers.from_term(term)
        self.energy = energy
    # Data

    @property
    def energy(self) -> pint.Quantity:
        return self._energy

    @energy.setter
    @default_units("Ry")
    def energy(self, value: pint.Quantity):
        self._energy = value

    # Identifiers

    @property
    def name(self) -> str:
        return f"{self.configuration} {self.term}"

    def __repr__(self) -> str:
        return f"State({self.name} {self.energy:0.4g~P})"

    @property
    def __key(self):
        return (self.configuration, self.quantum_numbers)

    def __eq__(self, other):
        if isinstance(other, State):
            return self.__key == other.__key
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, State):
            return self.energy < other.energy
        return NotImplemented

    def __hash__(self):
        return hash(self.__key)

    def match(self, pattern):
        return pattern.lower() in self.name.lower()

    # Derived properties

    @property
    def term(self) -> str:
        return self.quantum_numbers.term

    @property
    def coupling(self) -> Coupling:
        return get_coupling(self.quantum_numbers)

    @property
    def gJ(self) -> float:
        if self.coupling == Coupling.LS:
            L, S, J = (getattr(self.quantum_numbers, n) for n in ('L', 'S', 'J'))
            return g_lande_fine_LS(L, S, J)
        else:
            raise NotImplementedError(
                f"Lande factor calculation is implemented only for LS coupling, but state has {self.coupling}")

    @property
    def g(self) -> float:
        return self.gJ

    @property
    def sublevels(self) -> list[float]:
        return magnetic_sublevels(self.quantum_numbers.J)

    @property
    def Gamma(self) -> pint.Quantity:
        if self._atom is None:
            warnings.warn("State not attached to an Atom: no transitions available")
            return self._ureg("0 MHz")
        transitions_to = self._atom.transitions_to(self)
        if len(transitions_to) == 0:
            return self._ureg("0 MHz")
        return sum([tr.Gamma for tr in transitions_to])

    @property
    def lifetime(self) -> pint.Quantity:
        try:
            return (1 / (2 * pi * self.Gamma)).to('seconds')
        except ZeroDivisionError:
            return self._ureg("inf seconds")


class HyperfineState(State):
    def __init__(self, configuration: str, term: str, energy: pint.Quantity,
                 I: float, F: float,
                 atom=None, _ureg=None):
        super().__init__(configuration, term, energy, atom, _ureg)
        self._quantum_numbers = QuantumNumbers.from_term(term, I=I, F=F)

    @property
    def name(self):
        F = self.quantum_numbers.F
        return f"{self.configuration} {self.term} F={Fraction(F)}"

    def __repr__(self):
        return f"HyperfineState({self.name} {self.energy:0.4g~P})"

    @property
    def gF(self) -> float:
        J, I, F = (getattr(self.quantum_numbers, n) for n in ('J', 'I', 'F'))
        return g_lande_hyperfine(J, I, F, self.gJ)

    @property
    def g(self) -> float:
        return self.gF

    @property
    def sublevels(self) -> list[float]:
        return magnetic_sublevels(self.quantum_numbers.F)


def hyperfine_manifold(state: State, I: float, Ahf: pint.Quantity, Bhf: pint.Quantity = 0) -> list[HyperfineState]:
    J = state.quantum_numbers.J
    Fs = couple_angular_momenta(J, I)

    hstates = []
    for F in Fs:
        energy = state.energy + hyperfine_shift(J, I, F, Ahf, Bhf)
        hs = HyperfineState(state.configuration, state.term, energy, I, F, state._ureg, state._atom)
        hstates.append(hs)
    return hstates


# class StateRegistry(TypeRegistry):
#     def __init__(self, data=[], *args, **kwargs):
#         super().__init__(data=data, type=State, *args, **kwargs)

#     def __call__(self, key):
#         if isinstance(key, int):
#             return self[key]
#         elif isinstance(key, str):
#             try:
#                 return self(self._ureg.Quantity(key))
#             except (pint.errors.UndefinedUnitError, pint.errors.DimensionalityError):
#                 pass

#             try:
#                 return next(state for state in self if state.match(key))
#             except StopIteration:
#                 pass

#             raise KeyError(f"no state {key} found")
#         elif isinstance(key, float):
#             energy = self._ureg.Quantity(key, "E_h")
#             return min(self, key=lambda state: abs(state.energy - energy))
#         elif isinstance(key, self._ureg.Quantity):
#             return min(self, key=lambda state: abs(state.energy - key))
#         elif isinstance(key, Iterable):
#             return StateRegistry([self(item) for item in key], ureg=self._ureg)
#         else:
#             raise TypeError(
#                 "key must be integer index, term string, energy, or iterable"
#             )

#     def match(self, **kwargs):
#         return self.filter(
#             lambda state: all(
#                 getattr(state, key, lambda: None) == val
#                 for key, val in kwargs.items()
#                 if key not in ["energy", "En"]
#             )
#         )
