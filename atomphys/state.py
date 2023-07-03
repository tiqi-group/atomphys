#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 06/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>


import pint

from . import _ureg
from .term import Coupling, QuantumNumbers
from .util import default_units

from fractions import Fraction
from atomphys.calc.hyperfine import hyperfine_shift
from atomphys.calc.zeeman import g_lande_fine_LS, g_lande_hyperfine
from atomphys.calc.angular_momentum import couple_angular_momenta, magnetic_sublevels


class State:

    _configuration: str
    _energy: pint.Quantity
    _quantum_numbers: QuantumNumbers
    _ureg: pint.UnitRegistry
    _atom = None

    def __init__(self, configuration: str, term: str, energy: pint.Quantity, ureg=None, atom=None) -> None:
        self._configuration = configuration
        self._term = term
        self._energy = pint.Quantity(energy)
        self._quantum_numbers = QuantumNumbers.from_term(term)

        self._ureg = _ureg
        self._atom = atom
        if ureg is not None:
            self._ureg = ureg
        if atom is not None:
            self._ureg = atom._ureg

    def __repr__(self):
        return f"State({self.configuration} {self.term} {self.energy:0.4g~P})"

    def __eq__(self, other):
        return self.energy == other.energy and self.quantum_numbers == other.quantum_numbers

    def __lt__(self, other):
        return self.energy < other.energy

    def match(self, pattern):
        name = f"{self.configuration} {self.term}"
        return pattern.lower() in name.lower()

    @property
    def configuration(self) -> str:
        return self._configuration

    @property
    def term(self) -> str:
        return self.quantum_numbers.term

    @property
    def energy(self) -> pint.Quantity:
        return self._energy

    @energy.setter
    @default_units("E_h")
    def energy(self, energy: pint.Quantity):
        self._energy = energy

    @property
    def quantum_numbers(self) -> QuantumNumbers:
        return self._quantum_numbers

    @property
    def coupling(self):
        if self.quantum_numbers.L is not None and self.quantum_numbers.S is not None:
            return Coupling.LS
        if self.quantum_numbers.J1 is not None and self.quantum_numbers.J2 is not None:
            return Coupling.jj
        if self.quantum_numbers.S2 is not None and self.quantum_numbers.K is not None:
            return Coupling.LK

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
    def Gamma(self):
        # return sum([transition.Gamma for transition in self.transitions_down])
        raise NotImplementedError

    @property
    def lifetime(self):
        return 1 / self.Gamma


class HyperfineState(State):
    def __init__(self, configuration: str, term: str, energy: pint.Quantity,
                 I: float, F: float,
                 ureg=None, atom=None) -> None:
        super().__init__(configuration, term, energy, ureg, atom)
        self._quantum_numbers = QuantumNumbers.from_term(term, I=I, F=F)

    def __repr__(self):
        F = self.quantum_numbers.F
        return f"HyperfineState({self.configuration} {self.term} F={Fraction(F)} {self.energy:0.4g~P})"

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
