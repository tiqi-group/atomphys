#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 06/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>


import pint
from math import pi
import warnings
from fractions import Fraction
from typing import TYPE_CHECKING

from .quantum_numbers import QuantumNumbers
from .util import default_units

from .calc.hyperfine import hyperfine_shift
from .calc.zeeman import g_lande_fine_LS, g_lande_hyperfine, field_sensitivity
from .calc.angular_momentum import couple_angular_momenta, magnetic_sublevels
from .calc.coupling import get_coupling, Coupling

if TYPE_CHECKING:
    from .atom import Atom


class State:

    """
    This class represents a quantum state in atomic physics.

    The idea behind this class is that all properties of the state are re-derived from
    few basic properties, such as the electron configuration, the term symbol and the energy of the atom.

    Attributes:
        configuration (str): The electron configuration of the state.
        quantum_numbers (QuantumNumbers): Quantum numbers of the state.
        _energy (pint.Quantity): Energy of the state.
        _ureg (pint.UnitRegistry): Unit registry for the state.

    Methods:
        __init__(configuration, term, energy, atom=None, _ureg=None): Initializes a new instance of the State class.
        energy: Property that gets or sets the energy of the state.
        name: Property that gets the name of the state.
        __repr__(): Returns a string representation of the state.
        __eq__(other): Determines whether the current state is equal to another state.
        __lt__(other): Determines whether the current state is less than another state.
        __hash__(): Returns a hash value for the state.
        match(pattern): Determines whether the state name matches a pattern.
        term: Property that gets the term of the state.
        coupling: Property that gets the coupling of the state.
        atom: Property that gets or sets the atom of the state.
        gJ: Property that gets the Landé g-factor for the state.
        g: Property that gets the Landé g-factor for the state.
        sublevels: Property that gets the magnetic sublevels of the state.
        transitions_from: Property that gets the transitions from the state.
        transitions_to: Property that gets the transitions to the state.
        Gamma: Property that gets the total decay rate of the state.
        decay_branching_ratios: Property that gets the decay branching ratios of the state.
        lifetime: Property that gets the lifetime of the state.
    """

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

    def __repr__(self) -> str:
        return f"State({self.name} {self.energy:0.4g~P})"

    @property
    def energy(self) -> pint.Quantity:
        return self._energy

    @energy.setter
    @default_units("Ry")
    def energy(self, value: pint.Quantity):
        self._energy = value

    @property
    def frequency(self) -> pint.Quantity:
        return (self.energy.to('THz', 'sp'))

    @property
    def angular_frequency(self) -> pint.Quantity:
        return (self.energy.to('1/s', 'sp') * self._ureg('2*pi'))

    @property
    def name(self) -> str:
        return f"{self.configuration} {self.term}"

    @property
    def __key(self):
        return (self.configuration, self.quantum_numbers, self.energy)

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

    @property
    def atom(self) -> 'Atom':
        return self._atom

    @atom.setter
    def atom(self, value: 'Atom'):
        self._atom = value

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
    def sublevels_field_sensitivity(self) -> dict[float: pint.Quantity]:
        return {m: field_sensitivity(self.g, m, self._ureg) for m in self.sublevels}

    @property
    def transitions_from(self) -> list:
        if self._atom is None:
            warnings.warn("State not attached to an Atom: no transitions available")
            return list()
        return self._atom.transitions_from(self)

    @property
    def transitions_to(self) -> list:
        if self._atom is None:
            warnings.warn("State not attached to an Atom: no transitions available")
            return list()
        return self._atom.transitions_to(self)

    @property
    def transitions(self) -> list:
        return self.transitions_from + self.transitions_to

    @property
    def Gamma(self) -> pint.Quantity:
        if len(self.transitions_to) == 0:
            return self._ureg("0 MHz")
        return sum([tr.Gamma for tr in self.transitions_to])

    @property
    def decay_branching_ratios(self) -> dict:
        transitions_to = self.transitions_to.copy()
        transitions_to.sort(key=lambda tr: tr.Gamma.m, reverse=True)
        return {tr.state_i: (tr.Gamma / self.Gamma).m for tr in transitions_to}

    @property
    def lifetime(self) -> pint.Quantity:
        try:
            return (1 / (2 * pi * self.Gamma)).to('seconds')
        except ZeroDivisionError:
            return self._ureg("inf seconds")

    def to_json(self):
        return {
            'configuration': self.configuration,
            'term': self.term,
            'energy': f"{self.energy.to('Ry'):f~P}"
        }


class HyperfineState(State):
    """
    This class represents a hyperfine quantum state in atomic physics, and it extends the State class.

    Methods:
        __init__(configuration, term, energy, I, F, atom=None, _ureg=None): Initializes a new instance of the HyperfineState class.
        name: Property that gets the name of the state.
        __repr__(): Returns a string representation of the state.
        gF: Property that gets the Landé g-factor for the hyperfine state.
        g: Property that gets the Landé g-factor for the hyperfine state.
        sublevels: Property that gets the magnetic sublevels of the hyperfine state.
    """

    def __init__(self, configuration: str, term: str, energy: pint.Quantity,
                 I: float, F: float,  # noqa: E741
                 atom=None, _ureg=None):
        super().__init__(configuration, term, energy, atom, _ureg)
        self.quantum_numbers = QuantumNumbers.from_term(term, I=I, F=F)

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


def hyperfine_manifold(state: State, I: float, Ahf: pint.Quantity, Bhf: pint.Quantity = 0) -> list[HyperfineState]:  # noqa: E741
    """
    This function generates a list of hyperfine states for a given state.

    Arguments:
        state (State): The state for which to generate the hyperfine states.
        I (float): The nuclear spin quantum number.
        Ahf (pint.Quantity): The magnetic dipole hyperfine structure constant.
        Bhf (pint.Quantity, optional): The electric quadrupole hyperfine structure constant. Defaults to 0.

    Returns:
        list[HyperfineState]: A list of hyperfine states.
    """
    J = state.quantum_numbers.J
    Fs = couple_angular_momenta(J, I)

    hstates = []
    for F in Fs:
        energy = state.energy + hyperfine_shift(J, I, F, Ahf, Bhf)
        hs = HyperfineState(state.configuration, state.term, energy, I, F, state._atom, state._ureg)
        hstates.append(hs)
    return hstates
