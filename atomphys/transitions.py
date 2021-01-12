import csv
import io
import urllib.request
from .util import sanitize_energy
from .data import nist
from .states import State

from math import pi as π
from math import inf


try:
    from . import _ureg, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None


class TransitionRegistry(list):
    def __getitem__(self, key):
        if isinstance(key, int):
            return super().__getitem__(key)
        if isinstance(key, slice):
            return TransitionRegistry(super().__getitem__(key))
        else:
            raise TypeError('key must be integer, slice, or term string')

    def __repr__(self):
        repr = f'{len(self)} Transitions (\n'
        if self.__len__() <= 6:
            for transition in self:
                repr += (str(transition) + '\n')
        else:
            for transition in self[:3]:
                repr += (str(transition) + '\n')
            repr += '...\n'
            for transition in self[-3:]:
                repr += (str(transition) + '\n')
        repr = repr[:-1] + ')'
        return repr

    def __add__(self, other):
        assert isinstance(other, TransitionRegistry)
        return TransitionRegistry(list(self) + list(other))

    def up_from(self, state: State):
        return TransitionRegistry(transition for transition in self if transition.i == state)

    def down_from(self, state: State):
        return TransitionRegistry(transition for transition in self if transition.f == state)

    def to_dict(self):
        return [transition.to_dict() for transition in self]


class Transition(dict):

    _USE_UNITS = False
    _ureg = {}
    _state_i: State = None
    _state_f: State = None

    def __init__(self, USE_UNITS=False, ureg=None, **transition):
        self._USE_UNITS = USE_UNITS and _HAS_PINT
        if ureg and self._USE_UNITS:
            self._ureg = ureg
        elif self._USE_UNITS:
            self._ureg = _ureg
        else:
            self._ureg = {}

        if not self._USE_UNITS:
            self._ureg['hbar'] = 1
            self._ureg['h'] = 2*π
            self._ureg['ε_0'] = 1/(4*π)
            self._ureg['c'] = 137.03599908356244

        if 'Gamma' in transition:
            if self._USE_UNITS:
                Gamma = self._ureg.Quantity(transition['Gamma'])
            else:
                Gamma = float(transition['Gamma'])
        elif 'Aki(s^-1)' in transition:
            if self._USE_UNITS:
                Gamma = self._ureg.Quantity(
                    float(transition['Aki(s^-1)']), 's^-1').to('Eh/hbar')
            else:
                Gamma = 2.4188843265856806e-17 * float(transition['Aki(s^-1)'])
        else:
            Gamma = 0.0

        if 'Ei' in transition:
            if self._USE_UNITS:
                Ei = self._ureg.Quantity(transition['Ei'])
            else:
                Ei = float(transition['Ei'])
        elif 'Ei(Ry)' in transition:
            if self._USE_UNITS:
                Ei = self._ureg.Quantity(float(sanitize_energy(
                    transition['Ei(Ry)'])), 'Ry').to('Eh')
            else:
                Ei = 0.5 * float(sanitize_energy(transition['Ei(Ry)']))
        else:
            Ei = 0.0

        if 'Ef' in transition:
            if self._USE_UNITS:
                Ef = self._ureg.Quantity(transition['Ef'])
            else:
                Ef = float(transition['Ef'])
        elif 'Ek(Ry)' in transition:
            if self._USE_UNITS:
                Ef = self._ureg.Quantity(float(sanitize_energy(
                    transition['Ek(Ry)'])), 'Ry').to('Eh')
            else:
                Ef = 0.5 * float(sanitize_energy(transition['Ek(Ry)']))
        else:
            Ef = 0.0

        super(Transition, self).__init__({'Ei': Ei, 'Ef': Ef, 'Gamma': Gamma})

    def __repr__(self):
        fmt = '0.4g~P' if self._USE_UNITS else '0.4g'
        if self.i is not None:
            state_i = f'{self.i.valence} {self.i.term}'
        else:
            state_i = f'{self.Ei:{fmt}}'
        if self.f is not None:
            state_f = f'{self.f.valence} {self.f.term}'
        else:
            state_f = f'{self.Ef:{fmt}}'

        return (
            f'Transition({state_i} <---> {state_f}, '
            f'λ={self.λ_nm:{fmt}}, '
            f'Γ=2π×{self.Γ_MHz/(2*π):{fmt}})'
        )

    def to_dict(self):
        return {'Ei': str(self.Ei), 'Ef': str(self.Ef), 'Gamma': str(self.Gamma)}

    @property
    def Ei(self):
        return self['Ei']

    @property
    def Ef(self):
        return self['Ef']

    @property
    def Gamma(self):
        return self['Gamma']

    @property
    def Γ(self):
        return self['Gamma']

    @property
    def Gamma_MHz(self):
        return self.Γ_MHz

    @property
    def Γ_MHz(self):
        if self._USE_UNITS:
            return self.Γ.to('MHz')
        else:
            return self.Γ * 41341373335.18245  # E_h / hbar / MHz

    @property
    def i(self):
        return self._state_i

    @property
    def f(self):
        return self._state_f

    @property
    def ω(self):
        ℏ = self._ureg['hbar']
        return (self.Ef-self.Ei)/ℏ

    @property
    def angular_frequency(self):
        return self.ω

    @property
    def ν(self):
        return self.ω/(2*π)

    @property
    def frequency(self):
        return self.ν

    @property
    def λ(self):
        c = self._ureg['c']
        try:
            return c/self.ν
        except ZeroDivisionError:
            return inf

    @property
    def wavelength(self):
        return self.λ

    @property
    def λ_nm(self):
        if self._USE_UNITS:
            return self.λ.to('nm')
        else:
            return self.λ * 0.052917721090397746  # nm / a_0

    @property
    def wavelength_nm(self):
        return self.λ_nm

    @property
    def saturation_intensity(self):
        h = self._ureg['h']
        c = self._ureg['c']
        return π*h*c*self.Γ/(3*self.λ**3)

    @property
    def branching_ratio(self):
        r = self.Γ * self.f.τ
        if self._USE_UNITS and isinstance(r, self._ureg.Quantity):
            r = r.m
        return r
