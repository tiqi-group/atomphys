#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import pint
import numpy as np
from math import pi
from numpy.typing import ArrayLike

from .util import default_units, make_alias, make_alias_with_setter


class ElectricField:
    def __init__(self, frequency: pint.Quantity, _ureg=None) -> None:
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg
        self._frequency = frequency

    @property
    def frequency(self) -> pint.Quantity:
        return self._frequency.to('THz')
    
    @frequency.setter
    @default_units('THz')
    def frequency(self, value: pint.Quantity):
        self._frequency = value
    
    @property
    def angular_frequency(self) -> pint.Quantity:
        return (self.frequency.to('1/s'))*self._ureg('2*pi')
    
    @angular_frequency.setter
    @default_units('2pi/s')
    def angular_frequency(self, value: pint.Quantity):
        self.frequency = value / self._ureg('2*pi')
    


    nu = make_alias_with_setter('frequency')
    ν = make_alias_with_setter('frequency')
    ω = make_alias_with_setter('angular_frequency')
    omega = make_alias_with_setter('angular_frequency')


    def field(self, x, y, z):
        raise NotImplementedError

    def gradient(self, x, y, z):
        raise NotImplementedError


    def __add__(self, other):
        if not isinstance(other, ElectricField):
            raise TypeError('Both objects must be an instance of ElectricField')
        if not self.frequency == other.frequency:
            raise ValueError('Can only sum fields at the same frequency')
        return SumElectricField(self, other)

    @staticmethod
    def _ravel_coords(*args):
        args = np.broadcast_arrays(*args)
        shape = args[0].shape
        args = list(map(np.ravel, args))
        X = np.stack(args, axis=1).astype(float)
        return shape, X


class SumElectricField(ElectricField):
    def __init__(self, field_a: ElectricField, field_b: ElectricField):
        super().__init__(field_a.wavelength, field_a._ureg)
        self._field_a = field_a
        self._field_b = field_b

    def field(self, x, y, z):
        return self._field_a.field(x, y, z) + self._field_b.field(x, y, z)

    def gradient(self, x, y, z):
        return self._field_a.gradient(x, y, z) + self._field_b.gradient(x, y, z)


class PlaneWaveElectricField(ElectricField):
    def __init__(self, E0: float, polarization: ArrayLike, wavevector: ArrayLike,
                 frequency: pint.Quantity, _ureg=None) -> None:
        super().__init__(frequency, _ureg)
        assert np.dot(polarization, wavevector) == 0, "Polarization must be perpendicular to wavevector"
        self._epsilon = np.asarray(polarization) / np.linalg.norm(polarization)
        self._kappa = np.asarray(wavevector) / np.linalg.norm(wavevector)
        self._field_amplitude = E0

    def _phase(self, X):
        xk = (np.dot(X, self.k).to('')).magnitude[0]
        return np.exp(1j * xk)

    def phase(self, x, y, z):
        shape, X = self._ravel_coords(x, y, z)
        return self._phase(X).reshape(shape)

    def field(self, x, y, z):
        shape, X = self._ravel_coords(x, y, z)
        return self._epsilon.reshape((1,) * len(shape) + (-1,)) * self._field_amplitude * self._phase(X).reshape(shape + (1,))

    def gradient(self, x, y, z):
        # outer product
        return np.einsum('i,...j->...ij', 1j * self.wavevector, self.field(x, y, z))

    @property
    def wavelength(self):
        return (self._ureg('c') / self.frequency).to('nm')
    
    @wavelength.setter
    @default_units('nm')
    def wavelength(self, value):
        self.frequency = self._ureg('c') / value
    
    @property
    def k(self):
        return (self._kappa/self.wavelength).to('1/nm')* self._ureg('2*pi')

    @property
    def polarization(self):
        return self._epsilon
    
    @polarization.setter
    def polarization(self, value):
        self._epsilon = value / np.linalg.norm(value)
    
    
    λ = make_alias_with_setter('wavelength')
    E0 = make_alias_with_setter('_field_amplitude')
    wavevector = make_alias('k')
    epsilon = make_alias_with_setter('polarization')
    eps = make_alias_with_setter('polarization')
    ε = make_alias_with_setter('polarization')


