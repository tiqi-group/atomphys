#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import pint
import numpy as np
from math import pi
from numpy.typing import ArrayLike

from .util import default_units


class ElectricField:
    def __init__(self, wavelength: pint.Quantity, _ureg=None) -> None:
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg
        self.wavelength = wavelength

    def field(self, x, y, z):
        raise NotImplementedError

    def gradient(self, x, y, z):
        raise NotImplementedError

    @property
    def wavelength(self) -> pint.Quantity:
        return self._wavelength

    @wavelength.setter
    @default_units('nm')
    def wavelength(self, value: pint.Quantity):
        self._wavelength = value

    @property
    def frequency(self):
        return self.wavelength.to('THz', 'sp')

    @property
    def omega(self):
        return 2 * pi * self.frequency

    def __add__(self, other):
        if not isinstance(other, ElectricField):
            raise TypeError('Both objects must be an instance of ElectricField')
        if not self.wavelength == other.wavelength:
            raise ValueError('Can only sum fields at the same wavelength')
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
    def __init__(self, field: float, polarization: ArrayLike, wavevector: ArrayLike,
                 wavelength: pint.Quantity, _ureg=None) -> None:
        super().__init__(wavelength, _ureg)
        assert np.dot(polarization, wavevector) == 0, "Polarization must be perpendicular to wavevector"
        self._epsilon = np.asarray(polarization) / np.linalg.norm(polarization)
        self._kappa = np.asarray(wavevector) / np.linalg.norm(wavevector)
        self._field = field * self._ureg('V/m')

    @default_units('m')
    def _phase(self, X):
        return np.exp(1j * np.dot(X, self.k))

    def phase(self, x, y, z):
        shape, X = self._ravel_coords(x, y, z)
        return self._phase(X).reshape(shape)

    def field(self, x, y, z):
        shape, X = self._ravel_coords(x, y, z)
        return self._epsilon.reshape((1,) * len(shape) + (-1,)) * self._field * self._phase(X).reshape(shape + (1,))

    def gradient(self, x, y, z):
        # outer product
        return np.einsum('i,...j->...ij', 1j * self.wavevector, self.field(x, y, z))

    @property
    def k(self):
        return self._kappa * (2 * pi / self.wavelength.to('m'))

    @property
    def polarization(self):
        return self._epsilon

    @property
    def wavevector(self):
        return self.k  # alias! I should not
