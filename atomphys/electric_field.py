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

    def field(self):
        raise NotImplementedError

    def gradient(self):
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


class PlaneWaveElectricField(ElectricField):
    def __init__(self, field: float, polarization: ArrayLike, wavevector: ArrayLike,
                 wavelength: pint.Quantity, _ureg=None) -> None:
        super().__init__(wavelength, _ureg)
        assert np.dot(polarization, wavevector) == 0, "Polarization must be perpendicular to wavevector"
        self._polarization = np.asarray(polarization) / np.linalg.norm(polarization)
        self._wavevector = np.asarray(wavevector) / np.linalg.norm(wavevector)
        self._field = field * self._ureg('V/m')

    def field(self):
        return self._polarization * self._field

    def gradient(self):
        return np.outer(self._polarization, self._wavevector) * 1j * self._field * (2 * pi / self.wavelength.to('m'))

    @property
    def polarization(self):
        return self._polarization

    @property
    def wavevector(self):
        return self._wavevector * (2 * pi / self.wavelength.to('m'))
