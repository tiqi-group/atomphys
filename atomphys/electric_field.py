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
        return self._frequency.to("THz")

    @frequency.setter
    @default_units("THz")
    def frequency(self, value: pint.Quantity):
        self._frequency = value

    @property
    def angular_frequency(self) -> pint.Quantity:
        return (self.frequency.to("1/s")) * self._ureg("2*pi")

    @angular_frequency.setter
    @default_units("2pi/s")
    def angular_frequency(self, value: pint.Quantity):
        self.frequency = value / self._ureg("2*pi")

    nu = make_alias_with_setter("frequency")
    ν = make_alias_with_setter("frequency")
    ω = make_alias_with_setter("angular_frequency")
    omega = make_alias_with_setter("angular_frequency")

    def field(self, x, y, z):
        raise NotImplementedError

    def gradient(self, x, y, z):
        raise NotImplementedError

    def __add__(self, other):
        if not isinstance(other, ElectricField):
            raise TypeError("Both objects must be an instance of ElectricField")
        if not self.frequency == other.frequency:
            raise ValueError("Can only sum fields at the same frequency")
        return SumElectricField(self, other)

    @staticmethod
    def _ravel_coords(*args):
        args = np.broadcast_arrays(*args)
        shape = args[0].shape
        args = list(map(np.ravel, args))
        X = np.stack(args, axis=1).astype(float)
        return shape, X


class LaserField(ElectricField):
    def __init__(
        self,
        polarization,
        direction_of_propagation,
        frequency=None,
        wavelength=None,
        intensity=None,
        power=None,
        waist=None,
        detuning=None,
        _ureg=None,
    ):
        # Call the ElectricField constructor
        super().__init__(frequency, _ureg)

        # Make sure the user has specified either intensity or both power and waist
        if intensity is None and (power is None or waist is None):
            raise ValueError("Must specify either intensity or both power and waist")

        # If the user specified intensity, use it
        if intensity is not None:
            self.intensity = intensity

        # If the user specified power and waist, calculate the intensity
        if power is not None and waist is not None:
            self.intensity = self.calculate_intensity(power, waist)

        if frequency is not None:
            self.frequency = frequency
        elif wavelength is not None:
            self.wavelength = wavelength

        self._power = power
        self._waist = waist
        assert np.dot(polarization, direction_of_propagation) == 0, "Polarization must be perpendicular to wavevector"
        self._epsilon = np.asarray(polarization) / np.linalg.norm(polarization)
        self._kappa = np.asarray(direction_of_propagation) / np.linalg.norm(direction_of_propagation)

        if detuning is not None:
            self._detuning = detuning
        else:
            self._detuning = 0 * self._ureg("MHz")
        # Calculates electric field amplitude from intensity

    @staticmethod
    def calculate_intensity(power, waist):
        # Use the formula for the intensity of a Gaussian beam:
        # I = 2P/(pi*w^2)
        import numpy as np

        return 2 * power / (np.pi * waist**2)

    @property
    def detuning(self):
        return self._detuning

    @detuning.setter
    @default_units("MHz")
    def detuning(self, value):
        self._detuning = value

    @property
    def frequency(self) -> pint.Quantity:
        return (self._frequency - self._detuning).to("THz")

    @frequency.setter
    @default_units("THz")
    def frequency(self, value: pint.Quantity):
        self._frequency = value

    @property
    def wavelength(self):
        return self._ureg("c") / self._frequency

    @wavelength.setter
    def wavelength(self, value):
        self._frequency = self._ureg("c") / value

    @property
    def field_amplitude(self):
        return (self._field_amplitude).to("V/m")

    @field_amplitude.setter
    @default_units("V/m")
    def field_amplitude(self, value):
        self._field_amplitude = value

    @property
    def intensity(self):
        return (self._field_amplitude**2 * self._ureg("c*epsilon_0") / 2).to("mW/mm^2")

    @intensity.setter
    @default_units("W/cm^2")
    def intensity(self, value):
        self._field_amplitude = np.sqrt(2 * value / self._ureg("c*epsilon_0"))

    @property
    def power(self):
        return (self._power).to("mW")

    @power.setter
    @default_units("W")
    def power(self, value):
        self._power = value
        if self.waist is not None:
            self.intensity = self.calculate_intensity(value, self.waist)

    @property
    def waist(self):
        return self._waist

    @waist.setter
    @default_units("um")
    def waist(self, value):
        self._waist = value
        if self.power is not None:
            self.intensity = self.calculate_intensity(self.power, value)

    @property
    def wavevector(self):
        return self._kappa * self.angular_frequency / self._ureg("c")

    k = make_alias("wavevector")

    def field(self):
        return self._epsilon * self._field_amplitude

    def gradient(self):
        return np.einsum("i,...j->...ij", 1j * self.wavevector, self.field())


class PlaneWaveElectricField(ElectricField):
    def __init__(
        self, E0: float, polarization: ArrayLike, wavevector: ArrayLike, frequency: pint.Quantity, _ureg=None
    ) -> None:
        super().__init__(frequency, _ureg)
        assert np.dot(polarization, wavevector) == 0, "Polarization must be perpendicular to wavevector"
        self._epsilon = np.asarray(polarization) / np.linalg.norm(polarization)
        self._kappa = np.asarray(wavevector) / np.linalg.norm(wavevector)
        self._field_amplitude = E0

    def _phase(self, X):
        xk = (np.dot(X, self.k).to("")).magnitude[0]
        return np.exp(1j * xk)

    def phase(self, x, y, z):
        shape, X = self._ravel_coords(x, y, z)
        return self._phase(X).reshape(shape)

    def field(self, x, y, z):
        shape, X = self._ravel_coords(x, y, z)
        return (
            self._epsilon.reshape((1,) * len(shape) + (-1,))
            * self._field_amplitude
            * self._phase(X).reshape(shape + (1,))
        )

    def gradient(self, x, y, z):
        # outer product
        return np.einsum("i,...j->...ij", 1j * self.wavevector, self.field(x, y, z))

    @property
    def wavelength(self):
        return (self._ureg("c") / self.frequency).to("nm")

    @wavelength.setter
    @default_units("nm")
    def wavelength(self, value):
        self.frequency = self._ureg("c") / value

    @property
    def k(self):
        return (self._kappa / self.wavelength).to("1/nm") * self._ureg("2*pi")

    @property
    def polarization(self):
        return self._epsilon

    @polarization.setter
    def polarization(self, value):
        self._epsilon = value / np.linalg.norm(value)

    λ = make_alias_with_setter("wavelength")
    E0 = make_alias_with_setter("_field_amplitude")
    wavevector = make_alias("k")
    epsilon = make_alias_with_setter("polarization")
    eps = make_alias_with_setter("polarization")
    ε = make_alias_with_setter("polarization")


class SumElectricField(ElectricField):
    def __init__(self, field_a: ElectricField, field_b: ElectricField):
        super().__init__(field_a.wavelength, field_a._ureg)
        self._field_a = field_a
        self._field_b = field_b

    def field(self, x, y, z):
        return self._field_a.field(x, y, z) + self._field_b.field(x, y, z)

    def gradient(self, x, y, z):
        return self._field_a.gradient(x, y, z) + self._field_b.gradient(x, y, z)
