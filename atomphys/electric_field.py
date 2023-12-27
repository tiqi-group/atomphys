#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import pint
import numpy as np
from numpy.typing import ArrayLike

from .util import default_units, set_default_units, make_alias, make_alias_with_setter


class ElectricField:
    def __init__(self, frequency: pint.Quantity, E0: pint.Quantity, _ureg=None) -> None:
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg
        self._frequency = frequency
        self._field_amplitude = set_default_units(E0, "V/m", self._ureg)

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

    @property
    def field_amplitude(self):
        return (self._field_amplitude).to("V/m")

    @field_amplitude.setter
    @default_units("V/m")
    def field_amplitude(self, value):
        self._field_amplitude = value

    @property
    def intensity(self):
        return (self._field_amplitude**2 * self._ureg("c*epsilon_0") / 2).to("mW/cm^2")

    @intensity.setter
    @default_units("W/cm^2")
    def intensity(self, value):
        self._field_amplitude = np.sqrt(2 * value / self._ureg("c*epsilon_0"))

    E0 = make_alias_with_setter("_field_amplitude")
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

    def _ravel_coords(self, *args):
        args = tuple(map(lambda _x: set_default_units(_x, "m", self._ureg), args))
        args = np.broadcast_arrays(*args)
        shape = args[0].shape
        args = list(map(np.ravel, args))
        X = np.stack(args, axis=1).astype(float)
        return shape, X


class GaussianBeam(ElectricField):
    def __init__(
        self,
        polarization=None,
        direction_of_propagation=None,
        phi=None,
        gamma=None,
        alpha=None,
        frequency=None,
        wavelength=None,
        intensity=None,
        power=None,
        waist=None,
        detuning=None,
        _ureg=None,
    ):
        """
        Args:
            polarization (array-like): Polarization vector
            direction_of_propagation (array-like): Direction of propagation of the laser beam
            frequency (pint.Quantity): Frequency of the laser (NOT ANGULAR FREQUENCY) - One can provide instead wavelength
            phi (float): Angle between the laser beam and the magnetic field, where the B field is defined to be aligned with z (Radians)
            gamma (float): Angle between the polarization vector and the plane defined by n and B (Radians)
            alpha (float): Phase of the polarization vector (Radians) as defined in https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/111013/eth-48865-02.pdf?sequence=2&isAllowed=y
                            'For Left or Right circularily polarised state use +- pi/2'
            wavelength (pint.Quantity): Wavelength of the laser - One can provide instead frequency
            intensity (pint.Quantity): Intensity of the laser - One can provide instead power and waist
            power (pint.Quantity): Power of the laser - One can provide instead intensity
            waist (pint.Quantity): Waist of the laser (either single value or array of [waist_x, waist_y] - One can provide instead intensity
            detuning (pint.Quantity): Detuning of the laser: (its detuning in MHz not in 2*pi*MHz) - i.e. Δ = ν-ν0
            _ureg (pint.UnitRegistry): Unit registry
        """

        # Call the ElectricField constructor
        super().__init__(frequency, E0=0, _ureg=_ureg)

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
        if polarization is None and direction_of_propagation is None:
            direction_of_propagation = np.array([np.sin(phi), 0, np.cos(phi)])
            polarization = np.array([-np.cos(gamma)*np.cos(phi), np.exp(1j*alpha)*np.sin(gamma), np.cos(gamma)*np.sin(phi)])

        assert (
            np.abs(np.dot(polarization, direction_of_propagation)) < 1e-6
        ), "Polarization must be perpendicular to wavevector"
        self._epsilon = np.asarray(polarization) / np.linalg.norm(polarization)
        self._kappa = np.asarray(direction_of_propagation) / np.linalg.norm(
            direction_of_propagation
        )

        if detuning is not None:
            self._detuning = detuning
        else:
            self._detuning = 0 * self._ureg("MHz")

    @staticmethod
    def calculate_intensity(power, waist):
        # Use the formula for the intensity of a Gaussian beam:
        # I = 2P/A
        # P ... Power
        # A ... Area of the beam
        # Circular
        area = np.pi * waist**2
        return 2 * power / area

    @property
    def detuning(self):
        return self._detuning


    @detuning.setter
    @default_units("MHz")
    @default_units("MHz")
    def detuning(self, value):
        self._detuning = value

    @property
    def frequency(self) -> pint.Quantity:
        return (self._frequency - self._detuning).to("THz")

    @frequency.setter
    @default_units("THz")
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
        return (self._field_amplitude**2 * self._ureg("c*epsilon_0") / 2).to(
            "mW/mm^2"
        )

    @intensity.setter
    @default_units("W/cm^2")
    def intensity(self, value):
        self._field_amplitude = np.sqrt(2 * value / self._ureg("c*epsilon_0"))

    @property
    def power(self):
        return (self._power).to("mW")

    @power.setter
    @default_units("W")
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
    @default_units("um")
    def waist(self, value):
        self._waist = value
        if self.power is not None:
            self.intensity = self.calculate_intensity(self.power, value)

    @property
    def wavevector(self):
        wavevector = self._kappa * self.angular_frequency / self._ureg("c")
        return wavevector

    k = make_alias("wavevector")

    def field(self, x=0, y=0, z=0):
        return self._epsilon * self._field_amplitude

    def gradient(self, x=0, y=0, z=0):
        return np.einsum("i,...j->...ij", 1j * self.wavevector, self.field(x, y, z))


class PlaneWaveElectricField(ElectricField):
    def __init__(
        self,
        E0: float,
        polarization: ArrayLike,
        wavevector: ArrayLike,
        frequency: pint.Quantity,
        _ureg=None,
    ) -> None:
        super().__init__(frequency, _ureg)
        assert (
            np.dot(polarization, wavevector) == 0
        ), "Polarization must be perpendicular to wavevector"
        self._epsilon = np.asarray(polarization) / np.linalg.norm(polarization)
        self._kappa = np.asarray(wavevector) / np.linalg.norm(wavevector)
        self._field_amplitude = E0

    def _phase(self, X):
        xk = (np.dot(X, self.k).to("")).magnitude[0]
        xk = (np.dot(X, self.k).to("")).magnitude[0]
        return np.exp(1j * xk)

    def phase(self, x=0, y=0, z=0):
        shape, X = self._ravel_coords(x, y, z)
        return self._phase(X).reshape(shape)

    def field(self, x=0, y=0, z=0):
        shape, X = self._ravel_coords(x, y, z)
        return (
            self._epsilon.reshape((1,) * len(shape) + (-1,))
            * self._field_amplitude
            * self._phase(X).reshape(shape + (1,))
        )

    def gradient(self, x=0, y=0, z=0):
        # outer product
        return np.einsum("i,...j->...ij", 1j * self.wavevector, self.field(x, y, z))

    @property
    def wavelength(self):
        return (self._ureg("c") / self.frequency).to("nm")

    @wavelength.setter
    @default_units("nm")
    def wavelength(self, value):
        self.frequency = self._ureg("c") / value

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


class SumElectricField(ElectricField):
    def __init__(self, field_a: ElectricField, field_b: ElectricField):
        super().__init__(field_a.wavelength, field_a._ureg)
        self._field_a = field_a
        self._field_b = field_b

    def field(self, x=0, y=0, z=0):
        return self._field_a.field(x, y, z) + self._field_b.field(x, y, z)

    def gradient(self, x=0, y=0, z=0):
        return self._field_a.gradient(x, y, z) + self._field_b.gradient(x, y, z)
