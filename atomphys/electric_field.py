#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pint
import numpy as np
from .utils.utils import default_units, set_default_units


class ElectricField:
    def __init__(self, frequency: pint.Quantity, _ureg=None) -> None:
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg
        self.frequency = frequency

    @property
    def frequency(self) -> pint.Quantity:
        return (self._frequency).to("THz")

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

    # def field(self, x, y, z):
    #     raise NotImplementedError

    # def gradient(self, x, y, z):
    #     raise NotImplementedError

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
        frequency: pint.Quantity,
        waist: pint.Quantity,
        power: pint.Quantity,
        polarization,
        direction_of_propagation,
        _ureg: pint.UnitRegistry,
    ):
        # Call the ElectricField constructor
        super().__init__(frequency, _ureg=_ureg)

        if frequency is None:
            raise ValueError("Must specify frequency")
        if (power is None or waist is None):
            raise ValueError("Must specify both power and waist")
        if (polarization is None or direction_of_propagation is None):
            raise ValueError("Must specify polarization and direction of propagation")

        self.frequency = frequency
        self._power = power
        self._waist = waist
        self.power = power
        self.waist = waist
        if polarization is not None and direction_of_propagation is not None:
            self.polarization = np.asarray(polarization) 
            self.direction_of_propagation = np.asarray(direction_of_propagation) 
        assert (
            np.abs(np.dot(self._polarization, self._direction_of_propagation)) < 1e-6
        ), "Polarization must be perpendicular to wavevector"

    @staticmethod
    def calculate_intensity(power, waist):
        """ Calculate the peak intensity of a Gaussian beam given the power and waist """
        return 2 * power / (np.pi * waist**2)

    @property
    def frequency(self) -> pint.Quantity:
        return (self._frequency).to("THz")

    @frequency.setter
    @default_units("THz")
    def frequency(self, value: pint.Quantity):
        self._frequency = value

    @property
    def wavelength(self):
        """ Returns the wavelength (in vacuum) of the Gaussian beam in nm """
        return (self._ureg("c") / self.frequency).to("nm")

    @property
    def direction_of_propagation(self):
        """ Returns the direction of propagation of the Gaussian beam, which is effectively unitless, unit wavevector. """
        return self._direction_of_propagation

    @direction_of_propagation.setter
    def direction_of_propagation(self, value):
        self._direction_of_propagation = value / np.linalg.norm(value)

    @property
    def polarization(self):
        """ Returns the Jones polarization vector of the electric field in cartesian coordinates."""
        return self._polarization

    @polarization.setter
    def polarization(self, value):
        self._polarization = value / np.linalg.norm(value)

    @property
    def intensity(self):
        """ Returns the peak intensity of the Gaussian beam in W/cm^2 """
        return self._field_amplitude**2 * self._ureg("c*epsilon_0") / 2

    @property
    def power(self):
        """ Returns the power of the Gaussian beam in W, mW, uW, nW, or pW depending on the magnitude of the power. This would be the total power of the beam. The one that you would measure with a powermeter in a lab."""
        power_in_watts = self._power.to("W").magnitude

        if power_in_watts >= 1:
            return self._power.to("W")
        elif power_in_watts >= 1e-3:
            return self._power.to("mW")
        elif power_in_watts >= 1e-6:
            return self._power.to("uW")
        elif power_in_watts >= 1e-9:
            return self._power.to("nW")
        else:
            return self._power.to("pW")

    @power.setter
    def power(self, value):
        self._power = value
        intensity = self.calculate_intensity(self._power, self._waist)
        self._field_amplitude = np.sqrt(2 * intensity / self._ureg("c*epsilon_0"))

    @property
    def waist(self):
        """ Returns the waist of the Gaussian beam. """
        return self._waist

    @waist.setter
    @default_units("um")
    def waist(self, value):
        self._waist = value
        intensity = self.calculate_intensity(self._power, self._waist)
        self._field_amplitude = np.sqrt(2 * intensity / self._ureg("c*epsilon_0"))

    @property
    def wavevector(self):
        """ Returns the wavevector of the Gaussian beam."""
        wavevector = self.direction_of_propagation * self.angular_frequency / self._ureg("c")
        return wavevector

    def field(self):
        return self.polarization * self._field_amplitude

    def gradient(self):
        return np.einsum("i,...j->...ij", 1j * self.wavevector, self.field())


# class PlaneWaveElectricField(ElectricField):
#     def __init__(
#         self,
#         E0: float,
#         polarization: ArrayLike,
#         wavevector: ArrayLike,
#         frequency: pint.Quantity,
#         _ureg=None,
#     ) -> None:
#         super().__init__(frequency, _ureg)
#         assert (
#             np.dot(polarization, wavevector) == 0
#         ), "Polarization must be perpendicular to wavevector"
#         self._epsilon = np.asarray(polarization) / np.linalg.norm(polarization)
#         self._kappa = np.asarray(wavevector) / np.linalg.norm(wavevector)
#         self._field_amplitude = E0

#     def _phase(self, X):
#         xk = (np.dot(X, self.k).to("")).magnitude[0]
#         xk = (np.dot(X, self.k).to("")).magnitude[0]
#         return np.exp(1j * xk)

#     def phase(self, x=0, y=0, z=0):
#         shape, X = self._ravel_coords(x, y, z)
#         return self._phase(X).reshape(shape)

#     def field(self, x=0, y=0, z=0):
#         shape, X = self._ravel_coords(x, y, z)
#         return (
#             self._epsilon.reshape((1,) * len(shape) + (-1,))
#             * self._field_amplitude
#             * self._phase(X).reshape(shape + (1,))
#         )

#     def gradient(self, x=0, y=0, z=0):
#         # outer product
#         return np.einsum("i,...j->...ij", 1j * self.wavevector, self.field(x, y, z))

#     @property
#     def wavelength(self):
#         return (self._ureg("c") / self.frequency).to("nm")

#     @wavelength.setter
#     @default_units("nm")
#     def wavelength(self, value):
#         self.frequency = self._ureg("c") / value

#     @property
#     def k(self):
#         return (self._kappa / self.wavelength).to("1/nm") * self._ureg("2*pi")

#     @property
#     def polarization(self):
#         return self._epsilon


#     @polarization.setter
#     def polarization(self, value):
#         self._epsilon = value / np.linalg.norm(value)


class SumElectricField(ElectricField):
    def __init__(self, field_a: ElectricField, field_b: ElectricField):
        super().__init__(field_a.frequency, field_a._ureg)
        self._field_a = field_a
        self._field_b = field_b

    def field(self, x=0, y=0, z=0):
        return self._field_a.field(x, y, z) + self._field_b.field(x, y, z)

    def gradient(self, x=0, y=0, z=0):
        return self._field_a.gradient(x, y, z) + self._field_b.gradient(x, y, z)
