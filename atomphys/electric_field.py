#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>, Philip Leindecker <pleindecker@ethz.ch>

import pint
import numpy as np
from numpy.typing import ArrayLike
from .util import default_units, set_default_units, make_alias, make_alias_with_setter


class ElectricField:
    def __init__(self, frequency: pint.Quantity, E0: pint.Quantity, _ureg=None) -> None:
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
    def wavelength(self):
        return (self._ureg("c") / self.frequency).to("nm")

    @wavelength.setter
    @default_units("nm")
    def wavelength(self, value):
        self.frequency = self._ureg("c") / value

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
        frequency: pint.Quantity = None,
        wavelength: pint.Quantity = None,
        detuning: pint.Quantity = None,
        intensity: pint.Quantity = None,
        power: pint.Quantity = None,
        waist: float | list[float] = None,
        polarization = None,
        direction_of_propagation = None,
        phi: pint.Quantity = None,
        gamma: pint.Quantity = None,
        theta: pint.Quantity = None,
        alpha: pint.Quantity = None,
        _ureg: pint.UnitRegistry = None,
    ):
        """
        Gaussian Beam class is designed in such way to handle hierarchal redundancies. 
        The idea is that one can define laser in multiple ways. 
        For instance, one can define the polarization vector and the direction of propagation or 
        one can define the angle of the laser beam and the polarization vector. 
        The class will then calculate the missing parameters. 

        Construction of this hierarchal redundancies is following, 
        when missing a base property, the base property will be calculated from the remaining properties. Base property are listed below:
        e.g. frequency is the base property, and wavelength will be only used when frequency is missing

        (frequency) <- (wavelength)         (+) (detuning)
        (intensity) <- (power, waist)
        (polarization, direction_of_propagation) <- (phi, gamma, alpha)


        References:
            - 
            - Transformation of phi, gamma, alpha to propagation vector from Gillen Beck's Master Thesis (p. 26)
            - Transformation of phi, gamma, alpha to polarization vector from Roland's Matts PhD Thesis (p. 19)


        Args:
            frequency: Frequency of the laser (NOT ANGULAR FREQUENCY)
            wavelength: Wavelength of the laser
            detuning: Detuning of the laser: (its detuning in MHz not in 2*pi*MHz)

            intensity: Intensity of the laser [Pint Quantity]
            power: Power of the laser [Pint Quantity]
            waist: Waist of the laser [Pint Quantity]

            polarization: Polarization vector of the laser beam - doesn't have to be normalized, as it will be later normalized
            direction_of_propagation: Direction of the propagation of the laser beam - doesn't have to be normalized, as it will be later normalized
            phi: angle between the axis of the k-vector and the quantization axis
            gamma: angle between the axis of the polarization elipse and the quantization axis
            alpha: phase of the polarization vector [Fill out how to define right and left circularly polarized states]

            _ureg: Unit registry
        """

        # Call the ElectricField constructor
        super().__init__(frequency, E0=0, _ureg=_ureg)

        if detuning is not None:
            self._detuning = detuning
        else:
            self._detuning = 0 * self._ureg("MHz")

        if frequency is None and wavelength is None:
            raise ValueError("Must specify either frequency or wavelength")
        
        if intensity is None and (power is None or waist is None):
            raise ValueError("Must specify either intensity or both power and waist")
        
        if (polarization is None or direction_of_propagation is None) and (phi is None or gamma is None or alpha is None):
            raise ValueError("Must specify either polarization and direction of propagation or phi, gamma and alpha")

        if frequency is not None:
            self.frequency = frequency
        elif wavelength is not None:
            self.wavelength = wavelength

        # If the user specified intensity, use it
        if intensity is not None:
            self.intensity = intensity
        elif power is not None and waist is not None:
            self.intensity = self.calculate_intensity(power, waist)

        self._power = power
        self._waist = waist
        if polarization is not None and direction_of_propagation is not None:
            self._epsilon = np.asarray(polarization) / np.linalg.norm(polarization)
            self._n = np.asarray(direction_of_propagation) / np.linalg.norm(
                direction_of_propagation
            )
        elif phi is not None and gamma is not None and alpha is not None:
            self.phi = phi
            self.gamma = gamma
            self.alpha = alpha
            self.calculate_n_and_epsilon()
        
        assert (
            np.abs(np.dot(self._epsilon, self._n)) < 1e-6
        ), "Polarization must be perpendicular to wavevector"

    @staticmethod
    def calculate_n_and_epsilon(self):
        if self.phi is not None and self.gamma is not None and self.alpha is not None:
            phi = self.phi.m
            gamma = self.gamma.m
            alpha = self.alpha.m
            self.epsilon = np.array([-np.cos(gamma)*np.cos(phi), np.exp(1j*alpha)*np.sin(gamma), np.cos(gamma)*np.sin(phi)])
            self.n = np.array([np.sin(phi), 0, np.cos(phi)])
            assert (
                np.abs(np.dot(self._epsilon, self._n)) < 1e-6
            ), "Polarization must be perpendicular to wavevector"

    @staticmethod
    def calculate_intensity(power, waist):
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
        return (self._frequency + self._detuning).to("THz")

    @frequency.setter
    @default_units("THz")
    def frequency(self, value: pint.Quantity):
        self._frequency = value


    @property
    def n(self):
        return self._n
    
    @n.setter
    def n(self, value):
        self._n = value / np.linalg.norm(value)
    
    @property
    def epsilon(self):
        return self._epsilon
    
    @epsilon.setter
    def epsilon(self, value):
        self._epsilon = value / np.linalg.norm(value)

    @property
    def gamma(self):
        return self._gamma.to('rad')
    
    @gamma.setter
    @default_units("rad")
    def gamma(self, value):
        self.set_attribute('gamma', value)
        self.calculate_n_and_epsilon()

    @property
    def alpha(self):
        return self._alpha.to('rad')
    
    @alpha.setter
    @default_units("rad")
    def alpha(self, value):
        self.set_attribute('alpha', value)
        self.calculate_n_and_epsilon()

    @property
    def theta(self):
        return self._theta.to('rad')
    
    @theta.setter
    @default_units("rad")
    def theta(self, value):
        self.set_attribute('theta', value)
        self.calculate_n_and_epsilon()

    @property
    def intensity(self):
        return self._field_amplitude**2 * self._ureg("c*epsilon_0") / 2
    
    @intensity.setter
    @default_units("W/cm^2")
    def intensity(self, value):
        self._field_amplitude = np.sqrt(2 * value / self._ureg("c*epsilon_0"))

    @property
    def power(self):
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
    @default_units("W")
    def power(self, value):
        self._power = value

    @property
    def waist(self):
        return self._waist

    @waist.setter
    @default_units("um")
    def waist(self, value):
        self._waist = value

    @property
    def wavevector(self):
        wavevector = self.n * self.angular_frequency / self._ureg("c")
        return wavevector

    def field(self):
        return self._epsilon * self._field_amplitude

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
        super().__init__(field_a.wavelength, field_a._ureg)
        self._field_a = field_a
        self._field_b = field_b

    def field(self, x=0, y=0, z=0):
        return self._field_a.field(x, y, z) + self._field_b.field(x, y, z)

    def gradient(self, x=0, y=0, z=0):
        return self._field_a.gradient(x, y, z) + self._field_b.gradient(x, y, z)
