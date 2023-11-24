#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>, Philip Leindecker <pleindecker@ethz.ch>

import pint
import numpy as np
from numpy.typing import ArrayLike
from .util import default_units


class ElectricField:
    def __init__(self, frequency: pint.Quantity = None, wavelength: pint.Quantity = None, detuning: pint.Quantity = 0, _ureg=None) -> None:
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg

        # Initialize frequency from given frequency or calculate it from given wavelength and detuning
        self._detuning = detuning
        if frequency is not None:
            self._frequency = frequency
        elif wavelength is not None:
            self._frequency = self._ureg("c") / wavelength
        else:
            raise ValueError("Either frequency or wavelength must be provided.")
        
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


class GaussianBeam(ElectricField):

    def __init__(
        self,
        frequency: pint.Quantity=None,
        wavelength: pint.Quantity=None,
        phi: pint.Quantity=None,
        gamma: pint.Quantity=None,
        alpha: pint.Quantity=None,
        theta: pint.Quantity=None,
        power: pint.Quantity=None,
        waist: float | list[float]=None,
        detuning: pint.Quantity=None,
        _ureg: pint.UnitRegistry=None,
    ):
        """
        Args:
            frequency: Frequency of the laser (NOT ANGULAR FREQUENCY) - One can provide instead wavelength.
            wavelength: Wavelength of the laser - One can provide instead frequency.
            phi: Angle between the laser beam and the magnetic field, where the B field is defined to be aligned with z (Radians).
            gamma: Angle between the polarization vector and the plane defined by n and B (Radians).
            alpha: Phase of the polarization vector (Radians). For Left or Right circularily polarised state use +- pi/2.
            theta: Angle of the laser beam to fully define the direction in 3D.
            power: Power of the laser - One can provide instead intensity.
            waist: Waist of the laser (either single value or array of [waist_x, waist_y].
            detuning: Detuning of the laser: (its detuning in MHz not in 2*pi*MHz)
            _ureg: Unit registry
        """

        super().__init__(frequency, wavelength, detuning, _ureg)
        self._phi = phi
        self._gamma = gamma
        self._alpha = alpha
        self._theta = theta
        self._power = power
        self._waist = waist
        self._detuning = detuning


    # JSON METHODS

    @classmethod
    def from_json(cls, json_data: dict, _ureg: pint.UnitRegistry):

        def parse_unit_value(json_data: dict, value_name: str):
            return json_data[value_name]["value"] * _ureg(json_data[value_name]["units"])
        
        if 'frequency' in json_data:
            frequency = parse_unit_value(json_data, 'frequency')
        elif 'wavelength' in json_data:
            wavelength = parse_unit_value(json_data, 'wavelength')
        phi = parse_unit_value(json_data, 'phi')
        gamma = parse_unit_value(json_data, 'gamma')
        alpha = parse_unit_value(json_data, 'alpha')
        theta = parse_unit_value(json_data, 'theta')
        power = parse_unit_value(json_data, 'power')
        if 'value' in json_data['waist']:
            waist = parse_unit_value('waist')
        else:
            waist = [parse_unit_value(json_data['waist'], 'wx'), parse_unit_value(json_data['waist'], 'wy')]
        detuning = parse_unit_value('detuning')

        return cls(frequency, wavelength, phi, gamma, alpha, theta, power, waist, detuning)
    

    def to_json(self):
        def serialize_unit_value(quantity):
            return {'value': quantity.magnitude, 'units': str(quantity.units)}

        data = {
            'wavelength': serialize_unit_value(self.wavelength),
            'phi': serialize_unit_value(self.phi),
            'gamma': serialize_unit_value(self.gamma),
            'alpha': serialize_unit_value(self.alpha),
            'theta': serialize_unit_value(self.theta),
            'power': serialize_unit_value(self.power),
            'detuning': serialize_unit_value(self.detuning)
        }
        if isinstance(self.waist, list):
            data['waist'] = {
                'wx': serialize_unit_value(self.waist[0]),
                'wy': serialize_unit_value(self.waist[1])
            }
        else:
            data['waist'] = serialize_unit_value(self.waist)

        return data

    
    # PROPERTIES

    def set_attribute(self, attribute_name, value):
        """
        A method to safely set the attributes of the class instance,
        ensuring the beam's validity after each setting.

        Args:
        attribute_name (str): The name of the attribute to change.
        value: The value to set the attribute to.
        """
        try:
            setattr(self, '_' + attribute_name, value)
            self._beam_is_valid()
        except ValueError as e:
            print(f"An error occurred: {e}")

    @property
    def phi(self):
        return self._phi.to('rad')
    
    @phi.setter
    @default_units("rad")
    def phi(self, value):
        self.set_attribute('phi', value)

    @property
    def gamma(self):
        return self._gamma.to('rad')
    
    @gamma.setter
    @default_units("rad")
    def gamma(self, value):
        self.set_attribute('gamma', value)

    @property
    def alpha(self):
        return self._alpha.to('rad')
    
    @alpha.setter
    @default_units("rad")
    def alpha(self, value):
        self.set_attribute('alpha', value)

    @property
    def theta(self):
        return self._theta.to('rad')
    
    @theta.setter
    @default_units("rad")
    def theta(self, value):
        self.set_attribute('theta', value)

    @property
    def power(self):
        return (self._power).to("mW")
    
    @power.setter
    @default_units("W")
    def power(self, value):
        self._power = value

    @property
    def waist(self):
        return self._waist
    
    @waist.setter
    def waist(self, value):
        self._waist = value

    # CALCULATED PROPERTIES

    @property
    def area(self):
        # Elliptic Beam
        if isinstance(self.waist, list) and len(self.waist) == 2:
            return np.pi * self.waist[0] * self.waist[1]
        # Circular Beam
        else:
            return np.pi * self.waist**2

    @property
    def propagation_vector(self):
        # TODO: Use the theta vector
        phi = self.phi.to('rad').m
        # theta = self.theta.to('rad').m
        return np.array([np.sin(phi), 0, np.cos(phi)])
    
    @property
    def polarization_vector(self):
        phi = self.phi.to('rad').m
        gamma = self.gamma.to('rad').m
        alpha = self.alpha.to('rad').m
        return np.array([-np.cos(gamma)*np.cos(phi), np.exp(1j*alpha)*np.sin(gamma), np.cos(gamma)*np.sin(phi)])

    @property
    def intensity(self):
        return (2 * self.power / self.area).to("mW/mm^2")

    @property
    def wavevector(self):
        kappa = np.asarray(self.propagation_vector) / np.linalg.norm(self.propagation_vector)
        return kappa * self.angular_frequency / self._ureg("c")
    
    @property
    def field_amplitude(self):
        return (np.sqrt(2 * self.intensity / self._ureg("c*epsilon_0"))).to("V/m")
    
    # Methods
    # TODO: Change them to properties?
    
    def field(self):
        epsilon = np.asarray(self.polarization_vector) / np.linalg.norm(self.polarization_vector)
        return epsilon * self.field_amplitude
    
    def gradient(self):
        return np.einsum("i,...j->...ij", 1j * self.wavevector, self.field())
    
    # PRIVATE METHODS

    def _beam_is_valid(self):
        if np.abs(np.dot(self.polarization_vector, self.propagation_vector)) >= 1e-6:
            raise ValueError("Invalid beam: Polarization is not perpendicular to wavevector.")


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

    def field(self, x, y, z):
        return self._field_a.field(x, y, z) + self._field_b.field(x, y, z)

    def gradient(self, x, y, z):
        return self._field_a.gradient(x, y, z) + self._field_b.gradient(x, y, z)
