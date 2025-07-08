import pint
from pint import Quantity, UnitRegistry
import numpy as np
from .utils.utils import default_units, set_default_units
from abc import ABC, abstractmethod


class ElectricField(ABC):
    def __init__(self, frequency: Quantity, _ureg=None) -> None:
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg
        self.frequency = frequency

    @property
    def frequency(self) -> Quantity:
        return self._frequency.to("THz")

    @frequency.setter
    @default_units("THz")
    def frequency(self, value: Quantity):
        self._frequency = value

    @property
    def angular_frequency(self) -> Quantity:
        return self.frequency.to("1/s") * self._ureg("_2pi")

    @angular_frequency.setter
    @default_units("_2pi/s")
    def angular_frequency(self, value: Quantity):
        self.frequency = value / self._ureg("_2pi")

    @property
    def field_amplitude(self):
        return self._field_amplitude.to("V/m")

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

    @abstractmethod
    def field(self, x, y, z):
        raise NotImplementedError

    @abstractmethod
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
        frequency: Quantity,
        waist: Quantity,
        power: Quantity,
        polarization: np.ndarray | list,
        direction_of_propagation: np.ndarray | list,
        _ureg: UnitRegistry,
    ):
        super().__init__(frequency, _ureg=_ureg)

        self._waist = waist
        self._power = power
        self._update_field_amplitude()
        self._polarization = np.asarray(polarization)
        self._polarization = self._polarization / np.linalg.norm(self._polarization)
        self.direction_of_propagation = np.asarray(direction_of_propagation) / np.linalg.norm(np.asarray(direction_of_propagation))
        self._validate_beam()

    def _validate_beam(self):
        if np.abs(np.dot(self.polarization, self.direction_of_propagation)) >= 1e-6:
            raise ValueError("Polarization must be perpendicular to the wavevector")

    @classmethod
    def from_json(cls, json_data: dict, _ureg: UnitRegistry):
        def parse_unit_value(json_data: dict, key: str) -> Quantity:
            return _ureg.Quantity(json_data[key]["value"], json_data[key]["units"])

        return cls(
            frequency=parse_unit_value(json_data, "frequency"),
            waist=parse_unit_value(json_data, "waist"),
            power=parse_unit_value(json_data, "power"),
            polarization=np.array(json_data["polarization"]),
            direction_of_propagation=np.array(json_data["direction_of_propagation"]),
            _ureg=_ureg,
        )

    def to_json(self):
        def serialize_unit_value(quantity):
            return {"value": quantity.magnitude, "units": str(quantity.units)}

        data = {
            "frequency": serialize_unit_value(self.frequency),
            "waist": serialize_unit_value(self.waist),
            "power": serialize_unit_value(self.power),
            "polarization": self.polarization.tolist(),
            "direction_of_propagation": self.direction_of_propagation.tolist(),
        }
        return data

    @property
    def waist(self) -> Quantity:
        """Returns the waist of the Gaussian beam."""
        return self._waist

    @waist.setter
    @default_units("um")
    def waist(self, value):
        self._waist = value
        self._update_field_amplitude()

    @property
    def power(self) -> Quantity:
        """Returns the power of the Gaussian beam in W, mW, uW, nW, or pW depending on the magnitude of the power. This would be the total power of the beam. The one that you would measure with a powermeter in a lab."""
        return self._power.to_compact()

    @power.setter
    def power(self, value: Quantity):
        self._power = value
        self._update_field_amplitude()

    @property
    def polarization(self) -> np.ndarray:
        """Returns the Jones polarization vector of the electric field in cartesian coordinates."""
        return self._polarization

    @polarization.setter
    def polarization(self, value):
        self._polarization = value / np.linalg.norm(value)
        self._validate_beam()

    @property
    def direction_of_propagation(self):
        """Returns the direction of propagation of the Gaussian beam, which is effectively unitless, unit wavevector."""
        return self._direction_of_propagation

    @direction_of_propagation.setter
    def direction_of_propagation(self, value: np.ndarray | list):
        self._direction_of_propagation = np.asarray(value) / np.linalg.norm(value)
        self._validate_beam()

    @staticmethod
    def calculate_intensity(power: Quantity, waist: Quantity) -> Quantity:
        """Calculate the peak intensity of a Gaussian beam given the power and waist."""
        return 2 * power / (np.pi * waist**2)

    def _update_field_amplitude(self):
        intensity = self.calculate_intensity(self._power, self._waist)
        self._field_amplitude = np.sqrt(2 * intensity / self._ureg("c*epsilon_0"))

    @property
    def wavevector(self) -> np.ndarray:
        """Returns the wavevector of the Gaussian beam."""
        return self.direction_of_propagation * self.angular_frequency / self._ureg("c")

    def field(self):
        return self.polarization * self._field_amplitude

    def gradient(self):
        return np.einsum("i,...j->...ij", 1j * self.wavevector, self.field())


class SumElectricField(ElectricField):
    def __init__(self, field_a: ElectricField, field_b: ElectricField):
        super().__init__(field_a.frequency, field_a._ureg)
        self._field_a = field_a
        self._field_b = field_b

    def field(self, x=0, y=0, z=0):
        return self._field_a.field(x, y, z) + self._field_b.field(x, y, z)

    def gradient(self, x=0, y=0, z=0):
        return self._field_a.gradient(x, y, z) + self._field_b.gradient(x, y, z)
