from math import inf
from math import pi as π

import pint

from . import _ureg
from .util import default_units


class Laser:
    _ureg: pint.UnitRegistry
    __omega = pint.Quantity
    __linewidth = pint.Quantity
    __electric_field = pint.Quantity


    def __init__(self, ureg=None, laser=None, **new_laser):
        if ureg is not None:
            self._ureg = ureg
        else:
            self._ureg = _ureg

        self.omega = 0
        self.__linewidth = 0
        self.__electric_field = 0

        if laser is not None:
            self._ureg = laser._ureg
            self.__omega = laser.__omega
            self.__linewidth = laser.__linewidth
            self.__electric_field = laser.__electric_field


        for attr in new_laser:
            if attr in dir(self):
                self.__setattr__(attr, new_laser[attr])

    def __repr__(self):
        fmt = "0.4g~P"
        return f"Laser(λ={self.λ:{fmt}})"

    # ---------
    # Frequency
    # ---------

    @property
    def omega(self):
        return self.__omega

    @omega.setter
    @default_units("THz")
    def omega(self, ω):
        self.__omega = ω

    @property
    def angular_frequency(self):
        return self.omega

    @angular_frequency.setter
    def angular_frequency(self, ω):
        self.omega = ω

    @property
    def ω(self):
        return self.omega

    @ω.setter
    def ω(self, ω):
        self.omega = ω

    @property
    def ν(self):
        return self.omega / (2 * π)

    @ν.setter
    def ν(self, ν):
        self.omega = 2 * π * ν

    @property
    def nu(self):
        return self.ν

    @nu.setter
    def nu(self, ν):
        self.ν = ν

    @property
    def frequency(self):
        return self.ν

    @frequency.setter
    def frequency(self, ν):
        self.ν = ν

    # ----------
    # Wavelength
    # ----------

    @property
    def wavelength(self):
        c = self._ureg.c
        try:
            return (2 * π * c / self.omega).to("nm")
        except ZeroDivisionError:
            return inf * self._ureg("nm")

    @wavelength.setter
    @default_units("nm")
    def wavelength(self, λ):
        c = self._ureg.c
        self.omega = 2 * π * c / λ

    @property
    def λ(self):
        return self.wavelength

    @λ.setter
    def λ(self, λ):
        self.wavelength = λ

    # ---------
    # Linewidth
    # ---------

    @property
    def linewidth(self):
        return self.__linewidth

    @linewidth.setter
    @default_units("Hz")
    def linewidth(self, linewidth):
        self.__linewidth = linewidth

    # --------------
    # Electric Field
    # --------------

    @property
    def electric_field(self):
        return self.__electric_field

    @electric_field.setter
    @default_units("V/m")
    def electric_field(self, E):
        self.__electric_field = E

    @property
    def E(self):
        return self.electric_field

    @E.setter
    def E(self, E):
        self.electric_field = E

    @property
    def intensity(self):
        c = self._ureg.c
        ε_0 = self._ureg.ε_0
        return self.electric_field ** 2 * (c * ε_0 / 2)

    @intensity.setter
    @default_units("W/m^2")
    def intensity(self, I):
        c = self._ureg.c
        ε_0 = self._ureg.ε_0
        self.electric_field = (2 * I / (c * ε_0)) ** (1 / 2)

    @property
    def I(self):
        return self.intensity

    @I.setter
    def I(self, I):
        self.intensity = I

 
