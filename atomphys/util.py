#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>


from functools import wraps
from typing import Callable
import pint

_ureg = pint.get_application_registry()


def set_default_units(quantity, unit, _ureg=_ureg):
    if isinstance(quantity, str):
        quantity = _ureg(quantity)
    if not isinstance(quantity, _ureg.Quantity) or quantity.dimensionless:
        quantity = _ureg.Quantity(quantity, unit)
    if not quantity.check(unit):
        raise ValueError(f"must have units equivalent to {unit}")
    return quantity


def default_units(unit: str):
    def decorator(setter):
        @wraps(setter)
        def wrapper(self, quantity, *args, **kwargs):
            quantity = set_default_units(quantity, unit, self._ureg)
            return setter(self, quantity, *args, **kwargs)

        return wrapper

    return decorator


def fsolve(func: Callable, x0, x1=None, tol: float = 1.49012e-08, maxfev: int = 100):
    """
    Find the roots of a function using the secant method.

    Return the roots of the equation `func(x) = 0` given a
    starting estimate `x0`. A second starting estimate `x1`
    for the next iteration of the secant method can be supplied.

    Arguments:
        func: a function `f(x)` that takes a single argument `x`
        x0: The starting estimate for the roots of `func(x) = 0``
        x1: A second starting estimate for the next iteration of
        the secant method. Defaults to `1.0001 * x0`.
        tol: The calculation will terminate if the relative
            error between two consecutive iterates is at most `tol`
        maxfev: The maximum number of calls to the function.

    Returns:
        The root of the function `f(x)`.
    """
    if x1 is None:
        x1 = x0 * 1.001 if x0 else 0.001
    fx0, fx1 = func(x0), func(x1)
    i = 2
    while (abs(fx0) > 0) and (abs((fx1 - fx0) / fx0) > tol) and (i < maxfev + 1):
        x0, x1 = x1, x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        fx0, fx1 = fx1, func(x1)
        i += 1
    return x1


def make_alias(attr_name: str, get_unit: str = None, multiplication_factor: float = 1):
    @property
    def prop(self):
        if get_unit is None:
            return multiplication_factor * getattr(self, attr_name)
        else:
            return multiplication_factor * getattr(self, attr_name).to(get_unit)
    return prop



def make_alias_with_setter(attr_name: str, get_unit: str = None):
    @property
    def prop(self):
        if get_unit is None:
            return getattr(self, attr_name)
        else:
            return getattr(self, attr_name).to(get_unit)

    @prop.setter
    def prop(self, value):
        if get_unit is None:
            setattr(self, attr_name, value)
        else:
            setattr(self, attr_name, value.to(get_unit))

    return prop
