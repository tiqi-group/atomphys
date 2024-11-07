#!/usr/bin/env python
# -*- coding:utf-8 -*-

from functools import wraps
from typing import Callable
import pint
import numpy as np

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


def spherical_basis_vector(q):
    """
    Function to compute the normalized spherical basis vector c(q) in cartesian corrindates..

    Args:
    q (int): Quantum number, should be -1, 0, or 1.

    Returns:
    np.array: Corresponding basis vector.
    """
    if q == 1:
        return (-1) * np.array([1, -1j, 0]) / np.sqrt(2)
        return (-1) * np.array([1, -1j, 0]) / np.sqrt(2)
    elif q == 0:
        return np.array([0, 0, 1])
    elif q == -1:
        return np.array([1, 1j, 0]) / np.sqrt(2)
    else:
        raise ValueError("q should be -1, 0, or 1.")

        raise ValueError(f"Invalid value {q} for q. It must be one of (-1, 0, 1).")


def spherical_basis_second_rank_tensor(q):
    """
    Function to compute the second rank tensor c(q)_ij in cartesian coordinates.

    Args:
    q (int): Quantum number, should be -1, 0, or 1.

    Returns:
    np.array: Corresponding second rank tensor.
    """

    if q == 2:
        return (1 / np.sqrt(6)) * np.array([[1, -1j, 0], [-1j, -1, 0], [0, 0, 0]])
        return (1 / np.sqrt(6)) * np.array([[1, -1j, 0], [-1j, -1, 0], [0, 0, 0]])
    elif q == 1:
        return (1 / np.sqrt(6)) * np.array([[0, 0, -1], [0, 0, 1j], [-1, 1j, 0]])
        return (1 / np.sqrt(6)) * np.array([[0, 0, -1], [0, 0, 1j], [-1, 1j, 0]])
    elif q == 0:
        return (1 / 3) * np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 2]])
        return (1 / 3) * np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 2]])
    elif q == -1:
        return (1 / np.sqrt(6)) * np.array([[0, 0, 1], [0, 0, 1j], [1, 1j, 0]])
        return (1 / np.sqrt(6)) * np.array([[0, 0, 1], [0, 0, 1j], [1, 1j, 0]])
    elif q == -2:
        return (1 / np.sqrt(6)) * np.array([[1, 1j, 0], [1j, -1, 0], [0, 0, 0]])
        return (1 / np.sqrt(6)) * np.array([[1, 1j, 0], [1j, -1, 0], [0, 0, 0]])
    else:
        raise ValueError(
            f"Invalid value {q} for q. It must be one of (-2, -1, 0, 1, 2)."
        )


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


def inthalf(value):
    if abs(value - round(value)) < 1e-8:  # Close enough to an integer
        return int(round(value))  # Convert to an integer
    elif abs(value * 2 - round(value * 2)) < 1e-8:  # Close enough to a half-integer
        return round(value * 2) / 2  # Convert to a half-integer (like 0.5, 1.5, etc.)
    else:
        raise ValueError(f"Cannot convert {value} to integer or half-integer.")
