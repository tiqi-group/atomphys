#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>
#
# Reference for this calculation is James 1998: Quantum dynamics of cold trappedions, with application to quantum computation
# https://arxiv.org/abs/quant-ph/9702053
#
# Whilst reading it, be careful about the definitions of j' and j. It is not consistent throughout the paper.
# I took j' to be the lower level states, as it agrees with other sources and makes sense intuitively
# -> it increases spontaneous decay rate the more possibilities of decays to the lower states there are.

import pint
from sympy.physics.wigner import wigner_3j as w3j
from atomphys.utils.utils import (
    spherical_basis_second_rank_tensor,
    spherical_basis_vector,
    inthalf,
)


def reduced_dipole_matrix_element(
    A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry
):
    """
    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_f: Angular momentum of the upper state
        _ureg: Unit registry

    Returns:
        Reduced dipole matrix element [a0]
    """

    return (_ureg("3/(4*c*alpha)") * (2 * J_f + 1) / k**3 * A) ** (1 / 2)


def reduced_electric_dipole_matrix_element(
    A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry
):
    """
    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_f: Angular momentum of the upper state
        _ureg: Unit registry

    Returns:
        Reduced electric dipole matrix element [e a0]
    """

    return _ureg("e") * reduced_dipole_matrix_element(A, k, J_f, _ureg)


def reduced_quadrupole_matrix_element(
    A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry
):
    """
    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_f: Angular momentum of the upper state
        _ureg: Unit registry

    Returns:
        Reduced dipole matrix element [a0**2]
    """

    return (_ureg("15/(c*alpha)") * (2 * J_f + 1) / k**5 * A) ** (1 / 2)


def reduced_electric_quadrupole_matrix_element(
    A: pint.Quantity, k: pint.Quantity, J_f: float, _ureg: pint.UnitRegistry
):
    """
    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_f: Angular momentum of the upper state
        _ureg: Unit registry

    Returns:
        Reduced dipole matrix element [e a0**2]
    """

    return _ureg("e") * reduced_quadrupole_matrix_element(A, k, J_f, _ureg)


def dipole_matrix_element_basis(
    A: pint.Quantity,
    k: pint.Quantity,
    J_i: float,
    J_f: float,
    mJ_i: float,
    mJ_f: float,
    q: float,
    _ureg: pint.UnitRegistry,
):
    """
    Dipole Matrix Element is reduced dipole matrix element decorated with clebsh gordan coefficients (wigner 3j sybol),
    Dipole Matrix Element is reduced dipole matrix element decorated with clebsh gordan coefficients (wigner 3j sybol),
    hence expressing the coupling between different magnetic sublevels.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        q: Spherical basis vector
        _ureg: Unit registry

    Returns:
        Dipole matrix element [a0]
    """

    dme = 0 * _ureg("a0")
    d = reduced_dipole_matrix_element(A, k, J_f, _ureg)

    w3j_coeff = float(
        w3j(inthalf(J_f), 1, inthalf(J_i), inthalf(-mJ_f), q, inthalf(mJ_i))
    )
    dme += w3j_coeff * d
    
    return dme.to("a0").magnitude * _ureg("a0")

def dipole_matrix_element(
    A: pint.Quantity,
    k: pint.Quantity,
    J_i: float,
    J_f: float,
    mJ_i: float,
    mJ_f: float,
    _ureg: pint.UnitRegistry,
):
    """
    Dipole Matrix Element is reduced dipole matrix element decorated with clebsh gordan coefficients (wigner 3j sybol),
    Dipole Matrix Element is reduced dipole matrix element decorated with clebsh gordan coefficients (wigner 3j sybol),
    hence expressing the coupling between different magnetic sublevels.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry

    Returns:
        Dipole matrix element [a0]
    """

    dme = 0 * _ureg("a0")
    d = reduced_dipole_matrix_element(A, k, J_f, _ureg)

    for q in range(-1, 2):
        w3j_coeff = float(
            w3j(inthalf(J_f), 1, inthalf(J_i), inthalf(-mJ_f), q, inthalf(mJ_i))
        )
        dme += w3j_coeff * spherical_basis_vector(q) * d
    return dme.to("a0").magnitude * _ureg("a0")


def quadrupole_matrix_element(
    A: pint.Quantity,
    k: pint.Quantity,
    J_i: float,
    J_f: float,
    mJ_i: float,
    mJ_f: float,
    _ureg: pint.UnitRegistry,
):
    """
    Quadrupole Matrix Element is reduced quadrupole matrix element decorated with clebsh gordan coefficients (wigner 3j sybol),
    hence expressing the coupling between different magnetic sublevels.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry

    Returns:
        Quadrupole matrix element [a0**2]
    """

    qme = 0 * _ureg("a0**2")
    red_q = reduced_quadrupole_matrix_element(A, k, J_f, _ureg)
    for q in range(-2, 3):
        w3j_coeff = float(
            w3j(inthalf(J_f), 2, inthalf(J_i), inthalf(-mJ_f), q, inthalf(mJ_i))
        )
        sbt = spherical_basis_second_rank_tensor(q)
        qme += w3j_coeff * sbt * red_q
    return qme.to("a0**2").magnitude * _ureg("a0**2")


def electric_dipole_matrix_element(
    A: pint.Quantity,
    k: pint.Quantity,
    J_i: float,
    J_f: float,
    mJ_i: float,
    mJ_f: float,
    _ureg: pint.UnitRegistry,
):
    """
    See dipole matrix element for details how is it different from reduced matrix elements.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry

    Returns:
        Electric dipole matrix element [e a0]
    """

    return _ureg("e") * dipole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)


def electric_quadrupole_matrix_element(
    A: pint.Quantity,
    k: pint.Quantity,
    J_i: float,
    J_f: float,
    mJ_i: float,
    mJ_f: float,
    _ureg: pint.UnitRegistry,
):
    """
    See quadrupole matrix element for details how is it different from reduced matrix elements.

    Args:
        A: Einstein coefficient [1/s]
        k: Wavenumber of a transition [1/m]
        J_i: Angular momentum of the lower state
        J_f: Angular momentum of the upper state
        mJ_i: Magnetic quantum number of the lower state
        mJ_f: Magnetic quantum number of the upper state
        _ureg: Unit registry


    Returns:
        Electric quadrupole matrix element [e a0]
    """

    return (
        _ureg("e") * quadrupole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)
    ).to("e*a0**2")