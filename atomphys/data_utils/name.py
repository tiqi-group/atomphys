#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import re

re_atom_name = re.compile(r"^(?P<A>\d*)(?P<element>[A-Za-z]+)(?P<charge>\d*\+|\++)?$")


def parse_atom_name(name: str) -> tuple[int, str, int]:
    """
    Parse atom name

    Args:
        name (str): Atom name in the format (A)element(charge) with:

        - A (optional): the atomic mass number, composed by any number of digits (0 or more).
        - element: The element symbol, composed by any number of uppercase or lowercase letters (1 or more).
        - charge (optional): the ionization charge,
            either an integer number followed by a single "+" or any number of "+" characters.

        Examples: Yb, 23Na, Ca+, Sr++, 133Ba+, Ar13+

    Returns:
        a tuple containing the extracted information in the following order:

        - A (int): The value of "A" as an integer.
        - element (str): The content of "element" as a string.
        - charge (int): Either the number in the "charge" group if present, or the number of "+" characters.

        If no match is found, raise ValueError
    """
    # TODO move this docstring perhaps somewhere else
    match = re_atom_name.match(name)
    if match:
        A = int(match.group("A")) if match.group("A") else 0
        element = match.group("element")
        charge = match.group("charge")
        if charge:
            charge = charge.count("+") if set(charge) == {"+"} else int(charge[:-1])
        else:
            charge = 0
        return A, element, charge
    else:
        raise ValueError(f"Invalid atom name {name}")
