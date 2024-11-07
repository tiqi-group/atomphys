#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

from dataclasses import dataclass, asdict
from fractions import Fraction
from atomphys.utils.term import print_term, parse_term


@dataclass(frozen=True, kw_only=True)
class QuantumNumbers:
    J: float | None = None
    L: int | None = None
    S: float | None = None
    K: float | None = None
    S2: float | None = None
    J1: float | None = None
    J2: float | None = None
    I: float | None = None
    F: float | None = None
    parity: int | None = None
    ionization_limit: bool | None = None

    def __getitem__(self, key):
        return self.__dict__[key]

    @staticmethod
    def from_term(term: str, I=None, F=None):
        data = parse_term(term)
        return QuantumNumbers(**data, I=I, F=F)

    @property
    def term(self) -> str:
        return print_term(**asdict(self))

    def as_dict(self):
        return asdict(self)

    def __repr__(self) -> str:
        out = [f"{key}={Fraction(value)}" for key, value in asdict(self).items() if value is not None]
        return f"QuantumNumbers({', '.join(out)})"
