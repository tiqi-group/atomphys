import enum
import re
from fractions import Fraction
from dataclasses import dataclass, asdict


class Coupling(enum.Enum):
    LS = "LS"  # Russell-Saunders coupling
    jj = "jj"
    LK = "LK"  # pair coupling


LS_term = re.compile(r"^(?P<S>\d+)(?P<L>[A-Z])(?P<parity>\*)?(?P<J>\d+(/2)?)$")
r"""
Regular expression for LS terms.

Matches a pattern that consists of:
- ^: Start of the string.
- (?P<S>\d+): Named capturing group 'S' matching one or more digits.
- (?P<L>[A-Z]): Named capturing group 'L' matching a single uppercase letter.
- (?P<parity>\*)?: Named capturing group 'parity' matching an asterisk character '*' zero or one time (optional).
- (?P<J>\d+(/2)?): Named capturing group 'J' matching a sequence of digits, optionally followed by '/2'.
- $: End of the string.
"""

J1J2_term = re.compile(r"^\((?P<J1>\d+(/2)?),(?P<J2>\d+(/2)?)\)(?P<parity>\*)?(?P<J>\d+(/2)?)$")
r"""
Regular expression for J1J2 terms.

Matches a pattern that consists of:
- ^: Start of the string.
- \(: Literal parentheses '(' character.
- (?P<J1>\d+(/2)?): Named capturing group 'J1' matching a sequence of digits, optionally followed by '/2'.
- ,: Literal comma ',' character.
- (?P<J2>\d+(/2)?): Named capturing group 'J2' matching a sequence of digits, optionally followed by '/2'.
- \): Literal parentheses ')' character.
- (?P<parity>\*)?: Named capturing group 'parity' matching an asterisk character '*' zero or one time (optional).
- (?P<J>\d+(/2)?): Named capturing group 'J' matching a sequence of digits, optionally followed by '/2'.
- $: End of the string.
"""

LK_term = re.compile(r"^(?P<S2>\d+)\[(?P<K>\d+(/2)?)\](?P<parity>\*)?(?P<J>\d+(/2)?)$")
r"""
Regular expression for LK terms.

Matches a pattern that consists of:
- ^: Start of the string.
- (?P<S2>\d+): Named capturing group 'S2' matching one or more digits.
- \[: Literal square bracket '[' character.
- (?P<K>\d+(/2)?): Named capturing group 'K' matching a sequence of digits, optionally followed by '/2'.
- \]: Literal square bracket ']' character.
- (?P<parity>\*)?: Named capturing group 'parity' matching an asterisk character '*' zero or one time (optional).
- (?P<J>\d+(/2)?): Named capturing group 'J' matching a sequence of digits, optionally followed by '/2'.
- $: End of the string.
"""


L = {
    "S": 0, "P": 1, "D": 2, "F": 3,
    "G": 4, "H": 5, "I": 6, "K": 7,
    "L": 8, "M": 9, "N": 10, "O": 11,
    "Q": 12, "R": 13, "T": 14, "U": 15,
    "V": 16, "W": 17, "X": 18, "Y": 19
}

L_inv = {value: key for key, value in L.items()}


@dataclass(frozen=True, kw_only=True)
class QuantumNumbers:
    J: float | None = None
    S: float | None = None
    L: int | None = None
    J1: float | None = None
    J2: float | None = None
    S2: float | None = None
    K: float | None = None
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


def parse_term(term: str) -> dict:
    """
    parse term symbol string in NIST ASD
    """
    if "Limit" in term:
        return {"ionization_limit": True}

    match = LS_term.search(term)
    if match is None:
        match = J1J2_term.search(term)
    if match is None:
        match = LK_term.search(term)
    if match is None:
        raise ValueError(f"Invalid term {term}")

    def convert(key, value):
        if key in ["S", "S2"]:
            return float(Fraction((int(value) - 1) / 2))
        if key in ["J1", "J2", "J", "K", "I", "F"]:
            return float(Fraction(value))
        if key == "L":
            return L[value]
        if key == "parity":
            return -1 if value == "*" else 1

    term = {
        key: convert(key, value) for key, value in match.groupdict().items()
    }

    return term


def print_term(J=None, S=None, L=None,
               J1=None, J2=None, S2=None, K=None,
               I=None, F=None,
               parity=None, ionization_limit=None) -> str:

    if ionization_limit is not None and ionization_limit:
        return "Ionization Limit"

    P = "*" if parity == -1 else ""

    J = f"{Fraction(J)}"

    # Russell-Saunders coupling
    if S is not None and L is not None:
        return f"{2 * S + 1:g}{L_inv[L]}{P}{J}"

    # J1J2 coupling
    if J1 is not None and J2 is not None:
        return f"({Fraction(J1)},{Fraction(J2)}){P}{J}"

    # LK coupling
    if S2 is not None and K is not None:
        return f"{2 * S2 + 1:g}[{Fraction(K)}]{P}{J}"

    raise ValueError("Invalid arguments")
