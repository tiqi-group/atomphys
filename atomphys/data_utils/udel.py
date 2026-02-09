"""
Data loader for the UDel ATOM portal (https://www1.udel.edu/atom/).

The ATOM portal at the University of Delaware provides high-precision
atomic data calculated using linearized coupled-cluster methods,
including energies, transition rates, matrix elements, polarizabilities,
and hyperfine constants.

Data is accessed via a GraphQL API at https://atom.ece.udel.edu/graphql.
"""

import re
import requests
from typing import List

from atomphys.data_utils.name import parse_atom_name
from atomphys.data_utils.cache import disk_cache

GRAPHQL_ENDPOINT = "https://atom.ece.udel.edu/graphql"

# 1 Rydberg = 109737.31568160 cm^-1 (CODATA 2018)
_RY_TO_CM = 109737.31568160

# Regex for states with explicit term info (divalent, HCI)
# e.g.: "3d4s <sup>3</sup>D<sub>1</sub>" or "4s<sup>2</sup> <sup>1</sup>S<sub>0</sub>"
_re_explicit_term = re.compile(
    r"<sup>(\d+)</sup>([A-Z])<sub>([^<]+)</sub>$"
)

# Regex for states without explicit term info (monovalent)
# e.g.: "5s<sub>1/2</sub>" or "4f<sub>7/2</sub>"
_re_monovalent_state = re.compile(r"^(\d+[a-z])<sub>([^<]+)</sub>$")

_orbital_to_L = {
    "s": "S",
    "p": "P",
    "d": "D",
    "f": "F",
    "g": "G",
    "h": "H",
    "i": "I",
}


def _name_to_udel_title(name: str) -> str:
    """Convert atomphys atom name to UDel portal title format.

    Examples:
        "Rb" -> "Rb1", "Ca+" -> "Ca2", "Cs6+" -> "Cs7"
    """
    _A, element, charge = parse_atom_name(name)
    return f"{element}{charge + 1}"


def _graphql_query(query: str, variables: dict = None) -> dict:
    """Execute a GraphQL query against the UDel ATOM portal."""
    payload = {"query": query}
    if variables:
        payload["variables"] = variables
    resp = requests.post(GRAPHQL_ENDPOINT, json=payload, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    if "errors" in data:
        raise ValueError(f"UDel ATOM portal GraphQL error: {data['errors']}")
    return data["data"]


def _parse_state_html(html_state: str) -> tuple:
    """Parse UDel HTML state name into (configuration, term).

    Handles two formats:
    - Explicit term: "3d4s <sup>3</sup>D<sub>1</sub>" -> ("3d4s", "3D1")
    - Monovalent: "5s<sub>1/2</sub>" -> ("5s", "2S1/2")

    Returns:
        tuple: (configuration, term) where term includes J.
    """
    # Try explicit term format (divalent, HCI)
    m = _re_explicit_term.search(html_state)
    if m:
        config_html = html_state[: m.start()].strip()
        # Remove HTML tags from configuration part
        config = re.sub(r"<[^>]+>", "", config_html)
        multiplicity = m.group(1)
        L_letter = m.group(2)
        J = m.group(3)
        term = f"{multiplicity}{L_letter}{J}"
        return config, term

    # Try monovalent format
    m = _re_monovalent_state.match(html_state.strip())
    if m:
        config = m.group(1)
        J = m.group(2)
        l_char = config[-1]
        L_letter = _orbital_to_L.get(l_char, "?")
        # Monovalent atoms are doublets (S=1/2, so 2S+1=2)
        term = f"2{L_letter}{J}"
        return config, term

    raise ValueError(f"Cannot parse UDel state HTML: {html_state}")


_ENERGIES_QUERY = """
query GetElementEnergies($title: String) {
    element(title: $title) {
        id
        title
        titleDisplay
        NISTASDTitle
        energies {
            id
            state
            energy
            energyUncertainty
            isFromTheory
        }
    }
}
"""

_TRANSITIONS_QUERY = """
query GetElementTransitions($title: String) {
    element(title: $title) {
        transitionRates {
            stateOne
            stateOneConfiguration
            stateOneTerm
            stateOneJ
            stateTwo
            stateTwoConfiguration
            stateTwoTerm
            stateTwoJ
            transitionRate
            transitionRateUncertainty
        }
    }
}
"""


@disk_cache
def fetch_energies(title: str, refresh_cache=False) -> list:
    """Fetch energy levels from the UDel ATOM portal.

    Args:
        title: UDel element title (e.g. "Rb1", "Ca2")
        refresh_cache: Whether to refresh cached data

    Returns:
        List of energy level dictionaries from the GraphQL API.
    """
    data = _graphql_query(_ENERGIES_QUERY, {"title": title})
    element = data.get("element")
    if element is None:
        raise ValueError(
            f"Element '{title}' not found on UDel ATOM portal. "
            f"Available elements can be viewed at https://www1.udel.edu/atom/"
        )
    return element.get("energies", [])


@disk_cache
def fetch_transitions(title: str, refresh_cache=False) -> list:
    """Fetch transition rates from the UDel ATOM portal.

    Args:
        title: UDel element title (e.g. "Rb1", "Ca2")
        refresh_cache: Whether to refresh cached data

    Returns:
        List of transition rate dictionaries from the GraphQL API.
    """
    data = _graphql_query(_TRANSITIONS_QUERY, {"title": title})
    element = data.get("element")
    if element is None:
        raise ValueError(
            f"Element '{title}' not found on UDel ATOM portal. "
            f"Available elements can be viewed at https://www1.udel.edu/atom/"
        )
    return element.get("transitionRates", [])


def _build_state_info_lookup(transitions_data: list) -> dict:
    """Build a lookup from state HTML display name to (configuration, term).

    Uses transition data which has explicit configuration, term, and J fields.
    """
    lookup = {}
    for tr in transitions_data:
        state_one_html = tr["stateOne"]
        config_one = tr["stateOneConfiguration"]
        term_one = tr["stateOneTerm"] + tr["stateOneJ"]
        lookup[state_one_html] = (config_one, term_one)

        state_two_html = tr["stateTwo"]
        config_two = tr["stateTwoConfiguration"]
        term_two = tr["stateTwoTerm"] + tr["stateTwoJ"]
        lookup[state_two_html] = (config_two, term_two)

    return lookup


def parse_states(
    energies_data: list, transitions_data: list
) -> List[dict]:
    """Parse UDel energy and transition data into atomphys state format.

    Args:
        energies_data: Energy data from the UDel portal GraphQL API.
        transitions_data: Transition data from the UDel portal GraphQL API.

    Returns:
        List of state dictionaries compatible with load_from_database().
    """
    state_info = _build_state_info_lookup(transitions_data)

    states = []
    for e in energies_data:
        state_html = e["state"]
        energy_cm = e["energy"]

        # Get configuration and term from transitions lookup (preferred)
        # or fall back to parsing the HTML state name
        if state_html in state_info:
            config, term = state_info[state_html]
        else:
            try:
                config, term = _parse_state_html(state_html)
            except ValueError:
                continue  # Skip states we cannot parse

        # Convert energy from cm^-1 to Rydbergs
        energy_ry = energy_cm / _RY_TO_CM

        states.append(
            {
                "configuration": config,
                "term": term,
                "energy": f"{energy_ry} Ry",
            }
        )

    return states


def parse_transitions(
    transitions_data: list, energies_data: list
) -> List[dict]:
    """Parse UDel transition data into atomphys transition format.

    Args:
        transitions_data: Transition data from the UDel portal GraphQL API.
        energies_data: Energy data from the UDel portal GraphQL API.

    Returns:
        List of transition dictionaries compatible with load_from_database().
    """
    # Build energy lookup: state HTML display name -> energy in cm^-1
    energy_lookup = {e["state"]: e["energy"] for e in energies_data}

    parsed = []
    for tr in transitions_data:
        A = tr["transitionRate"]
        if A is None or A == 0:
            continue

        term_one = tr["stateOneTerm"] + tr["stateOneJ"]
        term_two = tr["stateTwoTerm"] + tr["stateTwoJ"]

        energy_one = energy_lookup.get(tr["stateOne"])
        energy_two = energy_lookup.get(tr["stateTwo"])

        if energy_one is None or energy_two is None:
            continue

        # Convert energies from cm^-1 to Rydbergs
        energy_one_ry = energy_one / _RY_TO_CM
        energy_two_ry = energy_two / _RY_TO_CM

        # stateOne is the upper (decaying) state, stateTwo is the lower state
        # In atomphys convention: state_i = lower, state_f = upper
        parsed.append(
            {
                "A": f"{A} s^-1",
                "state_i": {
                    "energy": f"{energy_two_ry} Ry",
                    "term": term_two,
                },
                "state_f": {
                    "energy": f"{energy_one_ry} Ry",
                    "term": term_one,
                },
            }
        )

    return parsed


def load_from_udel(name: str, refresh_cache: bool = False):
    """Load atomic data from the UDel ATOM portal.

    The UDel ATOM portal (https://www1.udel.edu/atom/) provides high-precision
    atomic data calculated using linearized coupled-cluster methods.

    Available elements include monovalent atoms (Li, Na, K, Rb, Cs, Fr),
    singly-charged ions (Be+, Mg+, Ca+, Sr+, Ba+, Ra+), divalent atoms
    (Mg, Ca, Sr), and several highly charged ions.

    Args:
        name: Atom name in atomphys format (e.g. "Rb", "Ca+", "Sr", "Ba+")
        refresh_cache: Whether to refresh cached data from the portal

    Returns:
        Tuple of (name, states_data, transitions_data) compatible with
        load_from_database().
    """
    title = _name_to_udel_title(name)

    energies_data = fetch_energies(title, refresh_cache)
    transitions_data = fetch_transitions(title, refresh_cache)

    states = parse_states(energies_data, transitions_data)
    transitions = parse_transitions(transitions_data, energies_data)

    return name, states, transitions
