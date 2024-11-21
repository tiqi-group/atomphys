import csv
import io
import re
import requests
from typing import List
from html.parser import HTMLParser

from atomphys.data_utils.name import parse_atom_name
from atomphys.data_utils.cache import disk_cache
from atomphys.utils.term import vaildate_term

re_monovalent = re.compile(r"^[a-z0-9]*p6\.(?P<n>\d+)[a-z]$")


def remove_annotations(s: str) -> str:
    """
    Removes specific characters used for annotations from a string.

    This function is designed to process strings representing energy levels or
    similar scientific data from the NIST Atomic Spectra Database (ASD). It
    removes any characters used for annotations and uncertainties, along with
    the specific "&dagger;" symbol.

    An alternative approach using regex was considered but discarded in favor of
    the current method due to the significant performance improvement of the latter.

    Parameters:
    s (str): The input string that needs annotation characters removed.

    Returns:
    str: The string with specific annotation characters removed.
    """

    # Strip the following characters from both ends of the string: '(', ')', '[', ']', 'a', 'l', 'u', 'x', 'y', 'z', ' ', '+', '?'.
    # These characters are often used for annotations or uncertainties in NIST ASD data.
    stripped_string = s.strip("()[]aluxyz +?")

    # Replace the "&dagger;" symbol with an empty string, effectively removing it.
    # The "&dagger;" symbol is often used to denote different states or variations in the data.
    cleaned_string = stripped_string.replace("&dagger;", "")

    return cleaned_string


def tokenize_name(name):
    A, element, charge = parse_atom_name(name)
    element = element.lower()
    ionization_state = "i" * (charge + 1)
    token = f"{element} {ionization_state}"
    return token


def load_from_nist(name, refresh_cache=False):
    # TODO: do something with the mass number
    token = tokenize_name(name)
    states_data = parse_states(fetch_states(token, refresh_cache))
    transitions_data = parse_transitions(fetch_transitions(token, refresh_cache))
    return name, states_data, transitions_data


class NISTHTMLParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_title = False
        self.in_error = False
        self.title = ""
        self.error_message = ""

    def handle_starttag(self, tag, attrs):
        if tag == "h2":
            self.in_title = True
        elif tag == "font" and ("color", "red") in attrs:
            self.in_error = True

    def handle_endtag(self, tag):
        if tag == "h2":
            self.in_title = False
        elif tag == "font":
            self.in_error = False

    def handle_data(self, data):
        if self.in_title:
            self.title = data
        elif self.in_error:
            self.error_message = data


def query_nist_database(url: str, params: dict):
    resp = requests.get(url, params=params)
    if resp.text.startswith("<"):
        parser = NISTHTMLParser()
        parser.feed(resp.text)
        raise ValueError(f"{parser.title} : {parser.error_message}")
    else:
        return resp.text


@disk_cache
def fetch_states(atom, refresh_cache=False):
    # option "off" is invalid, to not include a field just do not query for it
    url = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"
    values = {
        "spectrum": atom,
        "units": 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        "format": 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        "multiplet_ordered": 1,  # energy ordred
        "term_out": "on",  # output the term symbol string
        "conf_out": "on",  # output the configutation string
        "level_out": "on",  # output the energy level
        # "unc_out": "on",  # uncertainty on energy
        "j_out": "on",  # output the J level
        "g_out": "on",  # output the g-factor
        "lande_out": "on",  # output experimentally measured g-factor
    }

    resp_text = query_nist_database(url, values)
    data = list(
        csv.DictReader(io.StringIO(resp_text), dialect="excel-tab", restkey="None")
    )
    return data


def _parse_state(state: dict):
    term = state["Term"] + state["J"]
    if vaildate_term(term):
        _parsed_data = {
            "configuration": state["Configuration"],
            "term": term.replace("*", ""),
            "energy": remove_annotations(state["Level (Ry)"]) + " Ry",
            # "g": None if state["g"] == "" else float(state["g"]),
            # "Lande": None if "Lande" not in state or state["Lande"] == "" else float(state["Lande"])
        }
        # n_match = re_monovalent.match(state["Configuration"])
        # if n_match:
        #     _parsed_data['n'] = int(n_match['n'])
    else:
        _parsed_data = None
    return _parsed_data


def parse_states(data: List[dict]):
    states = []
    for state in data:
        parsed = _parse_state(state)
        if parsed is not None:
            states.append(parsed)
    return states


@disk_cache
def fetch_transitions(atom, refresh_cache=False):
    # the NIST url and GET options.
    url = "http://physics.nist.gov/cgi-bin/ASD/lines1.pl"
    values = {
        "spectra": atom,
        "format": 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        "en_unit": 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        "line_out": 2,  # only with {1: transition , 2: level classifications}
        #        "show_av": 5, - was causing some errors
        "allowed_out": 1,
        "forbid_out": 1,
        "enrg_out": "on",
        "term_out": "on",
        "J_out": "on",
        "no_spaces": "on",
    }

    resp_text = query_nist_database(url, values)
    data = list(csv.DictReader(io.StringIO(resp_text), dialect="excel-tab"))
    return data


def _parse_transition(transition: dict):
    A = transition["Aki(s^-1)"] + " s^-1"
    term_i = transition["term_i"] + transition["J_i"]
    term_k = transition["term_k"] + transition["J_k"]
    # data = {
    #     k.lower(): v for k, v in transition.items()
    #     if k in ['Acc', 'Type']
    # }

    if A and term_i and term_k:
        _parsed_data = {
            "A": A,
            "state_i": {
                "energy": remove_annotations(transition["Ei(Ry)"]) + " Ry",
                "term": term_i.replace("*", ""),
            },
            "state_f": {
                "energy": remove_annotations(transition["Ek(Ry)"]) + " Ry",
                "term": term_k.replace("*", ""),
            },
            # **data
        }
    return _parsed_data


def parse_transitions(data: list[dict]):
    return [_parse_transition(tr) for tr in data]
