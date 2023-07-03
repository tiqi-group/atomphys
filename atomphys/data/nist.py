import csv
import io
import re
import requests
from typing import List
from html.parser import HTMLParser

from .name import parse_atom_name
from atomphys.term import print_term
from atomphys.util import disk_cache

re_monovalent = re.compile(r"^[a-z0-9]*p6\.(?P<n>\d+)[a-z]$")


def remove_annotations(s: str) -> str:
    """remove annotations from energy strings in NIST ASD"""
    # re_energy = re.compile("-?\\d+\\.\\d*|$")
    # return re_energy.findall(s)[0]

    # this is about 3.5× faster than re.findall, but it's less flexible
    # overall this can make a several hundred ms difference when loading
    return s.strip("()[]aluxyz +?").replace("&dagger;", "")


def tokenize_name(name):
    A, element, charge = parse_atom_name(name)
    element = element.lower()
    ionization_state = 'i' * (charge + 1)
    token = f"{element} {ionization_state}"
    return token


def load_from_nist(name, refresh_cache=False):
    # TODO: do something with the mass number
    token = tokenize_name(name)
    states_data = parse_states(fetch_states(token, refresh_cache))
    transitions_data = parse_transitions(fetch_transitions(token, refresh_cache))
    return states_data, transitions_data


class NISTHTMLParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_title = False
        self.in_error = False
        self.title = ''
        self.error_message = ''

    def handle_starttag(self, tag, attrs):
        if tag == 'h2':
            self.in_title = True
        elif tag == 'font' and ('color', 'red') in attrs:
            self.in_error = True

    def handle_endtag(self, tag):
        if tag == 'h2':
            self.in_title = False
        elif tag == 'font':
            self.in_error = False

    def handle_data(self, data):
        if self.in_title:
            self.title = data
        elif self.in_error:
            self.error_message = data


def query_nist_database(url: str, params: dict):
    resp = requests.get(url, params=params)
    if resp.text.startswith('<'):
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
    data = list(csv.DictReader(io.StringIO(resp_text), dialect="excel-tab", restkey="None"))
    return data


def parse_states1(state: dict):
    term = state["Term"] + state["J"]
    if term:
        _parsed_data = {
            "energy": remove_annotations(state["Level (Ry)"]) + " Ry",
            "term": term,
            "configuration": state["Configuration"],
            "g": None if state["g"] == "" else float(state["g"]),
            "Lande": None if "Lande" not in state or state["Lande"] == "" else float(state["Lande"])
        }
        n_match = re_monovalent.match(state["Configuration"])
        if n_match:
            _parsed_data['n'] = int(n_match['n'])
    else:
        _parsed_data = {}
    return _parsed_data


def parse_states(data: List[dict]):
    return [parse_states1(state) for state in data]


@disk_cache
def fetch_transitions(atom, refresh_cache=False):
    # the NIST url and GET options.
    url = "http://physics.nist.gov/cgi-bin/ASD/lines1.pl"
    values = {
        "spectra": atom,
        "format": 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        "en_unit": 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        "line_out": 2,  # only with {1: transition , 2: level classifications}
        "show_av": 5,
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


def parse_transitions(data: List[dict]):
    parsed = []
    for transition in data:
        A = transition["Aki(s^-1)"] + "s^-1"
        term_i = print_term(term=transition["term_i"], include_parity=True, J=transition["J_i"])
        term_k = print_term(term=transition["term_k"], include_parity=True, J=transition["J_k"])

        print(transition['term_k'], transition['J_k'], term_k)
        if A and term_i and term_k:
            _parsed_data = {
                "A": A, "type": transition["Type"],
                "state_i": {
                    "energy": remove_annotations(transition["Ei(Ry)"]) + " Ry",
                    "term": term_i,
                },
                "state_f": {
                    "energy": remove_annotations(transition["Ek(Ry)"]) + " Ry",
                    "term": term_k,
                },
            }
            parsed.append(_parsed_data)
    return parsed
