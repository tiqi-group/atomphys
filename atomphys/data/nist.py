import csv
import io
import re
import requests
from typing import List
from html.parser import HTMLParser

from atomphys.term import print_term
from atomphys.util import disk_cache

re_monovalent = re.compile(r"^[a-z0-9]*p6\.(?P<n>\d+)[a-z]$")
# re_atom_name = re.compile(r"^(?P<A>\d*)(?P<element>[A-Za-z]*)(?P<charge>(?:\d+\+|\++)+)$")
re_atom_name = re.compile(r"^(?P<A>\d*)(?P<element>[A-Za-z]*)(?P<charge>\d*\+|\++)?$")


def remove_annotations(s: str) -> str:
    """remove annotations from energy strings in NIST ASD"""
    # re_energy = re.compile("-?\\d+\\.\\d*|$")
    # return re_energy.findall(s)[0]

    # this is about 3.5× faster than re.findall, but it's less flexible
    # overall this can make a several hundred ms difference when loading
    return s.strip("()[]aluxyz +?").replace("&dagger;", "")


class NISTHTMLParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_title = False
        self.in_error = False
        self.title = ''
        self.error_message = ''

    def handle_starttag(self, tag, attrs):
        if tag == 'title':
            self.in_title = True
        elif tag == 'font' and ('color', 'red') in attrs:
            self.in_error = True

    def handle_endtag(self, tag):
        if tag == 'title':
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
    if 'text/plain' in resp.headers['content-type']:
        return resp.text
    elif 'text/html' in resp.headers['content-type']:
        parser = NISTHTMLParser()
        parser.feed(resp.text)
        raise ValueError(f"{parser.title} : {parser.error_message}")
    else:
        raise ValueError("Invalid response")


def parse_atom_name(name):
    """
    Parse atom name

    Args:
        input_string (str): The input string containing the information.

    Returns:
        tuple or None: A tuple containing the extracted information or None if no match found.

    Description:
        The function extracts information from the input string based on the following pattern:

        - A: Any number of digits (0 or more).
        - element: Any number of uppercase or lowercase letters (0 or more).
        - charge: Either any number of digits followed by a single "+" or any number of "+" characters.

        The function returns a tuple containing the extracted information in the following order:

        - A (int): The value of "A" as an integer.
        - element (str): The content of "element" as a string.
        - num_charge (int): Either the number in the "charge" group if present, or the number of "+" characters.

        If no match is found, None is returned.
    """
    # TODO move this docstring perhaps somewhere else
    match = re_atom_name.match(name)
    if match:
        A = int(match.group('A')) if match.group('A') else 0
        element = match.group('element')
        charge = match.group('charge')
        if charge:
            num_charge = charge.count('+') if len(set(charge)) == 1 else int(charge[:-1])
        else:
            num_charge = 0
        return A, element, num_charge
    else:
        return None


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
        # "lande_out": "on",  # output experimentally measured g-factor
    }

    resp_text = query_nist_database(url, values)
    data = list(csv.DictReader(io.StringIO(resp_text), dialect="excel-tab", restkey="None"))
    return data


def parse_states(data: List[dict]):
    return [
        {
            **{
                "energy": remove_annotations(state["Level (Ry)"]) + " Ry",
                "term": print_term(state["Term"], include_parity=True, J=state["J"]),
                "configuration": state["Configuration"],
                "g": None if state["g"] == "" else float(state["g"]),
            },
            **(
                {"n": int(re_monovalent.match(state["Configuration"])["n"])}
                if re_monovalent.match(state["Configuration"])
                else {}
            ),
        }
        for state in data
        if print_term(state["Term"], J=state["J"])
    ]


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
    return [
        {
            "state_i": {
                "energy": remove_annotations(transition["Ei(Ry)"]) + " Ry",
                "term": print_term(term=transition["term_i"], J=transition["J_i"]),
            },
            "state_f": {
                "energy": remove_annotations(transition["Ek(Ry)"]) + " Ry",
                "term": print_term(term=transition["term_k"], J=transition["J_k"]),
            },
            "A": transition["Aki(s^-1)"] + "s^-1",
            "type": transition["Type"],
        }
        for transition in data
        if (
            transition["Aki(s^-1)"] and
            print_term(term=transition["term_i"], J=transition["J_i"]) and
            print_term(term=transition["term_k"], J=transition["J_k"])
        )
    ]


def load_from_nist(name, refresh_cache):
    # TODO: do something with the mass number
    A, element, num_charge = parse_atom_name(name)
    element = element.lower()
    ionization_state = 'i' * (num_charge + 1)
    token = f"{element} {ionization_state}"

    states_data = parse_states(fetch_states(token, refresh_cache))
    transitions_data = parse_transitions(fetch_transitions(token, refresh_cache))
    return states_data, transitions_data
