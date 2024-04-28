import json
from atomphys.data_utils.name import parse_atom_name


def load_from_json(filename):
    with open(filename, 'r') as fp:
        data = json.load(fp)
    name = data['name']
    parse_atom_name(name)  # only to validate here
    states_data = data['states']
    transitions_data = data['transitions']
    return name, states_data, transitions_data
