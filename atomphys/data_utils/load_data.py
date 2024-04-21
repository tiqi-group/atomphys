from ..atom import Atom
from ..utils.utils import set_default_units
from . import nist, json
from ..state import State
from ..transition import Transition


def load_from_database(
    name: str, states_data: list[dict], transitions_data: list[dict], energy_cutoff="inf", remove_isolated=False
) -> Atom:
    """
    Returns an atom object given the database input.

    This method returns an atom object given the database input. The database input is a list of dictionaries
    containing the data for each state and transition. The data is loaded into the atom object and the states and
    transitions are added to the graph. The states and transitions are filtered by the energy cutoff and isolated
    states are removed if the remove_isolated flag is set to True.

    Args:
        name (str): Name of the atom
        states_data (list[dict]): List of dictionaries containing the state data
        transitions_data (list[dict]): List of dictionaries containing the transition data
        energy_cutoff (str, optional): Energy cutoff for states. Defaults to "inf".
        remove_isolated (bool, optional): Remove isolated states. Defaults to False.

    Returns:
        Atom: Atom object
    """
    atom = Atom(name)
    print(f"Loading atom {name}")

    # add states
    energy_cutoff = set_default_units(energy_cutoff, "Ry", atom._ureg)
    states = [State(**d) for d in states_data if atom._ureg(d["energy"]) < energy_cutoff]
    atom.add_states(states)
    print(f"Added {len(states)} states")

    # add transitions
    transitions = []
    unmatched = 0
    for d in transitions_data:
        si = atom._match_term_and_energy(d["state_i"]["term"], d["state_i"]["energy"])
        sf = atom._match_term_and_energy(d["state_f"]["term"], d["state_f"]["energy"])
        if si is not None and sf is not None and sf.energy > si.energy:
            tr = Transition(si, sf, d["A"])
            transitions.append(tr)
        else:
            unmatched += 1
    if unmatched:
        print(f"Dropping {unmatched} unmatched transitions")
    atom.add_transitions(transitions)
    print(f"Added {len(transitions)} transitions")

    if remove_isolated:
        atom.remove_isolated(copy=False)

    return atom


def from_nist(name: str, energy_cutoff="inf", remove_isolated=False, refresh_cache=False) -> Atom:
    """
    Returns an atom from the NIST database

    This function is a wrapper for the `load_from_database` function that loads the data from the NIST database.


    Args:
        name (str): Name of the atom to load
        energy_cutoff (str, optional): Energy cutoff for states. Defaults to "inf".
        remove_isolated (bool, optional): Remove isolated states. Defaults to False.
        refresh_cache (bool, optional): Refresh the cache. Defaults to False.

    Returns:
        Atom: Atom object


    """
    name, states_data, transitions_data = nist.load_from_nist(name, refresh_cache)
    return load_from_database(name, states_data, transitions_data, energy_cutoff, remove_isolated)


def from_json(filename: str, energy_cutoff="inf", remove_isolated=False) -> Atom:
    """
    Returns an atom from data in JSON format

    Args:
        filename (str): Path of the data file. The filename without extension (stem) must be a valid atom name.
        energy_cutoff (str, optional): Energy cutoff for states. Defaults to "inf".
        remove_isolated (bool, optional): Remove isolated states. Defaults to False.
        refresh_cache (bool, optional): Refresh the cache. Defaults to False.

    Returns:
        Atom: Atom object


    """
    name, states_data, transitions_data = json.load_from_json(filename)
    return load_from_database(name, states_data, transitions_data, energy_cutoff, remove_isolated)
