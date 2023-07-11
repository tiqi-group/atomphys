import networkx as nx
import pint
from copy import deepcopy
from typing import Union, List, Dict

from .state import State
from .transition import Transition
from .util import default_units, set_default_units
from .data import nist


class Atom:
    def __init__(self, name: str, _ureg=None) -> None:
        self.name = name
        self._graph = nx.DiGraph()
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg

    def __repr__(self) -> str:
        return f"Atom({self.name} {len(self.states)} states {len(self.transitions)} transitions)"

    def copy(self):
        return deepcopy(self)


    """Graph modification functions"""

    def add_state(self, s: State):
        """
            Add state to the graph associated with the Atom

            Args:
                s (State): State object to add to the Atom
        """
        if s in self.states:
            raise ValueError(f"State {s} already in atom")
        
        s.atom = self
        s._ureg = self._ureg
        self._graph.add_node(s)

    def remove_state(self, s: State):
        """
            Remove state and all associated transitions from the graph associated with the Atom

            Args:
                s (State): State object to remove from the Atom
        """
        if s not in self.states:
            raise ValueError(f"State {s} not in atom")
        
        self._graph.remove_node(s)

    def add_states(self, states: list[State]):
        """
            Add states to the graph associated with the Atom

            Args:
                states (list[State]): List of State objects to add to the Atom
        """

        for s in states:
            self.add_state(s)
    
    def remove_states(self, states: list[State]):
        """
        Removes states and all associated transitions from the graph associated with the Atom

        Args:
            states (list[State]): List of State objects to remove from the Atom
                
        """
        for s in states:
            self.remove_state(s)
    

    def add_transition(self, tr: Transition):
        """
        Adds transition to the graph associated with the Atom

        Args:
            tr (Transition): Transition object to add to the Atom
                
        """

        if tr.state_i not in self._graph:
            raise ValueError(f"Transition initial state {tr.state_i} not in atom")
        if tr.state_f not in self._graph:
            raise ValueError(f"Transition final state {tr.state_f} not in atom")
        if tr in self.transitions:
            raise ValueError(f"Transition {tr} already in atom")
        
        tr._ureg = self._ureg
        self._graph.add_edge(tr.state_i, tr.state_f, transition=tr)

    def add_transitions(self, transitions: list[Transition]):
        """
        Adds transitions from a list to the graph associated with the Atom

        Args:
            transitions (list[Transition]): List of Transition objects to add to the Atom
                
        """

        for tr in transitions:
            self.add_transition(tr)

    def remove_transition(self, tr: Transition):
        """
        Removes transition from the graph associated with the Atom

        Args:
            tr (Transition): Transition object to remove from the Atom
        """
        if tr not in self.transitions:
            raise ValueError(f"Transition {tr} not in atom")

        self._graph.remove_edge(tr.state_i, tr.state_f)
    
    def remove_transitions(self, transitions: list[Transition]):
        """
        Removes transitions from a list from the graph associated with the Atom

        Args:
            transitions (list[Transition]): List of Transition objects to remove from the Atom
                
        """
        for tr in transitions:
            self.remove_transition(tr)

    def remove_isolated(self, copy=True):
        """
        Removes isolated states(states without transitions) from the Atom
        
        Args:
            copy (bool, optional): Return a copy of the Atom with isolated states removed. Defaults to True.
        
        Returns:
            Atom: Atom with isolated states removed
        """
        
        
        atom = self.copy() if copy else self
        g = atom._graph
        isolated = list(nx.isolates(g))
        g.remove_nodes_from(isolated)
        print(f"Removed {len(isolated)} states without transitions")
        return atom
    
    def remove_states_above_energy(self, energy: pint.Quantity, copy=True):
        """
        Removes states with energy above a given value from the Atom
        
        Args:
            energy (pint.Quantity): Energy above which states will be removed
            copy (bool, optional): Return a copy of the Atom with states removed. Defaults to True.
        
        Returns:
            Atom: Atom with states above energy removed
        """
        atom = self.copy() if copy else self
        states = atom.states
        for s in states:
            if s.energy > energy:
                atom.remove_state(s)
        print(f"Removed {len(states) - len(atom.states)} states above {energy}")
        return atom
    
    def remove_states_below_energy(self, energy: pint.Quantity, copy=True):
        """
        Removes states with energy below a given value from the Atom
        
        Args:
            energy (pint.Quantity): Energy below which states will be removed
            copy (bool, optional): Return a copy of the Atom with states removed. Defaults to True.
        
        Returns:
            Atom: Atom with states below energy removed
        """
        atom = self.copy() if copy else self
        states = atom.states
        for s in states:
            if s.energy < energy:
                atom.remove_state(s)
        print(f"Removed {len(states) - len(atom.states)} states below {energy}")
        return atom


    """Retrieval functions"""
    @property
    def states(self) -> list[State]:
        """
        Returns a list of states in the Atom
    
        Returns:
            - list[State]: List of states in the Atom
        """
        return list(self._graph.nodes)

    @property
    def isolated_states(self) -> list[State]:
        """
        Returns a list of states without transitions

        Returns:
            - list[State]: List of states in an Atom without associated transitions
        """
        return list(nx.isolates(self._graph))


    def get_state(self, key: str | pint.Quantity):
        """
        Retrieves a state from the Atom by parsing name or energy.

        This method iterates through the states of an Atom and returns the state
        that matches the given key. The key can be either the name of the state to
        retrieve, which follows the format (configuration, quantum_numbers), or the
        energy of the state to retrieve. If no matching state is found, it returns None.

        Args:
            key (str | pint.Quantity): The name or energy of the state to retrieve.

        Returns:
            State: A State object with a name or energy matching the key. Returns None if no match is found.
        
        Raises:
            StopIteration: If no matching state is found.
        """
        if isinstance(key, str):
            return self.get_state_by_name(key)
        if isinstance(key, pint.Quantity):
            return self.get_state_by_energy(key)

    def get_state_by_name(self, key: str):
        """
        Retrieves a state from the Atom by parsing name.

        This method iterates through the states of an Atom and returns the state
        that matches the given key. The key is the name of the state to retrieve, 
        which follows the format (configuration, quantum_numbers). If no matching 
        state is found, it returns None.

        Args:
            key (str): The name of the state to retrieve.

        Returns:
            State: A State object with a name matching the key. Returns None if no match is found.
            
        Raises:
            StopIteration: If no matching state is found.
        """
        try:
            return next(state for state in self.states if state.match(key))
        except StopIteration:
            return None

    @default_units('Ry')
    def get_state_by_energy(self, energy: pint.Quantity):
        """
        Retrieves a state from the Atom by parsing energy.

        This method iterates through the states of an Atom and returns the state that is the closest match to the given energy.

        Args:
            energy (pint.Quantity): The energy of the state to retrieve.
        
        Returns:
            State: A State object with an energy closest to the given energy.
        
        Raises:
            StopIteration: If no matching state is found.

        """
        return min(self.states, key=lambda state: abs(state.energy - energy))

    def _match_term(self, key: str):
        """
        Returns a list of states matching term-key

        This method returns a list of states matching the given term-key. The term-key is the name of the state to retrieve,
        which follows the format (configuration, quantum_numbers). If no matching state is found, it returns None.

        Args:
            key (str): The name of the state to retrieve.
        
        Returns:
            list[State]: A list of State objects with a name matching the key. Returns None if no match is found.
        
        Raises:
            StopIteration: If no matching state is found.
        """
        return [s for s in self.states if s.match(key)]

    def _match_term_and_energy(self, key: str, energy: pint.Quantity):
        """
        Returns a state matching term and energy

        This method returns a state matching the given term and energy. The term is the name of the state to retrieve,
        which follows the format (configuration, quantum_numbers). The energy is the energy of the state to retrieve.
        If no matching state is found, it returns None.

        Args:
            key (str): The name of the state to retrieve.
            energy (pint.Quantity): The energy of the state to retrieve. It should be in units of energy.
        
        Returns:
            State: A State object with a name matching the key and energy closest to the energy argument. Returns None if no match is found.
        """

        matching_states = self._match_term(key)
        if len(matching_states) == 0:
            return None
        energy = set_default_units(energy, 'Ry', self._ureg)
        return min(matching_states, key=lambda state: abs(state.energy - energy))


    @property
    def transitions(self) -> list[Transition]:
        """
        Retrieve all transitions in the Atom

        This method returns a list of all transitions in the Atom. The transitions are returned in the order they are
        stored in the Atom.

        Returns:
            list[Transition]: A list of all transitions in the Atom.
        """
        return list(nx.get_edge_attributes(self._graph, 'transition').values())
    

    def transitions_from(self, state: State) -> dict[State, Transition]:
        """
        Retrieve all transitions from a given state

        This method returns a dictionary of transitions from a given state. The keys of the dictionary are the states.
        The values of the dictionary are the transitions from the given state to the corresponding state key.

        Args:
            state (State): The state to retrieve transitions from.
        
        Returns:
            dict[State, Transition]: A dictionary of transitions from the given state.
        
        Raises:
            KeyError: If the given state is not in the Atom.
        """
        return {node: edge['transition'] for node, edge in self._graph.succ[state].items()}

    def transitions_to(self, state: State) -> dict[State, Transition]:
        """
        Retrieve all transitions to a given state

        This method returns a dictionary of transitions to a given state. The keys of the dictionary are the states.
        The values of the dictionary are the transitions from the corresponding state key to the given state.

        Args:
            state (State): The state to retrieve transitions to.

        Returns:
            dict[State, Transition]: A dictionary of transitions to the given state.
        
        Raises:
            KeyError: If the given state is not in the Atom.
        """

        return {node: edge['transition'] for node, edge in self._graph.pred[state].items()}

    @default_units('nm')
    def get_transition_by_wavelength(self, wavelength: str | float | pint.Quantity):
        """
        Retrieve a transition by parsing wavelength

        This method returns a transition from the Atom by parsing wavelength. The wavelength of the transition
        is compared to the wavelength argument, and the transition with the closest wavelength is returned.

        Args:
            wavelength (str | float | pint.Quantity): The wavelength of the transition to retrieve.
        
        Returns:
            Transition: A Transition object with a wavelength closest to the wavelength argument.
        
        Raises:
            StopIteration: If no matching transition is found.
        """
        return min(self.transitions, key=lambda tr: abs(tr.wavelength - wavelength))


"""Database loading functions"""
def load_from_database(name: str, states_data: list[dict], transitions_data: list[dict],
                    energy_cutoff="inf", remove_isolated=False) -> Atom:
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
    energy_cutoff = set_default_units(energy_cutoff, 'Ry', atom._ureg)
    states = [State(**d) for d in states_data if atom._ureg(d['energy']) < energy_cutoff]
    atom.add_states(states)
    print(f"Added {len(states)} states")

    # add transitions
    transitions = []
    unmatched = 0
    for d in transitions_data:
        si = atom._match_term_and_energy(d['state_i']['term'], d['state_i']['energy'])
        sf = atom._match_term_and_energy(d['state_f']['term'], d['state_f']['energy'])
        if si is not None and sf is not None:
            tr = Transition(si, sf, d['A'])
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


def from_nist(name: str, energy_cutoff="inf", remove_isolated=False, refresh_cache=False) -> "Atom":
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
    states_data, transitions_data = nist.load_from_nist(name, refresh_cache)
    return load_from_database(name, states_data, transitions_data, energy_cutoff, remove_isolated)


