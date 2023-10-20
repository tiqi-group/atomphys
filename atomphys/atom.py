from copy import deepcopy
import networkx as nx
import pint
from .transition import Transition
from .util import default_units, set_default_units
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .state import State


class Atom:
    def __init__(self, name: str, _ureg=None) -> None:
        self.name = name
        self._graph = nx.DiGraph()
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg
        self._ureg.define('_2pi = 2 * pi')

    def __repr__(self) -> str:
        return f"Atom({self.name} {len(self.states)} states {len(self.transitions)} transitions)"

    def copy(self):
        return deepcopy(self)

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

    def remove_all_but_states(self, states: list[State], copy: bool = False):
        """
        Removes all states and all associated transitions from the graph associated with the Atom except the ones in the list

        Args:
            states (list[State]): List of State objects to keep in the Atom
        Returns:
            Atom: Atom object with only the states in the list

        """
        atom = self.copy() if copy else self
        for s in self.states:
            if s not in states:
                atom.remove_state(s)
        return atom

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
        return atom

    def remove_states_above_energy(self, energy: pint.Quantity, copy=True, remove_isolated=False):
        """
        Removes states with energy above a given value from the Atom

        Args:
            energy (pint.Quantity): Energy above which states will be removed
            copy (bool, optional): Return a copy of the Atom with states removed. Defaults to True.
            remove_isolated (bool, optional): Remove isolated states after removing states above energy. Defaults to True.

        Returns:
            Atom: Atom with states above energy removed
        """
        atom = self.copy() if copy else self
        states = atom.states
        for s in states:
            if s.energy > energy:
                atom.remove_state(s)
        if remove_isolated:
            atom.remove_isolated(copy=False)
        return atom

    def remove_states_below_energy(self, energy: pint.Quantity, copy=True, remove_isolated=True):
        """
        Removes states with energy below a given value from the Atom

        Args:
            energy (pint.Quantity): Energy below which states will be removed
            copy (bool, optional): Return a copy of the Atom with states removed. Defaults to True.
            remove_isolated (bool, optional): Remove isolated states after removing states below energy. Defaults to True.

        Returns:
            Atom: Atom with states below energy removed
        """
        atom = self.copy() if copy else self
        states = atom.states
        for s in states:
            if s.energy < energy:
                atom.remove_state(s)
        if remove_isolated:
            atom.remove_isolated(copy=False)
        return atom
    
    def remove_transitions_above_wavelength(self, wavelength: pint.Quantity, copy=True):
        atom = self.copy() if copy else self
        for t in atom.transitions:
            if t.wavelength > wavelength:
                atom.remove_transition(t)
        return atom
    
    def remove_transitions_below_wavelength(self, wavelength: pint.Quantity, copy=True):
        atom = self.copy() if copy else self
        for t in atom.transitions:
            if t.wavelength < wavelength:
                atom.remove_transition(t)
        return atom

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
            return self.get_state_by_term(key)
        if isinstance(key, pint.Quantity):
            return self.get_state_by_energy(key)

    def get_state_by_term(self, key: str):
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

        This method iterates through the states of an Atom and returns 
        the state that is the closest match to the given energy.

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

        This method returns a list of states matching the given term-key 
        The term-key is the name of the state to retrieve, 
        which follows the format (configuration, quantum_numbers).
        If no matching state is found, it returns None.

        Args:
            key (str): The name of the state to retrieve.

        Returns:
            list[State]: A list of State objects with a name matching the key. 
            Returns None if no match is found.

        Raises:
            StopIteration: If no matching state is found.
        """
        return [s for s in self.states if s.match(key)]

    def _match_term_and_energy(self, key: str, energy: pint.Quantity):
        """
        Returns a state matching term and energy

        This method returns a state matching the given term and energy.
        The term is the name of the state to retrieve, which follows the format 
        (configuration, quantum_numbers). The energy is the energy of the state to retrieve.
        If no matching state is found, it returns None.

        Args:
            key (str):              The name of the state to retrieve.
            energy (pint.Quantity): The energy of the state to retrieve. 
                                    It should be in units of energy.

        Returns:
            State:  A State object with a name matching the key and energy closest to the energy argument.
                    Returns None if no match is found.
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

    def transitions_from(self, state: State) -> list[Transition]:
        """
        Retrieve all transitions from a given state

        This method returns a list of transitions from a given state. The transitions are returned in the order they are
        stored in the Atom.

        Args:
            state (State): The state to retrieve transitions from.

        Returns:
            list[Transition]: A list of transitions from the given state.

        Raises:
            KeyError: If the given state is not in the Atom.
        """
        return [edge['transition'] for node, edge in self._graph.succ[state].items()]

    def transitions_to(self, state: State) -> list[Transition]:
        """
        Retrieve all transitions to a given state

        This method returns a list of all transitions to a given state. 
        The transitions are returned in the order they are stored in the Atom.

        Args:
            state (State): The state to retrieve transitions to.

        Returns:
            list[Transition]: A list of transitions to the given state.

        Raises:
            KeyError: If the given state is not in the Atom.
        """

        return [edge['transition'] for node, edge in self._graph.pred[state].items()]

    def transition_between(self, state_i: State, state_f: State) -> Transition:
        """
        Retrieve a transition between two states

        This method returns a transition between two states.
        The transition is returned in the order it is stored in the Atom.

        Args:
            state_i (State): The initial state of the transition.
            state_f (State): The final state of the transition.

        Returns:
            Transition: A transition between the given states.

        Raises:
            KeyError: If the given state is not in the Atom.
        """
        if state_f not in self._graph.succ[state_i]:
            return None
        return self._graph.get_edge_data(state_i, state_f)['transition']

    @default_units('nm')
    def get_transition_by_wavelength(self, wavelength: str | float | pint.Quantity):
        """
        Args:
            wavelength (str | float | pint.Quantity):
                The wavelength of the transition to retrieve. [m]

        Returns:
            Transition: A Transition object with a wavelength closest to the given wavelength.
        """
        return min(self.transitions, key=lambda tr: abs(tr.wavelength - wavelength))
