# import json

# import pint

# from . import _ureg
# from .data import nist
# from .state import State, StateRegistry
# from .transition import Transition, TransitionRegistry

import pint
import networkx as nx
from copy import deepcopy

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

    # database loading
    @staticmethod
    def from_nist(name: str, energy_cutoff="inf", remove_isolated=False, refresh_cache=False) -> "Atom":
        states_data, transitions_data = nist.load_from_nist(name, refresh_cache)
        return load_from_database(name, states_data, transitions_data, energy_cutoff, remove_isolated)

    def copy(self):
        return deepcopy(self)

    def add_state(self, s: State):
        s._atom = self
        s._ureg = self._ureg
        self._graph.add_node(s)

    def remove_state(self, s: State):
        """Remove state and all associated transitions
        """
        self._graph.remove_node(s)

    def add_states(self, states: list[State]):
        for s in states:
            self.add_state(s)

    def add_transition(self, tr: Transition):
        if tr.state_i not in self._graph:
            raise ValueError(f"Transition initial state {tr.state_i} not in atom")
        if tr.state_f not in self._graph:
            raise ValueError(f"Transition final state {tr.state_f} not in atom")
        tr._ureg = self._ureg
        self._graph.add_edge(tr.state_i, tr.state_f, transition=tr)

    def add_transitions(self, transitions: list[Transition]):
        for tr in transitions:
            self.add_transition(tr)

        # states

    @property
    def states(self) -> list[State]:
        return list(self._graph.nodes)

    @property
    def isolated_states(self) -> list[State]:
        return list(nx.isolates(self._graph))

    # def get_state(self, key: str | pint.Quantity):
    #     if isinstance(key, str):
    #         ret = _get_state_by_name(self.states, key)
    #     if ret is None:
    #         if isinstance(key, pint.Quantity):
    #             return _get_state_by_energy(self.states, key)

    def remove_isolated(self, copy=True):
        atom = self.copy() if copy else self
        g = atom._graph
        isolated = list(nx.isolates(g))
        g.remove_nodes_from(isolated)
        print(f"Removed {len(isolated)} states without transitions")
        return atom

    def get_state_by_name(self, key: str):
        try:
            return next(state for state in self.states if state.match(key))
        except StopIteration:
            return None

    @default_units('Ry')
    def get_state_by_energy(self, energy: str | float | pint.Quantity):
        return min(self.states, key=lambda state: abs(state.energy - energy))

    def _match_term_and_energy(self, term: str, energy: str | float | pint.Quantity):
        matching_states = self._match_term(term)
        if len(matching_states) == 0:
            return None
        energy = set_default_units(energy, 'Ry', self._ureg)
        return min(matching_states, key=lambda state: abs(state.energy - energy))

    def _match_term(self, term: str):
        return [s for s in self.states if s.match(term)]

    # transitions

    @property
    def transitions(self) -> list[Transition]:
        return list(nx.get_edge_attributes(self._graph, 'transition').values())

    def transitions_from(self, state: State) -> list[Transition]:
        return [edge['transition'] for edge in self._graph.succ[state].values()]

    def transitions_to(self, state: State) -> list[Transition]:
        return [edge['transition'] for edge in self._graph.pred[state].values()]

    @default_units('nm')
    def get_transition_by_wavelength(self, wavelength: str | float | pint.Quantity):
        return min(self.transitions, key=lambda tr: abs(tr.wavelength - wavelength))


def load_from_database(name: str, states_data: list[dict], transitions_data: list[dict],
                       energy_cutoff="inf", remove_isolated=False) -> Atom:
    """Default loading from database
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
