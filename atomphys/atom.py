# import json

# import pint

# from . import _ureg
# from .data import nist
# from .state import State, StateRegistry
# from .transition import Transition, TransitionRegistry

import networkx as nx
from .state import State
from .transition import Transition

from .util import default_units
import pint


class Atom:
    def __init__(self, _ureg=None) -> None:
        self._graph = nx.DiGraph()
        self._ureg = pint.get_application_registry() if _ureg is None else _ureg

    def add_state(self, state: State):
        self._graph.add_node(state)

    def add_states(self, states: list[State]):
        self._graph.add_nodes_from(states)

    def add_transition(self, tr: Transition):
        if tr.state_i not in self._graph:
            raise ValueError(f"Transition initial state {tr.state_i} not in atom")
        if tr.state_f not in self._graph:
            raise ValueError(f"Transition final state {tr.state_f} not in atom")
        self._graph.add_edge(tr.state_i, tr.state_f, transition=tr)

    def add_transitions(self, transitions: list[Transition]):
        for tr in transitions:
            self.add_transition(tr)

    # states

    @property
    def states(self) -> list[State]:
        return list(self._graph.nodes)

    # def get_state(self, key: str | pint.Quantity):
    #     if isinstance(key, str):
    #         ret = _get_state_by_name(self.states, key)
    #     if ret is None:
    #         if isinstance(key, pint.Quantity):
    #             return _get_state_by_energy(self.states, key)

    def get_state_by_name(self, key: str):
        try:
            return next(state for state in self.states if state.match(key))
        except StopIteration:
            return None

    @default_units('Ry')
    def get_state_by_energy(self, energy: str | float | pint.Quantity):
        return min(self.states, key=lambda state: abs(state.energy - energy))

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
