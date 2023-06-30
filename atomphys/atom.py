import json

import pint

from . import _ureg
from .data import nist
from .state import State, StateRegistry
from .transition import Transition, TransitionRegistry


class Atom:
    """An atom object, containing states and transitions

    Attributes:
        name (str): The name of the atom
    """

    _ureg: pint.UnitRegistry
    name: str
    __states: StateRegistry
    __transitions: TransitionRegistry

    def __init__(self, atom=None, ureg=None, refresh_cache=False):
        self._ureg = ureg if ureg is not None else _ureg
        self.name = ""
        self.__states = StateRegistry(atom=self)
        self.__transitions = TransitionRegistry(atom=self)

        if atom is None:
            return

        try:
            self.load(atom)
        except FileNotFoundError:
            self.load_nist(atom, refresh_cache)

        self.__transitions.sort()
        for state in self.__states:
            state.transitions.sort()

    def __call__(self, key):
        try:
            return self.__states(key)
        except KeyError:
            pass
        try:
            return self.__transitions(key)
        except KeyError:
            pass
        raise KeyError(f"no property of atom {self.name} {key}")

    def __repr__(self):
        ground_state = (
            f"Ground State: {self.__states[0].term}\n" if self.__states else ""
        )
        return (
            f"{ground_state}{len(self.__states)} States\n"
            f"{len(self.__transitions)} Transitions"
        )

    def _load_states(self, states):
        self.__states = StateRegistry(
            sorted([State(**state, atom=self) for state in states]),
            atom=self,
        )

    def _load_transitions(self, transitions):
        self.__transitions = TransitionRegistry(
            [Transition(**transition, atom=self) for transition in transitions],
            atom=self,
        )

    def to_dict(self):
        return {
            "name": self.name,
            "states": self.__states.to_dict(),
            "transitions": self.__transitions.to_dict(),
        }

    def save(self, filename):
        with open(filename, "w") as file:
            json.dump(self.to_dict(), file, indent=4, ensure_ascii=False)

    def load(self, filename):
        with open(filename) as file:
            data = json.load(file)

        self.name = data["name"]
        self._load_states(data["states"])
        self._load_transitions(data["transitions"])

    def load_nist(self, name, refresh_cache=False):
        self.name = name
        states_data, transitions_data = nist.load_from_nist(name, refresh_cache)
        self._load_states(states_data)
        self._load_transitions(transitions_data)

    @property
    def states(self) -> StateRegistry:
        """StateRegistry: the atomic states."""
        return self.__states

    @property
    def transitions(self) -> TransitionRegistry:
        """TransitionRegistry: the atomic transitions."""
        return self.__transitions

    @property
    def units(self):
        """pint.UnitRegistry(): readonly access to the pint UnitRegistry used by the atom."""
        return self._ureg
