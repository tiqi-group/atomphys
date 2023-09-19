from atomphys import State, Transition, Atom
from pint import get_application_registry
import pytest


from math import pi


@pytest.fixture(scope="module")
def rubidium():
    return Atom("Rb")


@pytest.fixture(scope="module")
def two_levels_atom():
    ureg = get_application_registry()

    s_s = State(configuration='1s', term="1S0", energy=0)
    s_p = State('2p', '1P1', energy=ureg('400 nm').to('Ry'))
    tr_sp = Transition(state_i=s_s, state_f=s_p, A=1 / ureg('1 ns') / 2 / pi)

    a2 = Atom('Simple_two_levels')
    a2.add_states([s_s, s_p])
    a2.add_transition(tr_sp)

    return a2


@pytest.fixture(scope="module")
def three_levels_atom():
    ureg = get_application_registry()

    s_s = State(configuration='1s', term="1S0", energy=0)
    s_p = State('2p', '1P1', energy=ureg('400 nm').to('Ry'))
    s_d = State('3d', '1D2', energy=ureg('600 nm').to('Ry'))

    tr_sp = Transition(state_i=s_s, state_f=s_p, A=1 / ureg('1 ns') / 2 / pi)
    tr_sd = Transition(s_s, s_d, A=1 / ureg('1 s') / 2 / pi)
    tr_dp = Transition(s_d, s_p, A=1 / ureg('10 ns') / 2 / pi)

    a3 = Atom('Simple_three_levels')
    a3.add_states([s_s, s_p, s_d])
    a3.add_transitions([tr_sp, tr_sd, tr_dp])

    return a3
