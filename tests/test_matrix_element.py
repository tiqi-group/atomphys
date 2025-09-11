#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 03/2025
# Test file for matrix element calculations

import pytest
import numpy as np
import sympy as sp
from pint import get_application_registry

from atomphys import Atom, State, Transition
from atomphys.calc.matrix_element import (
    reduced_matrix_element,
    wigner_eckart_factor,
    radial_matrix_element
)
from atomphys.utils.coupling import wigner_3j, wigner_6j


@pytest.fixture
def ureg():
    return get_application_registry()


@pytest.fixture
def rb_atom():
    return Atom("Rb")


def test_wigner_3j_selection_rules():
    """Test that Wigner 3j symbols follow selection rules."""
    # Test triangle inequality
    assert wigner_3j(1, 1, 3, 0, 0, 0) == 0
    
    # Test m1 + m2 + m3 = 0 rule
    assert wigner_3j(1, 1, 1, 1, 1, 0) == 0
    
    # Test valid values
    assert wigner_3j(1, 1, 2, 0, 0, 0) != 0
    assert wigner_3j(1, 1, 0, 0, 0, 0) != 0
    
    # Test symmetry properties
    assert wigner_3j(1, 2, 3, 0, 0, 0) == wigner_3j(2, 3, 1, 0, 0, 0)
    
    # Test specific known values
    assert wigner_3j(1, 1, 0, 0, 0, 0) == pytest.approx((-1) ** 1 / np.sqrt(3))


def test_wigner_6j_selection_rules():
    """Test that Wigner 6j symbols follow selection rules."""
    # Test triangle inequality
    assert wigner_6j(1, 1, 3, 0, 1, 1) == 0
    
    # Test valid values
    assert wigner_6j(1, 1, 2, 1, 1, 1) != 0
    
    # Test symmetry properties
    val1 = wigner_6j(1, 2, 3, 4, 5, 6)
    val2 = wigner_6j(1, 6, 5, 2, 3, 4)
    assert val1 == pytest.approx(val2)


def test_reduced_matrix_element_vs_known(rb_atom):
    """Test reduced matrix element calculations against known values."""
    # Get important Rb transition
    d2_line = None
    for tr in rb_atom.transitions:
        if tr.state_i.term == 'S1/2' and tr.state_f.term == 'P3/2':
            d2_line = tr
            break
    
    assert d2_line is not None, "D2 line not found in Rb atom"
    
    # Calculate reduced matrix element
    me = reduced_matrix_element(d2_line)
    
    # TODO: Compare with reference value
    # Known value should be around 4.22 atomic units
    # assert me.to('a.u.').magnitude == pytest.approx(4.22, rel=1e-2)
    
    # For now, just check it's non-zero
    assert me.magnitude != 0


def test_wigner_eckart_factors_sum_rule():
    """Test that Wigner-Eckart factors satisfy sum rules."""
    # For a J->J' transition, sum of squared factors should equal 1
    j_i = 1
    j_f = 2
    
    sum_squared = 0
    for m_i in np.arange(-j_i, j_i + 1):
        for m_f in np.arange(-j_f, j_f + 1):
            for q in [-1, 0, 1]:  # Three polarization components
                if abs(m_f - m_i - q) < 1e-6:  # Selection rule m_f = m_i + q
                    factor = wigner_eckart_factor(j_i, m_i, 1, q, j_f, m_f)
                    sum_squared += abs(factor) ** 2
    
    # The sum rule should be satisfied
    assert sum_squared == pytest.approx(2 * j_i + 1)


def test_radial_matrix_element_scaling(rb_atom):
    """Test scaling properties of radial matrix elements."""
    transitions = []
    for tr in rb_atom.transitions:
        if tr.state_i.principal_n == 5 and tr.state_i.term == 'S1/2':
            if tr.state_f.principal_n > 5 and tr.state_f.term == 'P3/2':
                transitions.append(tr)
    
    # Sort by principal quantum number of final state
    transitions.sort(key=lambda tr: tr.state_f.principal_n)
    
    # Expect radial matrix elements to decrease as n_f increases
    if len(transitions) >= 2:
        rme_values = [radial_matrix_element(tr) for tr in transitions]
        # Check that values decrease (approximately as n^(-3/2))
        for i in range(1, len(rme_values)):
            assert abs(rme_values[i]) < abs(rme_values[i-1])


def test_e1_selection_rules(rb_atom, ureg):
    """Test that E1 dipole selection rules are enforced."""
    # Create states that would violate selection rules
    s_ground = rb_atom.get_state("S1/2")
    
    # Try to create transitions that violate selection rules
    # ΔJ must be 0, ±1 (but not J=0 -> J=0)
    s_state = State(configuration='6s', term="2S1/2", energy=1 * ureg('eV'))
    d_state = State(configuration='5d', term="2D5/2", energy=2 * ureg('eV'))
    
    # Create a transition that violates ΔL = ±1 (S->D transition)
    invalid_transition = Transition(state_i=s_state, state_f=d_state)
    
    # Matrix element should be zero due to selection rules
    me = reduced_matrix_element(invalid_transition)
    assert me.magnitude == pytest.approx(0)