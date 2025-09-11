#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 03/2025
# Test file for Zeeman effect calculations

import pytest
import numpy as np
from pint import get_application_registry

from atomphys import Atom, State
from atomphys.calc.zeeman import (
    zeeman_shift,
    zeeman_hamiltonian,
    lande_g_factor
)


@pytest.fixture
def ureg():
    return get_application_registry()


@pytest.fixture
def rb_atom():
    return Atom("Rb")


def test_lande_g_factor_calculation():
    """Test Landé g-factor calculation against analytical formula."""
    # Test with various J, L, S combinations
    # g_J = 1 + [J(J+1) + S(S+1) - L(L+1)] / [2J(J+1)]
    
    # Case 1: S1/2 state (L=0, S=1/2, J=1/2) -> g_J = 2
    g_s12 = lande_g_factor(J=1/2, L=0, S=1/2)
    assert g_s12 == pytest.approx(2.0)
    
    # Case 2: P1/2 state (L=1, S=1/2, J=1/2) -> g_J = 2/3
    g_p12 = lande_g_factor(J=1/2, L=1, S=1/2)
    assert g_p12 == pytest.approx(2/3)
    
    # Case 3: P3/2 state (L=1, S=1/2, J=3/2) -> g_J = 4/3
    g_p32 = lande_g_factor(J=3/2, L=1, S=1/2)
    assert g_p32 == pytest.approx(4/3)
    
    # Case 4: D5/2 state (L=2, S=1/2, J=5/2) -> g_J = 6/5
    g_d52 = lande_g_factor(J=5/2, L=2, S=1/2)
    assert g_d52 == pytest.approx(6/5)


def test_zeeman_shift_proportional_to_mj(rb_atom, ureg):
    """Test that Zeeman shifts are proportional to mJ."""
    # Get S1/2 state
    s_state = rb_atom.get_state("S1/2")
    
    # Calculate shifts for all mJ values at 1 Tesla
    B_field = 1.0 * ureg("tesla")
    shifts = {}
    
    for mJ in s_state.sublevels:
        shifts[mJ] = zeeman_shift(s_state, mJ, B_field)
    
    # Should be two values: mJ = -1/2, +1/2
    assert len(shifts) == 2
    
    # Shifts should be proportional to mJ with opposite signs
    mJ_values = sorted(list(shifts.keys()))
    assert mJ_values[0] == -0.5
    assert mJ_values[1] == 0.5
    
    # Check proportionality and sign
    assert shifts[0.5] == -shifts[-0.5]  # Equal magnitude, opposite sign
    
    # For S1/2 with g_J=2, shift should be μ_B * g_J * mJ * B
    # At 1T, μ_B = 1.4 MHz/G * 10^4 G = 14 GHz
    expected_shift = ureg("bohr_magneton") * 2 * 0.5 * B_field
    assert shifts[0.5].to("GHz") == pytest.approx(expected_shift.to("GHz"))


def test_zeeman_shift_proportional_to_B(rb_atom, ureg):
    """Test that Zeeman shifts are proportional to B field strength."""
    # Get S1/2 state
    s_state = rb_atom.get_state("S1/2")
    mJ = 0.5
    
    # Calculate shifts at different field strengths
    B_fields = [0.1, 0.5, 1.0, 2.0] * ureg("tesla")
    shifts = [zeeman_shift(s_state, mJ, B) for B in B_fields]
    
    # Check proportionality to B
    ratios = [shifts[i]/B_fields[i] for i in range(len(shifts))]
    
    # All ratios should be the same
    for i in range(1, len(ratios)):
        assert ratios[i] == pytest.approx(ratios[0], rel=1e-10)


def test_zeeman_hamiltonian_eigenvalues(rb_atom, ureg):
    """Test that Zeeman Hamiltonian eigenvalues match analytical shifts."""
    # Get P3/2 state
    p_state = rb_atom.get_state("P3/2")
    
    # Define B field
    B_field = 0.5 * ureg("tesla")
    B_vector = np.array([0, 0, 1]) * B_field.magnitude
    
    # Calculate Zeeman Hamiltonian
    H_zeeman = zeeman_hamiltonian(p_state, B_vector)
    
    # Get eigenvalues
    eigenvalues = np.linalg.eigvalsh(H_zeeman)
    
    # Calculate expected shifts analytically
    g_factor = lande_g_factor(J=3/2, L=1, S=1/2)
    expected_shifts = []
    for mJ in [-3/2, -1/2, 1/2, 3/2]:
        shift = ureg("bohr_magneton") * g_factor * mJ * B_field
        expected_shifts.append(shift.to("J").magnitude)
    
    # Sort both sets of values
    eigenvalues.sort()
    expected_shifts.sort()
    
    # Compare values
    assert eigenvalues == pytest.approx(expected_shifts)


def test_arbitrary_field_direction(rb_atom, ureg):
    """Test Zeeman effect with arbitrary field direction."""
    # Get P3/2 state
    p_state = rb_atom.get_state("P3/2")
    
    # Define B field with arbitrary direction
    B_magnitude = 0.5 * ureg("tesla")
    B_direction = np.array([1, 1, 1]) / np.sqrt(3)  # Unit vector along (1,1,1)
    B_vector = B_direction * B_magnitude.magnitude
    
    # Calculate Zeeman Hamiltonian for this field
    H_zeeman = zeeman_hamiltonian(p_state, B_vector)
    
    # Get eigenvalues
    eigenvalues = np.linalg.eigvalsh(H_zeeman)
    
    # For arbitrary field direction, we should still get the same eigenvalues
    # as for field along z-axis, as only the magnitude matters
    H_zeeman_z = zeeman_hamiltonian(p_state, np.array([0, 0, 1]) * B_magnitude.magnitude)
    eigenvalues_z = np.linalg.eigvalsh(H_zeeman_z)
    
    # Compare eigenvalues (should be the same)
    assert eigenvalues == pytest.approx(eigenvalues_z)
    
    # However, eigenvectors would be different (not tested here)


def test_weak_vs_strong_field(rb_atom, ureg):
    """Compare behavior in weak vs. strong field regimes."""
    # This is more of a demonstration than a test, as the behavior
    # depends on the relative strengths of Zeeman and hyperfine interactions
    
    # TODO: Implement this test once hyperfine + Zeeman interaction is fully supported
    pass