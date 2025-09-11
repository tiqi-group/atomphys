#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 03/2025
# Test file for hyperfine structure calculations

import pytest
import numpy as np
from pint import get_application_registry

from atomphys import Atom
from atomphys.calc.hyperfine import (
    hyperfine_hamiltonian,
    hyperfine_shift,
    calculate_hyperfine_levels
)
from atomphys.state import HyperfineState


@pytest.fixture
def ureg():
    return get_application_registry()


@pytest.fixture
def rb87_atom():
    """Return a Rb-87 atom with I=3/2."""
    # This will use the default isotope for NIST data
    return Atom("Rb")


def test_hyperfine_constants(rb87_atom):
    """Test that hyperfine A, B constants match literature values."""
    # Get ground state of Rb (5s)
    ground_state = rb87_atom.get_state("S1/2")
    
    # Rb-87 ground state A coefficient should be around 3.417 GHz
    A_coeff = hyperfine_shift(ground_state, A_only=True)
    expected_A = 3.417 * rb87_atom._ureg("GHz")
    
    # Check if approximately correct (adjust precision based on your model's accuracy)
    assert A_coeff.to("GHz") == pytest.approx(expected_A, rel=1e-1)


def test_hyperfine_level_spacing(rb87_atom):
    """Test that hyperfine level spacing follows F(F+1) rule."""
    # Get ground state and calculate hyperfine levels
    ground_state = rb87_atom.get_state("S1/2")
    I = 3/2  # Nuclear spin for Rb-87
    hyperfine_levels = calculate_hyperfine_levels(ground_state, I)
    
    # For J=1/2, I=3/2, we expect F=1,2
    # Energy spacing should follow E_F = A/2 * [F(F+1) - I(I+1) - J(J+1)]
    assert len(hyperfine_levels) == 2
    
    # Sort by increasing F
    hyperfine_levels.sort(key=lambda x: x[0])
    
    # Extract F values and energies
    F_values = [level[0] for level in hyperfine_levels]
    energies = [level[1] for level in hyperfine_levels]
    
    # Check F values are correct
    assert F_values == pytest.approx([1, 2])
    
    # Energy spacing should follow F(F+1) rule
    # For ground state with J=1/2, I=3/2:
    # E_F=1 = A/2 * [1*2 - 3/2*5/2 - 1/2*3/2] = A/2 * (-2)
    # E_F=2 = A/2 * [2*3 - 3/2*5/2 - 1/2*3/2] = A/2 * (+1)
    # Energy difference should be 3*A/2
    energy_diff = energies[1] - energies[0]
    expected_diff = 3 * hyperfine_shift(ground_state, A_only=True) / 2
    
    assert energy_diff == pytest.approx(expected_diff, rel=1e-2)


def test_hyperfine_hamiltonian_eigenvalues(rb87_atom, ureg):
    """Test that hyperfine Hamiltonian eigenvalues match expected level structure."""
    # Get P1/2 state (for which both A and B terms are relevant)
    p_state = rb87_atom.get_state("P1/2")
    I = 3/2  # Nuclear spin for Rb-87
    
    # Calculate hyperfine Hamiltonian
    H_hf = hyperfine_hamiltonian(p_state, I)
    
    # Get eigenvalues
    eigenvalues = np.linalg.eigvalsh(H_hf)
    
    # Sort eigenvalues
    eigenvalues.sort()
    
    # Convert to energies
    energies = eigenvalues * ureg("J")
    
    # Check we have correct number of levels
    # For J=1/2, I=3/2, we expect F=1,2 with total dimension (2*1+1) + (2*2+1) = 3+5 = 8
    # But each level is (2*F+1) degenerate
    # So we expect 2 distinct eigenvalues
    unique_energies = set([round(e.magnitude, 10) for e in energies])
    assert len(unique_energies) == 2
    
    # Check that energy splitting follows expected pattern
    # Calculate energies another way for comparison
    expected_levels = calculate_hyperfine_levels(p_state, I)
    expected_levels.sort(key=lambda x: x[0])  # Sort by F
    expected_energies = [level[1] for level in expected_levels]
    
    # Energy spacing should be consistent
    calculated_diff = energies[-1] - energies[0]
    expected_diff = expected_energies[1] - expected_energies[0]
    
    assert calculated_diff == pytest.approx(expected_diff, rel=1e-2)


def test_hyperfine_state_creation(rb87_atom):
    """Test creating HyperfineState objects and their properties."""
    # Get ground state
    ground_state = rb87_atom.get_state("S1/2")
    I = 3/2  # Nuclear spin for Rb-87
    
    # Create hyperfine states
    hfs_F1 = HyperfineState.from_fine_structure(ground_state, I=I, F=1)
    hfs_F2 = HyperfineState.from_fine_structure(ground_state, I=I, F=2)
    
    # Check quantum numbers are correctly set
    assert hfs_F1.F == 1
    assert hfs_F2.F == 2
    assert hfs_F1.I == 3/2
    assert hfs_F2.I == 3/2
    
    # Check mF ranges are correct
    assert set(hfs_F1.sublevels) == set([-1, 0, 1])
    assert set(hfs_F2.sublevels) == set([-2, -1, 0, 1, 2])
    
    # Check energy ordering
    assert hfs_F1.energy < hfs_F2.energy


# TODO: Add test for coupling to external fields
# def test_hyperfine_zeeman(rb87_atom, ureg):
#     """Test hyperfine Zeeman effect calculations."""
#     pass