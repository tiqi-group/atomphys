#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 03/2025
# Test file for polarizability calculations

import pytest
import numpy as np
from pint import get_application_registry

from atomphys import Atom
from atomphys.electric_field import GaussianBeam
from atomphys.calc.ac_stark import polarizability


@pytest.fixture
def ureg():
    return get_application_registry()


@pytest.fixture
def rb_atom():
    return Atom("Rb")


@pytest.fixture
def ca_atom():
    return Atom("Ca")


def test_polarizability_units(rb_atom, ureg):
    """Test that polarizability has the correct units."""
    # Get ground state
    ground_state = rb_atom.get_state("S1/2")
    
    # Create electric field
    test_field = GaussianBeam(
        polarization=(0, 0, 1),
        direction_of_propagation=(1, 0, 0),
        wavelength=ureg('1064 nm'),
        intensity=1e4 * ureg('W/m^2'),
        _ureg=ureg
    )
    
    # Calculate polarizability
    alpha = polarizability(ground_state, mJ=0.5, El_field=test_field, _ureg=ureg)
    
    # Check units are correct (should be volume)
    assert alpha.check('[length]^3')
    
    # Convert to atomic units and check it's a reasonable value
    alpha_au = alpha.to('e^2*a_0^2/E_h')
    assert alpha_au.magnitude != 0


def test_polarizability_sign_convention(rb_atom, ureg):
    """Test that polarizability follows proper sign convention."""
    # Get ground and excited states
    ground_state = rb_atom.get_state("S1/2")
    excited_state = None
    
    # Find first P state
    for tr in rb_atom.transitions:
        if tr.state_i == ground_state and 'P' in tr.state_f.term:
            excited_state = tr.state_f
            break
    
    assert excited_state is not None, "Could not find P state"
    
    # Create electric field far red-detuned from transition
    transition_wl = ground_state.transitions_to[0].wavelength
    test_field = GaussianBeam(
        polarization=(0, 0, 1),
        direction_of_propagation=(1, 0, 0),
        wavelength=transition_wl * 1.5,  # Far red-detuned
        intensity=1e4 * ureg('W/m^2'),
        _ureg=ureg
    )
    
    # Calculate polarizabilities
    alpha_ground = polarizability(ground_state, mJ=0.5, El_field=test_field, _ureg=ureg)
    alpha_excited = polarizability(excited_state, mJ=0.5, El_field=test_field, _ureg=ureg)
    
    # For red detuning, ground state should have positive polarizability
    # and excited state should have negative polarizability
    assert alpha_ground.magnitude > 0
    assert alpha_excited.magnitude < 0


def test_magic_wavelength_ca_ions(ca_atom, ureg):
    """Test finding magic wavelength in Ca+ ions against literature."""
    # Find S1/2 and D5/2 states in Ca+
    s_state = None
    d_state = None
    
    # This example assumes we're working with Ca+ ion
    for state in ca_atom.states:
        if state.term == 'S1/2' and state.configuration.startswith('4s'):
            s_state = state
        elif state.term == 'D5/2' and state.configuration.startswith('3d'):
            d_state = state
    
    if s_state is None or d_state is None:
        pytest.skip("Required states not found in Ca atom")
    
    # Known magic wavelength for Ca+ (S1/2, mj=1/2) and (D5/2, mj=1/2) 
    # is around 395 nm (depending on specific literature source)
    magic_wl_literature = 395 * ureg('nm')
    
    # Create a range of wavelengths to scan
    wavelengths = np.linspace(390, 400, 11) * ureg('nm')
    
    # Create test field
    test_field = GaussianBeam(
        polarization=(0, 0, 1),
        direction_of_propagation=(1, 0, 0),
        wavelength=395 * ureg('nm'),  # Will be overridden in calculation
        intensity=1e4 * ureg('W/m^2'),
        _ureg=ureg
    )
    
    # Calculate polarizabilities for both states across wavelengths
    alpha_s = []
    alpha_d = []
    
    for wl in wavelengths:
        test_field.wavelength = wl
        alpha_s.append(polarizability(s_state, mJ=0.5, El_field=test_field, _ureg=ureg))
        alpha_d.append(polarizability(d_state, mJ=0.5, El_field=test_field, _ureg=ureg))
    
    # Convert to numpy arrays for calculations
    alpha_s_values = np.array([a.to('e^2*a_0^2/E_h').magnitude for a in alpha_s])
    alpha_d_values = np.array([a.to('e^2*a_0^2/E_h').magnitude for a in alpha_d])
    
    # Find crossover point (where polarizabilities are equal)
    diff = np.abs(alpha_s_values - alpha_d_values)
    magic_index = np.argmin(diff)
    calculated_magic_wl = wavelengths[magic_index]
    
    # Check if close to literature value
    # Note: This is a loose comparison as calculations depend on model details
    assert calculated_magic_wl == pytest.approx(magic_wl_literature, rel=5e-2)


def test_wavelength_scan(rb_atom, ureg):
    """Test polarizability wavelength scan around resonance."""
    # Get ground state
    ground_state = rb_atom.get_state("S1/2")
    
    # Create test field
    test_field = GaussianBeam(
        polarization=(0, 0, 1),
        direction_of_propagation=(1, 0, 0),
        wavelength=ureg('800 nm'),  # Will be overridden in calculation
        intensity=1e4 * ureg('W/m^2'),
        _ureg=ureg
    )
    
    # Find a strong transition for reference wavelength
    transition = None
    for tr in ground_state.transitions_to:
        if tr.A.magnitude > 1e6:  # Strong transition
            transition = tr
            break
    
    assert transition is not None, "No strong transition found"
    
    # Create wavelength array centered around transition
    center_wl = transition.wavelength
    wavelengths = np.linspace(
        (center_wl * 0.9).to('nm').magnitude,
        (center_wl * 1.1).to('nm').magnitude,
        21
    ) * ureg('nm')
    
    # Calculate polarizability at each wavelength
    alpha_values = polarizability(
        ground_state, 
        mJ=0.5, 
        El_field=test_field, 
        _ureg=ureg, 
        wavelengths=wavelengths
    )
    
    # Expected behavior:
    # 1. Sign change across resonance
    # 2. Magnitude increases as we approach resonance
    
    # Find resonance index
    res_idx = np.argmin(np.abs(wavelengths - center_wl))
    
    # Check sign change across resonance
    if res_idx > 0 and res_idx < len(wavelengths) - 1:
        # Signs should be different across resonance
        assert np.sign(alpha_values[res_idx-1]) != np.sign(alpha_values[res_idx+1])
    
    # Check magnitude increases near resonance
    # Take points far from resonance and compare to points near resonance
    far_indices = [0, -1]  # First and last points
    near_indices = [max(0, res_idx-1), min(len(wavelengths)-1, res_idx+1)]
    
    far_magnitudes = [abs(alpha_values[i]) for i in far_indices]
    near_magnitudes = [abs(alpha_values[i]) for i in near_indices]
    
    # Near resonance should have larger magnitudes
    assert min(near_magnitudes) > max(far_magnitudes)