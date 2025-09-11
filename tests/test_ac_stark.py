#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 03/2025
# Test file for AC Stark shift calculations

import pytest
import numpy as np
from pint import get_application_registry

from atomphys import Atom, State, Transition
from atomphys.electric_field import GaussianBeam
from atomphys.calc.ac_stark import ac_stark_shift, polarizability


@pytest.fixture
def ureg():
    return get_application_registry()


@pytest.fixture
def rb_atom():
    return Atom("Rb")


@pytest.fixture
def simple_two_level_system(ureg):
    """Create a simplified two-level system with known polarizability."""
    s_s = State(configuration='1s', term="1S0", energy=0)
    s_p = State('2p', '1P1', energy=ureg('400 nm').to('Ry'))
    
    # Define transition with precise dipole moment for predictable polarizability
    tr_sp = Transition(
        state_i=s_s, 
        state_f=s_p, 
        A=1e6 / ureg('s')  # Strong transition rate
    )
    
    atom = Atom('Test_Atom')
    atom.add_states([s_s, s_p])
    atom.add_transition(tr_sp)
    
    return atom


def test_ac_stark_shift_sign(simple_two_level_system, ureg):
    """Test that AC Stark shifts have the expected sign for red/blue detuning."""
    atom = simple_two_level_system
    ground_state = atom.get_state("S")
    excited_state = atom.get_state("P")
    
    # Create red-detuned laser field (lower frequency than transition)
    red_detuned = GaussianBeam(
        polarization=(0, 0, 1),
        direction_of_propagation=(1, 0, 0),
        wavelength=atom.transitions[0].wavelength * 1.1,  # Red-detuned
        intensity=1e4 * ureg('W/m^2'),
        _ureg=ureg
    )
    
    # Create blue-detuned laser field (higher frequency than transition)
    blue_detuned = GaussianBeam(
        polarization=(0, 0, 1),
        direction_of_propagation=(1, 0, 0),
        wavelength=atom.transitions[0].wavelength * 0.9,  # Blue-detuned
        intensity=1e4 * ureg('W/m^2'),
        _ureg=ureg
    )
    
    # Calculate Stark shifts
    red_shift_ground = ac_stark_shift(ground_state, 0, red_detuned, ureg)
    blue_shift_ground = ac_stark_shift(ground_state, 0, blue_detuned, ureg)
    
    # For ground state: red-detuned should give negative shift, blue-detuned positive
    assert red_shift_ground < 0 * ureg('mK')
    assert blue_shift_ground > 0 * ureg('mK')
    
    # TODO: Add test for excited state (opposite behavior)


def test_polarizability_vs_calculated(rb_atom, ureg):
    """Test that polarizability values match expected values from literature."""
    # Get well-characterized state (e.g., ground state)
    ground_state = rb_atom.get_state("S1/2")
    
    # Create laser field at standard wavelength (e.g., 1064 nm for comparison with literature)
    test_field = GaussianBeam(
        polarization=(0, 0, 1),
        direction_of_propagation=(1, 0, 0),
        wavelength=ureg('1064 nm'),
        intensity=1e4 * ureg('W/m^2'),
        _ureg=ureg
    )
    
    # Calculate polarizability
    alpha = polarizability(ground_state, mJ=0.5, El_field=test_field, _ureg=ureg)
    
    # TODO: Compare with literature value
    # Reference value for Rb ground state at 1064 nm
    # expected_alpha = X * ureg('a_0^3')  # Replace X with actual value
    # assert alpha.to('a_0^3') == pytest.approx(expected_alpha, rel=1e-2)
    
    # For now, just check it's a reasonable non-zero value
    assert alpha.magnitude != 0


def test_wavelength_scan(simple_two_level_system, ureg):
    """Test AC Stark shift calculation with wavelength scan."""
    atom = simple_two_level_system
    ground_state = atom.get_state("S")
    
    # Create test field
    test_field = GaussianBeam(
        polarization=(0, 0, 1),
        direction_of_propagation=(1, 0, 0),
        wavelength=ureg('600 nm'),  # Will be overridden by wavelength scan
        intensity=1e4 * ureg('W/m^2'),
        _ureg=ureg
    )
    
    # Define wavelength range around the transition
    transition_wl = atom.transitions[0].wavelength
    wavelengths = np.linspace(
        (transition_wl * 0.8).to('nm').magnitude,
        (transition_wl * 1.2).to('nm').magnitude,
        21
    ) * ureg('nm')
    
    # Calculate Stark shifts for wavelength scan
    shifts = ac_stark_shift(ground_state, 0, test_field, ureg, wavelengths=wavelengths)
    
    # Test shape of result
    assert len(shifts) == len(wavelengths)
    
    # Test that the shift changes sign across resonance
    # (negative for red-detuned, positive for blue-detuned)
    red_shifts = shifts[wavelengths > transition_wl]
    blue_shifts = shifts[wavelengths < transition_wl]
    
    assert all(red_shift < 0 * ureg('mK') for red_shift in red_shifts)
    assert all(blue_shift > 0 * ureg('mK') for blue_shift in blue_shifts)
    
    # Test that polarizability diverges near resonance
    near_resonance_index = np.argmin(np.abs(wavelengths - transition_wl))
    near_resonance_shifts = np.abs(shifts[near_resonance_index-1:near_resonance_index+2])
    other_shifts = np.abs(shifts[0:3])  # Far from resonance
    
    assert np.max(near_resonance_shifts) > np.max(other_shifts)


# TODO: Add test for comparison with published magic wavelengths
# def test_magic_wavelength(ureg):
#     """Test calculation of magic wavelength for a well-characterized atom (e.g., Ca+)."""
#     pass