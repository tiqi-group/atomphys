#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 09/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import numpy as np
import qutip
from scipy.optimize import curve_fit

from atomphys.calc import Hamiltonians as hat
from atomphys.calc.rabi_frequency import Rabi_Frequency
from atomphys.electric_field import GaussianBeam

import pytest


def test_three_levels_decay(three_levels_atom):
    """test_three_levels_decay

    In a three leel atom, the ppulated excited state decays for dipole spontaneous emission.
    The decay time equals 1 / state.Gamma,
    and the final populations in the lower states are distributed
    according to the branching ratios between the decay channels.
    """
    atom = three_levels_atom

    n_states = sum(len(s.sublevels) for s in atom.states)

    kets = hat.kets(atom, atom.states)
    H0 = qutip.qzero(n_states)
    c_ops = hat.collapse_operators(atom, atom.states, atom._ureg)
    projectors = {}
    for s in atom.states:
        proj = sum([ket * ket.dag() for (_s, _), ket in kets.items() if _s == s])
        projectors[s] = proj

    e_ops = list(projectors.values())
    psi0 = kets[(atom.get_state("P"), -1)]

    t = np.arange(0, 0.1, 0.0005)
    sol = qutip.mesolve(H0, psi0, t, c_ops, e_ops)
    populations = np.asarray(sol.expect)

    expected_tau = (1 / atom.get_state("P").Gamma).to("us")
    expected_populations = list(atom.get_state("P").decay_branching_ratios.values())

    # ------- fit decay and check final populations

    def decay(t, tau):
        return np.exp(-t / tau)

    p, cov = curve_fit(decay, t, populations[1], (1,))
    fitted_tau = p[0] * atom._ureg("us")

    fitted_populations = populations[[0, 2], -1]

    assert expected_tau == pytest.approx(fitted_tau, rel=1e-4)
    assert expected_populations == pytest.approx(fitted_populations, rel=1e-4)


def test_two_levels_flop(two_levels_atom):
    """test_two_levels_flop

    In a two leel atom (J=0 <-> J=1) couple the m = 0 states
    with a pi polarized laser field at saturation intensity.
    The Rabi frequency equals transition.Gamma / sqrt(2).
    """
    atom = two_levels_atom

    eps = (0, 0, 1)
    e_z = (1, 0, 0)
    tr = atom.transitions[0]
    I0 = tr.saturation_intensity

    laser = GaussianBeam(
        polarization=eps,
        direction_of_propagation=e_z,
        wavelength=tr.wavelength,
        intensity=I0,
        detuning=0,
        _ureg=atom._ureg,
    )

    lasers = {laser: [tr]}
    states = atom.states
    kets = hat.kets(atom, states)
    H = hat.H0(atom, states, None) + hat.H_int(atom, states, lasers, None)
    projectors = {ss: ket * ket.dag() for ss, ket in kets.items()}
    e_ops = list(projectors.values())

    psi0 = kets[(atom.get_state("S"), 0)]
    t = np.arange(0, 0.2, 0.002)
    sol = qutip.mesolve(H, psi0, t, c_ops=None, e_ops=e_ops, progress_bar=True)

    populations = np.asarray(sol.expect).T

    def flop(t, rabi_freq):
        return 0.5 + 0.5 * np.cos(2 * np.pi * rabi_freq * t)

    p, cov = curve_fit(flop, t, populations[:, 0], (15,))

    expected_rabi_freq_at_isat = tr.Gamma.to("MHz") / np.sqrt(2) / 2 / np.pi
    fitted_rabi_freq = p[0] * atom._ureg("MHz")
    calculated_rabi_freq = (
        abs(float(Rabi_Frequency(laser, tr, mJ_i=0, mJ_f=0).to("MHz").m))
        * atom._ureg("MHz")
        / 2
        / np.pi
    )

    assert expected_rabi_freq_at_isat == pytest.approx(calculated_rabi_freq, rel=1e-7)
    assert expected_rabi_freq_at_isat == pytest.approx(fitted_rabi_freq, rel=1e-4)
