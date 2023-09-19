#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 09/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

from atomphys.calc import Hamiltonians as hat
import numpy as np
import qutip
from pint.testsuite.helpers import assert_quantity_almost_equal
import pytest


def test_three_levels_decay(three_levels_atom):
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

    expected_tau = (1 / atom.get_state('P').Gamma).to('us')
    expected_populations = list(atom.get_state('P').decay_branching_ratios.values())

    # ------- fit decay and check final populations
    from scipy.optimize import curve_fit

    def decay(t, tau):
        return np.exp(-t / tau)

    p, cov = curve_fit(decay, t, populations[1], (1,))
    fitted_tau = p[0] * atom._ureg('us')

    fitted_populations = populations[[0, 2], -1]

    assert_quantity_almost_equal(expected_tau, fitted_tau, rtol=1e-4)
    assert expected_populations == pytest.approx(fitted_populations, rel=1e-4)
