#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>


import matplotlib.pyplot as plt
from matplotlib.path import Path

import networkx as nx
from math import nan

from .state import State
from .atom import Atom


def plot_atom(atom: Atom, ax=None, energy_units='Ry'):
    if ax is None:
        fig, ax = plt.subplots()
    g = atom._graph
    pos = JE_graph_position(g.nodes, energy_units)
    # pos = nx.spring_layout(g)
    node_labels = {k: k.term for k in g.nodes}
    edge_labels = {edge: f"{tr.wavelength}" for edge, tr in nx.get_edge_attributes(g, 'transition').items()}

    # https://petercbsmith.github.io/marker-tutorial.html
    node_path = Path([(-1, 0), (1, 0)], [Path.MOVETO, Path.LINETO])  # just a simple horizontal line
    node_color = ['k' if s not in atom.isolated_states else 'C1' for s in atom.states]
    edge_color = [_wavelength_to_rgb(tr.wavelength.to('nm').m) for tr in atom.transitions]
    nx.draw_networkx_nodes(g, pos, node_shape=node_path, node_color=node_color, node_size=200, linewidths=2, ax=ax)
    nx.draw_networkx_labels(g, pos, node_labels, font_size=9, verticalalignment='bottom', ax=ax, clip_on=False)
    nx.draw_networkx_edges(g, pos, edge_color=edge_color, node_size=50, ax=ax)
    nx.draw_networkx_edge_labels(g, pos, edge_labels, font_size=9, clip_on=False, ax=ax)
    ax.tick_params(top=False, right=False, reset=True)
    ax.set(xlabel="Angular momentum [L]", ylabel=f"Energy [{energy_units}]", title=atom)


def plot_energy_histogram(atom: Atom, unit='Ry', ax=None, **hist_kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    energies = [s.energy.to(unit).m for s in atom.states]
    _kw = dict(bins=50)
    _kw.update(hist_kwargs)
    return ax.hist(energies, **_kw)


def x_pos_angular_momentum(state: State):
    qn = state.quantum_numbers
    if qn.ionization_limit is not None:
        x = 0
    elif qn.L is not None:
        x = qn.L
    elif qn.K is not None:
        x = qn.K
    elif qn.J is not None:
        x = qn.J
    else:
        x = nan
    return x


def JE_graph_position(states, energy_units: str):
    return {
        s: (x_pos_angular_momentum(s), s.energy.to(energy_units).m)
        for s in states
    }

# def spreaded_JE_graph_position(states, energy_units: str, threshold: float | None):
#     states = sorted(states, key=lambda s: s.energy)
#     ee = [s.energy.to(energy_units).m for s in states]
#     e1 = max(ee)
#     y = np.asarray(ee)
#     if threshold is not None:
#         threshold = threshold * e1
#         y_ab = abs(y[:, None] - y)
#         spread_factor = 1 / threshold / 10
#         for j in range(len(y)):
#             w = np.where(y_ab[j] < threshold)[0]
#             ym = y[w].mean()
#             d = y[w] - ym
#             y[w] += d * spread_factor / len(w)

#     return {s: (x_pos_angular_momentum(s), y[j]) for j, s in enumerate(states)}


def _wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an 
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range

    https://stackoverflow.com/a/44960748
    '''
    # TODO: make this a colormap
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 750:
        A = 1.
    else:
        A = 0.5
    if wavelength < 380:
        wavelength = 380.
    if wavelength > 750:
        wavelength = 750.
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R, G, B, A)
