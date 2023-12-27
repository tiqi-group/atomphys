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
from .util import set_default_units
from collections import defaultdict

# def plot_atom(atom: Atom, ax=None, energy_units='Ry', max_energy=nan, min_energy=nan):
#     if ax is None:
#         fig, ax = plt.subplots()
#     if max_energy is not nan:
#         atom = atom.remove_states_above_energy(max_energy)
#     if min_energy is not nan:
#         atom = atom.remove_states_below_energy(min_energy)

#     g = atom._graph
#     pos = JE_graph_position(g.nodes, energy_units)
#     # pos = nx.spring_layout(g)
#     node_labels = {k: k.term for k in g.nodes}
#     edge_labels = {edge: f"{tr.wavelength}" for edge, tr in nx.get_edge_attributes(g, 'transition').items()}

#     # https://petercbsmith.github.io/marker-tutorial.html
#     node_path = Path([(-1, 0), (1, 0)], [Path.MOVETO, Path.LINETO])  # just a simple horizontal line
#     node_color = ['k' if s not in atom.isolated_states else 'C1' for s in atom.states]
#     edge_color = [_wavelength_to_rgb(tr.wavelength.to('nm').m) for tr in atom.transitions]
#     nx.draw_networkx_nodes(g, pos, node_shape=node_path, node_color=node_color, node_size=200, linewidths=2, ax=ax)
#     nx.draw_networkx_labels(g, pos, node_labels, font_size=9, verticalalignment='bottom', ax=ax, clip_on=False)
#     nx.draw_networkx_edges(g, pos, edge_color=edge_color, node_size=50, ax=ax)
#     nx.draw_networkx_edge_labels(g, pos, edge_labels, font_size=9, clip_on=False, ax=ax)
#     ax.tick_params(top=False, right=False, reset=True)
#     ax.set(xlabel="Angular momentum [L]", ylabel=f"Energy [{energy_units}]", title=atom)


def plot_atom(
    atom: Atom,
    ax=None,
    energy_units="Ry",
    plot_transitions=True,
    max_energy=None,
    min_energy=None,
    max_J=None,
    state_label="name",
    introduce_offset=False,
    remove_isolated=False,
    energy_threshold=0.001,
    x_offset_factor=0.1,
    y_offset_factor=0.01,
):
    if ax is None:
        fig, ax = plt.subplots()

    if max_energy is not None:
        max_energy = set_default_units(max_energy, energy_units)
        atom = atom.remove_states_above_energy(max_energy, remove_isolated=remove_isolated, copy=True)
    if min_energy is not None:
        min_energy = set_default_units(min_energy, energy_units)
        atom = atom.remove_states_below_energy(min_energy, copy=True)
    if max_J is not None:
        for s in atom.states:
            if x_pos_angular_momentum(s) > max_J:
                atom.remove_state(s)

    g = atom._graph

    # Get original positions
    pos = JE_graph_position(g.nodes, energy_units, max_energy, min_energy)

    # Calculate offsets for nodes with similar energy and same angular momentum
    if introduce_offset:
        pos = offset_similar_nodes(pos, energy_threshold, x_offset_factor, y_offset_factor)

    node_labels = {k: getattr(k, state_label) for k in g.nodes}
    edge_labels = {edge: f"{tr.wavelength:~0.1fP}" for edge, tr in nx.get_edge_attributes(g, "transition").items()}

    # https://petercbsmith.github.io/marker-tutorial.html
    node_path = Path([(-1, 0), (1, 0)], [Path.MOVETO, Path.LINETO])  # just a simple horizontal line
    node_color = ["k" if s not in atom.isolated_states else "C1" for s in atom.states]
    edge_color = [_wavelength_to_rgb(tr.wavelength.to("nm").m) for tr in atom.transitions]
    nx.draw_networkx_nodes(g, pos, node_shape=node_path, node_color=node_color, node_size=200, linewidths=2, ax=ax)
    nx.draw_networkx_labels(g, pos, node_labels, font_size=9, verticalalignment="bottom", ax=ax, clip_on=False)
    if plot_transitions:
        nx.draw_networkx_edges(g, pos, edge_color=edge_color, node_size=50, ax=ax)
        nx.draw_networkx_edge_labels(
            g,
            pos,
            edge_labels,
            font_size=9,
            clip_on=False,
            ax=ax,
            bbox=dict(facecolor="none", edgecolor="none", boxstyle="round,pad=0.5"),
        )
    # nx.draw_networkx_edge_labels(g, pos, edge_labels, font_size=9, clip_on=False, ax=ax)
    ax.tick_params(top=False, right=False, reset=True)
    ax.set(xlabel="Angular momentum [L]", ylabel=f"Energy [{energy_units}]", title=atom)


def offset_similar_nodes(pos, energy_threshold, x_offset_factor, y_offset_factor):
    # Group nodes by x-coordinate
    x_dict = defaultdict(list)
    for node, (x, y) in pos.items():
        x_dict[x].append((node, y))

    # For each group, sort by y-coordinate and add offsets to nodes with similar y-coordinates
    new_pos = {}
    for x, node_list in x_dict.items():
        node_list.sort(key=lambda item: item[1])  # sort by y-coordinate
        for i, (node, y) in enumerate(node_list):
            x_offset = 0
            y_offset = 0
            for j in range(i):
                if i > 0:
                    prev_node, prev_y = node_list[i - j]
                    if abs(prev_y - y) < energy_threshold:
                        x_offset += x_offset_factor
                        y_offset += y_offset_factor
            new_pos[node] = (x + x_offset, y + y_offset)

    return new_pos


def plot_energy_histogram(atom: Atom, unit="Ry", ax=None, **hist_kwargs):
    """
    Plot a histogram of the energy levels of the atom.

    This function is a wrapper around matplotlib.pyplot.hist.

    Args:
        atom (Atom): The atom to plot.
        unit (str, optional): The unit of the energy. Defaults to 'Ry'.
        ax (matplotlib.axes.Axes, optional): The axes to plot on. Defaults to None.
        **hist_kwargs: Additional keyword arguments passed to matplotlib.pyplot.hist.

    """
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


def JE_graph_position(states, energy_units: str, max_energy=None, min_energy=None):
    # if max_energy is not None:
    #     states = [s for s in states if s.energy <= max_energy]
    # if min_energy is not None:
    #     states = [s for s in states if s.energy >= min_energy]
    return {s: (x_pos_angular_momentum(s), s.energy.to(energy_units).m) for s in states}


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
    """taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range

    https://stackoverflow.com/a/44960748
    """
    # TODO: make this a colormap
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 750:
        A = 1.0
    else:
        A = 0.5
    if wavelength < 380:
        wavelength = 380.0
    if wavelength > 750:
        wavelength = 750.0
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
