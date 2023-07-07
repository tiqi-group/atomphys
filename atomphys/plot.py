#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 07/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>


import numpy as np
import matplotlib.pyplot as plt
from math import nan

from .state import State
from .atom import Atom
import networkx as nx


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


def spread_energy_graph_position(states, threshold: float | None):
    states = sorted(states, key=lambda s: s.energy)
    ee = [s.energy.m for s in states]
    e1 = max(ee)
    y = np.asarray(ee) / e1
    if threshold is not None:
        y_ab = abs(y[:, None] - y)
        spread_factor = 1 / threshold / 10
        for j in range(len(y)):
            w = np.where(y_ab[j] < threshold)[0]
            ym = y[w].mean()
            d = y[w] - ym
            y[w] += d * spread_factor / len(w)

    return {s: (x_pos_angular_momentum(s), y[j]) for j, s in enumerate(states)}


def plot_atom(atom: Atom, ax=None, spread_threshold: float | None = None):

    if ax is None:
        fig, ax = plt.subplots()
    g = atom._graph
    pos = spread_energy_graph_position(g.nodes, spread_threshold)
    # pos = nx.spring_layout(g)
    labels = {k: k.term for k in g.nodes}
    nx.draw(g, pos, labels=labels, node_shape='o', node_size=100, ax=ax, clip_on=False)

    edge_labels = {edge: f"{tr.wavelength}" for edge, tr in nx.get_edge_attributes(g, 'transition').items()}
    _ = nx.draw_networkx_edge_labels(g, pos, edge_labels, clip_on=False, ax=ax)
