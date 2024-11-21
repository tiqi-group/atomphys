# atomphys <!-- omit from toc -->

<p align="left">
  <a href="https://mgrau.github.io/atomphys/"><img src="https://mgrau.github.io/atomphys/img/logo.svg" alt="atomphys logo"></a>
</p>

A Python package to help with atomic physics calculations.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14197195.svg)](https://doi.org/10.5281/zenodo.14197195)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://github.com/tiqi-group/atomphys/actions/workflows/tests.yml/badge.svg)](https://github.com/tiqi-group/atomphys/actions/workflows/tests.yml)
[![Codecov](https://img.shields.io/codecov/c/github/tiqi-group/atomphys)](https://app.codecov.io/gh/tiqi-group/atomphys)
[![PyPI](https://img.shields.io/pypi/v/atomphys)](https://pypi.org/project/atomphys/)

---

**Documentation**: [mgrau.github.io/atomphys/](https://mgrau.github.io/atomphys/)

**Source Code**: [github.com/tiqi-group/atomphys](https://github.com/mgrau/atomphys)

---

atomphys is meant to be a good starting off point for your atomic physics calculations. It can automate much of the frustrating process of searching for and compiling physical data and simple physical relations, and help you more quickly get to the good stuff.

It's designed with a natural interface and is easy to use.

## Installation

```bash
pip install atomphys
```

## Example

```python
>>> from atomphys import Atom
>>> Rb = Atom('Rb')

>>> print(Rb('S1/2').to('P1/2').λ.to('nm'))
795 nm

>>> print(Rb('P1/2').τ.to('ns'))
27.7 ns
```

## Features

- Integration with [Pint](https://pint.readthedocs.io/en/stable/) for robust handling of units
- Automatically fetch energy level and transition data from the [NIST Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database)
- Use transition data to calculate state lifetimes, polarizabilities, transition dipole moments, cross sections, and saturation intensities

## Requirements

Python 3.10+

atomphys makes extensive use of the excellent package [Pint](https://pint.readthedocs.io/en/stable/) to handle units.

## License

Atomphys is open source and released under the MIT license [(MIT)](https://opensource.org/license/mit/).

## Contributors

The project has been develoed in the [Trapped Ion Quantum Information](https://tiqi.ethz.ch/) (TIQI) group at [ETH Zurich](https://ethz.ch/), originally written by Matt Grau.

At a later stage, Carmelo Mordini and Wojtek Adamczyk contributed by extensively rewriting the package data structures, including access to different atomic databases, adding visualization tools and integration with `qutip`.
Philip Leindecker contributed with valuable advice of how to build the package such that its API is well integrateable to web-development.

Other contributors, listed alphabetically, are:

- Christoph Fischer
- Maria Radisch
- Will Cairncross
