# atomphys_tiqi

Atomphys is a python module scripted by Matt Grau that helps to make calculations regarding any atomic properties. This is forked version, which was then subsequently edited to match the requirements of tiqi. On top of it it allows us to make quick edits to our version. It was also changed, such it is easily implementable as a python module in our script - using poetry and submodules.

In case of any queries please feel free to contact me: <wadamczyk@phys.ethz.ch> - I should be able to answer some queries.

## Getting started

If you would like to install atomphys as a package:

```console
pip install git+https://gitlab.phys.ethz.ch/tiqi-projects/optical-trap/atomphys_tiqi.git
```

You can also use it as a package:

1) include it in your project as a submodule
2) poetry add module from a local location

## Functionality

Big part of documentation can be found here: <https://atomphys.org>

### Calculation of AC Stark shifts

What I found in the code of Matt, is that I dont fully understand where does the calculation of polarizability comes from. i.e. reference was not clearly given (I believe this is the reference <https://arxiv.org/pdf/1703.09950v1.pdf>). Therefore I spend a week looking into this, and thinking what could have gone wrong and what could have gone right, and decided to use different method of calculating polarizability. I didnt directly edit the polarizability file, but I added functionality to calcualte the AC Stark shifts.

This calculation follows from Thesis of Christoph Fisher Chapter 3.3.1.

Example of its use is

```python
from atomphys import Atom, Laser, State
import numpy as np
import matplotlib.pyplot as plt

Ca = Atom('Ca')

CaS = Ca('1S0')
CaP = Ca('1P1')
CaD = Ca('1D2')
CatP = Ca('3P1')

#Definition of Laser Tweezer
u = Ca._ureg   # Units Registry
eps  = (1,0,0) # Polarisation
e_z  = (0,0,1)
P    = 200*u('mW')
σ    = 10 * u('um')
x = np.linspace(-20, 20, 100)*u('um')
I0   = P/σ**2
laser_tweezer = Laser(λ=np.linspace(300, 1200, 600), ureg=u, I=I0)

CaS_0 = CaS.ACstark(mJ=0, laser=laser_tweezer, eps=eps, e_z=e_z, I=I0, u=u)
CaP_m1 = CaP.ACstark(mJ=-1, laser=laser_tweezer, eps=eps, e_z=e_z, I=I0, u=u)
CaP_0 = CaP.ACstark(mJ=0, laser=laser_tweezer, eps=eps, e_z=e_z, I=I0, u=u)
CaP_p1 = CaP.ACstark(mJ=1, laser=laser_tweezer, eps=eps, e_z=e_z, I=I0, u=u)
```

# atomphys

<p align="left">
  <a href="https://mgrau.github.io/atomphys/"><img src="https://mgrau.github.io/atomphys/img/logo.svg" alt="atomphys logo"></a>
</p>

<!--intro-start-->

A Python package to help with atomic physics calculations.

[![Tests](https://github.com/mgrau/atomphys/actions/workflows/tests.yml/badge.svg)](https://github.com/mgrau/atomphys/actions/workflows/tests.yml)
[![Codecov](https://img.shields.io/codecov/c/github/mgrau/atomphys)](https://app.codecov.io/gh/mgrau/atomphys)
[![GitHub](https://img.shields.io/github/license/mgrau/atomphys)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/atomphys)](https://pypi.org/project/atomphys/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mgrau/atomphys/main?urlpath=%2Ftree%2Fexamples)

---

**Documentation**: [mgrau.github.io/atomphys/](https://mgrau.github.io/atomphys/)

**Source Code**: [github.com/mgrau/atomphys](https://github.com/mgrau/atomphys)

---

atomphys is meant to be a good starting off point for your atomic physics calculations. It can automate much of the frustrating process of searching for and compiling physical data and simple physical relations, and help you more quickly get to the good stuff.

It's designed with a natural interface and is easy to use.

## Example

```python
>>> from atomphys import Atom
>>> Rb = Atom('Rb')

>>> print(Rb('S1/2').to('P1/2').λ.to('nm'))
795 nm

>>> print(Rb('P1/2').τ.to('ns'))
27.7 ns
```

## Installation

To install atomphys, simply use pip:

```console
pip install atomphys
```

## Features

- Integration with [Pint](https://pint.readthedocs.io/en/stable/) for robust handling of units
- Automatically fetch energy level and transition data from the [NIST Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database)
- Use transition data to calculate state lifetimes, polarizabilities, transition dipole moments, cross sections, and saturation intensities

## Requirements

Python 3.6+

atomphys makes extensive use of the excellent package [Pint](https://pint.readthedocs.io/en/stable/) to handle units.

<!--intro-end-->
