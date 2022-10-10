# atomphys_tiqi

Atomphys is a python module scripted by Matt Grau that helps to make calculations regarding any atomic properties. This is forked version, which was then subsequently edited to match the requirements of tiqi. On top of it it allows us to make quick edits to our version. It was also changed, such it is easily implementable as a python module in our script - using poetry and submodules.

In case of any queries please feel free to contact me: wadamczyk@phys.ethz.ch - I should be able to answer some queries.

## Getting started

If you would like to install atomphys as a package:
```
python -m pip install git+https://gitlab.phys.ethz.ch/tiqi-projects/optical-trap/atomphys_tiqi.git
```

You can also use it as a package:
1) include it in your project as a submodule
2) poetry add module from a local location


## Functionality: 

Big part of documentation can be found here: https://atomphys.org

### Calculation of AC Stark shifts:

What I found in the code of Matt, is that I dont fully understand where does the calculation of polarizability comes from. i.e. reference was not clearly given (I believe this is the reference https://arxiv.org/pdf/1703.09950v1.pdf). Therefore I spend a week looking into this, and thinking what could have gone wrong and what could have gone right, and decided to use different method of calculating polarizability. I didnt directly edit the polarizability file, but I added functionality to calcualte the AC Stark shifts.

This calculation follows from Thesis of Christoph Fisher Chapter 3.3.1.

Example of its use is 

```
from atomphys import Atom, Laser
import numpy as np
import matplotlib.pyplot as plt

Ca_ion = Atom('Ca+')
u = Ca_ion.units

#DEFINE PARAMETERS FOR THE BEAM
laser = Laser(λ=np.linspace(200, 1550, 10_000) * u.nm, A=0)
lam = laser.wavelength.to('nm').magnitude
I = 1*u.mW/(1*(u.micrometer)**2)
eps = (1,0,0) # Polarization vector of Electric wave
e_z = (0,0,1) # Quantization Axis

states_mJ_pairs = [('S1/2', 0.5),('P1/2', 0.5)] #pairs of statenames and mJ for which we want to calculate the things

AC_Stark_shifts = [(x[0], x[1], (Ca_ion(x[0]).ACstark(mJ=x[1], laser=laser, eps=eps, e_z=e_z, I=I, u=u)/u('k_B')).to('mK')) for x in states_mJ_pairs]

plt.figure(figsize=(9, 7))

for (state_name, mJ, ACshifts) in AC_Stark_shifts:
    plt.plot(lam, ACshifts, label=f"{state_name}: mJ = {mJ}")

plt.xlabel('Wavelength (nm)')
plt.ylabel('Trap depth (mK)')
plt.ylim(-1, 1)
plt.xlim(300,)
plt.legend()

plt.show()
```


