### Electric Field

Simultaneously, because majority of interest with atoms is the possibility to control them via electric fields, we introduce a second class which is electric field. Each electric field should have a frequency associated with it, as it is the property that doesn't depend on the medium. 

The class is electric field and not laser field, because we want to maintain generality. This electric field should have a function which tells us about the field and the gradient of the field at any location. These two properties are necessary to calculate the dipole rabi frequency and quadrupole rabi frequency. 

Some few examples of the fields that we included are:
- Plane Wave
- Gaussian Beam

User can define their own electric field class, which should be general enough to consider arbitrary electric fields, such as standing waves and other more complex fields.