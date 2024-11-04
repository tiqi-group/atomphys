### Electric Field Class

Because in atomic physics we want to control atoms via electric fields, we introduce electric field class. We decided to make it an electric field and not laser field, because we want to maintain high degree of generality. The general concept of the electric field should have associated frequency, field, and gradient of the field with it. These properties are necessary to calculate the dipole rabi frequency and quadrupole rabi frequency. If the user wants to consider more specific electric fields, it can be done by inheriting from the electric field class.

Some few examples of the fields that we included are:
- Plane Wave
- Gaussian Beam

User can define their own electric field class, such as standing waves and other more complex fields. Because of the general notion of the field, and gradient the user can analyse how the spatial structure of the field can influence the dynamics of the atom (Magnus Effect).