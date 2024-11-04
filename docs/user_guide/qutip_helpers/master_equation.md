# The Lindblad master equation

Representing atom as a density matrix, $\rho$, one can formulate the dynamics that captures coherent(laser-atom interaction) and incoherent(spontaneous decay) contributions to the dynamics of the electronic state of the atom. Such dynamics can be described via the Lindblad master equation:

$$
\frac{\partial \hat{\rho}}{\partial t}=-\frac{i}{\hbar}[\mathcal{H}, \hat{\rho}]+\mathcal{L}^{\mathrm{d}}(\hat{\rho})
$$

, where $\mathcal{H}$ is the hamiltonian of the system that captures the coherent dynamics and $\mathcal{L}^{\mathrm{d}}$ is Lindblad operator capturing incoherent contributions.

Such system can be solved for a steady state solution $\rho_{\infty}$, or analysed for a time dynamics.

Quite often translation of the situation in the lab into the form of the equations that consider all magnetic sublevels and polarizations is quite tideous and mistake prone. 

One of the goals of the package is to provide an easy interface such that the user can just define which lasers they have, which atom, which states they consider and it can spit out all relevant hamiltonians and lindblad operators to simply plug into this equations. 

Such collections of parameters can be then directly used in qutip to further calculate the density matrix dynamics.