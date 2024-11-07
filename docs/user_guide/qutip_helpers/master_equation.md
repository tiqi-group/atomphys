# The Lindblad master equation

Using density matrix, $\rho$, fromalism one can formulate the dynamics that captures coherent and incoherent contributions to the dynamics of the electronic state of the atom. Such dynamics can be described via the Lindblad master equation:

$$
\frac{\partial \hat{\rho}}{\partial t}=-\frac{i}{\hbar}[\mathcal{H}, \hat{\rho}]+\mathcal{L}^{\mathrm{d}}(\hat{\rho})
$$

, where $\mathcal{H}$ is the hamiltonian of the system that captures the coherent dynamics and $\mathcal{L}^{\mathrm{d}}$ is Lindblad operator capturing incoherent contributions.

Such system can be solved for a steady state solution $\rho_{\infty}$, or analysed for a time dynamics.

There are plenty of fantastic resources such as [qutip](https://qutip.org) to solve the master equation, however it is quite tideous and mistake prone to expand the Hamiltonian, $\mathcal{H}$, and the Lindblad operator, $\mathcal{L}^{\mathrm{d}}$, by hand.

One of the goals of the package is to provide an easy interface such that the user can just define which lasers they have, which atom, which states they consider and it can spit out all relevant hamiltonians and lindblad operators to simply plug into this equations. 

In this section we will show what kind of Hamiltonians and Lindblad operators can be generated by the package.