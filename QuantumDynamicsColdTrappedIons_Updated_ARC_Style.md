
# Atomphys Documentation: 

This document provides a guide to the Python module designed for atomic physics calculations. The idea behind this module is to make it easier for anyone to explore atomic physics, properties of atoms, and understand their electron dynamics in light fields. It aims to reduce the reinvention of the wheel. 

I hope this will be helpfull for many.

## Getting Started

### Installation

Install the necessary libraries via pip:

```sh
pip install git+https://gitlab.phys.ethz.ch/tiqi-projects/optical-trap/atomphys_tiqi.git
```
## Static Structure
### Atom
### State
### Transition
### Field

## Calculations

### Matrix Elements:

Before diving into the specific formulas and derivations used, it's worth pausing to consider the key concepts in atomic physics that we focus on and their significance. Atomic physics heavily relies on the calculation of matrix elements, denoted as $M_{if} = \left\langle f \mid \hat{M} \mid i \right\rangle$. This calculation involves the projection of an operator $\hat{M}$, which acts on an initial state $\mid i \rangle$, onto a final state $\mid f \rangle$. The value of $M_{if}$ is crucial because it measures the extent of interaction, or coupling, between the states $\mid i \rangle$ and $\mid f \rangle$ via the operator $\hat{M}$, when $\hat{M}$ is part of the Hamiltonian $\hat{H}$. Understanding this interaction provides insights into the dynamics of the overall state $\mid \psi(t) \rangle$.

However, computing these matrix elements is challenging due to their basis involved. Is there a way to simplify (break-down) this process using symmetry principles? Fortunately, the answer is yes, thanks to the Wigner-Eckart Theorem.

The Wigner-Eckart Theorem states that the matrix elements of spherical tensor operators, within the framework of angular momentum eigenstates, can be broken down into two components: one that is independent of the angular momentum orientation and another that is a Clebsch-Gordan coefficient.(https://en.wikipedia.org/wiki/Wigner%E2%80%93Eckart_theorem) 

$\left\langle j m\left|T_q^{(k)}\right| j^{\prime} m^{\prime}\right\rangle=\left\langle j^{\prime} m^{\prime} k q \mid j m\right\rangle\left\langle j\left\|T^{(k)}\right\| j^{\prime}\right\rangle = [j] \left(\begin{array}{ccc}j^{\prime} & k & j^* \\ m^{\prime} & q & m\end{array}\right)  \left\langle j\left\|T^{(k)}\right\| j^{\prime}\right\rangle$

Here, $[j]$ represents the dimension of the representation space of $j$, and $j^*$ indicates a complex conjugate representation to $j$. The parentheses contain the Wigner-3j symbol, and $\left\langle \mid\mid . \mid\mid \right\rangle$ denotes a reduced matrix element.

This theorem implies that to calculate any matrix element within an arbitrary $\mid n, j, m \rangle$ basis, we can initially compute the irreducible matrix elements between any $\mid n, j \rangle$ pairs. This foundational step allows us to determine any desired matrix element more clearly.

As the reduced matrix element, doesn't depend on basis this is the most basic block from which we will start the calculation.


#### Reduced Matrix Elements (Dipole and Quadrupole):
Following James (equation 5.9) and (equation 5.10) [1]

$A_{12}^{(E 1)} \equiv \sum_{m=-j}^j \bar{A}_{12}^{(E 1)}=\frac{4 c \alpha k_{12}^3}{3\left(2 j^{\prime}+1\right)}\left|\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle\right|^2$

$A_{12}^{(E 2)} \equiv \sum_{m=-j}^j \bar{A}_{12}^{(E 2)}=\frac{c \alpha k_{12}^5}{15\left(2 j^{\prime}+1\right)}\left|\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle\right|^2$

we can rewrite 

$\left|\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle\right| = \left( \frac{3\left(2 j^{\prime}+1\right)}{4 c \alpha k_{12}^3} A_{12}^{(E 1)} \right)^{\frac{1}{2}}$

$\left|\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle\right| = \left(\frac{15\left(2 j^{\prime}+1\right)}{c \alpha k_{12}^5} A_{12}^{(E 2)}\right)^{\frac{1}{2}}$

Equivalent result can be derived from Wigner-Weisskopf Decay. Example of such derivation can be found in Scully's Quantum Optics[2]


>`atomphys.calc.matrix_element.reduced_dipole_matrix_element(A, k, J_f, _ureg)`
>
>Calculates the reduced dipole matrix element.
>
>**Parameters:**
>- `A`: Einstein coefficient [1/s].
>- `k`: Wavenumber of the transition [1/m].
>- `J_f`: Angular momentum of the upper state.
>- `_ureg`: Unit registry.
>
>**Returns:** The reduced dipole matrix element [$a_0$].


>`atomphys.calc.matrix_element.reduced_quadrupole_matrix_element(A, k, J_f, _ureg)`
>
>Calculates the reduced dipole matrix element.
>
>**Parameters:**
>- `A`: Einstein coefficient [1/s].
>- `k`: Wavenumber of the transition [1/m].
>- `J_f`: Angular momentum of the upper state.
>- `_ureg`: Unit registry.
>
>**Returns:** The reduced quadrupole matrix element [$a_0^2$].


#### Matrix Elements

Now in order to transform it to actual quantities that we want we can use Wigner-Eckart theorem.

$\left\langle 1\left|\hat{r}_i\right| 2\right\rangle \epsilon_i=\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle \sum_{q=-1}^1\left(\begin{array}{ccc}j & 1 & j^{\prime} \\ -m_j & q & m_j^{\prime}\end{array}\right) c_i^{(q)} \epsilon_i$

$\left\langle 1\left|\hat{r}_i \hat{r}_j\right| 2\right\rangle \epsilon_i n_j=\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle \sum_{q=-2}^2\left(\begin{array}{ccc}j & 2 & j^{\prime} \\ -m & q & m^{\prime}\end{array}\right) c_{i j}^{(q)} \epsilon_i n_j$

keeping only tensorial parts of the equations we drop $\epsilon_i$ and $n_j$.

$\left\langle 1\left|\hat{r}_i\right| 2\right\rangle =\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle \sum_{q=-1}^1\left(\begin{array}{ccc}j & 1 & j^{\prime} \\ -m_j & q & m_j^{\prime}\end{array}\right) c_i^{(q)} $

$\left\langle 1\left|\hat{r}_i \hat{r}_j\right| 2\right\rangle=\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle \sum_{q=-2}^2\left(\begin{array}{ccc}j & 2 & j^{\prime} \\ -m & q & m^{\prime}\end{array}\right) c_{i j}^{(q)}$

Implementing this in code one can use it as

>`atomphys.calc.matrix_element.dipole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)`
> 
>Calculates the dipole matrix element
> 
>**Parameters**
>- `A`: Einstein coefficient [1/s].
>- `k`: Wavenumber of the transition [1/m].
>- `J_i`: Angular momentum of the lower state.
>- `J_f`: Angular momentum of the upper state.
>- `mJ_i`: Zeeman sublevel of the lower state.
>- `mJ_f`: Zeeman sublevel of the upper state.
>- `_ureg`: Unit registry.
>
>**Returns:** The dipole matrix element [$a_0$].

>`atomphys.calc.matrix_element.quadrupole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)`
> 
>Calculates the quadrupole matrix element
> 
>**Parameters**
>- `A`: Einstein coefficient [1/s].
>- `k`: Wavenumber of the transition [1/m].
>- `J_i`: Angular momentum of the lower state.
>- `J_f`: Angular momentum of the upper state.
>- `mJ_i`: Zeeman sublevel of the lower state.
>- `mJ_f`: Zeeman sublevel of the upper state.
>- `_ureg`: Unit registry.
>
>**Returns:** The quadrupole matrix element [$a_0^2$].


#### Electric Matrix Elements

Electric matrix elements are just the conventional dipole/quadrupole elements scaled by a funamental charge, $e$. 

$\left\langle 1\left|e \hat{r}_i\right| 2\right\rangle =e\left\langle 1\left|\hat{r}_i\right| 2\right\rangle$

$\left\langle 1\left|e \hat{r}_i \hat{r}_j\right| 2\right\rangle=e\left\langle 1\left|\hat{r}_i \hat{r}_j\right| 2\right\rangle$

>`atomphys.calc.matrix_element.electric_dipole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)`
> 
>Calculates the electric dipole matrix element
> 
>**Parameters**
>- `A`: Einstein coefficient [1/s].
>- `k`: Wavenumber of the transition [1/m].
>- `J_i`: Angular momentum of the lower state.
>- `J_f`: Angular momentum of the upper state.
>- `mJ_i`: Zeeman sublevel of the lower state.
>- `mJ_f`: Zeeman sublevel of the upper state.
>- `_ureg`: Unit registry.
>
>**Returns:** The electric dipole matrix element [$ea_0$].

>`atomphys.calc.matrix_element.electric_quadrupole_matrix_element(A, k, J_i, J_f, mJ_i, mJ_f, _ureg)`
> 
>Calculates the electric quadrupole matrix element
> 
>**Parameters**
>- `A`: Einstein coefficient [1/s].
>- `k`: Wavenumber of the transition [1/m].
>- `J_i`: Angular momentum of the lower state.
>- `J_f`: Angular momentum of the upper state.
>- `mJ_i`: Zeeman sublevel of the lower state.
>- `mJ_f`: Zeeman sublevel of the upper state.
>- `_ureg`: Unit registry.
>
>**Returns:** The electric quadrupole matrix element [$ea_0^2$].

### Rabi Frequencies:

Rabi Frequency is the fundamental quantity regarding interaction of atom with light. It tells us about how much coupling do we obtain between two eigenstates of our atom, given the presence of the oscilating electric field.

I decided in this documentation to re-derive the result to generalize it to arbitrary electric field and not just running wave. Analogous, but more specific derivations can be found in [Lindenfelser, Beck].

Let's consider Hamiltonian governing an electron wavefunction in an atom in the presence of electric field that oscilates with a single frequency $\omega$, described by a real vector field $\mathbf{A}_l(\mathbf{x}, t)$.

$\mathcal{H}=\mathcal{H}_0+\mathcal{H}_L$

This can be split into the bare hamiltonian, $\mathcal{H}_0$, and laser hamiltonian, $\mathcal{H}_L$, defined as below:

$\mathcal{H}_0=\sum_i \mathcal{E}_i|i\rangle\langle i|$

$\mathcal{H}_L=-\sum_l e \dot{\mathbf{x}}.\mathbf{A}_l=-\sum_{i, j, l}|i\rangle\left\langle i\left|e \dot{\mathbf{x}}.\mathbf{A}_l\right| j\right\rangle\langle j|$

, where $\mathcal{E}_i$ is eigenenergy of unperturbed eigenstate $|i\rangle$, $\mathbf{A}_l$ is a vector potential of the perturbing electric field. 

Given that $\mathbf{A}_l$ is a real vector field oscilating at single frequency $\omega$,  $\tilde{A}_l(\mathbf{x}, -\omega)$ =  $\tilde{A}_l^\dag(\mathbf{x}, \omega)$ and it can be re-written as [Oxford Optics notes (1.21)]: $\mathbf{A_l}=\frac{1}{2}\left(\mathbf{A}_{l}(\mathbf{x})e^{-i\omega_lt}+\mathbf{A}_{l}^{\dag}(\mathbf{x})e^{i\omega_lt}\right)$

Inserting $\mathbf{A}_l$ to the expression for the hamiltonian we obtain

$\mathcal{H}_L=-\sum_{i, j, l}|i\rangle \frac{1}{2}\left(\left\langle i\left|e \dot{\mathbf{x}}.\mathbf{A}_{l}(\mathbf{x})\right| j\right\rangle e^{-i \omega_l t}+\left\langle i\left|e \dot{\mathbf{x}}.\mathbf{A}_{l}^{\dag}(\mathbf{x}) \right| j\right\rangle e^{i \omega_l t}\right)\langle j|$

One can group terms in the bra-ket as Rabi Frequency $\Omega$ defined below:

$\Omega_{i j}= \begin{cases}\left\langle i\left|e \dot{\mathbf{x}}.\mathbf{A}_{l}(\mathbf{x}) / \hbar\right| j\right\rangle & (j \geq i) \\ \Omega_{j i}^{\dagger} & (j<i)\end{cases}$

Considering the motion of the electron relative to the nucleus and applying the Born-Oppenheimer approximation, the expression simplifies to:

$\Omega_{i j}= \begin{cases}\left\langle i\left|e \dot{\mathbf{r}}.\mathbf{A}_{l}(\mathbf{R}+\mathbf{r}) / \hbar\right| j\right\rangle & (j \geq i) \\ \Omega_{j i}^{\dagger} & (j<i)\end{cases}$

Utilizing the commutation relation $i \hbar \dot{\mathbf{r}}=\left[\mathbf{r}, \mathcal{H}_0\right]$, $\mathbf{E}=-\nabla \phi-\frac{\partial \mathbf{A}}{\partial t}$, and the radiation gauge in which $\nabla\phi=0$ we arrive at the final expression for the Rabi frequency:

$\Omega_{i j}= \begin{cases}-\left\langle i\left|e \mathbf{r}. \mathbf{E}_{l}(\mathbf{R}+\mathbf{r}) / \hbar\right| j\right\rangle & (j \geq i) \\ \Omega_{j i}^{\dagger} & (j<i)\end{cases}$

We can expand $\mathbf{r}.\mathbf{E}_l(\mathbf{R}+\mathbf{r})$ as $r_\mu E^\mu(R)+\frac{1}{2}\frac{\partial{E}}{}$

Taylor expanding $\mathbf{E}_l(\mathbf{R}+\mathbf{r})$ we get

$-\left\langle i\left|e \mathbf{r}. \mathbf{E}_{l}(\mathbf{R}) / \hbar\right| j\right\rangle +
-\left\langle i\left|e \mathbf{r}. \frac{\partial\mathbf{E}_{l}(\mathbf{R})}{\partial r_j}r_j(\mathbf{R}) / \hbar\right| j\right\rangle$

$\hbar\Omega_{i j}= \begin{cases} -\left\langle i\left|\mathbf{r} \right| j\right\rangle . \mathbf{E}_l & (\mathrm{DP}) \\ -i \frac{e E_{0, l}}{2 \hbar} \frac{2 \pi}{\lambda_l} e^{-i \mathbf{k}_l \mathbf{R}-i \varphi_l}\left\langle i\left|\left(\mathbf{r} \mathbf{n}_l\right)\left(\mathbf{r} \mathbf{\epsilon}_l\right)\right| j\right\rangle & (\mathrm{QP})\end{cases}$

> return (np.dot(E_field, d) * _ureg("e/hbar")).to("MHz")
>
> 
> qme = quadrupole_matrix_element(A=A, k=k, J_i=J_i, J_f=J_f, mJ_i=mJ_i, mJ_f=mJ_f, _ureg=_ureg)
>
>return (1 / 2 * np.sum(E_gradient * qme) * _ureg("e/hbar")).to("MHz")





$\Omega_{i j}= \begin{cases}-ie\omega_l \left\langle i\left|\mathbf{r} \right| j\right\rangle . \nabla{\mathbf{A}}& (\mathrm{DP}) \\ -i \frac{e E_{0, l}}{2 \hbar} \frac{2 \pi}{\lambda_l} e^{-i \mathbf{k}_l \mathbf{R}-i \varphi_l}\left\langle i\left|\left(\mathbf{r} \mathbf{n}_l\right)\left(\mathbf{r} \mathbf{\epsilon}_l\right)\right| j\right\rangle & (\mathrm{QP})\end{cases}$


## AC Stark Shifts:

In order to introduce the theory behind the AC Stark Shifts, I will closely follow Beloy's Thesis [Beloy 2009]. I intend to present the calculation of AC Stark Shift in its more primitive form, without expanding the formula into most basic units. This allows one to directly understand the simplicity of the AC Stark Shift calculations, not make any assumptions about the structure of the Electric Field or allowed transitions. Normal derivations of AC Stark Shift only consider dipole transitions driven by running wave field. 

Expression of AC Stark Shift in terms of Rabi Frequencies is possible, as we developed good apparatus to calculate Rabi Frequencies. 

Consider oscilating electric field 

$\mathbf{A_l}=\frac{1}{2}\left(\mathbf{A}_{l}(\mathbf{x})e^{-i\omega_lt}+\mathbf{A}_{l}^{\dag}(\mathbf{x})e^{i\omega_lt}\right)$

Such field interacts with an atom perturbing the hamiltonian by additional term(see section Rabi Frequencies):

$\mathcal{H}_L=-\sum_{i, j, l}|i\rangle \frac{1}{2}\left(\Omega_{ij}^{(l)} e^{-i \omega_l t}+ \Omega_{ij}^{\dag(l)} e^{i \omega_l t}\right)\langle j|$

It is noted that $\Omega_{ij}^{(l)}$ is a time independent operator. For small field amplitudes (condition), one can treat the $\mathcal{H}_L$ as small perturbation to $\mathcal{H}_0$. Using Floquet perturbation theory [Beloy 2009, Appendix D/Chapter 3], one can develop expression for the mean energy shifts $\delta \mathcal{E} = \delta \mathcal{E}^{(1)} + \delta \mathcal{E}^{(2)} + ...$ . Due to symmetry  $\delta \mathcal{E}^{(1)}$ is not present, and so we will start from evaluating $\delta \mathcal{E}^{(2)}$.

### The second order AC Stark Shift

Using derivation in [Beloy 6], we know that:

$\delta \mathcal{E}_i^{(2)}=\sum_{j} \frac{\Omega_{ij}^{(l)}\Omega_{ij}^{\dag(l)}}{\mathcal{E}_i-\mathcal{E}_{j}+\omega_l}+\sum_{j} \frac{\Omega_{ij}^{(l)}\Omega_{ij}^{\dag (l)}}{\mathcal{E}_i-\mathcal{E}_{j}-\omega_l}$

That's it. You would think that there is much more complicated maths in order to calculate AC Stark Shift, but acctually there isn't. All the long derivations presented in most thesis is to bring it to a simple form, which can only be done when according assumptions have been made. Because of numerical comfort we can keep above expression as final.

Thanks to such expression we didn't need to make any assumptions about the nature of the electric field, and about the nature of the electronic transitions. 

Taking second order AC Stark Shift as the largest contribution to overall stark shift we say:

$\delta \mathcal{E}_i \approx \delta \mathcal{E}_i^{(2)}$ 

>`atomphys.calc.ac_stark.ac_stark_shift(state, mJ, El_field, wavelengths, B, _ureg)`
> 
>Calculates the second order AC Stark Shift due to the Electric Dipole and Electric Quadrupole transitions present.
> 
>**Parameters**
>- `state`: Electronic state of an atom
>- `mJ`: Zeeman sublevel of a state of interest
>- `El_field`: Electric Field oscilating at frequency $\omega_l$
>- `wavelengths`: (optional) [nm] If an array of wavelengths is passed, then AC stark shift is calculated for all passed wavelengths simultaneously. If not, then the frequency of the electric field is taken. 
>- `B`: Additionally one can include the magnetic field that shifts the agnetic field sublevels
>- `_ureg`: Unit registry
>
>**Returns:** Energy Shift [mK $k_b$].

### Polarizability

Polarizability of an atom in state $|i \rangle$ , $\alpha_i$, is an AC Stark Shift, $\delta \mathcal{E}_i$  normalized by the electric field causing this stark shift, $E_0$.

$\alpha_i = \frac{-4 \left(\delta \mathcal{E}_i\right)}{E_0^2}$

>`atomphys.calc.ac_stark.polarizability(state, mJ, El_field, wavelengths, B, _ureg)`
> 
>Calculates the second order AC Stark Shift due to the Electric Dipole and Electric Quadrupole transitions present.
> 
>**Parameters**
>- `state`: Electronic state of an atom
>- `mJ`: Zeeman sublevel of a state of interest
>- `El_field`: Electric Field oscilating at frequency $\omega_l$
>- `wavelengths`: (optional) [nm] If an array of wavelengths is passed, then AC stark shift is calculated for all passed wavelengths simultaneously. If not, then the frequency of the electric field is taken. 
>- `B`: Additionally one can include the magnetic field that shifts the agnetic field sublevels
>- `_ureg`: Unit registry
>
>**Returns:** polarizability of the state [$e^2a_0^2E_h^{-1}$].

## References

James 1998: "Quantum dynamics of cold trapped ions, with application to quantum computation." [https://arxiv.org/abs/quant-ph/9702053](https://arxiv.org/abs/quant-ph/9702053)

[Quantum Optics Scully, Zubairy Page 206.]

[Jackson Chapter 5]

[Beloy 2009]

[Lindenfelser]

[Beck]

## Contributions:

This package is based on the package written by Matt Grau 'atomphys'. In order to change few data-structures and simplify overall package structure we (W.Adamczyk and C.Mordini) decided to rebuild it. This allowed us to add qutip helper functions, integrate other databases, and add visualising tools. States, Atoms and Transitions are now also not using registries, but graphs instead. P.Leindecker contributed valuable advice of how to build the package such that its API is well integrateable to web-development.

Wojciech Adamczyk <wadamczyk@phys.ethz.ch>,
Dr Carmelo Mordini <cmordini@phys.ehtz.ch>, Philipp Leindecker <pleindecker@phys.ethz.ch>,
Prof Matt Grau <mgrau@odu.edu>, Prof Jonathan Home <jhome@phys.ethz.ch>

## Comments:

> - include proper definition of wigner-3j symbol
> - understant [j] and include it in thedocumentation
> - understand whether replacing $\mathcal{E}^{(0,2)}$ with $\mathcal{E}^{(2)}$ is correct
> - understand assumptions of my calculation of the AC Stark Shift (Hyperfine Structure and so on...)
> - Check all absolute values, magnitudes and real and imaginary part of rabi frequencies and AC stark shifts
