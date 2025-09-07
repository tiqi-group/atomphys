# Atomphys Documentation: 

This document provides a guide to the Python module designed for atomic physics calculations. The idea behind this module is to explore atomic physics, properties of atoms, and understand their electron dynamics in light fields. 

The package should be primarily used for bosonic isotopes, where the hyperfine structure is not present. We decided to initially release the package with simplified structure such that to get a feeling for the usefullness of the tool to the community. In case of interest, we can extend this package to include calculations for the fermionic isotopes.

I hope this will be helpfull for many.

The documentation is divided into three subsections, which each discuss different potential use-cases of atomphys:
- (1) Core - Here we discuss how atomphys can be used as python API for NIST database and other databases
- (2) Calculations - Here we discuss the theory behind the calculations with which atomphys can help
- (3) Helper functions for qutip - Here we discuss how atomphys can be usef to build light-atom hamiltonians that can be directly parsed to qutip and can be solved

## Getting Started

### Installation

Install the necessary libraries via pip:

```sh
pip install git+https://github.com/tiqi-group/atomphys.git
```
## Static Structure
### Atom (Graph):

At the center of the atomphys package there is an atom. We decided to define atom as a graph data structure. This is because it is natural to think about an atom as a graph, where the states are the vertices, and the transitions are the edges of the graph. 

It allows us then to build upon great work that has been done on graphs in order to plot atoms, explore them, list all the properties, and their vertices, and edges. 


#### Loading of an atom:

But it would be extremaly tideous to build up such structure from scratch. Instead one can load an atom directly from NIST database, or from a custom json database. 

To load an atom from NIST, you can call from_nist() function. As a name you need to parse the element name of an atom you are interested. For instance if you are interested in Calcium atom/ion you need to call from_nist('Ca')/from_nist('Ca+').

> `atomphys.from_nist(name, _ureg)`
>
> Returns an atom with states and transitions found in NIST database. If NIST is not complete for your purposes, the missing transitions can be added manually. 

Alternatively one can use:

> `atomphys.from_json(name, _ureg)`
>
> Returns an atom with states and transitions found in the custom database. Database has to be a json file with a format given below

```
{"name": "Ca", 
"states": [{"configuration": "3p6.4s2", "term": "1S0", "energy": "0.0 Ry"}, {"configuration": "3p6.4s.4p", "term": "3P0", "energy": "0.1381289573728982 Ry"}, ... ], 
"transitions": [{"A": "2.74e+03 s^-1", "state_i": {"energy": "0.0 Ry", "term": "1S0"}, "state_f": {"energy": "0.138604292491823", "term": "3P1"}}, {"A": "4.89e+05 s^-1", "state_i": {"energy": "0.1381289573728982 Ry", "term": "3P0"}, "state_f": {"energy": "0.1853094352973106", "term": "3D1"}}...]
}
```
Above example 'Ca' database was taken from {Mills 2018}.

#### Basic Usage:
Upon importing your atom, it allows you to interact with it in a similar way you would interact with a databse / or a graph. It makes it easy to interact in an easy way with NIST database, or any other database that you would use. For the functionalities of the atom one should refer to the API reference, which should list all possible actions that can be performed. 

### States (Node) and Transitions (Edge):
Atom is formed by collection of states (nodes), connected via transitions (edges). Each State and Transition is a separate object with properties and possible functions associated with it. The representation of the atom as a graph could be used for the optimal repumping schemes, or optimal control schemes. 


## Calculations

### Matrix Elements:

Before diving into the specific formulas and derivations used, it's worth pausing to consider the key concepts in atomic physics that we focus on and their significance. Atomic physics heavily relies on the calculation of matrix elements, denoted as $M_{if} = \left\langle f \mid \hat{M} \mid i \right\rangle$. This calculation involves the projection of an operator $\hat{M}$, which acts on an initial state $\left| i \right>$, onto a final state $\left| f \right>$. The value of $M_{if}$ is crucial because it measures the extent of interaction, or coupling, between the states $\left| i \right>$ and $\left| f \right>$ via the operator $\hat{M}$, when $\hat{M}$ is part of the Hamiltonian $\hat{H}$. Understanding this interaction provides insights into the dynamics of the overall state $\left| \psi(t) \right>$.

However, computing these matrix elements is challenging due to their basis involved. Is there a way to simplify (break-down) this process using symmetry principles? Fortunately, the answer is yes, thanks to the Wigner-Eckart Theorem.

The Wigner-Eckart Theorem states that the matrix elements of spherical tensor operators, within the framework of angular momentum eigenstates, can be broken down into two components: one that is independent of the angular momentum orientation (so called reduced matrix element) and another that is a Clebsch-Gordan coefficient.(https://en.wikipedia.org/wiki/Wigner%E2%80%93Eckart_theorem) 

$\left\langle j m\left|T_q^{(k)}\right| j^{\prime} m^{\prime}\right\rangle=\left\langle j^{\prime} m^{\prime} k q \mid j m\right\rangle\left\langle j\left\|T^{(k)}\right\| j^{\prime}\right\rangle = [j] \left(\begin{array}{ccc}j^{\prime} & k & j^* \\ m^{\prime} & q & m\end{array}\right)  \left\langle j\left\|T^{(k)}\right\| j^{\prime}\right\rangle$

Here, $[j]$ represents the dimension of the representation space of $j$, and $j^*$ indicates a complex conjugate representation to $j$. The parentheses contain the Wigner-3j symbol, and $\left\langle \mid\mid . \mid\mid \right\rangle$ denotes a reduced matrix element.

As the reduced matrix element, doesn't depend on basis this is the most basic block from which we will start the calculation.


#### Reduced Matrix Elements (Dipole and Quadrupole):
Following James (equation 5.9) and (equation 5.10) [1]

$A_{12}^{(E 1)} \equiv \sum_{m=-j}^j \bar{A}_{12}^{(E 1)}=\frac{4 c \alpha k_{12}^3}{3\left(2 j^{\prime}+1\right)}\left|\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle\right|^2$

$A_{12}^{(E 2)} \equiv \sum_{m=-j}^j \bar{A}_{12}^{(E 2)}=\frac{c \alpha k_{12}^5}{15\left(2 j^{\prime}+1\right)}\left|\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle\right|^2$

we can rewrite 

$\left|\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle\right| = \left( \frac{3\left(2 j^{\prime}+1\right)}{4 c \alpha k_{12}^3} A_{12}^{(E 1)} \right)^{\frac{1}{2}}$

$\left|\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle\right| = \left(\frac{15\left(2 j^{\prime}+1\right)}{c \alpha k_{12}^5} A_{12}^{(E 2)}\right)^{\frac{1}{2}}$

Equivalent result can be derived from Wigner-Weisskopf Decay. Example of such derivation can be found in Scully's Quantum Optics[2]



#### Matrix Elements

Now in order to transform it to actual quantities that we want we can use Wigner-Eckart theorem.

$\left\langle 1\left|\hat{r}_i\right| 2\right\rangle \epsilon_i=\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle \sum_{q=-1}^1\left(\begin{array}{ccc}j & 1 & j^{\prime} \\ -m_j & q & m_j^{\prime}\end{array}\right) c_i^{(q)} \epsilon_i$

$\left\langle 1\left|\hat{r}_i \hat{r}_j\right| 2\right\rangle \epsilon_i n_j=\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle \sum_{q=-2}^2\left(\begin{array}{ccc}j & 2 & j^{\prime} \\ -m & q & m^{\prime}\end{array}\right) c_{i j}^{(q)} \epsilon_i n_j$

keeping only tensorial parts of the equations we drop $\epsilon_i$ and $n_j$.

$\left\langle 1\left|\hat{r}_i\right| 2\right\rangle =\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle \sum_{q=-1}^1\left(\begin{array}{ccc}j & 1 & j^{\prime} \\ -m_j & q & m_j^{\prime}\end{array}\right) c_i^{(q)} $

$\left\langle 1\left|\hat{r}_i \hat{r}_j\right| 2\right\rangle=\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle \sum_{q=-2}^2\left(\begin{array}{ccc}j & 2 & j^{\prime} \\ -m & q & m^{\prime}\end{array}\right) c_{i j}^{(q)}$


#### Electric Matrix Elements

Electric matrix elements are just the conventional dipole/quadrupole elements scaled by a funamental charge, $e$. 


### Rabi Frequencies:

Rabi Frequency is the fundamental quantity regarding interaction of atom with light. It tells us about how much coupling do we obtain between two eigenstates of our atom, given the presence of the oscilating electric field.

Given the complexities and the confusion surrounding the derivations of dipole and quadrupole Rabi frequencies, a re-derivation is presented here. This derivation draws from S Mavadia's thesis, Weissbluth book. Starting point is minimum-coupling-hamiltonian, for even more comprehensive form of the light-matter interaction refer to Weissbluth. This derivation operates within the Coulomb gauge and assumes a electromagnetic vacuum (no charged particles flying around). 

#### Minimum Coupling Hamiltonian

Consider an arbitrary electromagnetic vector field, $\mathbf{A}(\mathbf{r},t)$, within the Coulomb gauge ($\nabla \cdot \mathbf{A}(\mathbf{r}) = 0$) in vacuum ($\phi(\mathbf{r}, t)=0$).

The minimum-coupling Hamiltonian in Coulomb's gauge is expressed as:

$$ 
H=\sum_\alpha\left(\frac{\left[\mathbf{p}^{(\alpha)}-q^{(\alpha)} \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)\right]^2}{2 m^{(\alpha)}}\right)+H_F+V_{\text {Coul }}  +\frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \boldsymbol{\nabla} \times \mathbf{A}
$$

where $H_F$ represents the free field energy ($H_F=\frac{1}{2} \int \mathrm{d}^3 r\left(\epsilon_0 \mathbf{E}^2(\mathbf{r}, t)+\frac{1}{\mu_0} \mathbf{B}^2(\mathbf{r}, t)\right)$) and $V_{\text {Coul}}$ encompasses terms defining the atomic state, including Coulomb interactions and spin-orbit coupling.

Focusing on electron dynamics rather than absolute energy levels allows us to simplify our Hamiltonian. We can neglect constant energy terms $H_F$ and $\varepsilon_{Coul}^\alpha$, keeping only the terms that depend on $\mathbf{x}_\alpha$ or $\mathbf{p}_\alpha$. Coulomb's Gauge $\nabla \cdot \mathbf{A}(\mathbf{r}) = 0$ ensures that $\mathbf{p}_\alpha$ and $\mathbf{A}(\mathbf{x}_\alpha, t)$  commute and so we can rewrite ${H}$ as:

$$
\begin{aligned}
H &=\sum_\alpha\left(\frac{\left[\mathbf{p}^{(\alpha)}-q^{(\alpha)} \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)\right]^2}{2 m^{(\alpha)}}\right)+V_{\text {Coul }} + \frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \boldsymbol{\nabla} \times \mathbf{A}
\\&=\sum_\alpha \frac{\mathbf{p^{(\alpha)}}^2}{2 m^{(\alpha)}}+V_{\text {Coul }}+\sum_\alpha\left(-\frac{q^{(\alpha)} \mathbf{p}^{(\alpha)} \cdot \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)}{m^{(\alpha)}}+\frac{{q^{(\alpha)}}^2\mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)^2}{2 m^{(\alpha)}}\right) + 
\frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \boldsymbol{\nabla} \times \mathbf{A}
\end{aligned}
$$

Neglecting $\frac{q_\alpha^2\mathbf{A}\left(\mathbf{x}_\alpha, t\right)^2}{2 m_\alpha}$ and $\frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \boldsymbol{\nabla} \times \mathbf{A}$, and grouping terms together we get:



$$ 
\begin{aligned}
&H = H_0 + H_I \\
&H_0 = \sum_\alpha \frac{\mathbf{p^{(\alpha)}}^2}{2 m^{(\alpha)}}+V_{\text {Coul }} = \sum_i \mathcal{E}_i|i\rangle\langle i| \\
&H_I = \sum_\alpha-q^{(\alpha)} \mathbf{\dot{x}}^{(\alpha)} \cdot \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)
\end{aligned}
$$


$\mathbf{x} = \mathbf{R} + \mathbf{r}$, where $\mathbf{R}$ is a postion of the nucleus and $\mathbf{r}$ is position of an electron relative to the nucleus. Using Born-Oppenheimer approximation we can write $\mathbf{\dot{x}} = \mathbf{\dot{r}}$ neglecting $\mathbf{\dot{R}}$.


$$
H_I = \sum_\alpha-q_\alpha \mathbf{\dot{r}}_\alpha \cdot \mathbf{A}\left(\mathbf{R}+\mathbf{r}_\alpha, t\right)
$$
#### Multipole expansion
Taylor expanding $\mathbf{A}(\mathbf{R}+\mathbf{r}_\alpha, t)$, we get:

$$
H_I = \sum_\alpha-q_\alpha \dot{r}^{(\alpha)}_{\mu} \left(A^\mu\left(\mathbf{R}, t\right) + \partial^\nu A^\mu\left(\mathbf{R}, t\right)r^{(\alpha)}_{\nu} \right)
$$

From now on, Lets define $A^\mu = A^\mu(\mathbf{R}, t)$ 

$$
\begin{aligned}
H_I &= \sum_{\alpha, i, j} |i\rangle\left\langle i\left|e \dot{r}^{(\alpha)}_{\mu} \left(A^\mu + \partial^\nu A^\mu r^{(\alpha)}_{\nu} + ...\right) \right| j\right\rangle\langle j| \\
&= \sum_{\alpha, i, j} \left( |i\rangle\left\langle i\left|e \dot{r}^{(\alpha)}_{\mu} A^\mu\right| j\right\rangle\langle j| + |i\rangle\left\langle i\left|e \dot{r}^{(\alpha)}_{\mu} \partial^\nu A^\mu r^{(\alpha)}_{\nu} \right| j\right\rangle\langle j| + ... \right) \\
&= \sum_{\alpha, i, j} \left( e A^\mu|i\rangle\left\langle i\left| \dot{r}^{(\alpha)}_{\mu} \right| j\right\rangle\langle j| + e\partial^\nu A^\mu |i\rangle\left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} \right| j\right\rangle\langle j| + ... \right)
\end{aligned}
$$

As in the end we would like to see how the light field interacts with consecutive terms of the multipole expansion formed by the atom, we need to get rid of $\dot{r}_\mu^{(\alpha)}$.

As $\left[ r_\mu, p^2 \right] = 2i \hbar p_\mu$, then $i \hbar \dot{r}_\mu=\left[r_\mu, H_0\right]$

Lets then solve consecutive terms of the Taylor expansion

##### 0th Order term:
$$
\begin{aligned}
\left\langle i\left| \dot{r}^{(\alpha)}_{\mu} \right| j\right\rangle &=
-\frac{i}{\hbar} \left\langle i\left| [r^{(\alpha)}_{\mu}, H_0] \right| j\right\rangle  
\\&= -\frac{i}{\hbar} \left\langle i\left|r^{(\alpha)}_{\mu} H_0 - H_0 r^{(\alpha)}_{\mu} \right| j\right\rangle 
\\&= -i \left\langle i\left|r^{(\alpha)}_{\mu} \omega_j - \omega_i r^{(\alpha)}_{\mu} \right| j\right\rangle 
\\&= i \omega_{0} \left\langle i\left| r^{(\alpha)}_{\mu} \right| j\right\rangle
\end{aligned}
$$

, where $\omega_{0} = \omega_i-\omega_j$

##### 1st Order term:

$$
\begin{aligned}
\left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} \right| j\right\rangle 
&= -\frac{i}{\hbar}\left\langle i\left| [r^{(\alpha)}_{\mu}, H_0]  r^{(\alpha)}_{\nu} \right| j\right\rangle 
\\&=-\frac{i}{\hbar}\left\langle i\left| r^{(\alpha)}_{\mu}H_0r_\nu^{(\alpha)} - H_0r^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} \right| j\right\rangle
\end{aligned}
$$

This is more tricky, because now we need to to commute $H_0$ with $r_\nu^{(\alpha)}$ which as a result would give us $\dot{r}^{(\alpha)}_{\mu}$ again. Instead what we can do we can split the problem into symmetric and anti-symmetric part hoping it will get easier. 

$$
\left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} \right| j\right\rangle =
\frac{1}{2}\left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} +   r^{(\alpha)}_{\nu} \dot{r}^{(\alpha)}_{\mu} \right| j\right\rangle +
\frac{1}{2}\left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} -   r^{(\alpha)}_{\nu} \dot{r}^{(\alpha)}_{\mu} \right| j\right\rangle
$$

Lets solve the symmetric and antisymmetric part separately:

Symmetric part:

$$
\begin{aligned}
\frac{1}{2}\left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} +   r^{(\alpha)}_{\nu} \dot{r}^{(\alpha)}_{\mu} \right| j\right\rangle &= 
\frac{-i}{2\hbar}\left\langle i\left| [r^{(\alpha)}_{\mu}, H_0]  r^{(\alpha)}_{\nu} +   r^{(\alpha)}_{\nu} [r^{(\alpha)}_{\mu}, H_0] \right| j\right\rangle \\&= 
\frac{-i}{2\hbar}\left\langle i\left| -H_0r^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} +   r^{(\alpha)}_{\nu} r^{(\alpha)}_{\mu}H_0 \right| j\right\rangle \\&= 
\frac{1}{2}i\omega_{0}\left\langle i\left| r^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu}  \right| j\right\rangle
\end{aligned}
$$



Anti-symmetric part:

$$
\begin{aligned} \frac{1}{2} \partial^\nu A^\mu \left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} -   r^{(\alpha)}_{\nu} \dot{r}^{(\alpha)}_{\mu} \right| j\right\rangle &= \frac{1}{2} \varepsilon^{i \mu \nu}\varepsilon_{i}^{j k} \partial_j A_k \left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} \right| j\right\rangle \\&= 
\frac{1}{2} \varepsilon_{i}^{j k} \partial_j A_k \left\langle i\left| \varepsilon^{i \mu \nu}\dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} \right| j\right\rangle \\&= 
\frac{1}{2} \hbar \varepsilon_{i}^{j k} \partial_j A_k \left\langle i\left| L^i \right| j\right\rangle
\end{aligned}
$$

Therefore one can re-write 1st Order term as:

$$
\begin{aligned}
e \partial^\nu A^\mu\left\langle i\left| \dot{r}^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu} \right| j\right\rangle =
\frac{1}{2} ie\omega_0 \partial^\nu A^\mu \left\langle i\left| r^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu}  \right| j\right\rangle +
\frac{1}{2} \hbar e \varepsilon_{i}^{j k} \partial_j A_k \left\langle i\left| L^i \right| j\right\rangle
\end{aligned}
$$

, where the first term corresponds to the electric quadrupole coupling and the second term corresponds to magnetic dipole coupling

Collecting all the terms up to the 1st Order of Taylor expansion of $A_\mu$, we get:

$$

H_I =  \sum_{\alpha, i, j} |i\rangle\left( \underbrace{i e \omega_{0} A^\mu\left\langle i\left| r^{(\alpha)}_{\mu} \right| j\right\rangle }_\text{Electric Dipole} + 
\underbrace{\frac{1}{2} ie\omega_0 \partial^\nu A^\mu \left\langle i\left| r^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu}  \right| j\right\rangle}_\text{Electric Quadrupole} +
\underbrace{\frac{1}{2} \hbar e \varepsilon_{i}^{j k} \partial_j A_k \left\langle i\left| L^i \right| j\right\rangle}_\text{Magnetic Dipole} + ... \right) \langle j|

$$

In atomphys for now we don't support calculation of the magnetic dipole matrix elements, and so we will therefore only focus on the first two terms, assuming the light field interating with our atom doesn't have intrinsically large magnetic dipole moment. 



####  Constraining A-vector field:

Let us constrain our choice of $\mathbf{A}$ vector-field. In the end what we are interested in is an electric field $\mathbf{E}$ oscilating at single frequency $\omega$. As the electric field is an observable it must be real, and so $\tilde{E}(\mathbf{x}, -\omega)$ =  $\tilde{E}^\dag(\mathbf{x}, \omega)$, and so it can be written as [Oxford Optics notes (1.21)]: $\mathbf{E}=\frac{1}{2}\left(\mathbf{E}(\mathbf{x})e^{-i\omega_lt}+\mathbf{E}^{\dag}(\mathbf{x})e^{i\omega_lt}\right)$. 

Working in vacuum in Coulombs gauge we can write $\mathbf{E} = -\frac{\partial \mathbf{A}}{\partial t}$, hence

$$
\mathbf{A} =\frac{1}{2}\left(\left(\frac{1}{i\omega_l}\mathbf{E}(\mathbf{x})\right)e^{-i\omega_lt} + \left(\frac{1}{i\omega_l}\mathbf{E}(\mathbf{x})\right)^{\dag}e^{i\omega_lt}\right) = \frac{1}{2}\left(\mathbf{A}(\mathbf{x})e^{-i\omega_lt} + \mathbf{A}^{\dag}(\mathbf{x})e^{i\omega_lt}\right)
$$ 

The interaction then can be written as:

$$
\begin{aligned}
H_I &= \sum_{\alpha, i, j} |i\rangle\left\langle i\left|e \mathbf{\dot{r}}^{(\alpha)} \mathbf{A} \right| j\right\rangle\langle j| \\

&= \sum_{\alpha, i, j} |i\rangle \frac{1}{2} \left( \left\langle i\left|e \mathbf{\dot{r}}^{(\alpha)} \mathbf{A}(\mathbf{x}) \right| j\right\rangle e^{-i\omega_lt}
+ \left\langle i\left|e \mathbf{\dot{r}}^{(\alpha)} \mathbf{A}^{\dag}(\mathbf{x}) \right| j\right\rangle e^{i\omega_lt} \right) \langle j| \\

&= \sum_{\alpha, i, j} |i\rangle \frac{\hbar}{2} \left(\Omega_{ij} e^{-i\omega_lt}
+ \Omega^\dag_{ij} e^{i\omega_lt} \right) \langle j| 
\end{aligned}
$$

, where Rabi Frequency is defined as follows:

$$
\hbar\Omega_{i j}= \left\langle i\left|e \mathbf{\dot{r}}^{(\alpha)} \mathbf{A}(\mathbf{x}) \right| j\right\rangle
$$

As we saw, we can decompose it through the taylor expansion into the consecutive terms corresponding to different nature of the transition. Usually only one of the coupling types is dominant and the other can be neglected. The dominant type depends on the nature of the transition and the electric field structure.

$\Omega_{i j}= \begin{cases}
\frac{e\omega_{0}}{\hbar\omega_l}E^\mu\left\langle i\left| r^{(\alpha)}_{\mu} \right| j\right\rangle & (\mathrm{E1})
\\  \frac{e\omega_0}{2\hbar\omega_l} \partial^\nu E^\mu \left\langle i\left| r^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu}  \right| j\right\rangle & (\mathrm{E2})\\
\frac{1}{2} e \varepsilon_{\theta}^{\beta \gamma} \partial_\beta A_\gamma \left\langle i\left| L^\theta \right| j\right\rangle & (\mathrm{M1})\\
\end{cases}$

#### Final Remarks:
This is the final expression of the rabi-frequencies. As far as we are interested in only electric multipole expansion we took all required terms from the dirac equation. We worked in vacuum and in Coulomb gauge. Not switching the gauge allowed us to not to make any mistakes that arise from working in multiple gauges. 

Other common derivation is using PZW Gauge, which naturally has a form of multipole expansion. I, however, prefered not to work in it, as from what I have seen it wasnt a popular choice of understanding the problem. Coulomb's gauge was a preffered choice, however for many they missed a step to split the first Taylor expansion term into symmetric and antisymmetric part, which forced them to fudge a factor of 1/2. 


## AC Stark Shifts:

In order to introduce the theory behind the AC Stark Shifts, I will closely follow Beloy's Thesis [Beloy 2009]. I intend to present the calculation of AC Stark Shift in its more primitive form, without expanding the formula into most basic units. This allows one to directly understand the simplicity of the AC Stark Shift calculations, not make any assumptions about the structure of the Electric Field or allowed transitions. Normal derivations of AC Stark Shift only consider dipole transitions driven by running wave field. 

Expression of AC Stark Shift in terms of Rabi Frequencies is possible, as we developed good apparatus to calculate Rabi Frequencies. 

Consider oscilating electric field 

$\mathbf{A_l}=\frac{1}{2}\left(\mathbf{A}_{l}(\mathbf{x})e^{-i\omega_lt}+\mathbf{A}_{l}^{\dag}(\mathbf{x})e^{i\omega_lt}\right)$

Such field interacts with an atom perturbing the hamiltonian by additional term(see section Rabi Frequencies):

$\mathcal{H}_L=-\sum_{i, j, l}|i\rangle \frac{1}{2}\left(\Omega_{ij}^{(l)} e^{-i \omega_l t}+ \Omega_{ij}^{\dag(l)} e^{i \omega_l t}\right)\langle j|$

It is noted that $\Omega_{ij}^{(l)}$ is a time independent operator. For small field amplitudes (condition), one can treat the $\mathcal{H}_L$ as small perturbation to $\mathcal{H}_0$. Using Floquet perturbation theory [Beloy 2009, Appendix D/Chapter 3], one can develop expression for the mean energy shifts $\delta \mathcal{E} = \delta \mathcal{E}^{(1)} + \delta \mathcal{E}^{(2)} + ...$ . Due to symmetry  $\delta \mathcal{E}^{(1)}$ is not present, and so we will start from evaluating $\delta \mathcal{E}^{(2)}$.

### Floquet Perturbation Theory:

$
\left[H(\xi, t)-i \frac{\partial}{\partial t}\right] \Psi(\xi, t)=0 \\
H(\xi, t+\tau)=H(\xi, t)\\
\Psi(\xi, t)=\phi(\xi, t) e^{-i \mathcal{E} t},
$ such that $\phi(\xi, t+\tau)=\phi(\xi, t) \\
\mathcal{H}(\xi, t) \equiv H(\xi, t)-i \frac{\partial}{\partial t}\\
\mathcal{H}(\xi, t) \phi(\xi, t)=\mathcal{E} \phi(\xi, t)
$

Lets start our perturbation theory from here

$
H(\xi, t)=H_0(\xi)+V(\xi, t)\\
\mathcal{H}_0(\xi, t)=\mathcal{H}(\xi, t)-V(\xi, t)=H_0(\xi)-i \frac{\partial}{\partial t}\\
\mathcal{H}_0(\xi, t) \phi_{n q}^{(0)}(\xi, t)=\mathcal{E}_{n q}^{(0)} \phi_{n q}^{(0)}(\xi, t)\\
\phi_{n q}^{(0)}(\xi, t)=f_n(\xi) e^{i q \omega t}$ , where $\mathcal{E}_{n q}^{(0)}=E_n+q \omega$ and $H_0(\xi) f_n(\xi)=E_n f_n(\xi)$

### The second order AC Stark Shift

Using derivation in [Beloy 6], we know that:

$\delta \mathcal{E}_i^{(2)}=\sum_{j} \frac{\Omega_{ij}^{(l)}\Omega_{ij}^{\dag(l)}}{\mathcal{E}_i-\mathcal{E}_{j}+\omega_l}+\sum_{j} \frac{\Omega_{ij}^{(l)}\Omega_{ij}^{\dag (l)}}{\mathcal{E}_i-\mathcal{E}_{j}-\omega_l}$

That's it. You would think that there is much more complicated maths in order to calculate AC Stark Shift, but acctually there isn't. All the long derivations presented in most thesis is to bring it to a simple form, which can only be done when according assumptions have been made. Because of numerical comfort we can keep above expression as final.

Thanks to such expression we didn't need to make any assumptions about the nature of the electric field, and about the nature of the electronic transitions. 

Taking second order AC Stark Shift as the largest contribution to overall stark shift we say:

$\delta \mathcal{E}_i \approx \delta \mathcal{E}_i^{(2)}$ 


### Polarizability

Polarizability of an atom in state $|i \rangle$ , $\alpha_i$, is an AC Stark Shift, $\delta \mathcal{E}_i$  normalized by the electric field causing this stark shift, $E_0$.

$\alpha_i = \frac{-4 \left(\delta \mathcal{E}_i\right)}{E_0^2}$


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
