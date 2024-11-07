
### Rabi Frequencies

Rabi Frequency is the fundamental quantity regarding interaction of an atom with light. It tells us about how much coupling do we get between two eigenstates of our atom in the presence of the oscilating electric field.

Given the complexities and the confusion surrounding the derivations of dipole and quadrupole Rabi frequencies, a re-derivation is presented here. In this derivation we stay in a single gauge throughout. This derivation draws from [Weissbluth book](https://doi.org/10.1016/B978-0-12-744450-5.50026-8).

#### Minimum Coupling Hamiltonian

Consider an arbitrary electromagnetic vector field, $\mathbf{A}(\mathbf{r},t)$, within the Coulomb gauge ($\nabla \cdot \mathbf{A}(\mathbf{r}) = 0$) in vacuum ($\phi(\mathbf{r}, t)=0$). Throughout the derivation we will assume that the field is small enough so that we can treat it as a perturbation.

The minimum-coupling Hamiltonian in Coulomb's gauge is expressed as:

$$ 
H=\sum_\alpha\left(\frac{\left[\mathbf{p}^{(\alpha)}-q^{(\alpha)} \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)\right]^2}{2 m^{(\alpha)}}\right)+H_F+V_{\text{Coul}}  +\frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \boldsymbol{\nabla} \times \mathbf{A}
$$

where $H_F$ represents the free field energy ($H_F=\frac{1}{2} \int \mathrm{d}^3 r\left(\epsilon_0 \mathbf{E}^2(\mathbf{r}, t)+\frac{1}{\mu_0} \mathbf{B}^2(\mathbf{r}, t)\right)$) and $V_{\text{Coul}}$ contains terms defining the atomic state, including Coulomb interactions and spin-orbit coupling.

Focusing on electron dynamics rather than absolute energy levels allows us to neglect constant energy terms $H_F$ and $\varepsilon_{\text{Coul}}^\alpha$, keeping only the terms that depend on $\mathbf{x}_\alpha$ or $\mathbf{p}_\alpha$. Coulomb's Gauge $\nabla \cdot \mathbf{A}(\mathbf{r}) = 0$ ensures that $\mathbf{p}_\alpha$ and $\mathbf{A}(\mathbf{x}_\alpha, t)$ commute, so we can rewrite ${H}$ as:

$$
\begin{aligned}
H &=\sum_\alpha\left(\frac{\left[\mathbf{p}^{(\alpha)}-q^{(\alpha)} \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)\right]^2}{2 m^{(\alpha)}}\right)+V_{\text {Coul }} + \frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \boldsymbol{\nabla} \times \mathbf{A}
\\&=\sum_\alpha \frac{\mathbf{p^{(\alpha)}}^2}{2 m^{(\alpha)}}+V_{\text {Coul }}+\sum_\alpha\left(-\frac{q^{(\alpha)} \mathbf{p}^{(\alpha)} \cdot \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)}{m^{(\alpha)}}+\frac{{q^{(\alpha)}}^2\mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)^2}{2 m^{(\alpha)}}\right) + 
\frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \boldsymbol{\nabla} \times \mathbf{A}
\end{aligned}
$$

Neglecting less dominant terms of the interaction Hamiltonian $\frac{q_\alpha^2\mathbf{A}\left(\mathbf{x}_\alpha, t\right)^2}{2 m_\alpha}$, and grouping them together we get:

$$ 
\begin{aligned}
&H = H_0 + H_I + H_Z \\
&H_0 = \sum_\alpha \frac{\mathbf{p^{(\alpha)}}^2}{2 m^{(\alpha)}}+V_{\text {Coul }} = \sum_i \mathcal{E}_i\left|i\right>\left< i\right| \\
&H_I = \sum_\alpha-q^{(\alpha)} \mathbf{\dot{x}}^{(\alpha)} \cdot \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)
\\
&H_Z = \frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \mathbf{B}
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
\underbrace{\frac{1}{2} \hbar e \varepsilon_{i}^{j k} \partial_j A_k \left\langle i\left| L^i \right| j\right\rangle}_\text{Magnetic Dipole} + ... \right) \left< j\right|
$$

In atomphys for now we don't support calculation of the magnetic dipole matrix elements, and so we will therefore only focus on the first two terms, assuming the light field interating with our atom doesn't have intrinsically large magnetic dipole moment. 



####  Constraining A-vector field:

Let's constrain our choice of $\mathbf{A}$ vector-field. In the end what we are interested in is an electric field $\mathbf{E}$ oscilating at single frequency $\omega$. As the electric field is an observable, it must be real, and so $\tilde{E}(\mathbf{x}, -\omega)$ =  $\tilde{E}^\dagger(\mathbf{x}, \omega)$, and so it can be written as: $\mathbf{E}=\frac{1}{2}\left(\mathbf{E}(\mathbf{x})e^{-i\omega_lt}+\mathbf{E}^{\dagger}(\mathbf{x})e^{i\omega_lt}\right)$. 

Working in vacuum in Coulombs gauge we can write $\mathbf{E} = -\frac{\partial \mathbf{A}}{\partial t}$, hence

$$
\mathbf{A} =\frac{1}{2}\left(\left(\frac{1}{i\omega_l}\mathbf{E}(\mathbf{x})\right)e^{-i\omega_lt} + \left(\frac{1}{i\omega_l}\mathbf{E}(\mathbf{x})\right)^{\dagger}e^{i\omega_lt}\right) = \frac{1}{2}\left(\mathbf{A}(\mathbf{x})e^{-i\omega_lt} + \mathbf{A}^{\dagger}(\mathbf{x})e^{i\omega_lt}\right)
$$ 

The interaction then can be written as:

$$
\begin{aligned}
H_I &= \sum_{\alpha, i, j} \left|i\right> \left< i\left|e \mathbf{\dot{r}}^{(\alpha)} \mathbf{A} \right| j\right> \left< j\right| \\
&= \sum_{\alpha, i, j} \left|i\right> \frac{1}{2} \left( \left< i\left|e \mathbf{\dot{r}}^{(\alpha)} \mathbf{A}(\mathbf{x}) \right| j\right> e^{-i\omega_lt}
+ \left< i\left|e \mathbf{\dot{r}}^{(\alpha)} \mathbf{A}^{\dagger}(\mathbf{x}) \right| j\right> e^{i\omega_lt} \right) \left< j \right| \\
&= \sum_{\alpha, i, j} \left|i\right> \frac{\hbar}{2} \left(\Omega_{ij} e^{-i\omega_lt}
+ \Omega^\dagger_{ij} e^{i\omega_lt} \right) \left< j\right| 
\end{aligned}
$$

, where Rabi Frequency is defined as follows:

$$
\hbar\Omega_{i j}= \left< i\left|e \mathbf{\dot{r}}^{(\alpha)} \mathbf{A}(\mathbf{x}) \right| j\right>
$$

As we saw, we can decompose it through the Taylor expansion into the consecutive terms corresponding to different nature of the transition. Usually only one of the coupling types is dominant and the other can be neglected. The dominant type depends on the nature of the transition and the electric field structure.

$\Omega_{i j}= \begin{cases}
\frac{e\omega_{0}}{\hbar\omega_l}E^\mu\left\langle i\left| r^{(\alpha)}_{\mu} \right| j\right\rangle & (\mathrm{E1})
\\  \frac{e\omega_0}{2\hbar\omega_l} \partial^\nu E^\mu \left\langle i\left| r^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu}  \right| j\right\rangle & (\mathrm{E2})\\
\frac{1}{2} e \varepsilon_{\theta}^{\beta \gamma} \partial_\beta A_\gamma \left\langle i\left| L^\theta \right| j\right\rangle & (\mathrm{M1})\\
\end{cases}$

#### Final Remarks:
This is the final expression of the Rabi-frequencies. As far as we are interested in only electric multipole expansion we took all required terms from the dirac equation. We worked in vacuum and in Coulomb gauge. Not switching the gauge allowed us to not to make any mistakes that arise from working in multiple gauges. 

Other common derivation is using PZW Gauge, which naturally has a form of multipole expansion. I, however, prefered not to work in it, as from what I have seen it wasn't a popular choice of understanding the problem. Coulomb's gauge was a preffered choice, however for many they missed a step to split the first Taylor expansion term into symmetric and antisymmetric part, which forced them to fudge a factor of 1/2. 

It is important to note that the above derivation is valid for the cases when the interaction is weak, and so the perturbation theory is applicable. This is because $i \hbar \dot{r}_\mu=\left[r_\mu, H_0\right]$ uses this assumption. 