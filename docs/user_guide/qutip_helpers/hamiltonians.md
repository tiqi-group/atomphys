# Hamiltonians:

Let's work in Coulomb's gauge in vacuum, where $\nabla \cdot \mathbf{A}(\mathbf{r}) = 0$. In the presence of electromagnetic field the total Hamiltonian for an atom can be expressed as: 

$$
\begin{aligned}
H &=\underbrace{\frac{\mathbf{p}^2}{2 m}+V_{\text {Coul }}}_{H_0}
+\underbrace{\left(-\frac{q \mathbf{p} \cdot \mathbf{A}\left(\mathbf{x}, t\right)}{m}+\frac{q^2\mathbf{A}\left(\mathbf{x}, t\right)^2}{2 m}\right)}_{H_I} + 
\underbrace{\frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \mathbf{B}}_{H_Z}
\end{aligned}
$$



Neglecting less dominant term of the interaction Hamiltonian $\frac{q_\alpha^2\mathbf{A}\left(\mathbf{x}_\alpha, t\right)^2}{2 m_\alpha}$, and grouping them together we get:

$$ 
\begin{aligned}
&H = H_0 + H_I + H_Z \\
&H_0 = \sum_\alpha \frac{\mathbf{p^{(\alpha)}}^2}{2 m^{(\alpha)}}+V_{\text {Coul }} = \sum_i \mathcal{E}_i\left|i\right>\left< i\right| \\
&H_I = \sum_\alpha-q^{(\alpha)} \mathbf{\dot{x}}^{(\alpha)} \cdot \mathbf{A}\left(\mathbf{x}^{(\alpha)}, t\right)
\\
&H_Z = \frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \mathbf{B}
\end{aligned}
$$

, where $H_0$ is the Bare Hamiltonian, $H_I$ is the interaction Hamiltonian and $H_Z$ is the Zeeman Hamiltonian.

## $H_0$ - Bare Hamiltonian:

$$
H_0 = \sum_\alpha \frac{\mathbf{p^{(\alpha)}}^2}{2 m^{(\alpha)}}+V_{\text {Coul }} = \sum_i \mathcal{E}_i\left|i\right>\left< i\right|
$$

To calculate the Bare Hamiltonian the user can select an **atom**, and the **states**. Atomphys based on the atom database will then construct a Hamiltonian object (Matrix) with energies of the states given. The states will be primarily sorted in the increasing order of the energies, and secondarily sorted by their angular momentum projection along the $z$-axis - $m_j$.

## $H_Z$ - Zeeman Hamiltonian:

To calculate the zeeman hamiltonian the user can do it in two different ways. Firstly, one can do it with the B-field pointing in the same direction as the quantization axis. Secondly, one can do it with the B-field pointing in an aribtrary direction. We always keep the quantization axis in the $z$-direction as a default, and it is the magnetic field that can be changed. 

### Zeeman Hamiltonian with B-field pointing in the $z$-direction:

$$
H_Z = \frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \mathbf{B} = \sum_i \left( \mu_b g_J \mathbf{B}\cdot\hat{\mathbf{J}} \right) \left|i\right>\left< i\right| = \sum_i \left( \mu_b g_J B_z\cdot\hat{J_z} \right) \left|i\right>\left< i\right|
$$

### Zeeman Hamiltonian with arbitrary B-field direction:

$$
H_Z = \frac{e \hbar}{2 m c} \boldsymbol{\sigma} \cdot \mathbf{B} = \sum_i \left( \mu_b g_J \mathbf{B}\cdot\hat{\mathbf{J}} \right) \left|i\right>\left< i\right|
$$

$$
\begin{aligned}
\mathbf{B}\cdot\hat{\mathbf{J}} &= B_x\hat{J_x}+B_y\hat{J_y}+B_z\hat{J_z}\\
&= B_x\left(\frac{1}{2} \left(\hat{J_+} + \hat{J_-} \right) \right)+B_y\left(\frac{1}{2i} \left(\hat{J_+} - \hat{J_-} \right) \right)+B_z\hat{J_z}\\
&= \left(\frac{1}{2} \left( B_x-iB_y \right) \right)\hat{J_+}+\left(\frac{1}{2} \left(B_x + iB_y  \right) \right)\hat{J_-}+B_z\hat{J_z}
\end{aligned}
$$

As
> $$
> J_{+}|j, m\rangle=\hbar \sqrt{(j-m)(j+m+1)}|j, m+1\rangle=\hbar \sqrt{j(j+1)-m(m+1)}|j, m+1\rangle
> $$
> $$
> J_{-}|j, m\rangle=\hbar \sqrt{(j+m)(j-m+1)}|j, m-1\rangle=\hbar \sqrt{j(j+1)-m(m-1)}|j, m-1\rangle
> $$

$$
\begin{aligned}
\mathbf{B}\cdot\hat{\mathbf{J}} =\\&
\sum_{(n,l,s,j,m_j)\in \{S\}}\left(\frac{1}{2} \left( B_x-iB_y \right) \right)\hbar \sqrt{j(j+1)-m_j(m_j+1)}|n, l, s, j, m_j+1\rangle \langle n, l, s, j, m_j|
\\+&\sum_{(n,l,s,j,m_j)\in \{S\}}\left(\frac{1}{2} \left(B_x + iB_y  \right) \right)\hbar \sqrt{j(j+1)-m_j(m_j-1)}|n, l, s, j, m_j-1\rangle \langle n, l, s, j, m_j|
\\+&\sum_{(n,l,s,j,m_j)\in \{S\}}B_z \hbar m_j |n, l, s, j, m_j \rangle \langle n, l, s, j, m_j|
\end{aligned}
$$

Therefore 

$$
\begin{aligned}
H_B = 
-\mu_b g_J\mathbf{B}\cdot\hat{\mathbf{J}} =\\&
\sum_{\substack{(\zeta_1, \zeta_2)\in \{S\}|\\(n,l,s,j,m_j+1)^{\zeta_1}=(n,l,s,j,m_j)^{\zeta_2}}}
-\mu_b g_J (j^{(\zeta_1)}, l^{(\zeta_1)}, s^{(\zeta_1)} )\left(\frac{1}{2} \left( B_x-iB_y \right) \right) \sqrt{j^{(\zeta_1)}(j^{(\zeta_1)}+1)-m_j^{(\zeta_1)}(m_j^{(\zeta_1)}+1)}|\zeta_2 \rangle \langle \zeta_1 |
\\+&\sum_{\substack{(\zeta_1, \zeta_2)\in \{S\}|\\(n,l,s,j,m_j-1)^{\zeta_1}=(n,l,s,j,m_j)^{\zeta_2}}}
-\mu_b g_J (j^{(\zeta_1)}, l^{(\zeta_1)}, s^{(\zeta_1)} )\left(\frac{1}{2} \left(B_x + iB_y  \right) \right) \sqrt{j^{(\zeta_1)}(j^{(\zeta_1)}+1)-m_j^{(\zeta_1)}(m_j^{(\zeta_1)}-1)}|\zeta_2 \rangle \langle \zeta_1 |
\\+&\sum_{\zeta \in \{S\}}
-\mu_b g_J (j^{(\zeta)}, l^{(\zeta)}, s^{(\zeta)} )B_z m_j^{(\zeta)} |\zeta \rangle \langle \zeta |
\end{aligned}
$$


## $H_I$ - Off-diagonal part of the Interaction Hamiltonian:

The interaction Hamiltonian can be constructed by providing the **atom**, **states** and **laser**. Atomphys will then calculate the off-diagonal elements of the interaction Hamiltonian.

For the interaction Hamiltonian, quite standard procedure involves going into the interaction picture. We believe that the interaction picture is usually a more subtle step and as the users are physicists, we decided to allow the user to calculate the off-diagonal elements of the interaction hamiltonian stripping them of the time-dependent parts. 

It is at the discretion of the user to go into the interaction picture, by defining their own detuning as a diagonal Hamiltonian. 

$$
\begin{aligned}
H_I &= \sum_{i, j} \left|i\right> \left< i\left|e \mathbf{\dot{r}} \mathbf{A} \right| j\right> \left< j\right| \\
&= \sum_{i, j} \left|i\right> \frac{1}{2} \left( \left< i\left|e \mathbf{\dot{r}} \mathbf{A}(\mathbf{x}) \right| j\right> e^{-i\omega_lt}
+ \left< i\left|e \mathbf{\dot{r}} \mathbf{A}^{\dagger}(\mathbf{x}) \right| j\right> e^{i\omega_lt} \right) \left< j \right| \\
&= \sum_{, i, j} \left|i\right> \frac{\hbar}{2} \left(\Omega_{ij} e^{-i\omega_lt}
+ \Omega^\dagger_{ij} e^{i\omega_lt} \right) \left< j\right| 
\end{aligned}
$$

Taking out the time-dependent parts, we get:

$$
H_I = \sum_{i, j} \left|i\right> \frac{\hbar}{2} \left(\Omega_{ij}
+ \Omega^\dagger_{ij} \right) \left< j\right| 
$$

, where Rabi Frequency is defined as follows:

$$
\hbar\Omega_{i j}= \left< i\left|e \mathbf{\dot{r}} \mathbf{A}(\mathbf{x}) \right| j\right>
$$

Depending on the transition type, the Rabi-Frequencies will be calculated differently:

$\Omega_{i j}= \begin{cases}
\frac{e\omega_{0}}{\hbar\omega_l}E^\mu\left\langle i\left| r^{(\alpha)}_{\mu} \right| j\right\rangle & (\mathrm{E1})
\\  \frac{e\omega_0}{2\hbar\omega_l} \partial^\nu E^\mu \left\langle i\left| r^{(\alpha)}_{\mu}  r^{(\alpha)}_{\nu}  \right| j\right\rangle & (\mathrm{E2})\\
\frac{1}{2} e \varepsilon_{\theta}^{\beta \gamma} \partial_\beta A_\gamma \left\langle i\left| L^\theta \right| j\right\rangle & (\mathrm{M1})\\
\end{cases}$


## $H_I'$ - Diagonal part of the Interaction Hamiltonian:

For the diagonal part of the interaction Hamiltonian, the user can provide the **atom**, **states** and list of **detunings**. Atomphys will then create an adequate matrix with the detunings on the diagonal.


