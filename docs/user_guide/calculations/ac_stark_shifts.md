
## AC Stark Shifts:

In order to introduce the theory behind the AC Stark Shifts, We will closely follow [Beloy 2009](https://www.dereviankogroup.com/dereviankogroup/resources/Student_work/Beloy_Dissertation_2009.pdf). We intend to present the calculation of AC Stark Shift without splitting it into the different scalar, vector and tensor components. This allows us to directly understand the simplicity of the AC Stark Shift calculations, not make any assumptions about the structure of the Electric Field or allowed transitions. Normal derivations of AC Stark Shift only consider dipole transitions driven by running wave field. 

Consider oscilating electric field 

$\mathbf{A_l}=\frac{1}{2}\left(\mathbf{A}_{l}(\mathbf{x})e^{-i\omega_lt}+\mathbf{A}_{l}^\dagger(\mathbf{x})e^{i\omega_lt}\right)$

Such field interacts with an atom perturbing the hamiltonian by additional term (see section Rabi Frequencies):

$\mathcal{H}_L=-\sum_{i, j, l}\left|i\right> \frac{1}{2}\left(\Omega_{ij}^{(l)} e^{-i \omega_l t}+ \Omega_{ij}^{\dagger(l)} e^{i \omega_l t}\right)\left< j\right|$

It is noted that $\Omega_{ij}^{(l)}$ is a time independent operator. For small field amplitudes (condition), one can treat the $\mathcal{H}_L$ as small perturbation to $\mathcal{H}_0$. Using Floquet perturbation theory, one can develop expression for the mean energy shifts $\delta \mathcal{E} = \delta \mathcal{E}^{(1)} + \delta \mathcal{E}^{(2)} + ...$ . Due to symmetry  $\delta \mathcal{E}^{(1)}$ is not present, and so we will start from evaluating $\delta \mathcal{E}^{(2)}$.

Details of the Floquet Perturbation Theory can be found in [Beloy 2009](https://www.dereviankogroup.com/dereviankogroup/resources/Student_work/Beloy_Dissertation_2009.pdf).

### The second order AC Stark Shift

From [Beloy 2009](https://www.dereviankogroup.com/dereviankogroup/resources/Student_work/Beloy_Dissertation_2009.pdf) we know that:

$\delta \mathcal{E}_i^{(2)}=\frac{1}{4}\sum_{j} \frac{\Omega_{ij}^{(l)}\Omega_{ij}^{\dagger(l)}}{\mathcal{E}_i-\mathcal{E}_{j}+\omega_l}+\frac{1}{4}\sum_{j} \frac{\Omega_{ij}^{(l)}\Omega_{ij}^{\dagger (l)}}{\mathcal{E}_i-\mathcal{E}_{j}-\omega_l}$

In other sources you can find that people than take this expression and make assumption about the running-wave laser field, and consider only dipole transitions. This allows them to split the expression into scalar, vectorial and tensorial components. For the case of atomphys, it is very easy to calculate the Rabi Frequencies so we don't need to make any of those assumptions. This allows us to then treat Electric Fields of more complex vectorial structure, and see for instance Magnus effect, and equivalent. 

Taking second order AC Stark Shift as the largest contribution to overall stark shift we say:

$\delta \mathcal{E}_i \approx \delta \mathcal{E}_i^{(2)}$ 


### Polarizability

Polarizability of an atom in state $\left|i\right>$ , $\alpha_i$, is an AC Stark Shift, $\delta \mathcal{E}_i$  normalized by the electric field causing this stark shift, $E_0$.

$\alpha_i = \frac{-4 \left(\delta \mathcal{E}_i\right)}{E_0^2}$

In our package we decided to calculate overall polarizability without splitting it into scalar, vectorial and tensorial polarizability. This is because we wanted to generalize it to arbitrary types of transitions including light shifts coming from quadrupole transitions. 


### Off-resonant scattering




$$
\langle i| d_\mu\epsilon^\mu|j\rangle = -\Omega_{ij}\frac{\hbar \omega_l}{ \omega_{ki} E_0}
$$
$$
\begin{aligned}
\Gamma_{i \rightarrow f}&=\frac{I \omega^{\prime 3}}{\left(4 \pi \epsilon_0\right)^2 c^4 \hbar^3} \frac{8 \pi}{3} \sum_{q=-1}^1\left|
\sum_k\left(\langle f| d_q|k\rangle \frac{\langle k| \mathbf{d}\cdot\mathbf{\epsilon}_\ell|i\rangle}{\omega_{k i}-\omega}+\langle f| \mathbf{d}\cdot\mathbf{\epsilon}_\ell|k\rangle \frac{\langle k| d_q|i\rangle}{\omega_{k i}+\omega^{\prime}}\right)
\right|^2\\
&=\frac{I \omega^{\prime 3}}{\left(4 \pi \epsilon_0\right)^2 c^4 \hbar^3} \frac{8 \pi}{3} \sum_{q=-1}^1\left|
\sum_k\left(\langle f| d_q|k\rangle \frac{\langle k| \mathbf{d}\cdot\mathbf{\epsilon}_\ell|i\rangle}{\omega_{k i}-\omega}+\langle f| \mathbf{d}\cdot\mathbf{\epsilon}_\ell|k\rangle \frac{\langle k| d_q|i\rangle}{\omega_{k i}+\omega^{\prime}}\right)
\right|^2\\
&=\frac{I \omega^{\prime 3}}{\left(4 \pi \epsilon_0\right)^2 c^4 \hbar^3} \frac{8 \pi}{3} \sum_{q=-1}^1\left|
\sum_k\left(\langle f| d_q|k\rangle \frac{\Omega_{ki}}{\omega_{k i}-\omega}\frac{\hbar\omega}{\omega_{ki}E_0}+\Omega_{fk}\frac{\hbar\omega}{\omega_{fk}E_0}\frac{\langle k| d_q|i\rangle}{\omega_{k i}+\omega^{\prime}}\right)
\right|^2\\
&=\frac{c\epsilon_0E_0^2 \omega^2\omega^{\prime 3}}{2E_0^2\left(4 \pi \epsilon_0\right)^2 c^4 \hbar} \frac{8 \pi}{3} \sum_{q=-1}^1\left|
\sum_k\left(\langle f| d_q|k\rangle \frac{\Omega_{ki}}{\omega_{k i}-\omega}\frac{1}{\omega_{ki}}+\Omega_{fk}\frac{1}{\omega_{fk}}\frac{\langle k| d_q|i\rangle}{\omega_{k i}+\omega^{\prime}}\right)
\right|^2\\
&=\frac{\omega^2\omega^{\prime 3}}{12 \pi c \hbar} \sum_{q=-1}^1\left|
\sum_k\left(\langle f| d_q|k\rangle \frac{\Omega_{ki}}{\omega_{k i}-\omega}\frac{1}{\omega_{ki}}+\Omega_{fk}\frac{1}{\omega_{fk}}\frac{\langle k| d_q|i\rangle}{\omega_{k i}+\omega^{\prime}}\right)
\right|^2\\
&=\frac{\omega^2\left(\omega_\ell+\omega_{if}\right)^{3}}{12 \pi c \hbar} \sum_{q=-1}^1\left|
\sum_k\left(\langle f| d_q|k\rangle \frac{\Omega_{ki}}{\omega_{k i}-\omega}\frac{1}{\omega_{ki}}+\Omega_{fk}\frac{1}{\omega_{fk}}\frac{\langle k| d_q|i\rangle}{\omega_{k i}+\omega_\ell+\omega_{if}}\right)
\right|^2
\end{aligned}
$$

$$
\Gamma_{i \rightarrow f}=\frac{\omega^2\left(\omega_\ell+\omega_{if}\right)^{3}}{12 \pi c \hbar} \sum_{q=-1}^1\left|
\sum_k\left(\langle f| d_q|k\rangle \frac{\Omega_{ki}}{\omega_{k i}-\omega}\frac{1}{\omega_{ki}}+\Omega_{fk}\frac{1}{\omega_{fk}}\frac{\langle k| d_q|i\rangle}{\omega_{k i}+\omega_\ell+\omega_{if}}\right)
\right|^2
$$