
## AC Stark Shifts:

<span style="color:red">NOT FINISHED </span>

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