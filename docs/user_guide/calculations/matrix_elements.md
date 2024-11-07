
### Matrix Elements:

Before diving into the specific formulas and derivations used, it's worth pausing to consider the key concepts in atomic physics that we focus on and their significance. Atomic physics heavily relies on the calculation of matrix elements, denoted as $M_{if} = \left\langle f \mid \hat{M} \mid i \right\rangle$. This calculation involves the projection of an operator $\hat{M}$, which acts on an initial state $\left| i \right>$, onto a final state $\left| f \right>$. The value of $M_{if}$ is crucial because it measures the extent of interaction, or coupling, between the states $\left| i \right>$ and $\left| f \right>$ via the operator $\hat{M}$, when $\hat{M}$ is part of the Hamiltonian $\hat{H}$. Understanding this interaction provides insights into the dynamics of the overall state $\left| \psi(t) \right>$.

However, computing these matrix elements is challenging due to involved basis. One can simplify (break-down) this process using symmetry principles through Wigner-Eckart Theorem.

The [Wigner-Eckart Theorem](https://en.wikipedia.org/wiki/Wigner%E2%80%93Eckart_theorem) states that the matrix elements of spherical tensor operators, within the framework of angular momentum eigenstates, can be broken down into two components: one that is independent of the angular momentum orientation (so called reduced matrix element) and another that is a Clebsch-Gordan coefficient.

$\left\langle j m\left|T_q^{(k)}\right| j^{\prime} m^{\prime}\right\rangle=\left\langle j^{\prime} m^{\prime} k q \mid j m\right\rangle\left\langle j\left\|T^{(k)}\right\| j^{\prime}\right\rangle = [j] \left(\begin{array}{ccc}j^{\prime} & k & j^* \\ m^{\prime} & q & m\end{array}\right)  \left\langle j\left\|T^{(k)}\right\| j^{\prime}\right\rangle$

Here, $[j]$ represents the dimension of the representation space of $j$, and $j^*$ indicates a complex conjugate representation to $j$. The parentheses contain the Wigner-3j symbol, and $\left\langle \mid\mid . \mid\mid \right\rangle$ denotes a reduced matrix element.

As the reduced matrix element, doesn't depend on basis this is the most basic block from which we will start the calculation.


#### Reduced Matrix Elements (Dipole and Quadrupole):
Following [James 1998](https://doi.org/10.1007/s003400050373)

$A_{12}^{(E 1)} \equiv \sum_{m=-j}^j \bar{A}_{12}^{(E 1)}=\frac{4 c \alpha k_{12}^3}{3\left(2 j^{\prime}+1\right)}\left|\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle\right|^2$

$A_{12}^{(E 2)} \equiv \sum_{m=-j}^j \bar{A}_{12}^{(E 2)}=\frac{c \alpha k_{12}^5}{15\left(2 j^{\prime}+1\right)}\left|\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle\right|^2$

we can rewrite 

$\left|\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle\right| = \left( \frac{3\left(2 j^{\prime}+1\right)}{4 c \alpha k_{12}^3} A_{12}^{(E 1)} \right)^{\frac{1}{2}}$

$\left|\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle\right| = \left(\frac{15\left(2 j^{\prime}+1\right)}{c \alpha k_{12}^5} A_{12}^{(E 2)}\right)^{\frac{1}{2}}$

Equivalent result can be derived from Wigner-Weisskopf Decay.



#### Matrix Elements

Now in order to transform it to actual quantities that we want we can use Wigner-Eckart theorem.

$\left\langle 1\left|\hat{r}_i\right| 2\right\rangle \epsilon_i=\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle \sum_{q=-1}^1\left(\begin{array}{ccc}j & 1 & j^{\prime} \\ -m_j & q & m_j^{\prime}\end{array}\right) c_i^{(q)} \epsilon_i$

$\left\langle 1\left|\hat{r}_i \hat{r}_j\right| 2\right\rangle \epsilon_i n_j=\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle \sum_{q=-2}^2\left(\begin{array}{ccc}j & 2 & j^{\prime} \\ -m & q & m^{\prime}\end{array}\right) c_{i j}^{(q)} \epsilon_i n_j$

keeping only tensorial parts of the equations we drop $\epsilon_i$ and $n_j$.

$\left\langle 1\left|\hat{r}_i\right| 2\right\rangle =\left\langle 1\left\|r \mathbf{C}^{(1)}\right\| 2\right\rangle \sum_{q=-1}^1\left(\begin{array}{ccc}j & 1 & j^{\prime} \\ -m_j & q & m_j^{\prime}\end{array}\right) c_i^{(q)}$

$\left\langle 1\left|\hat{r}_i \hat{r}_j\right| 2\right\rangle=\left\langle 1\left\|r^2 \mathbf{C}^{(2)}\right\| 2\right\rangle \sum_{q=-2}^2\left(\begin{array}{ccc}j & 2 & j^{\prime} \\ -m & q & m^{\prime}\end{array}\right) c_{i j}^{(q)}$


#### Electric Matrix Elements

Electric matrix elements are just the conventional dipole/quadrupole elements scaled by a funamental charge, $e$. 