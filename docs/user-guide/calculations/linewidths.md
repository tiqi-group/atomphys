<span style="color:red">TD </span>

### Transition Specific A-Einstein Coefficients:


In atomphys we also allow to calculate transition specific linewidths. Each transition has an associated A-einstein-coefficient, which characterises the strength of this transition. Following (James, 1998) we can write this to be:

$$
A_{ij}^{(E 1)}=\frac{4 c \alpha k_{12}^3}{3} \sum_{q=-1}^1\left|\left\langle i\left|r C_q^{(1)}\right| j\right\rangle\right|^2
$$

$$
A_{ij}^{(E 2)}=\frac{c \alpha k_{12}^5}{15} \sum_{q=-2}^2\left|\left\langle i\left|r^2 C_q^{(2)}\right| j\right\rangle\right|^2
$$

This can be re-written as


$$
A_{i j}=\sum_{q=-m}^m\left(\begin{array}{ccc}
J_i & m & J_j \\
-m_{j_i} & q & m_{j_j}
\end{array}\right)^2 \left(2 J_j+1\right) A_{I J}
$$

with 
- $A_{ij}$ standing for a Zeeman sublevel specific transition  A-Einstein coefficient,
- $A_{IJ}$ standing for a Zeeman manifold transition  A-Einstein coefficient,


Similar result can be derived using combination of considerations outlined in the Rabi Frequency section and Wigner-Weisskopf formalism