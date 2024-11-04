# Hamiltonians:

There are multiple hamiltonians that can be explored. 










# Derivation of Zeeman Hamiltonian with arbitrary B-field

For the purpose of Pavel's project (Grating MOT simulations) we need to consider Zeeman Hamiltonian due to the field pointing in arbitrary direction. As in atomphys I express all states w.r.t. quantization axis pointing in direction $\hat{z}$, I need to now express the hamiltonian $\mathbf{B}\cdot\mathbf{J}$ in arbitrary basis.

This is rather simple, as we can simply expand all the operators in those basis.

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

Therefore

$$
\begin{aligned}
\mathbf{B}\cdot\hat{\mathbf{J}} =\\&
\sum_{(n,l,s,j,m_j)\in \{S\}}\left(\frac{1}{2} \left( B_x-iB_y \right) \right)\hbar \sqrt{j(j+1)-m_j(m_j+1)}|n, l, s, j, m_j+1\rangle \langle n, l, s, j, m_j|
\\+&\sum_{(n,l,s,j,m_j)\in \{S\}}\left(\frac{1}{2} \left(B_x + iB_y  \right) \right)\hbar \sqrt{j(j+1)-m_j(m_j-1)}|n, l, s, j, m_j-1\rangle \langle n, l, s, j, m_j|
\\+&\sum_{(n,l,s,j,m_j)\in \{S\}}B_z \hbar m_j |n, l, s, j, m_j \rangle \langle n, l, s, j, m_j|
\end{aligned}
$$

This translates to atomphys code of 

$$
\begin{aligned}
\mathbf{B}\cdot\hat{\mathbf{J}} =\\&
\sum_{\substack{(\zeta_1, \zeta_2)\in \{S\}|\\(n,l,s,j,m_j+1)^{\zeta_1}=(n,l,s,j,m_j)^{\zeta_2}}}\left(\frac{1}{2} \left( B_x-iB_y \right) \right)\hbar \sqrt{j^{(\zeta_1)}(j^{(\zeta_1)}+1)-m_j^{(\zeta_1)}(m_j^{(\zeta_1)}+1)}|\zeta_2 \rangle \langle \zeta_1 |
\\+&\sum_{\substack{(\zeta_1, \zeta_2)\in \{S\}|\\(n,l,s,j,m_j-1)^{\zeta_1}=(n,l,s,j,m_j)^{\zeta_2}}}\left(\frac{1}{2} \left(B_x + iB_y  \right) \right)\hbar \sqrt{j^{(\zeta_1)}(j^{(\zeta_1)}+1)-m_j^{(\zeta_1)}(m_j^{(\zeta_1)}-1)}|\zeta_2 \rangle \langle \zeta_1 |
\\+&\sum_{\zeta \in \{S\}}B_z \hbar m_j^{(\zeta)} |\zeta \rangle \langle \zeta |
\end{aligned}
$$

Therefore 

$$
\begin{aligned}
H_B = 
-\frac{\mu_b g_J}{\hbar}\mathbf{B}\cdot\hat{\mathbf{J}} =\\&
\sum_{\substack{(\zeta_1, \zeta_2)\in \{S\}|\\(n,l,s,j,m_j+1)^{\zeta_1}=(n,l,s,j,m_j)^{\zeta_2}}}
-\mu_b g_J (j^{(\zeta_1)}, l^{(\zeta_1)}, s^{(\zeta_1)} )\left(\frac{1}{2} \left( B_x-iB_y \right) \right) \sqrt{j^{(\zeta_1)}(j^{(\zeta_1)}+1)-m_j^{(\zeta_1)}(m_j^{(\zeta_1)}+1)}|\zeta_2 \rangle \langle \zeta_1 |
\\+&\sum_{\substack{(\zeta_1, \zeta_2)\in \{S\}|\\(n,l,s,j,m_j-1)^{\zeta_1}=(n,l,s,j,m_j)^{\zeta_2}}}
-\mu_b g_J (j^{(\zeta_1)}, l^{(\zeta_1)}, s^{(\zeta_1)} )\left(\frac{1}{2} \left(B_x + iB_y  \right) \right) \sqrt{j^{(\zeta_1)}(j^{(\zeta_1)}+1)-m_j^{(\zeta_1)}(m_j^{(\zeta_1)}-1)}|\zeta_2 \rangle \langle \zeta_1 |
\\+&\sum_{\zeta \in \{S\}}
-\mu_b g_J (j^{(\zeta)}, l^{(\zeta)}, s^{(\zeta)} )B_z m_j^{(\zeta)} |\zeta \rangle \langle \zeta |
\end{aligned}
$$