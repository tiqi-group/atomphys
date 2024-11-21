# Lindblad Jump Operators:

One can describe the incoherent dynamics of the system through the Lindblad master equation. For $k$ denoting different environment states the Lindblad operator $\mathcal{L}^{\mathrm{d}}$ can be written as:

$$
\mathcal{L}^{\mathrm{d}}(\hat{\rho}) = -\frac{1}{2}\sum_k \left(\hat{c}_k^\dagger \hat{c}_k \hat{\rho} + \hat{\rho} \hat{c}_k^\dagger \hat{c}_k - 2\hat{c}_k \hat{\rho} \hat{c}_k^\dagger \right)
$$

Depending on the cause of the decoherence, the jump operators $\hat{c}_k$ will be different. One can construct the jump operators for decoherence due to the laser dephasing, spontaneous emission, etc. 

In atomphys we will slowly build up the library of jump operators for different decoherence mechanisms.

## Jump Operators for Spontaneous Emission:
In Atomphys the user can automatically generate the spontaneous emission jump operators for their atom and states of interest.

For spontaneous emission, jump operators between levels $\left|f\right>$ and $\left|i\right>$ can be written as:

$$
\hat{c}_k = \sqrt{\Gamma_{fi}} \left|f\right>\left<i\right|
$$

, where $\Gamma_{fi}$ is the zeeman-sublevel transition specific linewidth.

The derivation of the master equation requires that for each final state of the environment, there can only be one jump operator. Therefore all decay channels leading to the same final state have to be grouped together.

$$
\hat{c}_k = \sum_{(i, f) \in k} \sqrt{\Gamma_{fi}} \left|f\right>\left<i\right|
$$

Grouping the decay channels together is not yet implemented in the atomphys package, but will be in the future.


