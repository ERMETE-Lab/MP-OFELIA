# Neutronics

The neutron Boltzmann equation, $\phi = \phi(\mathbf{x}, E,\mathbf{\Omega}, t)$

$$
\left\{
\begin{split}
\frac{1}{v(E)}\frac{\partial \phi}{\partial t}  =&-\mathbf{\Omega}\cdot \nabla \Phi - \Sigma_{t}(\mathbf{x}, E, t)\,\phi \\
&+\chi_g(E)\,(1-\beta)\int_{4\pi}\int_0^{+\infty}\bar{\nu}(E')\Sigma_{f}(\mathbf{x},E',t)\,\phi'\,dE'd\mathbf{\Omega}'\\
&+\int_{4\pi}\int_0^{+\infty}\Sigma(\mathbf{x}, E'\rightarrow E, \mathbf{\Omega}'\rightarrow \mathbf{\Omega}, t)\,\phi' \,dE'd\mathbf{\Omega}'+\sum_{j=1}^J\chi^{d_j}(E)\lambda_jC_j(\mathbf{x}, E, t)\\
\frac{\partial C_j(\mathbf{x}, E, t)}{\partial t} =& -\lambda_jC_j(\mathbf{x},E, t)+\beta_j\int_{4\pi}\int_0^{+\infty}\bar{\nu}(E')\Sigma_{f}(\mathbf{x},E',t)\,\phi'\,dE'd\mathbf{\Omega}'
\end{split}
\right.
$$

This system of $J+1$ equations are linear (if no thermal feedback are considered) and different discretisatio techniques can be used to solve it.

In this repository, the multigroup diffusion approximation is used 

$$
\left\{
\begin{aligned}
    \frac{1}{v_g}\frac{\partial\phi_g}{\partial t}=&+\nabla\cdot(D_g\nabla \phi_g)-\left(\Sigma_{a,g}+\sum_{g'\neq g}\Sigma_{s,g\rightarrow g'} + D_g B_{z,g}^2\right)\phi_g\\
    &+\sum_{g'\neq g}\Sigma_{s,g'\rightarrow g}\phi_{g'}+\chi_g\left(\frac{1-\beta}{k_{\text{eff}}}\sum_{g'}\nu_{g'}\Sigma_{f,g'}\phi_{g'}+\sum_{j=1}^J\lambda_jc_j\right)\\
    \frac{\partial c_j}{\partial t} =&  \frac{\beta_j}{k_{\text{eff}}}\sum_{g}\nu_g\Sigma_{f,g}\phi_{g} - \lambda_jc_j
\end{aligned}
\right.
$$
as well as a P1 and P3 implementation for 1D problems.

The governing equations are derived by performing mass balances and considering all possible ways neutrons can be born, be transferred and die. More details can be found in {cite:p}`DuderstadtHamilton`.

## Neutron Diffusion

The stationary version of the multi-group diffusion equation is an eigenvalue-eigenfunction problem for the effective multiplication factor $k_{\text{eff}}$

\begin{equation*}
-\nabla\cdot(D_g\nabla \phi_g)+\left(\Sigma_{a,g}+\sum_{g'\neq g}\Sigma_{s,g\rightarrow g'} + D_g B_{z,g}^2\right)\phi_g-\sum_{g'\neq g}\Sigma_{s,g'\rightarrow g}\phi_{g'}=\frac{1}{k_{\text{eff}}}\cdot \chi_g^p\sum_{g'}\nu_{g'}\Sigma_{f,g'}\phi_{g'}
\end{equation*}

for $g=1,\ dots, G$. Its algebraic formulation has the following form

\begin{equation*}
\mathbb{A}\boldsymbol{\phi} = \frac{1}{k_{\text{eff}}}
\mathbb{B}\boldsymbol{\phi}
\end{equation*}

Whereas, the transient problem is completed by another set of equations related to the precursors
\begin{equation*}
\left\{
\begin{split}
\frac{1}{v_g}\frac{\partial\phi_g}{\partial t}=&+\nabla\cdot(D_g\nabla \phi_g)-\left(\Sigma_{a,g}+\sum_{g'\neq g}\Sigma_{s,g\rightarrow g'} + D_g B_{z,g}^2\right)\phi_g\\
    &+\sum_{g'\neq g}\Sigma_{s,g'\rightarrow g}\phi_{g'}+\chi_g^p\frac{1-\beta}{k_{\text{eff}}}\sum_{g'}\nu_{g'}\Sigma_{f,g'}\phi_{g'}+\sum_{j=1}^J\chi_g^{d_j}\lambda_jc_j\\
\frac{\partial C_j}{\partial t}=&-\lambda_jC_j+\beta_j\sum_{g=1}^G\nu_{g'}\Sigma_{g',f}\phi_{g'}\qquad \text{ for }j= 1, \dots J
\end{split}
\right.
\end{equation*}

Adopting an implicit Euler scheme for time-advancement, the linear problem at each time step is

\begin{equation*}
\left\{
\begin{aligned}
    \mathbb{V}\frac{\boldsymbol{\Phi}_{i+1} - \boldsymbol{\Phi}_{i}}{\Delta t}&=-\mathbb{A}\boldsymbol{\phi}_{i+1} + \frac{1}{k_{\text{eff}}}
    \mathbb{B}\boldsymbol{\phi}_{i+1}+\mathbb{D}\boldsymbol{c}_{i+1}\\
    \frac{\boldsymbol{c}_{i+1} - \boldsymbol{c}_{i}}{\Delta t} &= \frac{1}{k_{\text{eff}}}\mathbb{B}^*\boldsymbol{\phi}_{i+1}-\mathbb{D}^*\boldsymbol{c}_{i+1}
\end{aligned}
\right.
\end{equation*}