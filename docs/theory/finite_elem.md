# Elements of Functional Analysis
This section collects basic knowledge of weak formulations and functional spaces to understand the math behind the finite element method {cite:p}`salsa, Gazzola_A3, quarteroni_2016`.

Let $\Omega\subset\mathbb{R}^{d}$ with $d=2,3$ and $\partial \Omega$ its boundary.

## $L^2$ and Sobolev spaces

Let $L^2(\Omega)$ be an Hilbert space of square integrable functions, i.e. 

$$
u\in L^2(\Omega)\Longleftrightarrow \int_\Omega |u|^2\,d\Omega <\infty
$$

Since $L^2$ is an Hilbert space, an inner product and an induced norm can be endowed

$$
(u, \, v) = \int_\Omega u\cdot v\,d\Omega \qquad \qquad
\|u\|_{L^2}^2 = \int_\Omega |u|^2\,d\Omega
$$

The same space can be easily extended to vector quantities and complex functions (the conjugate is introduced).

Let us introduced the notion of multi-index derivative  $D^{{\alpha}}$: given ${\alpha}\in\mathbb{N}^d$ of order $p=\sum_{i=1}^d \alpha_i$, the multi-index derivative of a function is defined as
\begin{equation*}
    D^{{\alpha}}u = \frac{\partial^pu}{\partial x_1^{\alpha_1} \dots \partial x_d^{\alpha_d}},
\end{equation*}
in which $d$ is usually 2 or 3 {cite}`Gazzola_A3`.

Let $\mathcal{H}^p(\Omega)$ be the Sobolev space of order $p$, i.e. 

$$
\mathcal{H}^p(\Omega) = \left\{u\in L^2(\Omega)\,:\,\int_\Omega \left|D^{{\alpha}}u\right|^2\,d\Omega <\infty \right\}
$$

All these derivatives are meant to be in the weak/distribution sense {cite}`Gazzola_A3`.

## Useful formulas

Let $\mathbf{u}\in[\mathcal{H}^1(\Omega)]^d$ a vector function, the Gauss (divergence) theorem states:

$$
\int_\Omega \nabla \cdot \mathbf{u}\,d\Omega = \int_{\partial\Omega}\mathbf{u}\cdot \mathbf{n}\,d\sigma
$$

From this theorem, the following formulas can be derived
\begin{equation*}
\int_\Omega (\nabla \cdot \mathbf{u}) \,p\,d\Omega = 
- \int_\Omega \mathbf{u} \cdot \nabla p\,d\Omega + \int_{\partial\Omega}(\mathbf{u}\cdot \mathbf{n})\,p\,d\sigma
\end{equation*}
given $p\in\mathcal{H}^1(\Omega)$.

````{prf:proof}
Recalling that $ \nabla\cdot (p\,\mathbf{u})=(\nabla \cdot \mathbf{u}) \,p + \mathbf{u} \cdot \nabla p$ we can rewrite the formula as

$$
\int_\Omega \left[(\nabla \cdot \mathbf{u}) \,p + \mathbf{u} \cdot \nabla p\right] \,d\Omega = 
\int_\Omega \nabla \cdot (p\mathbf{u})\,d\Omega = \int_{\partial\Omega}(\mathbf{u}\cdot \mathbf{n})\,p\,d\sigma
$$
from which thesis follows.
````

## Weak Formulations

Let $\Omega\subset\mathbb{R}^{d}$ with $d=2,3$ and $\partial \Omega = \Gamma_D\cup\Gamma_N$ its boundary, let us consider the following problem (in strong form):
\begin{equation*}
\left\{
    \begin{array}{ll}
        -\Delta u = f & \mbox{in } \Omega\\
        u = u_D & \mbox{on } \Gamma_D\\
        \displaystyle \frac{\partial u}{\partial \mathbf{n}} = g_N & \mbox{on } \Gamma_N
    \end{array}
\right.
\end{equation*}
We have a PDE to solve with its boundary condition (BC) associated. The strong solution must be differentiable twice in classical sense and, in general, it is not always possible to prove the well-posedness of these problems. Therefore, it comes the necessity to look for the solution into a broader functional space. Without entering into too much details on this topic, this sections aims to present the typical procedure to be followed when a weak formulation must be derived (its necessary to assemble the algebraic system coming from the Finite element method).

Let $\mathcal{V}\subset\mathcal{H}^1, \mathcal{V}_0\subset\mathcal{H}^1$ the trial and test space defined as

$$
\mathcal{V} = \left\{v\in\mathcal{H}^1:\;\left. v\right|_{\Gamma_D} = u_D\right\}\qquad 
\mathcal{V}_0 = \left\{v\in\mathcal{H}^1:\;\left. v\right|_{\Gamma_D} = 0\right\}
$$

Given $v\in\mathcal{V}_0$, the PDE is multipied (in $L^2$ sense) by $v$

$$
\int_\Omega -\Delta u\,v \,d\Omega=\int_\Omega f\,v\,d\Omega\qquad v\in\mathcal{V}_0
$$
and by applying the integration by parts formula and the effect of the BCs
\begin{equation*}
\int_\Omega \nabla u\cdot \nabla v \,d\Omega=\int_\Omega f\,v\,d\Omega+\int_{\Gamma_N}g\,v\,d\sigma\qquad v\in\mathcal{V}_0
\end{equation*}

### Galerkin approximation
Let us introduce a computational grid $\mathcal{T}_h$ (dim$\mathcal{T}_h = \mathcal{N}_h$) of the domain $\Omega$ and the finite dimensional spaces $\mathcal{V}_h\subset\mathcal{V}$ and $\mathcal{V}_{h,0}\subset\mathcal{V}_0$. The Galerkin approximation {cite}`quarteroni_2016, QuarteroniMateNumerica` reads: *find $u_h\in\mathcal{V}_h$ s.t.
\begin{equation*}
\int_\Omega \nabla u_h\cdot \nabla v_h \,d\Omega=\int_\Omega f_h\,v_h\,d\Omega+\int_{\Gamma_N}g\,v_h\,d\sigma\qquad v_h\in\mathcal{V}_{h,0}
\end{equation*}

### Finite element algebraic formulation
The finite dimensional formulation is necessary to assemble the associated linear system. 
Let us consider the homogeneous Dirichlet case, $u_D = 0$ (to which we can
always reduce through a lifting procedure), and let us introduce a basis $\left\{\phi_j\right\}_{j=1}^{\mathcal{N}_h}$ such that any function $u_h\in\mathcal{V}_h$ can be expressed as a linear combination of the basis functions

$$
u_h(\mathbf{x}) = \sum_{j=1}^{\mathcal{N}_h}\alpha_j \phi_j(\mathbf{x})
$$

```{note}

The choice of the basis functions results in different versions of the algebraic system. The most natural choice is given by the **hat functions** {cite:p}`QuarteroniMateNumerica`.
```

By taking $v_h=\phi_k$ as the test function, the Galerkin problem reduces to
\begin{equation*}
\int_\Omega \nabla u_h\cdot \nabla \phi_k \,d\Omega=\int_\Omega f_h\,\phi_k\,d\Omega+\int_{\Gamma_N}g\,\phi_k\,d\sigma\qquad k = 1, \dots, \mathcal{N}_h
\end{equation*}
then, the linear expansion can be introduced into the problem which results

$$
\sum_{j=1}^{\mathcal{N}_h}\alpha_j \int_\Omega \nabla \phi_j\cdot \nabla \phi_k \,d\Omega=\int_\Omega f_h\,\phi_k\,d\Omega+\int_{\Gamma_N}g\,\phi_k\,d\sigma\qquad k = 1, \dots, \mathcal{N}_h
$$
By defining the matrix $A_{kj} = \displaystyle\int_\Omega \nabla \phi_j\cdot \nabla \phi_k \,d\Omega$ and the vectors $f_k =  \displaystyle\int_\Omega f_h\,\phi_k\,d\Omega$ and $g_k =  \displaystyle\int_{\Gamma_N}g\,\phi_k\,d\sigma$, the problem into the following linear system in the unknowns $\alpha_j$
\begin{equation*}
\sum_{j=1}^{\mathcal{N}_h}A_{kj}\alpha_j = f_k + g_k\qquad k = 1, \dots, \mathcal{N}_h
\end{equation*}