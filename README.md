This repository contains the results and source code for the numerical experiments in

[**M. Vidal** & A. M. Aguilera. Novel whitening approaches in functional settings. Stat, 12( 1), e516.]( https://doi.org/10.1002/sta4.516)

Execute `Simulation.R` in R environments for further details.

See aslo the associated [Shiny app](https://mvidal.shinyapps.io/whitening/) and the [pfica](https://github.com/m-vidal/pfica) package v0.1.3.

:exclamation: Note:

Due to production errors, equation 3 in pp. 3 is written as

```math
\langle f,g\rangle_{\mathbb{M}}=\sum_{j=1}^{\infty}\lambda^{-1}_j\left\langle f,\gamma_{j}\right\rangle \left\langle g,\gamma_{j}\right\rangle =\left\langle \varGamma^{1/2\dagger}f,\varGamma^{1/2\dagger}g\right\rangle f, \quad g \in \mathbb{M},
```
while was originally written as

```math
\langle f,g\rangle_{\mathbb{M}}=\sum_{j=1}^{\infty}\lambda^{-1}_j\left\langle f,\gamma_{j}\right\rangle \left\langle g,\gamma_{j}\right\rangle =\left\langle \varGamma^{1/2\dagger}f,\varGamma^{1/2\dagger}g\right\rangle  \quad  f,g \in \mathbb{M}.
```

In section 4, the sentence "As 2tr $(\varGamma_{X\mathbb{X}})$ is the only *dependence* between the original and the whitened variable, the minimization problem can be reduced to the maximization of tr $(\varGamma_{X\mathbb{X}})$." reads also as "... is the only **dependent term**...".

In the Technical proofs (first paragraph), due to abuse of notation, in the sentence "Note that Condition 1 cannot be reached when $\langle X,\gamma_j\rangle^2=\lambda_j$, or for $c_j\rightarrow c>0$, $\langle X,\gamma_j\rangle^2=\lambda_jc_j$...", **$X$ stands for a deterministic function**.
