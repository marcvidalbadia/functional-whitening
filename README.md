This repository contains the results and source code for the numerical experiments in

[**Vidal, M**, Aguilera, AM. Novel whitening approaches in functional settings. Stat. 2022;e516.]( https://doi.org/10.1002/sta4.516)

Execute `Simulation.R` in R environments for further details.

See aslo the associated [Shiny app](https://mvidal.shinyapps.io/whitening/) and the [pfica](https://github.com/m-vidal/pfica) package v0.1.3.

Note:
Due to production errors, equation 3 in pp. 3 should be written as:

```math
\langle f,g\rangle_{\mathbb{M}}=\sum_{j=1}^{\infty}\lambda^{-1}_j\left\langle f,\gamma_{j}\right\rangle \left\langle g,\gamma_{j}\right\rangle =\left\langle \varGamma^{1/2\dagger}f,\varGamma^{1/2\dagger}g\right\rangle \quad \\\ f,g \in \mathbb{M}.
```
