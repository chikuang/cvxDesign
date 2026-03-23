# cvxDesign: Optimal experimental designs via convex optimization

2026-03-23

<!-- badges: start -->

<!-- [![R-CMD-check](https://github.com/chikuang/cvxDesign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chikuang/cvxDesign/actions/workflows/R-CMD-check.yaml) -->

[![](https://img.shields.io/github/languages/code-size/chikuang/cvxDesign.svg)](https://github.com/chikuang/cvxDesign)
<!-- badges: end -->

**Authors:**  
Chi-Kuang Yeh (Georgia State University)  
[![ORCID](https://img.shields.io/badge/ORCID-0000--0001--7057--2096-A6CE39?logo=orcid.png)](https://orcid.org/0000-0001-7057-2096)

Julie Zhou (University of Victoria)

``` r
# options(scipen = 999)
options(digits = 4)
```

## Description

`cvxDesign` is an R package for computing optimal experimental designs
using convex optimization.  
The package is intended for researchers and practitioners working on
approximate and exact design problems under a variety of criteria,
including classical criteria such as D-, A-, E-, and c-optimality, as
well as compound and multi-objective design problems.

The package is motivated by modern design settings where efficient
numerical optimization, equivalence-theorem-based verification, and
flexible model specification are all essential.

## Installation

You can install the development version of `cvxDesign` from GitHub by
running

``` r
# install.packages("pak")
pak::pak("chikuang/cvxDesign")
```

or

``` r
# install.packages("remotes")
remotes::install_github("chikuang/cvxDesign")
```

## Details

Consider a regression model of the form

$$
y_i = \eta(\mathbf{x}_i, \boldsymbol{\theta}) + \varepsilon_i, \quad i = 1, \ldots, n,
$$

where $y_i$ is the response observed at design point
$\mathbf{x}_i \in \mathcal{X}$,  
$\boldsymbol{\theta} \in \mathbb{R}^p$ is an unknown parameter vector,
and $\eta(\mathbf{x}, \boldsymbol{\theta})$ is the mean response
function, which may be linear or nonlinear in the parameters.

Optimal design seeks a design measure $\xi$ over the design space
$\mathcal{X}$ that optimizes a scalar criterion based on the information
matrix $M(\xi)$. Common criteria include:

- **D-optimality**, which maximizes $\log \det M(\xi)$;
- **A-optimality**, which minimizes $\mathrm{tr}\{M^{-1}(\xi)\}$;
- **E-optimality**, which maximizes the minimum eigenvalue of $M(\xi)$;
- **c-optimality**, which minimizes the variance of estimating a
  specified linear combination of parameters.

In many applications, the design problem can be formulated and solved
using convex optimization. This makes it possible to handle a wide range
of criteria and constraints in a unified computational framework.

The goal of `cvxDesign` is to provide a clean and extensible interface
for these computations, together with tools for:

- computing optimal approximate designs,
- studying exact designs on finite candidate sets,
- evaluating design efficiency,
- checking equivalence theorems,
- visualizing directional derivative or sensitivity functions.

## Examples

First load the package.

``` r
# load package
library(cvxDesign)
```

### D-optimal design for a quadratic regression model

Consider the quadratic regression model

$$
y_i = \beta_0 + \beta_1 x_i + \beta_2 x_i^2 + \varepsilon_i.
$$

We first define the regression function.

``` r
quad_reg <- function(x) {
  c(1, x, x^2)
}
```

Suppose the candidate design space is a finite grid on $[-1,1]$.

``` r
u <- seq(-1, 1, length.out = 101)
head(u)
```

    [1] -1.00 -0.98 -0.96 -0.94 -0.92 -0.90

A typical workflow in `cvxDesign` will look like the following.

``` r
dout <- calc_Dopt(
  u = u,
  f = quad_reg,
  drop_tol = 1e-4
)

dout$design |> round(3)
```

      point weight
    1    -1  0.333
    2     0  0.333
    3     1  0.333

``` r
dout$value
```

    [1] -1.91

``` r
dout$status
```

    [1] "optimal"

The returned object is expected to contain at least:

- the selected support points,
- the corresponding optimal weights,
- the value of the design criterion,
- and, when requested, information useful for equivalence-theorem
  checks.

### A-optimal design

Similarly, one may compute the A-optimal design by changing the
criterion.

``` r
aout <- calc_Aopt(
  u = u,
  f = quad_reg,
  drop_tol = 1e-4
)

aout$design |> round(3)
```

      point weight
    1    -1   0.25
    2     0   0.50
    3     1   0.25

``` r
aout$value
```

    [1] 8

``` r
aout$status
```

    [1] "optimal"

### c-optimal design

If the target is a specific linear combination $c^\top \theta$, then a
c-optimal design can be computed by supplying the vector `cvec`.

``` r
cout <- calc_copt(
  u = u,
  f = quad_reg,
  cVec = c(0.3, 0.4, 0.3),
  drop_tol = 1e-4
)

cout$design |> round(3)
```

      point weight
    1    -1  0.125
    2     1  0.875

``` r
cout$value
```

    [1] 0.16

``` r
cout$status
```

    [1] "optimal"

### Equivalence theorem

An important part of optimal design computation is verification.  
The package is designed to support equivalence-theorem-based diagnostics
and plotting.

``` r
quad_reg <- function(x) c(1, x, x^2)
u <- seq(-1, 1, length.out = 101)

dout <- calc_Dopt(u, quad_reg, drop_tol = 1e-4)
eq_d <- check_equivalence(dout, f = quad_reg)

print(eq_d)
```

    $candidate_points
      [1] -1.00 -0.98 -0.96 -0.94 -0.92 -0.90 -0.88 -0.86 -0.84 -0.82 -0.80 -0.78
     [13] -0.76 -0.74 -0.72 -0.70 -0.68 -0.66 -0.64 -0.62 -0.60 -0.58 -0.56 -0.54
     [25] -0.52 -0.50 -0.48 -0.46 -0.44 -0.42 -0.40 -0.38 -0.36 -0.34 -0.32 -0.30
     [37] -0.28 -0.26 -0.24 -0.22 -0.20 -0.18 -0.16 -0.14 -0.12 -0.10 -0.08 -0.06
     [49] -0.04 -0.02  0.00  0.02  0.04  0.06  0.08  0.10  0.12  0.14  0.16  0.18
     [61]  0.20  0.22  0.24  0.26  0.28  0.30  0.32  0.34  0.36  0.38  0.40  0.42
     [73]  0.44  0.46  0.48  0.50  0.52  0.54  0.56  0.58  0.60  0.62  0.64  0.66
     [85]  0.68  0.70  0.72  0.74  0.76  0.78  0.80  0.82  0.84  0.86  0.88  0.90
     [97]  0.92  0.94  0.96  0.98  1.00

    $directional_derivative
      [1] -1.564e-05 -1.712e-01 -3.252e-01 -4.628e-01 -5.850e-01 -6.926e-01
      [7] -7.862e-01 -8.667e-01 -9.348e-01 -9.913e-01 -1.037e+00 -1.072e+00
     [13] -1.098e+00 -1.115e+00 -1.123e+00 -1.125e+00 -1.119e+00 -1.106e+00
     [19] -1.088e+00 -1.065e+00 -1.037e+00 -1.005e+00 -9.686e-01 -9.295e-01
     [25] -8.878e-01 -8.437e-01 -7.979e-01 -7.507e-01 -7.025e-01 -6.538e-01
     [31] -6.048e-01 -5.559e-01 -5.076e-01 -4.600e-01 -4.136e-01 -3.685e-01
     [37] -3.251e-01 -2.836e-01 -2.442e-01 -2.072e-01 -1.728e-01 -1.410e-01
     [43] -1.122e-01 -8.644e-02 -6.384e-02 -4.452e-02 -2.858e-02 -1.611e-02
     [49] -7.157e-03 -1.768e-03  3.127e-05 -1.768e-03 -7.157e-03 -1.611e-02
     [55] -2.858e-02 -4.452e-02 -6.384e-02 -8.644e-02 -1.122e-01 -1.410e-01
     [61] -1.728e-01 -2.072e-01 -2.442e-01 -2.836e-01 -3.251e-01 -3.685e-01
     [67] -4.136e-01 -4.600e-01 -5.076e-01 -5.559e-01 -6.048e-01 -6.538e-01
     [73] -7.025e-01 -7.507e-01 -7.979e-01 -8.437e-01 -8.878e-01 -9.295e-01
     [79] -9.686e-01 -1.005e+00 -1.037e+00 -1.065e+00 -1.088e+00 -1.106e+00
     [85] -1.119e+00 -1.125e+00 -1.123e+00 -1.115e+00 -1.098e+00 -1.072e+00
     [91] -1.037e+00 -9.913e-01 -9.348e-01 -8.667e-01 -7.862e-01 -6.926e-01
     [97] -5.850e-01 -4.628e-01 -3.252e-01 -1.712e-01 -1.564e-05

    $support_points
    [1] -1  0  1

    $support_values
    [1] -1.564e-05  3.127e-05 -1.564e-05

    $max_violation
    [1] 3.127e-05

    $all_nonpositive
    [1] FALSE

    $support_equal_zero
    [1] FALSE

    $criterion
    [1] "D"

    $tol
    [1] 1e-06

    attr(,"class")
    [1] "cvx_equivalence"

``` r
plot_equivalence(eq_d)
```

![](README_files/figure-commonmark/equivalence%20for%20A-%20and%20D-optimalities-1.png)

``` r
# Aopt
aout <- calc_Aopt(u, quad_reg, drop_tol = 1e-4)
eq_a <- check_equivalence(aout, f = quad_reg)
plot_equivalence(eq_a)
```

![](README_files/figure-commonmark/equivalence%20for%20A-%20and%20D-optimalities-2.png)

``` r
# copt
cout <- calc_copt(u, quad_reg, cVec = c(1, 0.5, 0),
                  drop_tol = 1e-6)
eq_c <- check_equivalence(cout, f = quad_reg)
plot_equivalence(eq_c)
```

![](README_files/figure-commonmark/copt%20equivalence-1.png)

## For non-linear models

For non-linear models, the information matrix depends on the unknown
parameters. In this case, the package supports local optimal design
computation at a specified nominal parameter value. For instance, the
peleg model is given as $$
y_i = y_0 + \frac{x_i}{\theta_1 + \theta_2 x_i } + \varepsilon_i.
$$

One needs to calculate the partial derivatives of the mean function with
respect to the parameters, which will be used in the information matrix
calculation. In this case, it is

$$
\frac{\partial g(x,\theta)}{ \partial \theta_1} = -\frac{x}{(\theta_1+\theta_2 x)^2}, \quad
\frac{\partial g(x,\theta)}{\partial \theta_2} = -\frac{x^2}{(\theta_1 + \theta_2 x)^2},
$$ that is

$$
\frac{\partial g(x,\theta)}{\partial \theta} = \frac{-1}{(\theta_1+\theta_2x)^2} \begin{pmatrix} x \\ x^2
\end{pmatrix}^\top.
$$

``` r
peleg_grad <- function(t, theta = c(0.5, 0.05)) {
  k1 <- theta[1]
  k2 <- theta[2]
  d2 <- (k1 + k2 * t)^2

  c(
    -t / d2,
    -t^2 / d2
  )
}

u <- 0:180

res <- calc_Dopt(
  u = u,
  f = peleg_grad,
  solver = "CLARABEL",
  verbose = TRUE,
  drop_tol = 1e-4
)

res$design |> round(3)
```

      point weight
    1     9    0.5
    2   180    0.5

``` r
res$value
```

    [1] 14.88

``` r
res$status
```

    [1] "optimal"

## For multiple objective design

### Maximin design

The maximin design is a robust design that maximizes the minimum
efficiency across a set of candidate parameter values or different
design criteria. For instance, one may want to find a design that
performs well across multiple nominal parameter values in a non-linear
model. This can be achieved by specifying a set of candidate parameter
values and using the `calc_maximin` function.

``` r
library(cvxDesign)

## regression function
quad_reg <- function(x) c(1, x, x^2)

## candidate set
u <- seq(-1, 1, length.out = 41)
```

### Step 1: single-objective designs

``` r
res_D <- compute_design_SO(
  u = u,
  f = quad_reg,
  criterion = "D"
)

res_A <- compute_design_SO(
  u = u,
  f = quad_reg,
  criterion = "A"
)

## reference losses
loss_ref <- list(
  D = res_D$loss,
  A = res_A$loss
)

print(loss_ref)
```

    $D
    [1] 1.91

    $A
    [1] 8

### Step 2: maximin design

``` r
res_DA <- compute_maximin_design(
  u = u,
  f = quad_reg,
  loss_ref = loss_ref,
  criteria = c("D", "A")
)
```

### Step 3: inspect result

``` r
print(res_DA$design)
```

    # A tibble: 3 × 2
      point weight
      <dbl>  <dbl>
    1    -1  0.285
    2     0  0.430
    3     1  0.285

``` r
print(res_DA$loss)
```

    $D
    [1] 1.969

    $A
    [1] 8.159

``` r
print(res_DA$efficiency)
```

         D      A 
    0.9805 0.9805 

``` r
cat("tstar =", res_DA$tstar, "\n")
```

    tstar = 1.02 

``` r
cat("1 / tstar =", 1 / res_DA$tstar, "\n")
```

    1 / tstar = 0.9805 

``` r
cat("min efficiency =", min(res_DA$efficiency), "\n")
```

    min efficiency = 0.9805 

### Step 4: numerical checks

``` r
tol <- 1e-4

eq1 <- abs(res_DA$efficiency["D"] - res_DA$efficiency["A"])
eq2 <- abs(min(res_DA$efficiency) - 1 / res_DA$tstar)

cat("|eff_D - eff_A| =", eq1, "\n")
```

    |eff_D - eff_A| = 8.238e-08 

``` r
cat("|min(eff) - 1/tstar| =", eq2, "\n")
```

    |min(eff) - 1/tstar| = 4.024e-08 

``` r
stopifnot(eq1 < tol)
stopifnot(eq2 < tol)

cat("All checks passed.\n")
```

    All checks passed.

### Equivalence theorem

``` r
## =========================================================
## Maximin multi-criterion design: D + A + c
## =========================================================

quad_reg <- function(x) c(1, x, x^2)
u <- seq(-1, 1, length.out = 41)

## target contrast for c-optimality
cvec <- c(0, 1, 0)

## ---------------------------------------------------------
## Step 1: single-objective reference designs
## ---------------------------------------------------------

res_D <- compute_design_SO(
  u = u,
  f = quad_reg,
  criterion = "D"
)

res_A <- compute_design_SO(
  u = u,
  f = quad_reg,
  criterion = "A"
)

res_c <- compute_design_SO(
  u = u,
  f = quad_reg,
  criterion = "c",
  opts = list(cVec_c = cvec)
)

## reference losses
loss_ref <- list(
  D = res_D$loss,
  A = res_A$loss,
  c = res_c$loss
)

print(loss_ref)
```

    $D
    [1] 1.91

    $A
    [1] 8

    $c
    [1] 1

``` r
## ---------------------------------------------------------
## Step 2: maximin design
## ---------------------------------------------------------

res_DAc <- compute_maximin_design(
  u = u,
  f = quad_reg,
  loss_ref = loss_ref,
  criteria = c("D", "A", "c"),
  opts = list(cVec_c = cvec)
)

print(res_DAc$design)
```

    # A tibble: 3 × 2
      point weight
      <dbl>  <dbl>
    1    -1  0.375
    2     0  0.250
    3     1  0.375

``` r
print(res_DAc$loss)
```

    $D
    [1] 1.962

    $A
    [1] 10.67

    $c
    [1] 1.333

``` r
print(res_DAc$efficiency)
```

         D      A      c 
    0.9828 0.7500 0.7500 

``` r
cat("tstar =", res_DAc$tstar, "\n")
```

    tstar = 1.333 

``` r
cat("1 / tstar =", 1 / res_DAc$tstar, "\n")
```

    1 / tstar = 0.75 

``` r
cat("min efficiency =", min(res_DAc$efficiency), "\n")
```

    min efficiency = 0.75 

### Step 3: directional derivatives

``` r
dd_DAc <- calc_directional_derivatives(
  u = u,
  M = res_DAc$info_matrix,
  f = quad_reg,
  criteria = c("D", "A", "c"),
  cVec = cvec
)
```

### Step 4: eta weights for maximin equivalence theorem

``` r
eta_DAc <- calc_eta_weights_maximin(
  tstar = res_DAc$tstar,
  loss_ref = loss_ref,
  loss_model = res_DAc$loss,
  directional_derivatives = dd_DAc,
  criteria = c("D", "A", "c"),
  q = length(quad_reg(0))
)

print(eta_DAc)
```

            D         A         c 
    1.231e-06 4.167e-02 6.667e-01 

### Step 5: equivalence theorem check

``` r
eq_DAc <- check_equivalence_maximin(
  design_obj = res_DAc,
  directional_derivatives = dd_DAc,
  eta = eta_DAc
)

print(eq_DAc)
```

    $candidate_points
     [1] -1.00 -0.95 -0.90 -0.85 -0.80 -0.75 -0.70 -0.65 -0.60 -0.55 -0.50 -0.45
    [13] -0.40 -0.35 -0.30 -0.25 -0.20 -0.15 -0.10 -0.05  0.00  0.05  0.10  0.15
    [25]  0.20  0.25  0.30  0.35  0.40  0.45  0.50  0.55  0.60  0.65  0.70  0.75
    [37]  0.80  0.85  0.90  0.95  1.00

    $support_points
    [1] -1  0  1

    $eta
            D         A         c 
    1.231e-06 4.167e-02 6.667e-01 

    $combined_directional_derivative
     [1] -3.330e-07 -1.630e-01 -2.850e-01 -3.713e-01 -4.267e-01 -4.557e-01
     [7] -4.628e-01 -4.518e-01 -4.267e-01 -3.907e-01 -3.472e-01 -2.991e-01
    [13] -2.489e-01 -1.991e-01 -1.517e-01 -1.085e-01 -7.111e-02 -4.073e-02
    [19] -1.833e-02 -4.617e-03  9.990e-07 -4.617e-03 -1.833e-02 -4.073e-02
    [25] -7.111e-02 -1.085e-01 -1.517e-01 -1.991e-01 -2.489e-01 -2.991e-01
    [31] -3.472e-01 -3.907e-01 -4.267e-01 -4.518e-01 -4.628e-01 -4.557e-01
    [37] -4.267e-01 -3.713e-01 -2.850e-01 -1.630e-01 -3.330e-07

    $support_values
    [1] -3.33e-07  9.99e-07 -3.33e-07

    $max_violation
    [1] 9.99e-07

    $min_value
    [1] -0.4628

    $all_nonpositive
    [1] TRUE

    $support_equal_zero
    [1] TRUE

    $tol
    [1] 1e-06

    attr(,"class")
    [1] "cvx_equivalence_maximin"

### Step 6: plot

``` r
plot_equivalence_maximin(
  design_obj = res_DAc,
  directional_derivatives = dd_DAc,
  eta = eta_DAc,
  criteria = c("D", "A", "c"),
  cex_lab = 1.5,
  cex_axis = 1.3,
  cex_main = 1.2,
  mar = c(5, 7, 3, 1)
)
```

![](README_files/figure-commonmark/step%206-1.png)

## Planned features

- [x] Basic package infrastructure
- [x] Core support for convex-optimization-based design computation
- [x] D-optimality
- [x] A-optimality
- [ ] E-optimality
- [x] c-optimality
- [x] Compound and multi-objective criteria
- [x] Equivalence theorem (sensitivity functions) diagnostics
- [ ] Exact design support on finite candidate sets
- [ ] Benchmark examples and vignettes
- [ ] Shiny app
- [ ] JSS paper companion materials

## Package vision

The long-term goal of `cvxDesign` is to provide a general computational
framework for modern optimal design problems, including:

- classical regression design,
- heteroscedastic models,
- constrained design problems,
- compound and robust criteria,
- group testing design,
- and other application-driven problems arising in statistics,
  biostatistics, and machine learning.

## References

1.  Atkinson, Anthony C., Donev, Aleksandar N., and Tobias, Randall D.
    (2007). *Optimum Experimental Designs, with SAS*. Oxford University
    Press.
2.  Fedorov, Valerii V. (1972). *Theory of Optimal Experiments*.
    Academic Press.
3.  Pukelsheim, Friedrich. (2006). *Optimal Design of Experiments*.
    SIAM.
4.  Silvey, Samuel D. (1980). *Optimal Design*. Chapman and Hall.
5.  Yeh, Chi-Kuang and collaborators. Ongoing work on
    convex-optimization-based methods for compound and
    application-driven optimal design problems.
