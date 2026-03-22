# cvxDesign: Optimal experimental designs via convex optimization

2026-03-22

<!-- badges: start -->

<!-- [![R-CMD-check](https://github.com/chikuang/cvxDesign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chikuang/cvxDesign/actions/workflows/R-CMD-check.yaml) -->

[![](https://img.shields.io/github/languages/code-size/chikuang/cvxDesign.svg)](https://github.com/chikuang/cvxDesign)
<!-- badges: end -->

**Authors:**  
Chi-Kuang Yeh (Georgia State University)  
[![ORCID](https://img.shields.io/badge/ORCID-0000--0001--7057--2096-A6CE39?logo=orcid.png)](https://orcid.org/0000-0001-7057-2096)

Julie Zhou (University of Victoria)

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

    [1] -1.909543

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

    Equivalence theorem check
    Criterion          : D 
    Tolerance          : 1e-06 
    Max violation      : 3.94291e-05 
    All nonpositive    : FALSE 
    Support equal zero : FALSE 

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
y_i = \frac{\theta_1 x_i}{\theta_2 + x_i } + \varepsilon_i.
$$

One needs to calculate the partial derivatives of the mean function with
respect to the parameters, which will be used in the information matrix
calculation. In this case, it is $$
\frac{\partial \eta}{\partial \theta_1} = \frac{x}{\theta_2 + x}, \quad
\frac{\partial \eta}{\partial \theta_2} = -\frac{\theta_1 x}{(\theta_2 + x)^2}.
$$

``` r
peleg_reg <- function(x, theta) {
    d1 <- x / (theta[2] + x)
    d2 <- -theta[1] * x / (theta[2] + x)^2
    c(d1, d2)
}
u <- seq(0.1, 100, length.out = 1001)
theta_nom <- c(0.5, 0.05)
dout <- calc_Dopt(u, 
                  function(x) peleg_reg(x, theta = theta_nom), 
                  drop_tol = 1e-4)
dout$design |>
  dplyr::filter(.data$weight > 1e-3) |>
  round(3)
```

      point weight
    1   0.1    0.5
    2 100.0    0.5

``` r
dout$value
```

    [1] 0.2067205

``` r
dout$status
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

## -------------------------------
## Step 1: single-objective designs
## -------------------------------
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
    [1] 1.909543

    $A
    [1] 8

``` r
## -------------------------------
## Step 2: maximin design
## -------------------------------
res_DA <- compute_maximin_design(
  u = u,
  f = quad_reg,
  loss_ref = loss_ref,
  criteria = c("D", "A")
)

## -------------------------------
## Step 3: inspect result
## -------------------------------
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
    [1] 1.968502

    $A
    [1] 8.158781

``` r
print(res_DA$efficiency)
```

            D         A 
    0.9805388 0.9805387 

``` r
cat("tstar =", res_DA$tstar, "\n")
```

    tstar = 1.019848 

``` r
cat("1 / tstar =", 1 / res_DA$tstar, "\n")
```

    1 / tstar = 0.9805387 

``` r
cat("min efficiency =", min(res_DA$efficiency), "\n")
```

    min efficiency = 0.9805387 

``` r
## -------------------------------
## Step 4: numerical checks
## -------------------------------
tol <- 1e-4

eq1 <- abs(res_DA$efficiency["D"] - res_DA$efficiency["A"])
eq2 <- abs(min(res_DA$efficiency) - 1 / res_DA$tstar)

cat("|eff_D - eff_A| =", eq1, "\n")
```

    |eff_D - eff_A| = 8.237729e-08 

``` r
cat("|min(eff) - 1/tstar| =", eq2, "\n")
```

    |min(eff) - 1/tstar| = 4.024065e-08 

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
    [1] 1.909543

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
    [1] 1.961659

    $A
    [1] 10.66667

    $c
    [1] 1.333333

``` r
print(res_DAc$efficiency)
```

            D         A         c 
    0.9827780 0.7499999 0.7500000 

``` r
cat("tstar =", res_DAc$tstar, "\n")
```

    tstar = 1.333333 

``` r
cat("1 / tstar =", 1 / res_DAc$tstar, "\n")
```

    1 / tstar = 0.75 

``` r
cat("min efficiency =", min(res_DAc$efficiency), "\n")
```

    min efficiency = 0.7499999 

``` r
## ---------------------------------------------------------
## Step 3: directional derivatives
## ---------------------------------------------------------

dd_DAc <- calc_directional_derivatives(
  u = u,
  M = res_DAc$info_matrix,
  f = quad_reg,
  criteria = c("D", "A", "c"),
  cVec = cvec
)

## ---------------------------------------------------------
## Step 4: eta weights for maximin equivalence theorem
## ---------------------------------------------------------

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

               D            A            c 
    1.231246e-06 4.166653e-02 6.666650e-01 

``` r
## ---------------------------------------------------------
## Step 5: equivalence theorem check
## ---------------------------------------------------------

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
               D            A            c 
    1.231246e-06 4.166653e-02 6.666650e-01 

    $combined_directional_derivative
     [1] -3.329950e-07 -1.629517e-01 -2.850002e-01 -3.712849e-01 -4.266667e-01
     [6] -4.557292e-01 -4.627777e-01 -4.518401e-01 -4.266664e-01 -3.907288e-01
    [11] -3.472217e-01 -2.990619e-01 -2.488882e-01 -1.990618e-01 -1.516659e-01
    [16] -1.085061e-01 -7.111020e-02 -4.072822e-02 -1.833236e-02 -4.617062e-03
    [21]  9.989852e-07 -4.617062e-03 -1.833236e-02 -4.072822e-02 -7.111020e-02
    [26] -1.085061e-01 -1.516659e-01 -1.990618e-01 -2.488882e-01 -2.990619e-01
    [31] -3.472217e-01 -3.907288e-01 -4.266664e-01 -4.518401e-01 -4.627777e-01
    [36] -4.557292e-01 -4.266667e-01 -3.712849e-01 -2.850002e-01 -1.629517e-01
    [41] -3.329950e-07

    $support_values
    [1] -3.329950e-07  9.989852e-07 -3.329950e-07

    $max_violation
    [1] 9.989852e-07

    $min_value
    [1] -0.4627777

    $all_nonpositive
    [1] TRUE

    $support_equal_zero
    [1] TRUE

    $tol
    [1] 1e-06

    attr(,"class")
    [1] "cvx_equivalence_maximin"

``` r
## ---------------------------------------------------------
## Step 6: plot
## ---------------------------------------------------------

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

![](README_files/figure-commonmark/multi-equivalence%20theorem-1.png)

## Planned features

- [x] Basic package infrastructure
- [x] Core support for convex-optimization-based design computation
- [x] D-optimality
- [x] A-optimality
- [ ] E-optimality
- [ ] c-optimality
- [ ] Compound and multi-objective criteria
- [ ] Equivalence theorem (sensitivity functions) diagnostics
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
