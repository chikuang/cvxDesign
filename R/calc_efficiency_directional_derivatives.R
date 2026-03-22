#' Directional derivatives of D-, A-, and c-efficiencies (Kiefer--Wolfowitz kernel)
#'
#' Transforms the usual loss-based directional derivatives from
#' [calc_directional_derivatives()] into derivatives of **efficiencies**
#' \eqn{e_D, e_A, e_c} with respect to the same one-point perturbation of the
#' design. These are the objects that enter subgradient conditions for
#' **maximin efficiency** designs (maximize \eqn{\min_j e_j}), whereas raw
#' \eqn{d_D, d_A, d_c} are appropriate for **linear compound** loss criteria.
#'
#' Let \eqn{d_D^{\mathrm{kw}}}{dD_kw}, \eqn{d_A^{\mathrm{kw}}}{dA_kw},
#' \eqn{d_c^{\mathrm{kw}}}{dc_kw} denote the outputs of
#' [calc_directional_derivatives()]. With \eqn{p} parameters,
#' \eqn{T = \mathrm{tr}(M^{-1})}{T = tr(M^-1)}, \eqn{V = c'M^{-1}c}, reference
#' values from single-objective optima on the same candidate set, and current
#' efficiencies \eqn{e_D, e_A, e_c}:
#' \deqn{\frac{\mathrm{d} e_D}{\mathrm{d}\xi} = \frac{e_D}{p} \, d_D^{\mathrm{kw}}, \quad
#'       \frac{\mathrm{d} e_A}{\mathrm{d}\xi} = -\frac{e_A}{T} \, d_A^{\mathrm{kw}}, \quad
#'       \frac{\mathrm{d} e_c}{\mathrm{d}\xi} = -\frac{e_c}{V} \, d_c^{\mathrm{kw}}.}
#'
#' @param u Candidate design points.
#' @param M Information matrix at the design.
#' @param f Regression function.
#' @param criteria Character vector among `"D"`, `"A"`, `"c"`.
#' @param cVec Contrast for `"c"`.
#' @param efficiencies Named numeric vector with elements `D`, `A`, and/or `c`
#'   giving efficiencies at `M` (same normalization as [calc_maximin()]).
#' @param use_ginv Passed to [calc_directional_derivatives()].
#'
#' @return Named list with components `dD`, `dA`, and/or `dc` containing
#'   efficiency derivatives on `u` (same names as KW list for use with
#'   [calc_multi_directional_derivative()]).
#'
#' @export
calc_efficiency_directional_derivatives <- function(u,
                                                   M,
                                                   f,
                                                   criteria,
                                                   efficiencies,
                                                   cVec = NULL,
                                                   use_ginv = TRUE) {
  criteria <- unique(toupper(criteria))
  dd_kw <- calc_directional_derivatives(
    u, M, f,
    criteria = criteria,
    cVec = cVec,
    use_ginv = use_ginv
  )

  p <- ncol(M)
  tr_inv <- tryCatch(
    sum(diag(solve(M))),
    error = function(e) NA_real_
  )
  if (!is.finite(tr_inv) || tr_inv <= 0) {
    stop(
      "`M` must be nonsingular to form efficiency directional derivatives.",
      call. = FALSE
    )
  }

  Minv <- tryCatch(
    solve(M),
    error = function(e) {
      if (!use_ginv) {
        stop("`M` is singular and `use_ginv = FALSE`.", call. = FALSE)
      }
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package `MASS` is required for singular `M`.", call. = FALSE)
      }
      MASS::ginv(M)
    }
  )

  V <- NULL
  if ("C" %in% criteria) {
    if (is.null(cVec)) {
      stop("`cVec` is required when `criteria` includes \"c\".", call. = FALSE)
    }
    cv <- matrix(as.numeric(cVec), ncol = 1)
    V <- as.numeric(t(cv) %*% Minv %*% cv)
    if (!is.finite(V) || V <= 0) {
      stop("Variance c'M^{-1}c must be positive for c-efficiency derivatives.",
           call. = FALSE)
    }
  }

  out <- list()

  if ("D" %in% criteria) {
    if (is.null(efficiencies[["D"]])) {
      stop("`efficiencies` must contain element `D`.", call. = FALSE)
    }
    out$dD <- (efficiencies[["D"]] / p) * dd_kw$dD
  }
  if ("A" %in% criteria) {
    if (is.null(efficiencies[["A"]])) {
      stop("`efficiencies` must contain element `A`.", call. = FALSE)
    }
    out$dA <- -(efficiencies[["A"]] / tr_inv) * dd_kw$dA
  }
  if ("C" %in% criteria) {
    if (is.null(efficiencies[["c"]])) {
      stop("`efficiencies` must contain element `c`.", call. = FALSE)
    }
    out$dc <- -(efficiencies[["c"]] / V) * dd_kw$dc
  }

  out
}
