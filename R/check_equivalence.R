#' Check the equivalence theorem on a finite candidate set
#'
#' @param design_obj Output from calc_Dopt(), calc_Aopt(), or calc_copt().
#' @param f Regression function returning a numeric vector.
#' @param u Candidate design points. If omitted, uses design_obj$candidates.
#' @param tol Numerical tolerance for checking nonpositivity and equality.
#'
#' @return A list with directional derivative values and theorem checks.
#' @export
check_equivalence <- function(design_obj,
                              f,
                              u = NULL,
                              tol = 1e-6) {
  if (missing(design_obj) || is.null(design_obj$criterion)) {
    stop("`design_obj` must be a valid design object from cvxDesign.", call. = FALSE)
  }

  if (missing(f) || !is.function(f)) {
    stop("`f` must be a regression function.", call. = FALSE)
  }

  if (is.null(u)) {
    u <- design_obj$candidates
  }

  M <- design_obj$info_matrix
  Minv <- solve(M)
  criterion <- tolower(design_obj$criterion)

  deriv_vals <- vapply(u, function(x) {
    fx <- eval_regvec(x, f)

    if (criterion == "d") {
      p <- length(fx)
      as.numeric(t(fx) %*% Minv %*% fx - p)

    } else if (criterion == "a") {
      as.numeric(t(fx) %*% Minv %*% Minv %*% fx - sum(diag(Minv)))

    } else if (criterion == "c") {
      cVec <- design_obj$cVec
      as.numeric((t(fx) %*% Minv %*% cVec)^2 - as.numeric(t(cVec) %*% Minv %*% cVec))

    } else {
      stop("Equivalence theorem currently supports criteria 'D', 'A', and 'c'.",
           call. = FALSE)
    }
  }, numeric(1))

  support_pts <- design_obj$design$point
  support_idx <- match(support_pts, u)

  max_violation <- max(deriv_vals, na.rm = TRUE)
  support_vals <- deriv_vals[support_idx]

  out <- list(
    candidate_points = u,
    directional_derivative = deriv_vals,
    support_points = support_pts,
    support_values = support_vals,
    max_violation = max_violation,
    all_nonpositive = all(deriv_vals <= tol),
    support_equal_zero = all(abs(support_vals) <= max(10 * tol, tol)),
    criterion = design_obj$criterion,
    tol = tol
  )

  class(out) <- "cvx_equivalence"
  out
}
