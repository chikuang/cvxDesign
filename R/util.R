#' Build the candidate regression matrix
#'
#' @param u Candidate design points.
#' @param f Regression function returning a numeric vector.
#'
#' @return A numeric matrix of dimension length(u) x p.
#' @noRd
build_model_matrix <- function(u, f) {
  validate_design_inputs(u, f)

  Fmat <- t(vapply(
    X = u,
    FUN = function(x) {
      out <- f(x)
      if (is.matrix(out) && ncol(out) == 1) {
        out <- as.vector(out)
      }
      as.numeric(out)
    },
    FUN.VALUE = as.numeric(f(u[1]))
  ))

  colnames(Fmat) <- paste0("beta", seq_len(ncol(Fmat)) - 1)
  Fmat
}


#' Validate inputs for optimal design problems
#'
#' @param u Candidate design points.
#' @param f Regression function returning a numeric vector.
#'
#' @return Invisibly returns TRUE if inputs are valid.
#' @noRd
validate_design_inputs <- function(u, f) {
  if (missing(u) || length(u) == 0) {
    stop("`u` must be a non-empty candidate set.", call. = FALSE)
  }

  if (!is.function(f)) {
    stop("`f` must be a function.", call. = FALSE)
  }

  test_val <- f(u[1])

  if (!is.numeric(test_val)) {
    stop("`f(u[i])` must return a numeric vector.", call. = FALSE)
  }

  if (is.matrix(test_val) && ncol(test_val) == 1) {
    test_val <- as.vector(test_val)
  }

  if (!is.vector(test_val)) {
    stop("`f(u[i])` must return a numeric vector.", call. = FALSE)
  }

  invisible(TRUE)
}



#' Compute the information matrix for a finite approximate design
#'
#' @param weights A numeric vector of weights.
#' @param Fmat Candidate regression matrix.
#'
#' @return Information matrix.
#' @noRd
info_matrix <- function(weights, Fmat) {
  if (length(weights) != nrow(Fmat)) {
    stop("Length of `weights` must match the number of candidate points.", call. = FALSE)
  }

  crossprod(Fmat, weights * Fmat)
}



#' Evaluate regression vector at one design point
#'
#' @param x A single design point.
#' @param f Regression function.
#'
#' @return Numeric column vector.
#' @noRd
eval_regvec <- function(x, f) {
  out <- f(x)
  if (is.matrix(out) && ncol(out) == 1) {
    out <- as.vector(out)
  }
  as.numeric(out)
}
