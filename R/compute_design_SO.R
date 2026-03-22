#' Compute the weighted information matrix
#'
#' @param u Candidate design points.
#' @param w Numeric vector of design weights.
#' @param f Regression function returning a numeric vector.
#' @param info_weight Optional function returning a nonnegative scalar weight.
#'
#' @return A numeric information matrix.
#' @noRd
fim_matrix <- function(u, w, f, info_weight = NULL) {
  if (is.null(info_weight)) {
    info_weight <- function(x) 1
  }

  terms <- lapply(seq_along(u), function(i) {
    fi <- matrix(eval_regvec(u[i], f), ncol = 1L)
    ai <- as.numeric(info_weight(u[i]))

    if (length(ai) != 1L || !is.finite(ai) || ai < 0) {
      stop("`info_weight(x)` must return a finite nonnegative scalar.", call. = FALSE)
    }

    as.numeric(w[i]) * ai * tcrossprod(fi)
  })

  Reduce(`+`, terms)
}


#' Build the weighted information matrix expression for CVXR
#'
#' @param u Candidate design points.
#' @param w CVXR variable of weights.
#' @param f Regression function returning a numeric vector.
#' @param info_weight Optional function returning a nonnegative scalar weight.
#'
#' @return A CVXR matrix expression.
#' @noRd
fim_matrix_expr <- function(u, w, f, info_weight = NULL) {
  if (is.null(info_weight)) {
    info_weight <- function(x) 1
  }

  terms <- lapply(seq_along(u), function(i) {
    fi <- matrix(eval_regvec(u[i], f), ncol = 1L)
    ai <- as.numeric(info_weight(u[i]))

    if (length(ai) != 1L || !is.finite(ai) || ai < 0) {
      stop("`info_weight(x)` must return a finite nonnegative scalar.", call. = FALSE)
    }

    w[i] * ai * tcrossprod(fi)
  })

  Reduce(`+`, terms)
}


#' Extract the contrast vector needed by a criterion
#'
#' @param criterion Criterion name.
#' @param opts Named list of options.
#' @param p Dimension of the regression vector.
#'
#' @return A numeric column vector or NULL.
#' @noRd
get_contrast_vec <- function(criterion, opts, p) {
  criterion <- toupper(criterion)

  if (criterion == "DS") {
    cvec <- opts$cVec_Ds
    nm <- "cVec_Ds"
  } else if (criterion == "C") {
    cvec <- opts$cVec_c
    nm <- "cVec_c"
  } else {
    return(NULL)
  }

  if (is.null(cvec)) {
    stop("For criterion `", criterion, "`, supply `opts$", nm, "`.", call. = FALSE)
  }

  if (!is.numeric(cvec) || length(cvec) != p) {
    stop("`opts$", nm, "` must be a numeric vector of length ", p, ".", call. = FALSE)
  }

  matrix(as.numeric(cvec), ncol = 1L)
}


#' Compute the scalar loss from an information matrix
#'
#' @param M Information matrix.
#' @param criterion One of "D", "A", "Ds", "c", "E".
#' @param opts Named list of options.
#'
#' @return A numeric scalar loss.
#' @noRd
scalar_loss_from_M <- function(M, criterion, opts = list()) {
  cr <- toupper(as.character(criterion))

  safe_solve <- function(A, b = NULL) {
    tryCatch(
      {
        if (is.null(b)) solve(A) else solve(A, b)
      },
      error = function(e) {
        if (!requireNamespace("MASS", quietly = TRUE)) {
          stop("Matrix inversion failed and package `MASS` is required for `ginv()`.",
               call. = FALSE)
        }
        G <- MASS::ginv(A)
        if (is.null(b)) G else G %*% b
      }
    )
  }

  if (cr == "D") {
    return(-as.numeric(determinant(M, logarithm = TRUE)$modulus))
  }

  if (cr == "A") {
    Minv <- safe_solve(M)
    return(sum(diag(Minv)))
  }

  if (cr == "DS") {
    cvec <- get_contrast_vec("Ds", opts, ncol(M))
    return(as.numeric(t(cvec) %*% safe_solve(M, cvec)))
  }

  if (cr == "C") {
    cvec <- get_contrast_vec("c", opts, ncol(M))
    return(as.numeric(t(cvec) %*% safe_solve(M, cvec)))
  }

  if (cr == "E") {
    evals <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    return(-min(evals))
  }

  stop("Unknown criterion: ", criterion, call. = FALSE)
}


#' Compute a single-objective optimal approximate design
#'
#' Solves a convex optimization problem over a finite candidate set, following
#' the same template as the MATLAB function `compute_design_SO()`.
#'
#' @param u Candidate design points.
#' @param f A function returning the regression vector at a single design point.
#' @param criterion One of `"D"`, `"A"`, `"Ds"`, `"c"`, or `"E"`.
#' @param opts Named list of extra options. For `"Ds"` use `opts$cVec_Ds`; for
#'   `"c"` use `opts$cVec_c`.
#' @param info_weight Optional function returning a nonnegative scalar multiplier
#'   for each rank-one information contribution.
#' @param solver Solver passed to [CVXR::psolve()].
#' @param ... Additional arguments passed to [CVXR::psolve()].
#' @param support_tol Weights smaller than this are dropped from the reported
#'   support.
#' @param drop_tol Numerical tolerance for tiny solver noise before support
#'   cleanup.
#'
#' @return A list with components:
#'   \item{criterion}{Criterion name.}
#'   \item{design}{A tibble of support points and weights.}
#'   \item{weights}{Full weight vector over the candidate set.}
#'   \item{candidates}{Candidate points.}
#'   \item{info_matrix}{Optimal information matrix.}
#'   \item{loss}{Scalar loss evaluated from the final information matrix.}
#'   \item{value}{Same as `loss`.}
#'   \item{optval}{Raw solver objective value.}
#'   \item{status}{Solver status.}
#'   \item{solver}{Solver name.}
#'
#' @export
compute_design_SO <- function(u,
                              f,
                              criterion = c("D", "A", "Ds", "c", "E"),
                              opts = list(),
                              info_weight = NULL,
                              solver = "CLARABEL",
                              ...,
                              support_tol = 1e-4,
                              drop_tol = 1e-8) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package `tibble` is required.", call. = FALSE)
  }

  criterion <- match.arg(criterion)
  cr <- toupper(criterion)

  validate_design_inputs(u, f)

  p <- length(eval_regvec(u[1L], f))
  n <- length(u)

  w <- CVXR::Variable(n, nonneg = TRUE)
  del <- CVXR::Variable(1)

  M_expr <- fim_matrix_expr(u, w, f, info_weight)

  constraints <- list(
    sum(w) == 1,
    w >= 0
  )

  if (cr == "D") {
    constraints <- c(constraints, list(-CVXR::log_det(M_expr) <= del))
  } else if (cr == "A") {
    constraints <- c(constraints, list(CVXR::tr_inv(M_expr) <= del))
  } else if (cr == "E") {
    constraints <- c(constraints, list(-CVXR::lambda_min(M_expr) <= del))
  } else if (cr %in% c("DS", "C")) {
    cvec <- get_contrast_vec(criterion, opts, p)
    constraints <- c(constraints, list(CVXR::matrix_frac(cvec, M_expr) <= del))
  } else {
    stop("Unsupported criterion: ", criterion, call. = FALSE)
  }

  problem <- CVXR::Problem(
    objective = CVXR::Minimize(del),
    constraints = constraints
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  status <- CVXR::status(problem)

  if (!status %in% c("optimal", "optimal_inaccurate")) {
    stop("Solver did not return an optimal solution. Status: ", status, call. = FALSE)
  }

  w_opt <- as.numeric(CVXR::value(w))

  if (any(!is.finite(w_opt))) {
    stop("Solver returned non-finite weights.", call. = FALSE)
  }
  if (sum(w_opt) <= 0) {
    stop("Solver returned invalid weights.", call. = FALSE)
  }

  # Step 1: remove tiny numerical noise
  w_opt[abs(w_opt) < drop_tol] <- 0

  # Step 2: remove practically negligible support weights
  w_opt[w_opt < support_tol] <- 0

  if (sum(w_opt) <= 0) {
    stop("All weights were removed by `drop_tol` and `support_tol`.", call. = FALSE)
  }

  # Step 3: renormalize
  w_opt <- w_opt / sum(w_opt)

  M_hat <- fim_matrix(u, w_opt, f, info_weight)
  loss <- scalar_loss_from_M(M_hat, criterion, opts)

  design <- tibble::tibble(
    point = u,
    weight = w_opt
  )

  design <- design[design$weight > 0, ]

  out <- list(
    criterion = criterion,
    design = design,
    weights = w_opt,
    candidates = u,
    info_matrix = M_hat,
    loss = loss,
    value = loss,
    optval = as.numeric(optval),
    status = status,
    solver = solver,
    support_tol = support_tol,
    drop_tol = drop_tol
  )

  class(out) <- c("cvx_so_design", "cvx_design")
  out
}
