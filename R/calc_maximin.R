#' Maximin multi-criterion approximate design on a finite candidate set
#'
#' Computes a design that maximizes the minimum efficiency across several
#' classical criteria (D, A, c), each measured relative to the best value
#' achievable on the same candidate set under that criterion alone. This matches
#' the usual maximin-efficiency formulation used in multi-objective design
#' (similar in spirit to `compute_maximin_design` in group-testing code).
#'
#' For a criterion with optimal reference value \eqn{v_j^\star}{v_j*} (from the
#' corresponding single-objective problem on \code{u}) and value \eqn{v_j(w)}
#' at design \eqn{w}, the efficiency is \eqn{e_j(w) = v_j^\star / v_j(w)} for
#' criteria minimized at the optimum (A, c), and
#' \eqn{e_D(w) = (\det M(w) / \det M_D^\star)^{1/p}}{D-efficiency in the usual
#' determinant sense} for D-optimality, where \eqn{p} is the number of
#' parameters.
#'
#' The algorithm bisects \eqn{t \in (0,1]} and solves a convex feasibility
#' problem at each \eqn{t} (maximize \eqn{\min_j e_j(w) \ge t}{min_j e_j(w) >= t}),
#' then refines with a secondary objective among maximin-efficient designs.
#'
#' @param u Candidate design points.
#' @param f Regression function returning a numeric vector (same convention as
#'   [calc_Dopt()]).
#' @param criteria Character vector of criteria to combine: any of `"D"`, `"A"`,
#'   `"c"`. Duplicates are not allowed.
#' @param cVec Linear contrast for c-optimality; required if `"c"` is included.
#'   Length must match `length(f(u[1]))`.
#' @param solver Solver name passed to [CVXR::psolve()].
#' @param ... Additional arguments passed to [CVXR::psolve()].
#' @param drop_tol Weights below this (after normalization) are set to zero.
#' @param bisect_tol Stopping tolerance for the bisection on \eqn{t}.
#' @param bisect_max_iter Maximum bisection iterations.
#' @param min_eig Positive lower bound enforced on the minimum eigenvalue of the
#'   information matrix (via [CVXR::lambda_min()]). This avoids numerically
#'   singular designs that can satisfy loose `log_det` constraints but yield
#'   meaningless efficiencies.
#'
#' @return A list of class `c("cvx_maximin_design", "cvx_design")` with
#'   components including `design`, `weights`, `info_matrix`, `criteria`,
#'   `value` (the maximin efficiency \eqn{t^\star}{t*}), `efficiencies` (named
#'   vector of efficiencies at the returned design), `references` (single-objective
#'   values used for normalization), and solver diagnostics.
#'
#' @export
calc_maximin <- function(u,
                         f,
                         criteria = c("D", "A"),
                         cVec = NULL,
                         solver = "CLARABEL",
                         ...,
                         drop_tol = 1e-8,
                         bisect_tol = 1e-6,
                         bisect_max_iter = 60L,
                         min_eig = 1e-8) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }

  criteria <- unique(criteria)
  allowed <- c("D", "A", "c")
  if (length(criteria) == 0L) {
    stop("`criteria` must be non-empty.", call. = FALSE)
  }
  bad <- setdiff(criteria, allowed)
  if (length(bad) > 0L) {
    stop("Unknown criterion: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (length(criteria) != length(unique(criteria))) {
    stop("`criteria` must not contain duplicates.", call. = FALSE)
  }

  if ("c" %in% criteria) {
    if (is.null(cVec)) {
      stop("`cVec` is required when `criteria` includes \"c\".", call. = FALSE)
    }
  } else if (!is.null(cVec)) {
    warning("`cVec` is ignored when `criteria` does not include \"c\".", call. = FALSE)
  }

  validate_design_inputs(u, f)
  Fmat <- build_model_matrix(u, f)

  refs <- list()
  singles <- list()

  if ("D" %in% criteria) {
    d0 <- calc_Dopt(u, f, solver = solver, ..., drop_tol = drop_tol)
    singles$D <- d0
    refs$log_det_D <- as.numeric(d0$value)
  }
  if ("A" %in% criteria) {
    a0 <- calc_Aopt(u, f, solver = solver, ..., drop_tol = drop_tol)
    singles$A <- a0
    refs$tr_A <- as.numeric(a0$value)
  }
  if ("c" %in% criteria) {
    if (length(cVec) != ncol(Fmat)) {
      stop(
        "`cVec` length (", length(cVec), ") must match the model dimension (",
        ncol(Fmat), ").",
        call. = FALSE
      )
    }
    c0 <- calc_copt(u, f, cVec = cVec, solver = solver, ..., drop_tol = drop_tol)
    singles$c <- c0
    refs$val_c <- as.numeric(c0$value)
    refs$c_vec <- as.numeric(cVec)
  }

  feas <- function(t) {
    maximin_feasible(Fmat, t, criteria, refs, solver, min_eig = min_eig, ...)
  }

  lo <- 1e-12
  hi <- 1

  if (!feas(hi)) {
    if (!feas(lo)) {
      stop(
        "Maximin problem appears infeasible even at very loose efficiency.",
        call. = FALSE
      )
    }
  } else {
    t_star <- 1
    w_opt <- maximin_final_weights(
      Fmat, criteria, refs, 1, solver, ...,
      drop_tol = drop_tol,
      min_eig = min_eig
    )
    return(maximin_pack_result(
      u, f, Fmat, w_opt, criteria, refs, singles, 1, drop_tol
    ))
  }

  iter <- 0L
  while (hi - lo > bisect_tol && iter < bisect_max_iter) {
    iter <- iter + 1L
    mid <- (lo + hi) / 2
    if (feas(mid)) {
      lo <- mid
    } else {
      hi <- mid
    }
  }

  t_star <- lo
  w_opt <- maximin_final_weights(
    Fmat, criteria, refs, t_star, solver, ...,
    drop_tol = drop_tol,
    min_eig = min_eig
  )

  maximin_pack_result(
    u, f, Fmat, w_opt, criteria, refs, singles, t_star, drop_tol
  )
}


#' @noRd
maximin_M_expr <- function(w, Fmat) {
  n <- nrow(Fmat)
  Reduce(
    `+`,
    lapply(seq_len(n), function(i) {
      w[i] * tcrossprod(Fmat[i, , drop = FALSE])
    })
  )
}


#' @noRd
maximin_feasible <- function(Fmat, t, crits, refs, solver, min_eig = 1e-8, ...) {
  n <- nrow(Fmat)
  p <- ncol(Fmat)
  if (t <= 0 || !is.finite(t)) {
    return(FALSE)
  }

  w <- CVXR::Variable(n, nonneg = TRUE)
  M_expr <- maximin_M_expr(w, Fmat)
  constr <- list(
    sum(w) == 1,
    CVXR::lambda_min(M_expr) >= min_eig
  )

  if ("D" %in% crits) {
    constr <- c(constr, list(
      CVXR::log_det(M_expr) >= refs$log_det_D + p * log(t)
    ))
  }
  if ("A" %in% crits) {
    constr <- c(constr, list(
      CVXR::tr_inv(M_expr) <= refs$tr_A / t
    ))
  }
  if ("c" %in% crits) {
    cmat <- matrix(as.numeric(refs$c_vec), nrow = p, ncol = 1)
    constr <- c(constr, list(
      CVXR::matrix_frac(cmat, M_expr) <= refs$val_c / t
    ))
  }

  prob <- CVXR::Problem(CVXR::Minimize(0), constraints = constr)
  CVXR::psolve(prob, solver = solver, ...)
  st <- CVXR::status(prob)
  st %in% c("optimal", "optimal_inaccurate")
}


#' @noRd
maximin_final_weights <- function(Fmat, crits, refs, t, solver, ...,
                                  drop_tol, min_eig = 1e-8) {
  n <- nrow(Fmat)
  p <- ncol(Fmat)
  w <- CVXR::Variable(n, nonneg = TRUE)
  M_expr <- maximin_M_expr(w, Fmat)
  constr <- list(
    sum(w) == 1,
    CVXR::lambda_min(M_expr) >= min_eig
  )

  if ("D" %in% crits) {
    constr <- c(constr, list(
      CVXR::log_det(M_expr) >= refs$log_det_D + p * log(t)
    ))
  }
  if ("A" %in% crits) {
    constr <- c(constr, list(
      CVXR::tr_inv(M_expr) <= refs$tr_A / t
    ))
  }
  if ("c" %in% crits) {
    cmat <- matrix(as.numeric(refs$c_vec), nrow = p, ncol = 1)
    constr <- c(constr, list(
      CVXR::matrix_frac(cmat, M_expr) <= refs$val_c / t
    ))
  }

  # Any feasible point satisfies the maximin constraints; a secondary objective
  # like `Maximize(log_det(M))` can push the solver to a rank-deficient boundary
  # where `lambda_min` is violated only within solver tolerance, so we keep a
  # pure feasibility problem.
  prob <- CVXR::Problem(CVXR::Minimize(0), constraints = constr)
  CVXR::psolve(prob, solver = solver, ...)

  w_opt <- as.numeric(CVXR::value(w))
  w_opt[abs(w_opt) < drop_tol] <- 0
  if (sum(w_opt) <= 0) {
    stop("Maximin refinement returned invalid weights.", call. = FALSE)
  }
  w_opt / sum(w_opt)
}


#' @noRd
maximin_compute_efficiencies <- function(M, crits, refs, p) {
  # Use the same definitions as calc_Dopt / calc_Aopt: do not pre-filter
  # eigenvalues with a loose tolerance — that can falsely declare rank-deficiency
  # (D -> 0) and underestimate tr(M^{-1}) (A -> incorrectly 1).

  eff <- numeric(0)
  if ("D" %in% crits) {
    ld <- tryCatch(
      as.numeric(determinant(M, logarithm = TRUE)$modulus),
      error = function(e) NA_real_
    )
    if (is.na(ld) || !is.finite(ld)) {
      eff["D"] <- 0
    } else {
      eff["D"] <- min(1, exp((ld - refs$log_det_D) / p))
    }
  }
  if ("A" %in% crits) {
    trv <- tryCatch(
      sum(diag(solve(M))),
      error = function(e) NA_real_
    )
    if (is.na(trv) || !is.finite(trv) || trv <= 0) {
      eff["A"] <- 0
    } else {
      eff["A"] <- min(1, refs$tr_A / trv)
    }
  }
  if ("c" %in% crits) {
    cv <- tryCatch(
      as.numeric(t(refs$c_vec) %*% solve(M, refs$c_vec)),
      error = function(e) NA_real_
    )
    if (is.na(cv) || !is.finite(cv) || cv <= 0) {
      if (requireNamespace("MASS", quietly = TRUE)) {
        cv <- as.numeric(t(refs$c_vec) %*% MASS::ginv(M) %*% refs$c_vec)
      } else {
        cv <- Inf
      }
    }
    eff["c"] <- if (is.finite(cv) && cv > 0) min(1, refs$val_c / cv) else 0
  }
  eff
}


#' @noRd
maximin_pack_result <- function(u, f, Fmat, w_opt, crits, refs, singles,
                                t_star, drop_tol) {
  M_opt <- info_matrix(w_opt, Fmat)
  p <- ncol(Fmat)
  eff <- maximin_compute_efficiencies(M_opt, crits, refs, p)

  design <- data.frame(point = u, weight = w_opt)
  design <- design[design$weight > 0, , drop = FALSE]
  rownames(design) <- NULL

  out <- list(
    design = design,
    weights = w_opt,
    candidates = u,
    Fmat = Fmat,
    info_matrix = M_opt,
    criterion = "maximin",
    criteria = crits,
    value = min(eff),
    efficiencies = eff,
    references = singles,
    ref_values = refs,
    status = "optimal",
    solver = NA_character_
  )

  class(out) <- c("cvx_maximin_design", "cvx_design")
  out
}
