#' Standardize a criterion name
#'
#' @param cr A criterion name.
#'
#' @return A standardized criterion key.
#' @noRd
canon_crit_key <- function(cr) {
  key <- toupper(as.character(cr))

  if (key == "DS") {
    return("Ds")
  }
  if (key == "C") {
    return("c")
  }

  key
}


#' Standardize a vector of criterion names
#'
#' @param criteria A character vector of criterion names.
#'
#' @return A character vector of unique standardized criterion names.
#' @noRd
standardize_criteria <- function(criteria) {
  unique(vapply(criteria, canon_crit_key, character(1)))
}


#' Validate reference losses for maximin design
#'
#' @param loss_ref A named list of reference losses.
#' @param criteria A character vector of criteria.
#'
#' @return Invisibly returns TRUE.
#' @noRd
validate_loss_ref <- function(loss_ref, criteria) {
  if (!is.list(loss_ref)) {
    stop("`loss_ref` must be a named list.", call. = FALSE)
  }

  for (cr in criteria) {
    key <- canon_crit_key(cr)
    val <- loss_ref[[key]]

    if (is.null(val)) {
      stop(
        "`loss_ref$", key, "` is required when including criterion `", cr, "`.",
        call. = FALSE
      )
    }

    if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0) {
      stop(
        "`loss_ref$", key, "` must be a positive finite scalar.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


#' Compute criterion-specific losses at a given information matrix
#'
#' @param M An information matrix.
#' @param criteria A character vector of criteria.
#' @param opts A named list of criterion-specific options.
#'
#' @return A named list of losses.
#' @noRd
compute_losses_at_M <- function(M, criteria, opts = list()) {
  criteria <- standardize_criteria(criteria)

  out <- list()
  for (cr in criteria) {
    key <- canon_crit_key(cr)
    out[[key]] <- scalar_loss_from_M(M, cr, opts)
  }

  out
}


#' Compute maximin efficiencies from reference and achieved losses
#'
#' @param loss_ref A named list of reference losses.
#' @param loss A named list of achieved losses.
#' @param criteria A character vector of criteria.
#' @param q Parameter dimension.
#'
#' @return A named numeric vector of efficiencies.
#' @noRd
compute_efficiencies_maximin <- function(loss_ref, loss, criteria, q) {
  criteria <- standardize_criteria(criteria)

  eff <- numeric(0)

  for (cr in criteria) {
    key <- canon_crit_key(cr)
    ref_val <- loss_ref[[key]]
    cur_val <- loss[[key]]

    if (is.null(ref_val) || is.null(cur_val)) {
      next
    }

    if (key == "D") {
      eff["D"] <- exp(((-cur_val) - (-ref_val)) / q)
    } else {
      eff[key] <- ref_val / cur_val
    }
  }

  eff
}


#' Compute a maximin multi-criterion optimal approximate design
#'
#' Solves a joint convex optimization problem over a finite candidate set using
#' a maximin-efficiency formulation. The model is specified by `f(x)` and,
#' optionally, by a scalar weighting function `info_weight(x)` for each
#' rank-one contribution to the information matrix.
#'
#' Reference losses in `loss_ref` should be computed on the same scale as the
#' internal loss definition used by `compute_design_SO()`. For example, the
#' D-optimal loss is `-log det(M)`, the A-optimal loss is `tr(M^{-1})`, and the
#' c- and Ds-optimal losses are quadratic forms of the inverse information
#' matrix.
#'
#' @param u Candidate design points.
#' @param f A function returning the regression vector at a single design point.
#' @param loss_ref A named list of reference losses from single-objective
#'   designs.
#' @param criteria A character vector containing any of `"D"`, `"A"`, `"Ds"`,
#'   `"c"`, or `"E"`.
#' @param opts A named list of criterion-specific options. For `"Ds"`, supply
#'   `opts$cVec_Ds`; for `"c"`, supply `opts$cVec_c`.
#' @param info_weight An optional function returning a nonnegative scalar
#'   multiplier for each rank-one information contribution.
#' @param solver Solver passed to [CVXR::psolve()].
#' @param ... Additional arguments passed to [CVXR::psolve()].
#' @param support_tol Weights smaller than this are removed from the reported
#'   support.
#' @param drop_tol Numerical tolerance used to remove tiny solver noise before
#'   support cleanup.
#'
#' @return A list with components:
#' \item{criterion}{A character vector of included criteria.}
#' \item{design}{A tibble of support points and weights.}
#' \item{weights}{The full weight vector over the candidate set.}
#' \item{candidates}{The candidate design points.}
#' \item{loss}{A named list of achieved losses.}
#' \item{efficiency}{A named numeric vector of efficiencies.}
#' \item{value}{The minimum efficiency across the included criteria.}
#' \item{tstar}{The optimal value of the maximin variable.}
#' \item{M}{The optimal information matrix.}
#' \item{info_matrix}{Same as `M`.}
#' \item{optval}{The raw solver objective value.}
#' \item{status}{The solver status.}
#' \item{solver}{The solver name.}
#'
#' @export
compute_maximin_design <- function(u,
                                   f,
                                   loss_ref,
                                   criteria,
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

  validate_design_inputs(u, f)

  criteria <- standardize_criteria(criteria)
  validate_loss_ref(loss_ref, criteria)

  p <- length(eval_regvec(u[1L], f))
  n <- length(u)

  w <- CVXR::Variable(n, nonneg = TRUE)
  tstar <- CVXR::Variable(1)

  M_expr <- fim_matrix_expr(u, w, f, info_weight)

  constraints <- list(
    tstar >= 1e-8,
    w >= 0,
    sum(w) == 1
  )

  for (cr in criteria) {
    if (cr == "D") {
      constraints <- c(constraints, list(
        -CVXR::log_det(M_expr) <= loss_ref$D + p * log(tstar)
      ))
    } else if (cr == "A") {
      constraints <- c(constraints, list(
        CVXR::tr_inv(M_expr) <= tstar * loss_ref$A
      ))
    } else if (cr == "Ds") {
      c_ds <- get_contrast_vec("Ds", opts, p)
      constraints <- c(constraints, list(
        CVXR::matrix_frac(c_ds, M_expr) <= tstar * loss_ref$Ds
      ))
    } else if (cr == "c") {
      c_c <- get_contrast_vec("c", opts, p)
      constraints <- c(constraints, list(
        CVXR::matrix_frac(c_c, M_expr) <= tstar * loss_ref$c
      ))
    } else if (cr == "E") {
      constraints <- c(constraints, list(
        -CVXR::lambda_min(M_expr) <= loss_ref$E * CVXR::inv_pos(tstar)
      ))
    } else {
      stop("Unsupported criterion: ", cr, call. = FALSE)
    }
  }

  problem <- CVXR::Problem(
    objective = CVXR::Minimize(tstar),
    constraints = constraints
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  status <- CVXR::status(problem)

  if (!status %in% c("optimal", "optimal_inaccurate")) {
    stop(
      "Solver did not return an optimal solution. Status: ",
      status,
      call. = FALSE
    )
  }

  w_opt <- as.numeric(CVXR::value(w))
  t_opt <- as.numeric(CVXR::value(tstar))

  if (any(!is.finite(w_opt))) {
    stop("Solver returned non-finite weights.", call. = FALSE)
  }
  if (!is.finite(t_opt) || t_opt <= 0) {
    stop("Solver returned an invalid `tstar` value.", call. = FALSE)
  }
  if (sum(w_opt) <= 0) {
    stop("Solver returned invalid weights.", call. = FALSE)
  }

  w_opt[abs(w_opt) < drop_tol] <- 0
  w_opt[w_opt < support_tol] <- 0

  if (sum(w_opt) <= 0) {
    stop(
      "All weights were removed by `drop_tol` and `support_tol`.",
      call. = FALSE
    )
  }

  w_opt <- w_opt / sum(w_opt)

  M_hat <- fim_matrix(u, w_opt, f, info_weight)
  loss <- compute_losses_at_M(M_hat, criteria, opts)
  efficiency <- compute_efficiencies_maximin(loss_ref, loss, criteria, q = p)

  design <- tibble::tibble(
    point = u,
    weight = w_opt
  )
  design <- design[design$weight > 0, ]

  out <- list(
    criterion = criteria,
    design = design,
    weights = w_opt,
    candidates = u,
    loss = loss,
    efficiency = efficiency,
    value = min(efficiency, na.rm = TRUE),
    tstar = t_opt,
    M = M_hat,
    info_matrix = M_hat,
    optval = as.numeric(optval),
    status = status,
    solver = solver,
    support_tol = support_tol,
    drop_tol = drop_tol
  )

  class(out) <- c("cvx_maximin_design", "cvx_design")
  out
}
