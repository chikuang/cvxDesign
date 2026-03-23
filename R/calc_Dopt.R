calc_Dopt <- function(u,
                      f,
                      solver = "CLARABEL",
                      ...,
                      drop_tol = 1e-8,
                      ridge = 1e-8,
                      use_matrix_form = FALSE) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }

  validate_design_inputs(u, f)
  Fmat <- build_model_matrix(u, f)

  n <- nrow(Fmat)
  p <- ncol(Fmat)

  w <- CVXR::Variable(n, nonneg = TRUE)

  if (use_matrix_form) {
    F_cvx <- CVXR::Constant(Fmat)
    M_expr <- t(F_cvx) %*% CVXR::DiagVec(w) %*% F_cvx
  } else {
    M_expr <- Reduce(
      `+`,
      lapply(seq_len(n), function(i) {
        fi <- matrix(Fmat[i, ], ncol = 1)
        w[i] * (fi %*% t(fi))
      })
    )
  }

  problem <- CVXR::Problem(
    CVXR::Maximize(CVXR::log_det(M_expr + ridge * diag(p))),
    constraints = list(sum(w) == 1)
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  solver_status <- CVXR::status(problem)

  w_opt <- as.numeric(CVXR::value(w))
  w_opt[abs(w_opt) < drop_tol] <- 0
  w_opt <- w_opt / sum(w_opt)

  M_opt <- info_matrix(w_opt, Fmat)
  crit_val <- as.numeric(determinant(M_opt, logarithm = TRUE)$modulus)

  design <- data.frame(point = u, weight = w_opt)
  design <- design[design$weight > 0, , drop = FALSE]
  rownames(design) <- NULL

  out <- list(
    design = design,
    weights = w_opt,
    candidates = u,
    Fmat = Fmat,
    info_matrix = M_opt,
    criterion = "D",
    value = crit_val,
    optval = optval,
    status = solver_status,
    solver = solver,
    ridge = ridge
  )

  class(out) <- c("cvx_d_design", "cvx_design")
  out
}
