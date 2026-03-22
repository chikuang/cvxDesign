#' Equivalence theorem check for general criteria (single or maximin)
#'
#' For a **single** criterion (D, A, or c), checks that the directional derivative
#' is nonpositive on the candidate set and zero on the support of the optimal
#' design (standard equivalence theorem).
#'
#' For a **maximin** design ([calc_maximin()]), by default (`derivative =
#' "efficiency"`) the function uses [calc_efficiency_directional_derivatives()],
#' i.e. derivatives of the **efficiencies** \eqn{e_D, e_A, e_c} in the same
#' Kiefer--Wolfowitz sense, then forms \eqn{\sum_j \eta_j \, \mathrm{d} e_j /
#' \mathrm{d}\xi}{sum eta * de_j}. That matches the subgradient of the maximin
#' efficiency objective far better than weighting raw loss derivatives.
#'
#' Use `derivative = "kw_loss"` only if you want the old behaviour (derivatives
#' of \eqn{-\log\det M}{-log det M}, \eqn{\mathrm{tr}(M^{-1})}{tr(M^-1)}, etc.),
#' e.g. for comparison to linear compound designs.
#'
#' Weights \eqn{\eta} (see `eta_rule`):
#' * `"binding"`: uniform over criteria whose **efficiency** equals \eqn{t^\star}
#'   within `eta_tol`, zero elsewhere; if none match, uniform over all.
#' * `"uniform"`: \eqn{\eta_j = 1/k} for all \eqn{k} criteria (often clearer for
#'   comparing contributions in plots).
#' * `"optimize"`: choose \eqn{\eta} on the simplex to minimize
#'   \eqn{\max_u \sum_j \eta_j \, de_j(u)}{max sum eta de_j} over the candidate
#'   grid ([optimize_maximin_equiv_eta()]); best finite check when peaks remain
#'   slightly above zero with other rules.
#'
#' @param design_obj Output from [calc_Dopt()], [calc_Aopt()], [calc_copt()], or
#'   [calc_maximin()].
#' @param f Regression function (same as used for the design).
#' @param u Candidate design points; defaults to `design_obj$candidates`.
#' @param criteria Character vector of criteria to check. For single-objective
#'   designs, defaults to `design_obj$criterion`. For maximin, defaults to
#'   `design_obj$criteria`.
#' @param cVec Contrast for c-optimality; required when `"c"` is included and not
#'   stored on `design_obj` (for maximin, `ref_values$c_vec` is used if present).
#' @param eta Optional named numeric vector of weights (names like `"d"`, `"a"`,
#'   `"c"`). If `NULL`, maximin uses `eta_rule`; single-criterion uses weight 1.
#' @param eta_rule For maximin when `eta` is `NULL`: `"binding"`, `"uniform"`,
#'   or `"optimize"` (see Details).
#' @param tstar Maximin level; for maximin designs defaults to `design_obj$value`.
#' @param eta_tol Tolerance for matching each efficiency to `tstar` when forming
#'   automatic `eta` under `eta_rule = "binding"`.
#' @param tol Tolerance for nonpositivity and zero-on-support checks.
#' @param use_ginv Passed to [calc_directional_derivatives()].
#' @param derivative For **maximin** designs only: `"efficiency"` (default)
#'   uses [calc_efficiency_directional_derivatives()] so the combined curve
#'   reflects derivatives of \eqn{e_D, e_A, e_c}; `"kw_loss"` uses raw
#'   Kiefer--Wolfowitz derivatives of the underlying losses (diagnostic only for
#'   maximin efficiency).
#'
#' @return An object of class `c("cvx_equivalence_general", "cvx_equivalence")`
#'   with components including `directional_derivative` (the curve to plot:
#'   single criterion or combined), `directional_derivatives` (named list per
#'   criterion for maximin), `derivative` (`"efficiency"` or `"kw_loss"` for
#'   maximin), `eta` (weights, maximin only), `eta_optimization` (when
#'   `eta_rule = "optimize"`, output of [optimize_maximin_equiv_eta()]),
#'   `max_violation`, `all_nonpositive`, `support_equal_zero`, and `mode`
#'   (`"single"` or `"maximin"`).
#'
#' @seealso [check_equivalence()] for the original single-object API,
#'   [plot_equivalence()], [plot_multi_equivalence()].
#'
#' @examples
#' \donttest{
#' if (requireNamespace("CVXR", quietly = TRUE)) {
#'   ## Quadratic regression: f(x) = (1, x, x^2)'
#'   f <- function(x) c(1, x, x^2)
#'   u <- seq(-1, 1, length.out = 51)
#'
#'   ## --- Single-objective D-optimal design ---
#'   des <- calc_Dopt(u, f)
#'   eq <- check_equivalence_general(des, f)
#'   eq$max_violation
#'   eq$all_nonpositive
#'   eq$support_equal_zero
#'   plot_equivalence(eq)
#'
#'   ## Same model via compute_design_SO (MATLAB-style scalar del)
#'   des2 <- compute_design_SO(u, f, criterion = "D")
#'   eq2 <- check_equivalence_general(des2, f, criteria = "D")
#'
#'   ## --- Maximin (D + A): use efficiency derivatives (default) ---
#'   mm <- calc_maximin(u, f, criteria = c("D", "A"))
#'   eqm <- check_equivalence_general(mm, f, eta_rule = "uniform")
#'   plot_multi_equivalence(
#'     u,
#'     eqm$directional_derivatives,
#'     eqm$directional_derivative,
#'     criteria = mm$criteria
#'   )
#' }
#' }
#'
#' @section Calling the function:
#' If you see `unused argument (derivative = ...)` while `args(cvxDesign::check_equivalence_general)`
#' lists `derivative`, another object is masking the package function (e.g. an old
#' copy from `source()` or another session). Use the namespace-qualified call:
#' `cvxDesign::check_equivalence_general(...)` or run `find("check_equivalence_general")`
#' and remove the masking binding.
#'
#' @export
check_equivalence_general <- function(design_obj,
                                      f,
                                      u = NULL,
                                      criteria = NULL,
                                      cVec = NULL,
                                      eta = NULL,
                                      eta_rule = c("binding", "uniform", "optimize"),
                                      tstar = NULL,
                                      eta_tol = 1e-5,
                                      tol = 1e-6,
                                      use_ginv = TRUE,
                                      derivative = c("efficiency", "kw_loss")) {
  if (missing(design_obj) || is.null(design_obj$info_matrix)) {
    stop("`design_obj` must contain `info_matrix`.", call. = FALSE)
  }
  if (missing(f) || !is.function(f)) {
    stop("`f` must be a regression function.", call. = FALSE)
  }

  eta_opt_out <- NULL

  if (is.null(u)) {
    u <- design_obj$candidates
  }

  M <- design_obj$info_matrix
  maximin <- inherits(design_obj, "cvx_maximin_design") ||
    identical(design_obj$criterion, "maximin")

  if (maximin) {
    crits <- if (is.null(criteria)) design_obj$criteria else criteria
    if (is.null(crits) || length(crits) == 0L) {
      stop("Maximin `design_obj` must have `criteria`.", call. = FALSE)
    }
    crits <- unique(crits)

    if ("c" %in% crits && is.null(cVec)) {
      rv <- design_obj$ref_values
      if (!is.null(rv) && !is.null(rv$c_vec)) {
        cVec <- rv$c_vec
      } else if (!is.null(design_obj$references$c$cVec)) {
        cVec <- design_obj$references$c$cVec
      }
    }
    if ("c" %in% crits && is.null(cVec)) {
      stop("`cVec` is required for maximin designs with c-optimality.", call. = FALSE)
    }

    if (is.null(tstar)) {
      tstar <- design_obj$value
    }

    derivative <- match.arg(derivative)
    if (derivative == "efficiency") {
      dd_list <- calc_efficiency_directional_derivatives(
        u, M, f,
        criteria = crits,
        cVec = cVec,
        efficiencies = design_obj$efficiencies,
        use_ginv = use_ginv
      )
    } else {
      dd_list <- calc_directional_derivatives(
        u, M, f,
        criteria = crits,
        cVec = cVec,
        use_ginv = use_ginv
      )
    }

    if (is.null(eta)) {
      eta_rule <- match.arg(eta_rule)
      eta <- if (eta_rule == "uniform") {
        maximin_uniform_eta(crits)
      } else if (eta_rule == "optimize") {
        eta_opt_out <- optimize_maximin_equiv_eta(dd_list, crits)
        eta_opt_out$eta
      } else {
        maximin_default_eta(design_obj$efficiencies, crits, tstar, eta_tol)
      }
    } else {
      validate_eta_names(eta, crits)
    }

    combined <- calc_multi_directional_derivative(dd_list, eta)
    mode <- "maximin"
    crit_label <- paste0("maximin(", paste(crits, collapse = "+"), ")")
  } else {
    crits <- if (is.null(criteria)) design_obj$criterion else criteria
    if (length(crits) != 1L) {
      stop(
        "Single-objective `design_obj` needs exactly one `criteria` (or use `design_obj$criterion`).",
        call. = FALSE
      )
    }
    crits <- toupper(as.character(crits))

    if (crits == "C" && is.null(cVec)) {
      cVec <- design_obj$cVec
    }
    if (crits == "C" && is.null(cVec)) {
      stop("`cVec` is required for c-optimality.", call. = FALSE)
    }

    dd_list <- calc_directional_derivatives(
      u, M, f,
      criteria = crits,
      cVec = cVec,
      use_ginv = use_ginv
    )

    nm <- switch(
      tolower(crits),
      d = "dD",
      a = "dA",
      c = "dc",
      stop("Unknown criterion: ", crits, call. = FALSE)
    )
    combined <- dd_list[[nm]]
    eta <- NULL
    mode <- "single"
    crit_label <- crits
    derivative <- NA_character_
  }

  support_pts <- design_obj$design$point
  support_idx <- match(support_pts, u)
  if (anyNA(support_idx)) {
    support_idx <- match(support_pts, u, nomatch = 0L)
    support_idx <- support_idx[support_idx > 0L]
  }

  max_violation <- max(combined, na.rm = TRUE)
  support_vals <- combined[support_idx]
  all_ok <- all(combined <= tol, na.rm = TRUE)
  support_ok <- all(abs(support_vals) <= max(10 * tol, tol), na.rm = TRUE)

  out <- list(
    mode = mode,
    criterion = crit_label,
    criteria = if (maximin) crits else NULL,
    candidate_points = u,
    directional_derivatives = dd_list,
    eta = eta,
    directional_derivative = combined,
    derivative = if (maximin) derivative else NA_character_,
    support_points = support_pts,
    support_values = support_vals,
    max_violation = max_violation,
    all_nonpositive = all_ok,
    support_equal_zero = support_ok,
    tol = tol,
    tstar = if (maximin) tstar else NULL,
    eta_optimization = if (maximin && !is.null(eta_opt_out)) eta_opt_out else NULL
  )

  if (maximin && identical(derivative, "kw_loss") && max_violation > tol) {
    message(
      "Maximin design with kw_loss derivatives: ",
      "weighted curve may exceed 0; prefer derivative = \"efficiency\" for ",
      "maximin efficiency (see ?check_equivalence_general)."
    )
  }

  class(out) <- c("cvx_equivalence_general", "cvx_equivalence")
  out
}


#' @noRd
maximin_uniform_eta <- function(criteria) {
  n <- length(criteria)
  eta <- rep(1 / n, n)
  names(eta) <- tolower(criteria)
  eta
}


maximin_default_eta <- function(efficiencies, criteria, tstar, eta_tol) {
  eff <- efficiencies
  cr_low <- tolower(criteria)
  active <- vapply(criteria, function(ct) {
    key <- toupper(ct)
    if (!key %in% names(eff)) {
      stop("Efficiency `", key, "` not found in `efficiencies`.", call. = FALSE)
    }
    abs(eff[[key]] - tstar) <= eta_tol
  }, logical(1))

  n <- length(criteria)
  eta <- numeric(n)
  if (!any(active)) {
    eta[] <- 1 / n
  } else {
    eta[active] <- 1 / sum(active)
  }
  names(eta) <- cr_low
  eta
}


#' @noRd
validate_eta_names <- function(eta, criteria) {
  exp_n <- tolower(criteria)
  got <- names(eta)
  if (is.null(got)) {
    stop("`eta` must be named (d, a, c).", call. = FALSE)
  }
  if (!all(exp_n %in% tolower(got)) && !all(tolower(got) %in% exp_n)) {
    warning("`eta` names may not match `criteria`; verify weights.", call. = FALSE)
  }
  invisible(TRUE)
}
