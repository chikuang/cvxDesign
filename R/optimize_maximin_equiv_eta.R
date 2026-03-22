#' Weights that minimize the worst-case combined derivative on a grid
#'
#' Given efficiency-based directional derivatives \eqn{de_j(u_i)} (from
#' [calc_efficiency_directional_derivatives()]) on a candidate grid \eqn{u_i}, a
#' necessary finite-dimensional check for a **maximin efficiency** design is
#' that there exist convex weights \eqn{\eta} with \eqn{\sum_j \eta_j de_j(u_i)
#' \le 0}{sum eta de_j <= 0} for all \eqn{i}. Heuristic rules (binding
#' efficiencies, uniform) need not achieve this even when such an \eqn{\eta}
#' exists.
#'
#' This function solves
#' \deqn{\min_{\eta \in \Delta} \max_i \sum_j \eta_j \, de_j(u_i)}
#' over the simplex \eqn{\Delta} (grid search for \eqn{k \le 3}{k <= 3}).
#'
#' @param dd_list Output of [calc_efficiency_directional_derivatives()] (same
#'   names `dD`, `dA`, `dc`).
#' @param criteria Character vector of criteria in the same order as used for
#'   `dd_list` (e.g. `c("D", "A")`).
#' @param step Grid step on the simplex for \eqn{k = 3} (ignored for \eqn{k = 2}).
#'
#' @return A list with `eta` (named `d`, `a`, `c` style), `min_max` (the
#'   minimized maximum), `combined` (vector on the grid), and `convex_combination`
#'   (same as `combined`).
#'
#' @export
optimize_maximin_equiv_eta <- function(dd_list,
                                       criteria,
                                       step = 0.02) {
  criteria <- unique(toupper(criteria))
  k <- length(criteria)
  if (k < 2L) {
    stop("Need at least two criteria.", call. = FALSE)
  }

  get_col <- function(ct) {
    nm <- switch(ct,
      D = "dD",
      A = "dA",
      C = "dc",
      stop("Unknown criterion: ", ct, call. = FALSE)
    )
    if (is.null(dd_list[[nm]])) {
      stop("Missing `dd_list$", nm, "`.", call. = FALSE)
    }
    dd_list[[nm]]
  }

  D <- sapply(criteria, get_col)
  if (!is.matrix(D)) {
    D <- matrix(D, ncol = 1L)
  }
  n <- nrow(D)

  if (k == 2L) {
    best <- Inf
    best_al <- 0.5
    for (al in seq(0, 1, length.out = 801L)) {
      comb <- al * D[, 1L] + (1 - al) * D[, 2L]
      mx <- max(comb, na.rm = TRUE)
      if (mx < best) {
        best <- mx
        best_al <- al
      }
    }
    eta <- c(best_al, 1 - best_al)
    names(eta) <- tolower(criteria)
    comb <- best_al * D[, 1L] + (1 - best_al) * D[, 2L]
  } else if (k == 3L) {
    best <- Inf
    best_w <- rep(1 / 3, 3)
    for (a in seq(0, 1, by = step)) {
      for (b in seq(0, 1 - a, by = step)) {
        w <- c(a, b, 1 - a - b)
        comb <- as.vector(D %*% w)
        mx <- max(comb, na.rm = TRUE)
        if (mx < best) {
          best <- mx
          best_w <- w
        }
      }
    }
    eta <- best_w
    names(eta) <- tolower(criteria)
    comb <- as.vector(D %*% matrix(eta, ncol = 1))
  } else {
    stop("Only k = 2 or k = 3 criteria are supported for grid search.", call. = FALSE)
  }

  list(
    eta = eta,
    min_max = max(comb, na.rm = TRUE),
    combined = comb,
    convex_combination = comb
  )
}
