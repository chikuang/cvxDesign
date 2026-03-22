#' Compute eta weights for active criteria in a maximin design
#'
#' @param tstar Scalar maximin value.
#' @param loss_single Named list of single-objective losses.
#' @param loss_multi Named list or numeric vector of multi-objective losses.
#' @param criteria Character vector of criteria.
#' @param tol Tolerance for detecting active criteria.
#'
#' @return Numeric vector of eta weights summing to 1.
#' @export
calc_eta_weights <- function(tstar,
                             loss_single,
                             loss_multi,
                             criteria,
                             tol = 1e-6) {
  criteria_low <- tolower(criteria)

  get_loss <- function(x, nm) {
    if (is.list(x)) {
      x[[nm]]
    } else if (!is.null(names(x))) {
      x[[nm]]
    } else {
      NA_real_
    }
  }

  ratio_vec <- vapply(criteria_low, function(cr) {
    single_val <- get_loss(loss_single, toupper(cr))
    if (is.null(single_val) || is.na(single_val)) {
      single_val <- get_loss(loss_single, cr)
    }

    multi_val <- get_loss(loss_multi, toupper(cr))
    if (is.null(multi_val) || is.na(multi_val)) {
      multi_val <- get_loss(loss_multi, cr)
    }

    if (is.null(single_val) || is.na(single_val) ||
        is.null(multi_val) || is.na(multi_val)) {
      NA_real_
    } else {
      multi_val / single_val
    }
  }, numeric(1))

  active <- abs(ratio_vec - tstar) <= tol

  eta <- numeric(length(criteria_low))
  if (!any(active)) {
    eta[] <- 1 / length(criteria_low)
  } else {
    eta[active] <- 1 / sum(active)
  }

  names(eta) <- criteria_low
  eta
}
