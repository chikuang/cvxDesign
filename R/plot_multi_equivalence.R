#' Plot directional derivatives for multi-objective design
#'
#' @param u Candidate design points.
#' @param dd_list Named list of directional derivative vectors.
#' @param d_multi Weighted combined directional derivative.
#' @param criteria Character vector of criteria.
#' @param main_prefix Optional prefix for subplot titles.
#'
#' @return Invisibly returns NULL.
#' @export
plot_multi_equivalence <- function(u,
                                   dd_list,
                                   d_multi,
                                   criteria,
                                   main_prefix = NULL) {
  criteria_low <- tolower(criteria)
  k <- length(criteria_low)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (k == 2) {
    par(mfrow = c(1, 3))
  } else if (k == 3) {
    par(mfrow = c(2, 2))
  } else {
    par(mfrow = c(k, 1))
  }

  label_map <- list(
    d = expression(d[D](u[i], wstar)),
    a = expression(d[A](u[i], wstar)),
    c = expression(d[c](u[i], wstar)),
    ds = expression(d[D[s]](u[i], wstar))
  )

  title_tag <- letters[seq_len(k + 1)]

  for (j in seq_along(criteria_low)) {
    cr <- criteria_low[j]
    dd_name <- switch(
      cr,
      d = "dD",
      a = "dA",
      c = "dc",
      ds = "dDs"
    )

    plot(u, dd_list[[dd_name]], type = "l", lwd = 2,
         xlab = expression(u[i]),
         ylab = label_map[[cr]],
         main = paste0("(", title_tag[j], ")"))
    abline(h = 0, lty = 2)
    grid()
  }

  plot(u, d_multi, type = "l", lwd = 2,
       xlab = expression(u[i]),
       ylab = expression(sum(eta[j] * d[j](u[i], w^"*"), j == 1, m)),
       main = paste0("(", title_tag[k + 1], ")"))
  abline(h = 0, lty = 2)
  grid()

  invisible(NULL)
}
