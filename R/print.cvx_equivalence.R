#' @export
print.cvx_equivalence <- function(x, ...) {
  cat("Equivalence theorem check\n")
  if (inherits(x, "cvx_equivalence_general")) {
    cat("Mode               :", x$mode, "\n")
    if (!is.null(x$criteria)) {
      cat("Criteria           :", paste(x$criteria, collapse = ", "), "\n")
    }
    if (!is.null(x$eta)) {
      cat("Eta weights        :", paste(
        format(x$eta, digits = 4),
        collapse = ", "
      ), "\n")
    }
    if (!is.null(x$tstar)) {
      cat("t* (maximin)       :", signif(x$tstar, 6), "\n")
    }
  }
  cat("Criterion          :", x$criterion, "\n")
  cat("Tolerance          :", x$tol, "\n")
  cat("Max violation      :", signif(x$max_violation, 6), "\n")
  cat("All nonpositive    :", x$all_nonpositive, "\n")
  cat("Support equal zero :", x$support_equal_zero, "\n")
  invisible(x)
}
