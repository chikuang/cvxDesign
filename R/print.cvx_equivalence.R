#' @export
print.cvx_equivalence <- function(x, ...) {
  cat("Equivalence theorem check\n")
  cat("Criterion          :", x$criterion, "\n")
  cat("Tolerance          :", x$tol, "\n")
  cat("Max violation      :", signif(x$max_violation, 6), "\n")
  cat("All nonpositive    :", x$all_nonpositive, "\n")
  cat("Support equal zero :", x$support_equal_zero, "\n")
  invisible(x)
}
