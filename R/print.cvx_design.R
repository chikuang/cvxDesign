#' @export
print.cvx_so_design <- function(x, ...) {
  cat("Single-objective optimal design (compute_design_SO)\n")
  cat("Criterion          :", x$criterion, "\n")
  cat("Solver status      :", x$status, "\n")
  cat("Loss (value)       :", signif(x$loss, 6), "\n")
  cat("Support points     :", nrow(x$design), "\n")
  print(x$design, row.names = FALSE, ...)
  invisible(x)
}


#' @export
print.cvx_maximin_design <- function(x, ...) {
  if (!is.null(x$tstar) && !is.null(x$loss) && is.list(x$loss)) {
    cat("Joint maximin design (compute_maximin_design)\n")
    cat("Criteria           :", paste(x$criterion, collapse = ", "), "\n")
    cat("Solver status      :", x$status, "\n")
    cat("t*                 :", signif(as.numeric(x$tstar), 6), "\n")
    cat("Min efficiency     :", signif(x$value, 6), "\n")
    eff <- x$efficiency
    cat("Efficiencies       :", paste(
      names(eff), "=", signif(eff, 4),
      collapse = "; "
    ), "\n")
  } else {
    cat("Maximin approximate design (calc_maximin)\n")
    cat("Criteria           :", paste(x$criteria, collapse = ", "), "\n")
    cat("Min efficiency     :", signif(x$value, 6), "\n")
    eff <- x$efficiencies
    cat("Efficiencies       :", paste(
      names(eff), "=", signif(eff, 4),
      collapse = "; "
    ), "\n")
  }
  cat("Support points     :", nrow(x$design), "\n")
  print(x$design, row.names = FALSE, ...)
  invisible(x)
}


#' @export
print.cvx_d_design <- function(x, ...) {
  cat("D-optimal approximate design\n")
  cat("Solver status      :", x$status, "\n")
  cat("log det(M)         :", signif(x$value, 6), "\n")
  cat("Support points     :", nrow(x$design), "\n")
  print(x$design, row.names = FALSE, ...)
  invisible(x)
}


#' @export
print.cvx_a_design <- function(x, ...) {
  cat("A-optimal approximate design\n")
  cat("Solver status      :", x$status, "\n")
  cat("tr(M^{-1})         :", signif(x$value, 6), "\n")
  cat("Support points     :", nrow(x$design), "\n")
  print(x$design, row.names = FALSE, ...)
  invisible(x)
}


#' @export
print.cvx_c_design <- function(x, ...) {
  cat("c-optimal approximate design\n")
  cat("Solver status      :", x$status, "\n")
  cat("c'M^{-1}c          :", signif(x$value, 6), "\n")
  cat("Support points     :", nrow(x$design), "\n")
  print(x$design, row.names = FALSE, ...)
  invisible(x)
}
