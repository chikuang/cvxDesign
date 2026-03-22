test_that("compute_design_SO D-criterion returns valid design", {
  skip_if_not_installed("CVXR")

  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 21)

  res <- compute_design_SO(u, quad_reg, criterion = "D", solver = "CLARABEL")

  expect_s3_class(res, "cvx_so_design")
  expect_s3_class(res, "cvx_design")
  expect_equal(sum(res$weights), 1, tolerance = 1e-5)
  expect_true(all(res$weights >= -1e-10))
  expect_true(is.finite(res$loss))
  expect_equal(dim(res$info_matrix), c(3L, 3L))
})
