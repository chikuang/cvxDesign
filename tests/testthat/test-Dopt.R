test_that("D-optimal design returns valid weights", {

  quad_reg <- function(x) c(1, x, x^2)

  u <- seq(-1, 1, length.out = 21)

  res <- calc_Dopt(u, quad_reg)

  expect_true(inherits(res, "cvx_design"))

  expect_equal(sum(res$weights), 1, tolerance = 1e-6)

  expect_true(all(res$weights >= 0))
})
