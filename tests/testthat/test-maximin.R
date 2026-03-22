test_that("maximin design returns valid weights and efficiencies", {
  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 21)

  res <- calc_maximin(u, quad_reg, criteria = c("D", "A"))

  expect_true(inherits(res, "cvx_design"))
  expect_equal(sum(res$weights), 1, tolerance = 1e-5)
  expect_true(all(res$weights >= -1e-10))
  expect_true(res$value <= 1 + 1e-4)
  expect_true(all(res$efficiencies <= 1 + 1e-4))
  expect_equal(min(res$efficiencies), res$value, tolerance = 1e-4)
})

test_that("maximin with D, A, c requires cVec", {
  f <- function(x) c(1, x)
  u <- c(0, 1, 2)
  expect_error(calc_maximin(u, f, criteria = "c"), "cVec")
})
