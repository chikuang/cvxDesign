test_that("check_equivalence matches general for single D-opt design", {
  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 21)
  d <- calc_Dopt(u, quad_reg, drop_tol = 1e-4)

  eq1 <- check_equivalence(d, quad_reg, tol = 1e-3)
  eq2 <- check_equivalence_general(d, quad_reg, tol = 1e-3)

  expect_equal(eq1$directional_derivative, eq2$directional_derivative, tolerance = 1e-6)
  expect_true(eq2$all_nonpositive)
  expect_true(inherits(eq2, "cvx_equivalence_general"))
  expect_equal(eq2$mode, "single")
})

test_that("check_equivalence_general works for maximin design", {
  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 21)
  mm <- calc_maximin(u, quad_reg, criteria = c("D", "A"))

  eq <- suppressMessages(check_equivalence_general(mm, quad_reg, tol = 0.05))
  expect_true(eq$mode == "maximin")
  expect_true(length(eq$directional_derivatives$dD) == length(u))
  expect_true(length(eq$directional_derivative) == length(u))
  expect_equal(sum(eq$eta), 1, tolerance = 1e-8)
})

test_that("check_equivalence rejects maximin object", {
  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 11)
  mm <- calc_maximin(u, quad_reg, criteria = c("D", "A"))
  expect_error(check_equivalence(mm, quad_reg), "check_equivalence_general")
})

test_that("efficiency directional derivatives are smaller scale than KW for maximin", {
  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 51)
  mm <- calc_maximin(u, quad_reg, criteria = c("D", "A"))
  eq_kw <- suppressMessages(
    check_equivalence_general(mm, quad_reg, derivative = "kw_loss", tol = 1e6)
  )
  eq_ef <- check_equivalence_general(mm, quad_reg, derivative = "efficiency", tol = 1e6)
  expect_true(max(eq_ef$directional_derivative) < max(eq_kw$directional_derivative))
  expect_equal(eq_ef$derivative, "efficiency")
})

test_that("eta_rule optimize reduces max combined derivative vs binding", {
  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 41)
  mm <- calc_maximin(u, quad_reg, criteria = c("D", "A"))
  eq_b <- check_equivalence_general(mm, quad_reg, derivative = "efficiency", eta_rule = "binding")
  eq_o <- check_equivalence_general(mm, quad_reg, derivative = "efficiency", eta_rule = "optimize")
  expect_true(max(eq_o$directional_derivative) <= max(eq_b$directional_derivative) + 1e-8)
  expect_true(!is.null(eq_o$eta_optimization))
})
