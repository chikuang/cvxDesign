test_that("compute_maximin_design runs for D and A with reference losses", {
  skip_if_not_installed("CVXR")

  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 15)

  so_d <- compute_design_SO(u, quad_reg, criterion = "D", solver = "CLARABEL")
  so_a <- compute_design_SO(u, quad_reg, criterion = "A", solver = "CLARABEL")
  loss_ref <- list(D = so_d$loss, A = so_a$loss)

  res <- compute_maximin_design(
    u, quad_reg,
    loss_ref = loss_ref,
    criteria = c("D", "A"),
    solver = "CLARABEL"
  )

  expect_s3_class(res, "cvx_maximin_design")
  expect_s3_class(res, "cvx_design")
  expect_true(is.finite(res$tstar))
  expect_true(res$value <= 1 + 1e-3)
  expect_true(all(res$efficiency <= 1 + 1e-3))
  expect_equal(sum(res$weights), 1, tolerance = 1e-4)
})

test_that("print methods for SO and joint maximin run without error", {
  skip_if_not_installed("CVXR")

  f <- function(x) c(1, x)
  u <- c(-1, 0, 1)
  so <- compute_design_SO(u, f, criterion = "D", solver = "CLARABEL")
  expect_output(print(so), "Single-objective")

  so_d <- compute_design_SO(u, f, "D", solver = "CLARABEL")
  so_a <- compute_design_SO(u, f, "A", solver = "CLARABEL")
  mm <- compute_maximin_design(u, f, list(D = so_d$loss, A = so_a$loss), c("D", "A"), solver = "CLARABEL")
  expect_output(print(mm), "Joint maximin")
})
