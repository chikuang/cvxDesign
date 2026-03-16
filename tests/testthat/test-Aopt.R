test_that("A-optimal design works", {

  quad_reg <- function(x) c(1, x, x^2)

  u <- seq(-1, 1, length.out = 21)

  res <- calc_Aopt(u, quad_reg, drop_tol = 1e-4)

  expect_true(inherits(res, "cvx_design"))

  expect_true(res$value > 0)
})



# res$design |> round(3)
