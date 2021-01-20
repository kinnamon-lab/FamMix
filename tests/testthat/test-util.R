# Tests of exported utility functions

test_that("Bad estimates and covariance matrix arguments are caught", {

  # Bad theta_hat
  expect_error(
    print_ests(NULL, NULL),
    regexp = "Argument theta_hat must be a named numeric vector"
  )
  expect_error(
    print_ests(seq(1, 5), NULL),
    regexp = "Argument theta_hat must be a named numeric vector"
  )
  expect_error(
    print_ests(matrix(seq(1, 5)), NULL),
    regexp = "Argument theta_hat must be a named numeric vector"
  )
  expect_error(
    print_ests(letters[1:5], NULL),
    regexp = "Argument theta_hat must be a named numeric vector"
  )

  # Good theta_hat, bad V_theta_hat
  theta_hat <- c(20, -2, 1, 5)
  names(theta_hat) <- c("beta0", "beta1", "beta2", "beta3")
  expect_error(
    print_ests(theta_hat, NULL),
    regexp = "Argument V_theta_hat must be a numeric covariance matrix"
  )
  expect_error(
    print_ests(theta_hat, matrix(letters[1:16], nrow = 4, ncol = 4)),
    regexp = "Argument V_theta_hat must be a numeric covariance matrix"
  )
  expect_error(
    print_ests(theta_hat, seq(1, 16)),
    regexp = "Argument V_theta_hat must be a numeric covariance matrix"
  )
  expect_error(
    print_ests(theta_hat, matrix(seq(1, 16), nrow = 4, ncol = 4)),
    regexp = "Argument V_theta_hat must be a numeric covariance matrix"
  )
  expect_error(
    print_ests(theta_hat, diag(3)),
    regexp = "Argument V_theta_hat must be a numeric covariance matrix"
  )
  expect_error(
    print_ests(theta_hat, diag(5)),
    regexp = "Argument V_theta_hat must be a numeric covariance matrix"
  )

  # Warn if V_theta_hat has different names than theta_hat
  V_theta_hat <- diag(4)
  rownames(V_theta_hat) <- rev(names(theta_hat))
  expect_warning(
    capture_output(print_ests(theta_hat, V_theta_hat)),
    regexp = "Row/column names of V_theta_hat differ from those of theta_hat"
  )
  rownames(V_theta_hat) <- NULL
  colnames(V_theta_hat) <- rev(names(theta_hat))
  expect_warning(
    capture_output(print_ests(theta_hat, V_theta_hat)),
    regexp = "Row/column names of V_theta_hat differ from those of theta_hat"
  )
})

test_that("Bad trans argument is caught", {
  theta_hat <- c(20, -2, 1, 5)
  names(theta_hat) <- c("beta0", "beta1", "beta2", "beta3")
  expect_error(
    print_ests(theta_hat, diag(4), trans = 1),
    regexp = "Argument trans must be a function with exactly one argument"
  )
  expect_error(
    print_ests(theta_hat, diag(4), trans = pnorm),
    regexp = "Argument trans must be a function with exactly one argument"
  )
})

test_that("Printed output is correct", {
  ps <- c(0.90, 0.95, 0.975, 0.995)
  sds <- seq(1, 4)
  lbls <- paste("n", sds, ps, sep = "_")
  theta_hat <- qnorm(ps, sd = sds)
  names(theta_hat) <- lbls
  V_theta_hat <- diag(sds ^ 2)
  dimnames(V_theta_hat) <- list(lbls, lbls)
  expect_snapshot(print_ests(theta_hat, V_theta_hat))
  expect_snapshot(print_ests(theta_hat, V_theta_hat, trans = exp))
})
