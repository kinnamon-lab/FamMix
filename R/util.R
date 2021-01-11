# util.R - Utility functions for use with FamModel package

#' Make coefficient matrix
#'
#' Make coefficient matrix containing estimates, standard errors,
#' \eqn{100(1 - \alpha)}% CIs, Z statistics, and p-values given an input
#' vector of parameter estimates and their standard errors.
#'
#' @param theta_hat A named vector of parameter estimates.
#' @param V_theta_hat The covariance matrix of the parameter estimates in
#'   `theta_hat`.
#' @param trans A function accepting a single argument (e.g., `exp`) that
#'   is used to transform parameter estimates and confidence intervals. All
#'   tests and standard errors are still on the untransformed scale.
#' @param alpha The alpha level for the confidence intervals; default 0.05.
#' @param SEs `TRUE` (default) if standard errors should be reported.
#' @param CIs `TRUE` (default) if \eqn{100(1 - \alpha)}% CIs should be
#'   reported.
#' @param tests `TRUE`(default) if Wald test statistics and p-values based
#'   on the normal distribution should be reported.
#'
#' @return A `matrix` with row names given by the names of `theta_hat`
#'   containing the requested columns.
#'
#' @export
make_coef_mat <- function(theta_hat, V_theta_hat, trans, alpha = 0.05,
  SEs = TRUE, CIs = TRUE, tests = TRUE) {
  cl <- 100 * (1 - alpha)
  se_theta_hat <- sqrt(diag(V_theta_hat))
  trans_suff <- ""
  if (tests) Z <- theta_hat / se_theta_hat
  if (CIs) {
    lcl_theta <- theta_hat - se_theta_hat * stats::qnorm(1 - alpha / 2)
    ucl_theta <- theta_hat + se_theta_hat * stats::qnorm(1 - alpha / 2)
  }
  if (!missing(trans) && is.function(trans) &&
    length(formals(trans)) == 1) {
    trans_suff <- " (tr)"
    theta_hat <- trans(theta_hat)
    if (CIs) {
      lcl_theta <- trans(lcl_theta)
      ucl_theta <- trans(ucl_theta)
    }
  }
  coef_mat <- theta_hat
  header <- paste0("Estimate", trans_suff)
  if (SEs) {
    coef_mat <- cbind(coef_mat, se_theta_hat)
    header <- c(header, "SE")
  }
  if (CIs) {
    coef_mat <- cbind(coef_mat, lcl_theta, ucl_theta)
    header <- c(header, paste0(cl, "% LCL", trans_suff),
      paste0(cl, "% UCL", trans_suff))
  }
  if (tests) {
    coef_mat <- cbind(
      coef_mat, Z, stats::pchisq(Z ^ 2, df = 1, lower.tail = FALSE)
    )
    header <- c(header, "Z value", "Pr(>|Z|)")
  }
  colnames(coef_mat) <- header

  coef_mat
}

#' Print parameter estimates
#'
#' Prints estimates, standard errors, \eqn{100(1 - \alpha)}% CIs, test
#' statistics, and p-values to the console given a vector of estimates and their
#' covariance matrix. May also specify a function to transform the estimates and
#' 95% CIs.  Uses [make_coef_mat()] for underlying calculations.
#'
#' @param theta_hat A named vector of parameter estimates.
#' @param V_theta_hat The covariance matrix of the parameter estimates in
#'   `theta_hat`.
#' @param trans A function accepting a single argument (e.g., `exp`) that
#'   is used to transform parameter estimates and confidence intervals. All
#'   tests and standard errors are still on the untransformed scale.
#' @inheritDotParams make_coef_mat
#' @inheritDotParams stats::printCoefmat -x -cs.ind -tst.ind
#'
#' @export
print_ests <- function(theta_hat, V_theta_hat, trans, ...) {
  make_coef_mat_args <- list(
    quote(theta_hat), quote(V_theta_hat), quote(trans)
  )
  dots <- list(...)
  make_coef_mat_args <- c(
    make_coef_mat_args, dots[names(dots) %in% names(formals(make_coef_mat))]
  )
  coef_mat <- do.call(make_coef_mat, make_coef_mat_args)
  printCoefmat_args <- list(quote(coef_mat))
  if (!"Z value" %in% colnames(coef_mat)) {
    printCoefmat_args$cs.ind <- seq(1, ncol(coef_mat))
    printCoefmat_args$tst.ind <- integer(0)
  }
  printCoefmat_args <- c(
    printCoefmat_args,
    dots[names(dots) %in% names(formals(stats::printCoefmat))]
  )
  cat("\nParameter Estimates\n")
  cat("-------------------\n")
  if (!missing(trans) && is.function(trans) && length(formals(trans)) == 1) {
    cat("Transformation:", deparse(trans)[-1], "\n")
  }
  do.call(stats::printCoefmat, printCoefmat_args)
}

#' Executes LMM optimization
#'
#' The objective function made by [TMB::MakeADFun()] is a `list` containing an
#' environment and a series of functions enclosed by it. Thus, regardless of
#' whether `objfun` is copied, the function calls will still update the
#' environment in the object passed through `objfun` in the calling frame,
#' meaning that this function has the same side effect as executing the code in
#' the calling frame.
#'
#' @param objfun An objective function `list` from [TMB::MakeADFun()].
#' @param ... Additional parameters to pass to the `control` list for [optim()]
#'   with `method = "L-BGFS-B"`. Note that `parscale` and `fnscale` cannot
#'   be modified.
#'
#' @return A `list` with optimization results. See [optim()].
#'
#' @noRd
lmm_optim <- function(objfun, ...) {
  parm_lower <- rep(-Inf, length(objfun[["par"]]))
  parm_upper <- rep(Inf, length(objfun[["par"]]))
  names(parm_lower) <- names(parm_upper) <- names(objfun[["par"]])
  parm_lower[names(parm_lower) == "h2_a"] <- 0
  parm_upper[names(parm_upper) == "h2_a"] <- 1
  parm_lower[names(parm_lower) == "sigma"] <-
    sqrt(.Machine[["double.eps"]])
  # Transform problem by multiplying parameters by square root of their Hessian
  # diagonal elements at the initial estimates, which should improve
  # conditioning of the Hessian approximation and therefore numerical
  # convergence.  Note that optim converts the parameters by *dividing* by
  # parscale, so we need to set parscale to 1/d
  d <- sqrt(abs(diag(objfun[["he"]](objfun[["par"]]))))
  if (!all(is.finite(d))) {
    # If there are any -Inf, Inf, NA, or NaN values, do not scale
    d <- rep_len(1, length(d))
  } else {
    # Otherwise, set any scaling factors that are numerically <= 0 to the
    # minimum of those that are numerically > 0
    d[d <= sqrt(.Machine[["double.eps"]])] <-
      min(d[d > sqrt(.Machine[["double.eps"]])])
  }
  control <- list(...)
  control[["parscale"]] <- 1 / d
  # Scale objective function by inverse of initial value
  control[["fnscale"]] <- objfun[["fn"]](objfun[["par"]])
  if (is.null(control[["factr"]])) {
    control[["factr"]] <- 1e-10 / .Machine[["double.eps"]]
  }
  if (is.null(control[["pgtol"]])) {
    control[["pgtol"]] <- sqrt(.Machine[["double.eps"]])
  }
  optres <- optim(
    objfun[["par"]],
    objfun[["fn"]],
    objfun[["gr"]],
    method = "L-BFGS-B",
    lower = parm_lower,
    upper = parm_upper,
    control = control
  )

  optres
}
