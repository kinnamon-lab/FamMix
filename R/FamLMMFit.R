# FamLMMFit.R - Definition for FamLMMFit R6 subclass and methods

#' A \code{FamModelFit} subclass specialized for linear mixed models
#'
#' @description A specialized subclass of the
#'   \code{\link[FamModel]{FamModelFit}} abstract class for linear mixed model
#'   results. Methods of this subclass can be used to print model results in a
#'   nice format, access model objects, perform inference on general contrasts,
#'   and perform likelihood ratio tests on variance parameters.
#'
#' @export
FamLMMFit <- R6Class(
  "FamLMMFit",
  parent_env = getNamespace("FamModel"),
  lock_class = TRUE,
  inherit = FamModelFit,

  # Public members ============================================================

  public = list(

    # Constructor =============================================================

    #' @description Constructs a new instance of this class
    #'
    #' @param data The \code{\link{FamData}} object that produced this model
    #'   fit.
    #' @param formula The \code{\link[Formula]{Formula}} object that describes
    #'   the model.
    #' @param objfun A \code{\link[TMB]{MakeADFun}} object that has already been
    #'   run through an optimizer. This object is expected to have a "type"
    #'   attribute of "TMB".
    #' @param optres A \code{list} from the optimizer function detailing
    #'   results. This list is expected to have an attribute "opttype" that
    #'   contains a string identifying the function used for optimization.
    initialize = function(data, formula, objfun, optres) {
      if (class(data)[1] == "FamData") {
        private$data <- data
      } else {
        stop("Argument data must be a FamData object")
      }
      if (Formula::is.Formula(formula)) {
        private$formula <- formula
      } else {
        stop("Argument formula must be Formula object")
      }
      if (is.list(objfun) && attr(objfun, "type") == "TMB") {
        private$objfun <- objfun
        private$sd_report <- TMB::sdreport(objfun)
        par_list_names <- lapply(objfun[["env"]]$parameters, names)
        private$objfun_par_names <- names(objfun[["par"]])
        for (par_name in names(par_list_names)) {
          if (is.null(par_list_names[[par_name]])) next;
          private$objfun_par_names[private$objfun_par_names == par_name] <-
            if (par_name == "betas") {
              par_list_names[[par_name]]
            } else {
              paste(par_name, par_list_names[[par_name]], sep = ".")
            }
        }
        private$theta_hat <- private$sd_report[["par.fixed"]]
        names(private$theta_hat) <- private$objfun_par_names
        private$V_theta_hat <- private$sd_report[["cov.fixed"]]
        rownames(private$V_theta_hat) <- colnames(private$V_theta_hat) <-
          private$objfun_par_names
      } else {
        stop(
          "Argument objfun must be a list object with a type attribute ",
          "of TMB"
        )
      }
      if (is.list(optres) && attr(optres, "optfun") == "nlminb") {
        private$optres <- optres
      } else {
        stop(
          "Argument optres must be a list object with optfun attribute ",
          "of nlminb"
        )
      }
    },

    # Accessors ===============================================================

    #' @description Returns the \code{\link[Formula]{Formula}}
    #'   \code{\link[Formula]{Formula}} object used to construct the instance.
    get_formula = function() private$formula,

    #' @description Returns \code{\link[TMB]{MakeADFun}} object used to
    #'   construct the instance.
    get_objfun = function() private$objfun,

    #' @description Returns a vector of names derived from the original
    #'   formula that map to each element of \code{objfun$par}.
    get_objfun_par_names = function() private$objfun_par_names,

    #' @description Returns the \emph{conditional} residuals scaled by the
    #'   the inverse Cholesky factor of the covariance matrix of the
    #'   conditional distribution of the family members given the proband.
    #'   These are defined for family members only. If the model is correct,
    #'   these should approximate independent standard normal variables in
    #'   large samples.
    get_r_cholesky = function() private$objfun$report()[["r_cholesky"]],

    #' @description Returns the \emph{conditional} residuals, defined as
    #'   \code{y - X * b - d_pr}, where \code{X * b + d_pr} is the mean in
    #'   the family member conditional upon the data in the proband. These
    #'   are defined for family members only.
    get_r_cond = function() private$objfun$report()[["r_cond"]],

    #' @description Returns the \emph{population} residuals, defined as
    #'   \code{y - X * b}, for everyone, including the proband.
    get_r_pop = function() private$objfun$report()[["r_pop"]],

    #' @description Returns the \code{\link[TMB]{sdreport}} object produced
    #'   from the \code{\link[TMB]{MakeADFun}} object used to construct the
    #'   instance.
    get_sdreport = function() private$sd_report,

    # Print method ============================================================

    #' @description Prints contents of \code{FamLMMFit} objects with nice
    #'   formatting.
    #'
    #' @param ... Arguments passed on to \code{\link{print_ests}}.
    print = function(...) {
      theta_gr <- private$sd_report[["gradient.fixed"]]
      theta_hat <- private$theta_hat
      V_theta_hat <- private$V_theta_hat
      cat("\n===LINEAR MIXED MODEL RESULTS===\n")
      cat("DATA: ", private$data$get_data_name(), "\n", sep = "")
      cat("MEAN MODEL: ")
      print(formula(private$formula, rhs = 1))
      cat("VARIANCE PARAMETER GROUPS: ")
      print(formula(private$formula, lhs = 0, rhs = 2))
      cat(
        "\nFAMILIES USED: ",
        length(private$objfun[["env"]]$data[["f_sizes"]]), "\n", sep = ""
      )
      n_subj <- sum(private$objfun[["env"]]$data[["f_sizes"]])
      cat("SUBJECTS USED: ", n_subj, "\n", sep = "")
      cat(
        "PROBANDS: ", length(private$objfun[["env"]]$data[["f_pr_idxs"]]), "\n",
        sep = ""
      )
      cat("FAMILY SIZE DISTRIBUTION:\n")
      print(ftable(private$objfun[["env"]]$data[["f_sizes"]]))
      if (attr(private$optres, "optfun") == "nlminb") {
        cat(
          "CONVERGENCE ",
          ifelse(
            private$optres[["convergence"]] == 0,
            "ACHIEVED",
            "NOT ACHIEVED"
          ),
          " AT -2 LL = ", 2 * private$optres[["objective"]], " AFTER ",
          private$optres[["iterations"]], " ITERATIONS\n", sep = ""
        )
        cat("CRITERION: ", private$optres[["message"]], "\n", sep = "")
        cat(
          "MAX ABSOLUTE GRADIENT AT SOLUTION: ",
          format(max(abs(theta_gr)), scientific = TRUE),
          "\n", sep = ""
        )
      }
      v_parms <- grep("^(h2_a|sigma2)\\.*", names(theta_hat), value = TRUE)
      cat("\nVARIANCE PARAMETERS\n")
      print_ests(
        theta_hat[v_parms], V_theta_hat[v_parms, v_parms],
        CIs = FALSE, tests = FALSE
      )
      self$h2_a_lrts()
      cat("\nMEAN MODEL\n")
      mean_parms <- names(theta_hat)[!names(theta_hat) %in% v_parms]
      print_ests(
        theta_hat[mean_parms], V_theta_hat[mean_parms, mean_parms],
        ...
      )
      cat("\nNOTE: Wald tests and CIs are displayed in the above output\n")
    },

    # Inferential methods =====================================================

    #' @description Perform likelihood ratio test(s) for no polygenic effect(s)
    #'
    #' @details Tests that each narrow-sense heritability parameter is zero
    #'   individually using a likelihood ratio test. As these null hypotheses
    #'   are on the boundary of the parameter space, a 50:50 mixture of
    #'   chi-square(0) and chi-square(1) is used for the p-value (see Self and
    #'   Liang).
    #'
    #' @param print If \code{TRUE} (default), prints the likelihood ratio test
    #'   results in a nice format.
    #' @param ... Additional parameters to pass to the optimization function
    #'   (currently \code{\link{nlminb}}).
    #'
    #' @return A \code{data.table} containing likelihood ratio test results,
    #'   invisibly.
    #'
    #' @references
    #' Self SG, Liang K-Y. 1987. Asymptotic Properties of Maximum Likelihood
    #'   Estimators and Likelihood Ratio Tests Under Nonstandard Conditions.
    #'   \emph{J Am Stat Assoc} 82 (398): 605-10.
    #'   \url{https://doi.org/10.2307/2289471}.
    h2_a_lrts = function(print = TRUE, ...) {
      mod_data <- private$objfun[["env"]]$data
      parameters <- private$objfun[["env"]]$parameters
      h2_a_parms <- if (is.null(names(parameters[["h2_a"]]))) {
        "h2_a"
      } else {
        names(parameters[["h2_a"]])
      }
      parameters[["h2_a"]] <- rep(0, length(h2_a_parms))
      res <- rbindlist(lapply(h2_a_parms, function (x) {
        map_list <- list(
          h2_a = factor(seq(1, length(h2_a_parms)))
        )
        names(map_list[["h2_a"]]) <- h2_a_parms
        map_list[["h2_a"]][x] <- NA
        map_list[["h2_a"]] <- factor(map_list[["h2_a"]])
        objfun <- TMB::MakeADFun(
          mod_data,
          parameters,
          map = map_list,
          DLL = "sing_asc_lmm",
          method = NULL,
          silent = TRUE
        )
        if (attr(private$optres, "optfun") == "nlminb") {
          parm_lower <- rep(-Inf, length(objfun[["par"]]))
          parm_upper <- rep(Inf, length(objfun[["par"]]))
          names(parm_lower) <- names(parm_upper) <- names(objfun[["par"]])
          parm_lower[names(parm_lower) == "h2_a"] <- 0
          parm_upper[names(parm_upper) == "h2_a"] <- 1
          parm_lower[names(parm_lower) == "sigma2"] <-
            sqrt(.Machine[["double.eps"]])
          # Transform problem by scaling parameters by square root of their
          # Hessian diagonal elements at the initial estimates, which should
          # improve conditioning of the Hessian and therefore numerical
          # convergence. This is a non-adaptive version of the adaptive scaling
          # for unconstrained optimization in section 4c of the PORT
          # documentation included in the references for nlminb.
          d <- sqrt(abs(diag(objfun$he())))
          if (!all(is.finite(d))) {
            # If there are any -Inf, Inf, NA, or NaN values, do not scale
            d <- rep_len(1, length(d))
          } else {
            # Otherwise, set any scaling factors that are numerically <= 0 to
            # the minimum of those that are numerically > 0
            d[d < sqrt(.Machine[["double.eps"]])] <-
              min(d[d >= sqrt(.Machine[["double.eps"]])])
          }
          optres <- nlminb(
            objfun$par,
            objfun$fn,
            objfun$gr,
            objfun$he,
            lower = parm_lower,
            upper = parm_upper,
            scale = d,
            ...
          )
          converge <- (optres[["convergence"]] == 0)
        }

        # Gather LRT results
        max_abs_grad <- max(abs(TMB::sdreport(objfun)[["gradient.fixed"]]))
        lr_X2 <- 2 * (-private$optres[["objective"]] + optres[["objective"]])
        df_X2 <- 1
        p_X2 <- ifelse(
          abs(lr_X2) < sqrt(.Machine[["double.eps"]]),
          1,
          0.5 * pchisq(lr_X2, df_X2, lower.tail = FALSE)
        )

        data.table(
          h_0 = if (length(h2_a_parms) == 1) {
            "h2_a = 0"
          } else {
            paste0("h2_a.", x, " = 0")
          },
          converge,
          max_abs_grad,
          lr_X2,
          df_X2,
          p_X2
        )
      }))

      # In cases where the optimization failed to converge, set results to
      # missing
      res[
        converge != TRUE,
        `:=`(
          lr_X2 = NA,
          df_X2 = NA,
          p_X2 = NA
        )
      ]

      # Optional pretty printing
      if (print) {
        cat("\nLikelihood Ratio Tests\n")
        cat("----------------------\n")
        print(
          res[,
            .(
              Ho = h_0,
              `Max |g|` = format(max_abs_grad, scientific = TRUE),
              `LR X^2` = format(lr_X2, digits = 5),
              DF = df_X2,
              `Pr(> X^2)` = format.pval(p_X2)
            )
          ],
          row.names = FALSE
        )
        cat(
          "\nNOTE: P-values are calculated from a 50:50 mixture of ",
          "chi-square(0)\n      and chi-square(1) per Self and Liang (1987)\n",
          sep = ""
        )
      }

      invisible(res)
    }
  ),

  # Private members ===========================================================

  private = list(

    # Data members ============================================================

    formula = NULL,
    objfun = NULL,
    optfun = NULL,
    objfun_par_names = NULL,
    sd_report = NULL
  )
)
