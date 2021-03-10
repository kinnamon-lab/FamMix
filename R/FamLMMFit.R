# FamLMMFit.R - Definition for FamLMMFit R6 subclass and methods

#' A `FamModelFit` subclass specialized for linear mixed models
#'
#' @description A specialized subclass of the [`FamModelFit`] abstract class for
#'   linear mixed model results. Methods of this subclass can be used to print
#'   model results in a nice format, access model objects, perform inference on
#'   general contrasts, and perform likelihood ratio tests on variance
#'   parameters. Due to reliance on the internal structures of other objects,
#'   the object generator function is not exported, and instances of the
#'   `FamLMMFit` class can be initizalized *only from within functions in the
#'   `FamModel` namespace*.
#'
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
    #' @details Due to their reliance on the internal structures of other
    #' objects, instances of the `FamLMMFit` class can be initizalized *only
    #' from within functions in the `FamModel` namespace*.
    #'
    #' @param data The [`FamData`] object that produced this model
    #'   fit.
    #' @param formula The [`Formula::Formula`] object that describes
    #'   the model.
    #' @param objfun A [TMB::MakeADFun()] object that has already been run
    #'   through an optimizer. This object is expected to have a "type"
    #'   attribute of "LMM-SA".
    #' @param optres A `list` from the optimizer function detailing results.
    #'   This list is expected to have an attribute "opttype" that contains a
    #'   string identifying the function used for optimization.
    initialize = function(data, formula, objfun, optres) {
      # See documentation for packageName() regarding use of topenv(). Need to
      # move up 2 frames to get constructor calling environment because
      # initialize() is called by new()
      ctor_calling_env <- parent.frame(2)
      if (!identical(topenv(ctor_calling_env), getNamespace("FamModel"))) {
        stop(
          "FamLMMFit objects can only be constructed from within FamModel ",
          "package methods and functions"
        )
      }
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
      if (is.list(objfun) && attr(objfun, "type") == "LMM-SA") {
        private$objfun <- objfun
        private$theta_hat <- objfun[["env"]]$last.par.best
        # Numerical derivatives used for Hessian of objective function (negative
        # loglikelihood) by default; must supply Hessian calculated with
        # analytical derivatives
        hess <- Matrix::Matrix(objfun[["he"]](private$theta_hat))
        private$V_theta_hat <- as.matrix(Matrix::solve(hess))
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
        names(private$theta_hat) <- private$objfun_par_names
        rownames(private$V_theta_hat) <- colnames(private$V_theta_hat) <-
          private$objfun_par_names
        private$lmm_data <- lmm_data <- private$objfun[["env"]]$data

        # Verify that lmm_data rows map to specified IDs in original data member
        # after trip through MakeADFun(). Need to check:
        # 1) Row indices i corresponding to probands have incl_ids[i] mapping
        #    only to probands in the original FamData instance; likewise for
        #    non-probands
        # 2) For individual i in lmm_data model input, y[i], X[i, ], phi[i, ],
        #    and phi[, i] in lmm_data belong to individual incl_ids[i] in the
        #    original data set
        incl_ids <- as.data.table(lmm_data[["incl_ids"]])
        pr_mf_idxs <- lmm_data[["f_pr_idxs"]]
        non_pr_mf_idxs <- private$get_non_pr_mf_idxs()
        pr_ids <- incl_ids[pr_mf_idxs]
        non_pr_ids <- incl_ids[non_pr_mf_idxs]
        # Verify (1) above
        orig_data <- private$data$get_data()
        orig_phi <- private$data$get_phi()
        if (
          any(orig_data[pr_ids, pr != 1, on = .(fmid, id)]) ||
            any(orig_data[non_pr_ids, pr != 0, on = .(fmid, id)])
        ) {
          stop(
            "Proband and non-proband indexes or IDs from objective function ",
            "are incorrect"
          )
        }
        # Verify (2) above
        mm_objfun <- data.table(
          incl_ids,
          y = lmm_data[["y"]],
          lmm_data[["X"]]
        )
        mf_data <- model.frame(
          private$formula,
          orig_data[incl_ids, on = .(fmid, id)][order(fmid, id)],
          fmid = fmid,
          id = id
        )
        mm_data <- data.table(
          fmid = mf_data[["(fmid)"]],
          id = mf_data[["(id)"]],
          y = model.response(mf_data),
          model.matrix(private$formula, mf_data)
        )
        if (!isTRUE(all.equal(mm_objfun, mm_data))) {
          stop(
            "y[i] and X[i, ] in objfun data do not belong to individual ",
            "incl_ids[i] in the original data set for all i"
          )
        }
        incl_phi_ids <-
          if (!anyDuplicated(orig_data[, id])) {
            incl_ids[, as.character(id)]
          } else {
            incl_ids[, paste(fmid, id, sep = "/")]
          }
        if (!isTRUE(all.equal(
          as(as(lmm_data[["phi"]], "symmetricMatrix"), "dsCMatrix"),
          orig_phi[incl_phi_ids, incl_phi_ids]
        ))) {
          stop(
            "phi[i, ] and  phi[, i] in objfun data do not belong to ",
            "individual incl_ids[i] in the original phi matrix for all i"
          )
        }
      } else {
        stop(
          "Argument objfun must be a list object with a type attribute ",
          "of LMM-SA"
        )
      }
      if (is.list(optres)) {
        private$optres <- optres
      } else {
        stop("Argument optres must be a list object")
      }
    },

    # Accessors ===============================================================

    #' @description Returns the [`Formula::Formula`] object used to construct
    #'   the instance.
    get_formula = function() private$formula,

    #' @description Returns the original data with fitted values, residuals,
    #'   and diagnostics.
    #'
    #' @details See `vignette("linear_mixed_models")` for a detailed
    #'   description of each residual and diagnostic, suggestions for its use,
    #'   and relevant references.
    #'
    #' @param all If `FALSE` (default), include only observations that
    #'   contributed to the model fit. If `TRUE`, include all observations in
    #'   the [`FamData`] object `data` member. This can be useful for
    #'   generating a data set to refit the model without certain outliers.
    #'   Fitted values, residuals, and diagnostics columns will have `NA`
    #'   values for observations that did not contribute to the model fit.
    #'
    #' @return A [`data.table`] containing the `data` member of the
    #'   [`FamData`] object used to fit this model with the following
    #'   additional columns (notation defined in
    #'   `vignette("linear_mixed_models")`):
    #'   * `mu_hat`: The fitted population mean, \eqn{\hat{\mu}_{ij}}.
    #'   * `eta_hat`: The mean conditional upon the data in the proband,
    #'     \eqn{\hat{\eta}_{ij}}. These are defined for family members
    #'     (\eqn{j \neq j_i}) only and are `NA` in the proband.
    #'   * `r_c_hat`: The estimated Cholesky residuals,
    #'     \eqn{\hat{r}_{\mathrm{c}, ij}}. These are defined for family members
    #'     (\eqn{j \neq j_i}) only and are `NA` in the proband.
    #'   * `r_s_hat`: The estimated Pearson-type conditional residuals,
    #'     \eqn{\hat{r}_{\mathrm{s}, ij}}. These are defined for family members
    #'     (\eqn{j \neq j_i}) only and are `NA` in the proband.
    #'   * `c_star_hat`: Family-level chi-square goodness of fit statistic,
    #'     \eqn{\hat{c}^{*}_{i}}. The same value is provided for each family
    #'     member.
    #'   * `c_star_hat_df`: Degrees of freedom for \eqn{\hat{c}^{*}_{i}}, which
    #'     is equal to the number of non-probands in the family. The same value
    #'     is provided for each family member.
    #'   * `p_c_star_hat`: Probability of obtaining a deviate as or more extreme
    #'     than  \eqn{\hat{c}^{*}_{i}}, \eqn{\hat{p}_{\hat{c}^{*}_{i}}}. The
    #'     same value is provided for each family member.
    #'   * `r_star_hat`: Individual-level goodness of fit statistic,
    #'     \eqn{\hat{r}^{*}_{ij}}. These are defined for family members
    #'     (\eqn{j \neq j_i}) only and are `NA` in the proband.
    get_model_res = function(all = FALSE) {
      if (!all %in% c(TRUE, FALSE)) {
        stop("Argument all must be either TRUE or FALSE")
      }
      # REPORTed values at final parameter estimates calculated in C++ template
      # for checking
      objfun_report <- private$objfun[["report"]](private$theta_hat)
      lmm_data <- private$lmm_data
      incl_ids <- as.data.table(lmm_data[["incl_ids"]])
      pr_mf_idxs <- lmm_data[["f_pr_idxs"]]
      non_pr_mf_idxs <- private$get_non_pr_mf_idxs()
      parameters <- private$objfun[["env"]]$parameters
      beta_parms <- names(parameters[["betas"]])
      h2_a_parms <- if (is.null(names(parameters[["h2_a"]]))) {
        "h2_a"
      } else {
        paste0("h2_a.", names(parameters[["h2_a"]]))
      }
      sigma_parms <- if (is.null(names(parameters[["sigma"]]))) {
        "sigma"
      } else {
        paste0("sigma.", names(parameters[["sigma"]]))
      }

      # Make population covariance matrix
      h2_a <- Matrix::Diagonal(
        x = rep(
          lmm_data[["f_pops"]] %*% private$theta_hat[h2_a_parms],
          times = lmm_data[["f_sizes"]]
        )
      )
      sigma2 <- Matrix::Diagonal(
        x = rep(
          lmm_data[["f_pops"]] %*% private$theta_hat[sigma_parms] ^ 2,
          times = lmm_data[["f_sizes"]]
        )
      )
      sigma_mat <- as(
        2 * (sigma2 * h2_a) %*% lmm_data[["phi"]] + sigma2 - sigma2 * h2_a,
        "symmetricMatrix"
      )

      # All-subject diagnostics
      y <- lmm_data[["y"]]
      mu_hat <- as.numeric(lmm_data[["X"]] %*% private$theta_hat[beta_parms])
      if (!isTRUE(all.equal(mu_hat, objfun_report[["X_beta"]]))) {
        stop("Mismatch in mu_hat between C++ and R code")
      }
      dx_all <- data.table(incl_ids, mu_hat)

      # Sanity check of population residuals
      r_pop <- y - mu_hat
      if (!isTRUE(all.equal(r_pop, objfun_report[["r_pop"]]))) {
        stop("Mismatch in r_pop between C++ and R code")
      }

      # Family-member-only diagnostics
      sigma_mat_pr_pr <- sigma_mat[pr_mf_idxs, pr_mf_idxs]
      sigma_mat_non_pr_pr <- sigma_mat[non_pr_mf_idxs, pr_mf_idxs]
      eta_hat <- as.numeric(
        mu_hat[non_pr_mf_idxs] +
          sigma_mat_non_pr_pr %*%
          Matrix::solve(
            sigma_mat_pr_pr,
            y[pr_mf_idxs] - mu_hat[pr_mf_idxs]
          )
      )
      omega_mat <- as(
        sigma_mat[non_pr_mf_idxs, non_pr_mf_idxs] -
          sigma_mat_non_pr_pr %*%
          Matrix::solve(
            sigma_mat_pr_pr,
            Matrix::t(sigma_mat_non_pr_pr)
          ),
        "symmetricMatrix"
      )

      # Conditional and Cholesky residuals
      r_hat <- y[non_pr_mf_idxs] - eta_hat
      r_s_hat <- r_hat / Matrix::diag(omega_mat)
      r_c_hat <- as.numeric(
        Matrix::solve(Matrix::t(Matrix::chol(omega_mat)), r_hat)
      )

      # Adaptation of Q* from Hopper and Mattews (1982), referred to here as
      # r*
      omega_mat_inv <- as(Matrix::solve(omega_mat), "symmetricMatrix")
      T_hat_inv <- 1 / sqrt(Matrix::diag(omega_mat_inv))
      r_star_hat <- as.numeric(
          Matrix::Diagonal(x = T_hat_inv) %*% Matrix::solve(omega_mat, r_hat)
      )
      dx_non_pr <- data.table(
        incl_ids[non_pr_mf_idxs], eta_hat, r_hat, r_c_hat, r_s_hat, r_star_hat
      )

      # Add family-level diagnostics
      dx_non_pr[,
        `:=`(
          c_star_hat = sum(r_c_hat ^ 2),
          c_star_hat_df = .N
        ),
        by = .(fmid)
      ][,
        p_c_star_hat := pchisq(c_star_hat, c_star_hat_df, lower.tail = FALSE)
      ]
      # Beaty, Liang, and Rao (1987) note that the sum of the family chi-square
      # statistics, which is the sum of the squared Cholesky residuals, should
      # equal the sample size algebraically. An relative comparison within 1e-5
      # is performed to account for floating point error
      if (!isTRUE(all.equal(
        dx_non_pr[, sum(r_c_hat ^ 2)],
        nrow(dx_non_pr),
        tolerance = 1e-4
      ))) {
        stop("Squared Cholesky residuals do not sum to total sample size")
      }

      # Package results
      res <- private$data$get_data()
      if (!all) res <- res[incl_ids, on = .(fmid, id)]
      res <- merge(res, dx_all, by = c("fmid", "id"), all.x = TRUE)
      res <- merge(res, dx_non_pr, by = c("fmid", "id"), all.x = TRUE)
      # Ensure that original sort order inherited from FamData data member
      # is restored after merge
      setkeyv(res, key(private$data$get_data()))
      res
    },

    #' @description Returns [TMB::MakeADFun()] object used to construct the
    #'   instance.
    get_objfun = function() private$objfun,

    #' @description Returns a vector of names derived from the original formula
    #'   that map to each element of `objfun$par`.
    get_objfun_par_names = function() private$objfun_par_names,

    # Print method ============================================================

    #' @description Prints contents of `FamLMMFit` objects with nice formatting.
    #'
    #' @param ... Arguments passed on to [print_ests()].
    print = function(...) {
      lmm_data <- private$lmm_data
      theta_gr <- private$objfun[["gr"]](private$theta_hat)
      theta_hat <- private$theta_hat
      V_theta_hat <- private$V_theta_hat
      cat("\n===LINEAR MIXED MODEL RESULTS===\n")
      cat("DATA: ", private$data$get_data_name(), "\n", sep = "")
      cat("MEAN MODEL: ")
      print(formula(private$formula, rhs = 1), showEnv = FALSE)
      cat("VARIANCE PARAMETER GROUPS: ")
      print(formula(private$formula, lhs = 0, rhs = 2), showEnv = FALSE)
      cat(
        "\nFAMILIES USED: ",
        length(lmm_data[["f_sizes"]]), "\n", sep = ""
      )
      n_subj <- sum(lmm_data[["f_sizes"]])
      cat("SUBJECTS USED: ", n_subj, "\n", sep = "")
      cat(
        "PROBANDS: ", length(lmm_data[["f_pr_idxs"]]), "\n",
        sep = ""
      )
      cat("FAMILY SIZE DISTRIBUTION:\n")
      print(ftable(lmm_data[["f_sizes"]]))
      cat(
        "CONVERGENCE ",
        ifelse(
          private$optres[["convergence"]] == 0,
          "ACHIEVED",
          "NOT ACHIEVED"
        ),
        " AT -2 LL = ", 2 * private$optres[["value"]], "\n",
        "EVALUATIONS:\n", sep = ""
      )
      print(private$optres[["counts"]], sep = "")
      cat("MESSAGE: ", private$optres[["message"]], "\n", sep = "")
      cat(
        "MAX ABSOLUTE ELEMENT OF LL GRADIENT (g) AT SOLUTION: ",
        format(max(abs(theta_gr)), scientific = TRUE),
        "\n", sep = ""
      )
      # Hessian of objective function, which is negative Hessian of logliklihood
      # (-H)
      neg_hess_ll <- Matrix::Matrix(private$objfun[["he"]](private$theta_hat))
      min_ev_neg_hess_ll <-
        min(eigen(neg_hess_ll, only.values = TRUE)[["values"]])
      sv_neg_hess_ll <- svd(neg_hess_ll)[["d"]]
      rcond_neg_hess_ll <- min(sv_neg_hess_ll) / max(sv_neg_hess_ll)
      cat("NEGATIVE LL HESSIAN (-H) CHARACTERISTICS AT SOLUTION:\n")
      cat(
        "   SMALLEST EIGENVALUE: ",
        format(min_ev_neg_hess_ll, scientific = TRUE),
        "\n", sep = ""
      )
      cat(
        "   RECIPROCAL CONDITION NUMBER: ",
        format(rcond_neg_hess_ll, scientific = TRUE),
        "\n", sep = ""
      )
      # Note that the objective function is the negative loglikelihood, so the
      # gradient of this is negative the gradient of the
      # loglikelihood. Nonetheless, the two negatives cancel out
      cat(
        "SCALED LL GRADIENT (-g' * H^-1 * g) CRITERION AT SOLUTION: ",
        format(theta_gr %*% V_theta_hat %*% t(theta_gr), scientific = TRUE),
        "\n", sep = ""
      )
      v_parms <- grep("^(h2_a|sigma)(\\..+)?$", names(theta_hat), value = TRUE)
      cat("\nVARIANCE PARAMETERS\n")
      print_ests(
        theta_hat[v_parms], V_theta_hat[v_parms, v_parms],
        CIs = FALSE, tests = FALSE
      )
      self$get_h2_a_lrts()
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
    #'   individually using a likelihood ratio test. See
    #'   `vignette("linear_mixed_models")` for a detailed description of the
    #'   methods used. By default, results are cached in the [`FamLMMFit`]
    #'   object after the first call to avoid redundant optimizations on
    #'   subsequent calls.
    #'
    #' @param print If `TRUE` (default), prints the likelihood ratio test
    #'   results in a nice format.
    #' @param use_cached If `TRUE` (default), use cached likelihood ratio
    #'   test results if any exist. Otherwise, new likelihood ratio tests will
    #'   be run, and these results will replace any in the cache.
    #' @param ... Additional parameters to pass to the `control` list for
    #'   [optim()] with `method = "L-BGFS-B"`. Note that `parscale` and
    #'   `fnscale` cannot be modified.
    #'
    #' @return A [`data.table`] containing likelihood ratio test results,
    #'   invisibly.
    get_h2_a_lrts = function(print = TRUE, use_cached = TRUE, ...) {
      if (!use_cached || is.null(private$h2_a_lrts)) {
        mod_data <- private$objfun[["env"]]$data
        parameters <- private$objfun[["env"]]$parameters
        h2_a_parms <- if (is.null(names(parameters[["h2_a"]]))) {
          "h2_a"
        } else {
          names(parameters[["h2_a"]])
        }
        parameters[["h2_a"]] <- rep(0, length(h2_a_parms))
        names(parameters[["h2_a"]]) <- h2_a_parms
        res <- rbindlist(lapply(h2_a_parms, function(x) {
          map_list <- list(
            h2_a = seq(1, length(h2_a_parms))
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
          optres <- lmm_optim(objfun, ...)
          converge <- (optres[["convergence"]] == 0)

          # Gather LRT results
          neg_hess_ll <-
            Matrix::Matrix(objfun[["he"]](objfun[["env"]]$last.par.best))
          min_ev_neg_hess_ll <-
            min(eigen(neg_hess_ll, only.values = TRUE)[["values"]])
          sv_neg_hess_ll <- svd(neg_hess_ll)[["d"]]
          rcond_neg_hess_ll <- min(sv_neg_hess_ll) / max(sv_neg_hess_ll)
          theta_gr <- objfun[["gr"]](objfun[["env"]]$last.par.best)
          max_abs_grad <- max(abs(theta_gr))
          scaled_grad <- as.numeric(
            theta_gr %*% Matrix::solve(neg_hess_ll) %*% t(theta_gr)
          )
          if (converge) {
            cur_par <- if (length(h2_a_parms) == 1) {
              "h2_a"
            } else {
              paste0("h2_a.", x)
            }
            # Account for numerical imprecision in loglikelihood optimization
            # causing non-zero lr_X2 when estimated h2_a is exactly 0 and lr_X2
            # should be as well
            if (
              isTRUE(all.equal(optres[["value"]], private$optres[["value"]])) &&
                private$theta_hat[cur_par] == 0
            ) {
              lr_X2 <- 0
            } else {
              lr_X2 <- 2 * (optres[["value"]] - private$optres[["value"]])
            }
            if (lr_X2 < 0) {
              warning(
                "A negative LR statistic (", lr_X2, ") was obtained for ",
                "an estmiated h2_a = ", private$theta_hat[cur_par], " and ",
                "set to zero"
              )
              lr_X <- 0
            }
            p_X2 <- sum(0.5 * pchisq(lr_X2, c(0, 1), lower.tail = FALSE))
          } else {
            # In cases where the optimization failed to converge, set LRT X2 and
            # p-value to missing
            lr_X2 <- NA
            p_X2 <- NA
          }

          data.table(
            h_0 = if (length(h2_a_parms) == 1) {
              "h2_a = 0"
            } else {
              paste0("h2_a.", x, " = 0")
            },
            converge,
            max_abs_grad,
            scaled_grad,
            min_ev_neg_hess_ll,
            rcond_neg_hess_ll,
            lr_X2,
            p_X2
          )
        }))

        private$h2_a_lrts <- res
      } else {
        res <- private$h2_a_lrts
      }

      # Optional pretty printing
      if (print) {
        cat("\nLikelihood Ratio Tests\n")
        cat("----------------------\n")
        print(
          res[,
            .(
              Ho = h_0,
              `Max |g|` = format(max_abs_grad, scientific = TRUE),
              `-g' * H^-1 * g` = format(scaled_grad, scientific = TRUE),
              `Min lambda(-H)` = format(min_ev_neg_hess_ll, scientific = TRUE),
              `1 / kappa(-H)` = format(rcond_neg_hess_ll, scientific = TRUE),
              `LR X^2` = format(lr_X2, digits = 5),
              `Pr(> X^2)` = format.pval(p_X2)
            )
          ],
          row.names = FALSE
        )
      }

      invisible(res)
    }
  ),

  # Private members ===========================================================

  private = list(

    # Data members ============================================================

    formula = NULL,
    h2_a_lrts = NULL,
    lmm_data = NULL,
    objfun = NULL,
    objfun_par_names = NULL,
    sd_report = NULL,

    # Methods =================================================================

    get_non_pr_mf_idxs = function() {
      # This construct for non_pr_mf_idxs ensures that we are effectively
      # "striking out" the rows belonging to probands while otherwise
      # maintaining the order of the model frame
      pr_mf_idxs <- private$lmm_data[["f_pr_idxs"]]
      sort(setdiff(seq_along(private$lmm_data[["y"]]), pr_mf_idxs))
    }
  )
)
