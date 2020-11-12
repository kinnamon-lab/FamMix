# Contrast.R - Definition for Contrast R6 class and methods

#' Contrast estimated from an object inheriting from `FamModelFit`
#'
#' @description An R6 class that stores results for a contrast of the form
#'   \eqn{L \theta - m} estimated using a specified object inheriting from
#'   [`FamModelFit`]. Wald tests and confidence intervals for individual
#'   rows of the contrast matrix as well as an overall Wald chi-square test of
#'   the null hypothesis \eqn{L \theta - m = 0} are provided.
#'
#' @export
Contrast <- R6Class(
  "Contrast",
  parent_env = getNamespace("FamModel"),
  lock_class = TRUE,

  # Public members ============================================================

  public = list(

    # Constructor =============================================================

    #' @description Constructs a new instance of this class.
    #' @param model_fit An object inheriting from [`FamModelFit`].
    #' @param L_mat A contrast vector (1 df) or `matrix` (>1 df) containing
    #'   one contrast in each row.
    #' @param m An optional vector containing the null value for each contrast.
    #'   Will be set to the zero vector of length `nrow(L_mat)` if not
    #'   specified.
    initialize = function(model_fit, L_mat, m) {
      if (any(class(model_fit) == "FamModelFit")) {
        private$model_fit <- model_fit
      } else {
        stop("Argument model_fit must be an object inheriting from FamModelFit")
      }
      if (is.vector(L_mat, mode = "numeric")) {
        L_mat <- matrix(L_mat, nrow = 1, ncol = length(L_mat), byrow = TRUE)
      } else if (!is.matrix(L_mat)) {
        stop("Argument L_mat must be a numeric vector or matrix")
      }
      if (missing(m)) {
        m <- rep(0, nrow(L_mat))
      } else if (!is.vector(m, mode = "numeric") || length(m) != nrow(L_mat)) {
        stop(
          "Argument m must be a numeric vector with as many elements as ",
          "there are rows of L_mat"
        )
      }
      rownames(L_mat) <- paste("Row", seq(1, nrow(L_mat)))
      colnames(L_mat) <- names(model_fit$get_theta_hat())
      L_theta_hat_m <- drop(L_mat %*% model_fit$get_theta_hat()) - m
      V_L_theta_hat <- tcrossprod(
        L_mat %*% model_fit$get_V_theta_hat(),
        L_mat
      )
      X2 <- drop(
        crossprod(L_theta_hat - m, solve(V_L_theta_hat, L_theta_hat - m))
      )
      df_X2 <- Matrix::rankMatrix(L_mat)
      p_X2 <- pchisq(X2, df = df_X2, lower.tail = FALSE)
      private$L_mat <- L_mat
      private$m <- m
      private$L_theta_hat_m <- L_theta_hat_m
      private$V_L_theta_hat <- V_L_theta_hat
      private$X2 <- X2
      private$df_X2 <- df_X2
      private$p_X2 <- p_X2
    },

    # Accessors ===============================================================

    #' @description Returns `model_fit` object.
    get_model_fit = function() private$model_fit,

    #' @description Returns `L_mat`.
    get_L_mat = function() private$L_mat,

    #' @description Returns `m`.
    get_m = function() private$m,

    #' @description Returns \eqn{L \hat{\theta} - m}.
    get_L_theta_hat_m = function() private$L_theta_hat_m,

    #' @description Returns \eqn{L \hat{V}(\hat{\theta}) L^{'}}.
    get_V_L_theta_hat = function() private$V_L_theta_hat,

    #' @description Returns Wald chi-square statistic for null hypothesis
    #'   \eqn{L \theta - m = 0}.
    get_X2 = function() private$X2,

    #' @description Returns degrees of freedom of Wald chi-square statistic.
    get_df_X2 = function() private$df_X2,

    #' @description Returns p-value for Wald chi-square statistic.
    get_p_X2 = function() private$p_X2,

    # Print method ============================================================

    #' @description Formatted printing of the `Contrast` object.
    #' @param ... Arguments passed on to [print_ests()].
    print = function(...) {
      cat("\n=== CONTRAST ===\n")
      cat("\nContrast Matrix and Null Values (L_mat|m)'\n")
      cat("------------------------------------------\n")
      print(
        rbind(
          t(private$L_mat),
          `-` = rep("-", nrow(private$L_mat)),
          m = private$m
        ),
        quote = FALSE
      )
      print_ests(private$L_theta_hat_m, private$V_L_theta_hat, ...)
      cat("\nWald Chi-square Test of Ho: L_mat %*% theta_hat = m\n")
      cat("---------------------------------------------------\n")
      cat(
        "X2 = ", private$X2, ", df = ", private$df_X2,
        ", P(> X2) = ", format.pval(private$p_X2), "\n", sep = ""
      )
    }
  ),

  # Private members ===========================================================

  private = list(

    # Data members ============================================================

    model_fit = NULL,
    L_mat = NULL,
    m = NULL,
    L_theta_hat_m = NULL,
    V_L_theta_hat = NULL,
    X2 = NULL,
    df_X2 = NULL,
    p_X2 = NULL
  )
)
