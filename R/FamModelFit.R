# FamModelFit.R - Definition for FamModelFit R6 abstract class and methods

#' Abstract class for results of a fitted `FamModel` model
#'
#' @description An abstract R6 class defining the pattern for all `FamModelFit`
#'   subclasses created by fitting a model to a [`FamData`] object using one
#'   of its model-fitting methods. This class provides only a default
#'   constructor without any initializer.
#'
#' @export
FamModelFit <- R6Class(
  "FamModelFit",
  parent_env = getNamespace("FamModel"),
  lock_class = TRUE,

  # Public members ============================================================

  public = list(

    # Accessors ===============================================================

    #' @description Returns the [`FamData`] object that produced this model
    #'   fit.
    get_data = function() private$data,

    #' @description Returns `list` of optimization results.
    get_optres = function() private$optres,

    #' @description Returns \eqn{\hat{\theta}}, the vector of parameter
    #'   estimates.
    get_theta_hat = function() private$theta_hat,

    #' @description Returns \eqn{\hat{V}(\hat{\theta})}, the estimated
    #'   covariance matrix of the parameter estimates.
    get_V_theta_hat = function() private$V_theta_hat,

    # Print method ============================================================

    #' @description Formatted printing of the `FamModelFit` object.
    #'
    #' @param ... Arguments passed on to [print_ests()].
    print = function(...) {
      print_ests(private$theta_hat, private$V_theta_hat, ...)
    },

    # Inferential methods =====================================================

    #' @description Create a new `Contrast` object.
    #' @param L_mat A contrast vector (1 df) or `matrix` (>1 df) containing
    #'   one contrast in each row.
    #' @param m An optional vector containing the null value for each contrast.
    #'   Will be set to the zero vector of length `nrow(L_mat)` if not
    #'   specified.
    #' @return A [`Contrast`] object for the specified arguments and this
    #'   model fit.
    contrast = function(L_mat, m) Contrast$new(self, L_mat, m)
  ),

  # Private members ===========================================================

  private = list(

    # Data members ============================================================

    data = NULL,
    optres = NULL,
    theta_hat = NULL,
    V_theta_hat = NULL
  )
)
