# FamModelFit.R - Definition for FamModelFit R6 abstract class and methods

#' Abstract class for results of a fitted \code{FamModel} model
#'
#' @description An abstract R6 class defining the pattern for all
#'   \code{\link{FamModelFit}} subclasses created by fitting a model to a
#'   \code{\link{FamData}} object using one of its model-fitting methods.
#'
#' @export
FamModelFit <- R6Class(
  "FamModelFit",
  parent_env = getNamespace("FamModel"),
  lock_class = TRUE,

  # Public members ============================================================

  public = list(

    # Constructor =============================================================

    #' @description Stub constructor to prevent instantiation of this abstract
    #'   class.
    initialize = function() {
      stop(paste(
        class(self)[1],
        "is an abstract class that cannot be initialized."
      ))
    },

    # Accessors ===============================================================

    #' @description Returns the \code{\link{FamData}} object that produced this
    #'   model fit.
    get_data = function() private$data,

    #' @description Returns \code{list} of optimization results.
    get_optres = function() private$optres,

    #' @description Returns \code{theta_hat}, the vector of parameter
    #'   estimates.
    get_theta_hat = function() private$theta_hat,

    #' @description Returns \code{V_theta_hat}, the covariance matrix of the
    #'   parameter estimates.
    get_V_theta_hat = function() private$V_theta_hat,

    # Print method ============================================================

    #' @description Formatted printing of the \code{FamModelFit} object.
    #'
    #' @param ... Arguments passed on to \code{\link{print_ests}}.
    print = function(...) {
      print_ests(private$theta_hat, private$V_theta_hat, ...)
    },

    # Inferential methods =====================================================

    #' @description Create a new \code{\link{Contrast}} object.
    #' @param L_mat A contrast vector (1 df) or \code{matrix} (>1 df) containing
    #'   one contrast in each row.
    #' @param m An optional vector containing the null value for each contrast.
    #'   Will be set to the zero vector of length \code{nrow(L_mat)} if not
    #'   specified.
    #' @return A \code{\link{Contrast}} object for the specified arguments and
    #'   this model fit.
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
