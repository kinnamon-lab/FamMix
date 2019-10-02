# FamData.R - Definition for FamData R6 class and methods

#' Family data for mixed modeling
#'
#' A \code{R6} class that stores family data and a kinship matrix. Methods of
#' this class can then be used to fit genetic mixed models.
#'
#' @section Usage:
#' \strong{Constructors}
#' \preformatted{fdata <- FamData$new(data, family_id, indiv_id, proband, sex,
#'   maternal_id, paternal_id, mzgrp, dzgrp)}
#' \preformatted{fdata <- FamData$new(data, family_id, indiv_id, proband, phi)}
#' \strong{Accessors}
#' \preformatted{fdata$get_data()}
#' \preformatted{fdata$get_data_name()}
#' \preformatted{fdata$get_phi()}
#' \preformatted{fdata$get_var_map()}
#' \strong{Model Fitting}
#' \preformatted{fdata$lmm(formula, ...)}
#' \strong{Plotting}
#' \preformatted{fdata$plot_pedigree(famid)}
#'
#' @section Arguments:
#' \strong{Constructors}
#'
#' All constructor arguments other than \code{data} \emph{must be named}.
#' \describe{
#'   \item{\code{data}}{A \code{data.frame} or \code{data.table}
#'     containing all phenotypes and genetic variables of interest.}
#'   \item{\code{family_id}}{Character string containing the name of a numeric
#'     or character column in \code{data} containing a family identifier.}
#'   \item{\code{indiv_id}}{Character string containing the name of a numeric or
#'       character column in \code{data} containing a individual identifier.}
#'   \item{\code{proband}}{Character string containing the name of a numeric
#'     column in \code{data} containing 1 for probands and 0 otherwise.
#'     \code{NA} is not permitted.}
#'   \item{\code{sex}}{Character string containing the name of a numeric column
#'     in \code{data} containing the individual's sex (1 = "male", 2 = "female")
#'     or \code{NA} if unknown.}
#'   \item{\code{maternal_id}}{Character string containing the name of a column
#'     in \code{data} of the same type as \code{indiv_id} containing the
#'     \code{indiv_id} of the individual's mother or 0 (numeric),
#'     "" (character), or \code{NA} if a founder. See documentation for
#'     \code{\link[kinship2]{pedigree}}.}
#'   \item{\code{paternal_id}}{Character string containing the name of a column
#'     in \code{data} of the same type as \code{indiv_id} containing the
#'     \code{indiv_id} of the individual's father or or 0 (numeric),
#'     "" (character), or \code{NA} if a founder. See documentation for
#'     \code{\link[kinship2]{pedigree}}.}
#'   \item{\code{mzgrp}}{Character string containing the name of a numeric,
#'     character, or factor column in \code{data} containing the same value for
#'     all members of a monozygotic twin group. \code{NA} should be used for all
#'     individuals who are not monozygotic twins.}
#'   \item{\code{dzgrp}}{Character string containing the name of a numeric,
#'     character, or factor column in \code{data} containing the same value for
#'     all members of a dizygotic twin group. \code{NA} should be used for all
#'     individuals who are not dizygotic twins.}
#'   \item{\code{phi}}{A matrix coercible to a
#'     \code{\link[Matrix:dsCMatrix-class]{dsCMatrix}} containing the
#'     pairwise kinship coefficients for individuals. Note that
#'     \code{phi[i, j]} must contain the kinship coefficient between the
#'     individuals in \code{data[i, ]} and \code{data[j, ]}.}
#' }
#' \strong{Model Fitting}
#' \describe{
#'   \item{\code{formula}}{A \code{formula} object describing the mean model
#'     for a phenotype.}
#'   \item{\code{...}}{Additional arguments for \code{\link[stats]{nlminb}}.}
#' }
#' \strong{Plotting}
#' \describe{
#'   \item{\code{famid}}{A value of \code{family_id} from the original data
#'     set.}
#' }
#'
#' @section Methods:
#' \strong{Constructors}
#' \describe{
#'   \item{\code{$new()}}{Constructs a new instance of this class. See Usage
#'     for possible constructor signatures.}
#' }
#' \strong{Accessors}
#' \describe{
#'   \item{\code{$get_data()}}{Returns data member.}
#'   \item{\code{$get_data_name()}}{Returns name of input data set.}
#'   \item{\code{$get_phi()}}{Returns kinship matrix member.}
#'   \item{\code{$get_var_map()}}{Returns list mapping variable names in input
#'     data set to those in data member.}
#' }
#' \strong{Model Fitting}
#' \describe{
#'   \item{\code{$lmm()}}{Fits a linear mixed model. Returns a \code{list}
#'   object with the following contents:
#'   \describe{
#'     \item{\code{converge}}{Convergence code from
#'       \code{\link[stats]{nlminb}}. 0 is successful convergence.}
#'     \item{\code{message}}{Convergence message from
#'       \code{\link[stats]{nlminb}}.}
#'     \item{\code{max_grad}}{Maximum gradient component at solution.}
#'     \item{\code{deviance}}{-2 loglikelihood at solution.}
#'     \item{\code{parm_ests}}{\code{matrix} of parameter estimates, standard
#'       errors, and Wald tests/p-values.}
#'     \item{\code{vc_pl_cis}}{\code{matrix} of profile likelihood confidence
#'       intervals for variance components. \code{NA} if confidence level is
#'       not attainable or something went wrong in likelihood profiling.}
#'     \item{\code{lrts}}{\code{matrix} of likelihood ratio test results. Row
#'       names denote the null hypotheses, and columns indicate whether the
#'       null maximization converged, the X2 statistic, and the p-value
#'       calculated according to the approach of Self and Liang (1987).}
#'   }}
#' }
#' \strong{Plotting}
#' \describe{
#'   \item{\code{$plot_pedigree()}}{For objects initialized with pedigree data,
#'     allows plotting of individual pedigrees using \pkg{kinship2}.}
#' }
#'
#' @references
#'
#' Self SG, Liang K-Y. 1987. Asymptotic Properties of Maximum Likelihood
#'   Estimators and Likelihood Ratio Tests Under Nonstandard Conditions.
#'   \emph{J Am Stat Assoc} 82 (398): 605-10.
#'   \url{https://doi.org/10.2307/2289471}.
#'
#' @docType class
#' @name FamData-class
#' @aliases FamData
#' @importFrom kinship2 pedigree
#' @importFrom TMB MakeADFun
NULL

#' @export
FamData <- R6Class(
  classname = "FamData",
  parent_env = getNamespace("FamMix"),
  lock_class = TRUE,
  # Public members ============================================================
  public = list(
    # Constructor =============================================================
    initialize = function(data, ...) {
      data_name <- deparse(substitute(data))
      args <- list(...)
      # Pedigree data constructor ---------------------------------------------
      if (setequal(
        names(args),
        c(
          "family_id",
          "indiv_id",
          "maternal_id",
          "paternal_id",
          "sex",
          "proband",
          "mzgrp",
          "dzgrp"
        )
      )) {
        # Create variable name map
        private$var_map <- c(
          fmid = args$family_id,
          id = args$indiv_id,
          mid = args$maternal_id,
          pid = args$paternal_id,
          sex = args$sex,
          pr = args$proband,
          mzid = args$mzgrp,
          dzid = args$dzgrp
        )
        # Populate data member, rename to standard variable names for use
        # within FamData object, and sort by fmid, then id
        private$data <- data.table(data)
        setnames(
          private$data,
          old = private$var_map,
          new = names(private$var_map)
        )
        setkey(private$data, fmid, id)
        # Check variables
        if (
          !is.numeric(private$data$fmid) &
            !is.character(private$data$fmid)
        ) {
          stop(
            "family_id column contains improperly formatted data. ",
            "See ?FamData."
          )
        }
        if (!is.numeric(private$data$id) & !is.character(private$data$id)) {
          stop(
            "indiv_id column contains improperly formatted data. See ?FamData."
          )
        }
        if (class(private$data$mid) != class(private$data$id)) {
          stop(
            "maternal_id column contains improperly formatted data. ",
            "See ?FamData."
          )
        }
        if (class(private$data$pid) != class(private$data$id)) {
          stop(
            "paternal_id column contains improperly formatted data. ",
            "See ?FamData."
          )
        }
        if (
          !is.numeric(private$data$pr) |
            !setequal(unique(private$data$pr), c(0, 1))
        ) {
          stop(
            "proband column contains improperly formatted data. See ?FamData."
          )
        }
        if (
          !is.numeric(private$data$sex) |
            length(setdiff(unique(private$data$sex), c(1, 2, NA))) > 0
        ) {
          stop(
            "sex column contains improperly formatted data. See ?FamData."
          )
        }
        # Create relation data.table for twins
        rel <- NULL
        if (private$data[, any(!is.na(mzid))]) {
          rel <- private$data[
            !is.na(mzid),
            cbind(
              expand.grid(
                id1 = id,
                id2 = id,
                stringsAsFactors = FALSE
              ),
              code = 1
            ),
            keyby = .(fmid, mzid)
          ][id2 > id1][, mzid := NULL]
        }
        if (private$data[, any(!is.na(dzid))]) {
          rel <- rbind(
            rel,
            private$data[
              !is.na(dzid),
              cbind(
                expand.grid(
                  id1 = id,
                  id2 = id,
                  stringsAsFactors = FALSE
                ),
                code = 2
              ),
              keyby = .(fmid, dzid)
            ][id2 > id1][, dzid := NULL]
          )
        }
        # Create pedigreeList
        ped_args <- with(
          private$data,
          list(
            id = id,
            dadid = pid,
            momid = mid,
            sex = ifelse(is.na(sex), 3, sex),
            famid = fmid
          )
        )
        if (!is.null(rel)) {
          setnames(rel, old = "fmid", new = "famid")
          ped_args$relation <- rel
        }
        private$ped_list <-
          do.call(get("pedigree", getNamespace("kinship2")), ped_args)
        # Calculate autosomal kinship matrix
        private$phi <-
          kinship2::kinship(private$ped_list, chrtype = "autosome")
        # Data set and kinship matrix constructor -----------------------------
      } else if (
          setequal(names(args),
          c("family_id", "indiv_id", "proband", "phi"))
        ) {
        # Create variable name map
        private$var_map <- c(
          fmid = args$family_id,
          id = args$indiv_id,
          pr = args$proband
        )
        # Populate data member and rename to standard variable names for use
        # within FamData object. DO NOT SORT YET!
        private$data <- data.table(data)
        setnames(
          private$data,
          old = private$var_map,
          new = names(private$var_map)
        )
        # Check variables
        if (
          !is.numeric(private$data$fmid) &
            !is.character(private$data$fmid)
        ) {
          stop(
            "family_id column contains improperly formatted data. ",
            "See ?FamData."
          )
        }
        if (!is.numeric(private$data$id) & !is.character(private$data$id)) {
          stop(
            "indiv_id column contains improperly formatted data. See ?FamData."
          )
        }
        if (!is.numeric(private$data$pr) |
          !setequal(unique(private$data$pr), c(0, 1))) {
          stop(
            "proband column contains improperly formatted data. See ?FamData."
          )
        }
        # Cast kinship matrix argument to appropriate type
        phi <- as(args$phi, "symmetricMatrix")
        # Get "fmid/id" in original sort order from data member
        ids <- private$data[, paste(fmid, id, sep = "/")]
        # Make "fmid/id" the row and column names of kinship matrix
        dimnames(phi) <- list(ids, ids)
        # Sort data member member by fmid, then id
        setkey(private$data, fmid, id)
        # Get "fmid/id" from data member in new sort order
        sorted_ids <- private$data[, paste(fmid, id, sep = "/")]
        # Permute rows and columns of kinship matrix to be in this order and
        # assign to kinship matrix member
        private$phi <- phi[sorted_ids, sorted_ids]
        # Uncrecognized constructor -------------------------------------------
      } else {
        stop("Unrecognized constructor signature. Please see documentation.")
      }
    },
    # Accessors ===============================================================
    get_data = function() private$data,
    get_data_name = function() private$data_name,
    get_phi = function() private$phi,
    get_var_map = function() private$var_map,
    # Model fitting methods ===================================================
    lmm = function(formula, ...) {
      mod_data <- private$make_model_mats(formula)
      init_fits <- lm(formula, data = private$data[mod_data$incl])
      parameters <- list(
        betas = coef(init_fits),
        h2_g = 0,
        sigma2 = summary(init_fits)$sigma ^ 2
      )
      objfun <- TMB::MakeADFun(
        mod_data,
        parameters,
        DLL = "sing_asc_lmm",
        method = NULL,
        silent = TRUE
      )
      names(objfun$par) <-
        gsub("betas.", "", names(unlist(parameters)))
      parm_lower <- rep(-Inf, length(objfun$par))
      parm_upper <- rep(Inf, length(objfun$par))
      names(parm_lower) <- names(parm_upper) <- names(objfun$par)
      parm_lower["h2_g"] <- 0
      parm_upper["h2_g"] <- 1
      parm_lower["sigma2"] <- sqrt(.Machine$double.eps)
      # Transform problem by scaling parameters by square root of their Hessian
      # diagonal elements at the initial estimates, which should improve
      # conditioning of the Hessian and therefore numerical convergence. This
      # is a non-adaptive version of the adaptive scaling for unconstrained
      # optimization in section 4c of the PORT documentation included in the
      # references for nlminb.
      d <- sqrt(abs(diag(objfun$he())))
      if (!all(is.finite(d))) {
        # If there are any -Inf, Inf, NA, or NaN values, do not scale
        d <- rep_len(1, length(d))
      } else {
        # Otherwise, set any scaling factors that are numerically <= 0 to the
        # minimum of those that are numerically > 0
        d[d < sqrt(.Machine$double.eps)] <-
          min(d[d >= sqrt(.Machine$double.eps)])
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
      objfun_h2_g_zero <- TMB::MakeADFun(
        mod_data,
        parameters,
        map = list(h2_g = factor(NA)),
        DLL = "sing_asc_lmm",
        method = NULL,
        silent = TRUE
      )
      d_h2_g_zero <- sqrt(abs(diag(objfun_h2_g_zero$he())))
      if (!all(is.finite(d_h2_g_zero))) {
        d_h2_g_zero <- rep_len(1, length(d_h2_g_zero))
      } else {
        d_h2_g_zero[d_h2_g_zero < sqrt(.Machine$double.eps)] <-
          min(d_h2_g_zero[d_h2_g_zero >= sqrt(.Machine$double.eps)])
      }
      optres_h2_g_zero <- nlminb(
        objfun_h2_g_zero$par,
        objfun_h2_g_zero$fn,
        objfun_h2_g_zero$gr,
        objfun_h2_g_zero$he,
        lower = parm_lower[!names(parm_lower) %in% "h2_g"],
        upper = parm_upper[!names(parm_upper) %in% "h2_g"],
        scale = d_h2_g_zero,
        ...
      )
      report <- TMB::sdreport(objfun)
      pl_ci_handler <- function(parm) {
        res <- matrix(NA, nrow = 1, ncol = 2)
        dimnames(res) <- list(parm, c("lower", "upper"))
        res
      }
      vc_pl_cis <- rbind(
        tryCatch(
          confint(TMB::tmbprofile(
            objfun,
            "h2_g",
            parm.range = c(parm_lower["h2_g"], parm_upper["h2_g"]),
            trace = FALSE
          )),
          error = function(e) pl_ci_handler("h2_g"),
          warning = function(w) pl_ci_handler("h2_g")
        ),
        tryCatch(
          confint(TMB::tmbprofile(
            objfun,
            "sigma2",
            parm.range = c(parm_lower["sigma2"], parm_upper["sigma2"]),
            trace = FALSE
          )),
          error = function(e) pl_ci_handler("sigma2"),
          warning = function(w) pl_ci_handler("sigma2")
        )
      )
      lr_converge <- optres_h2_g_zero$convergence
      lr_X2 <- -2 * (-optres_h2_g_zero$objective - (-optres$objective))
      lr_p <- ifelse(abs(lr_X2 - 0) < .Machine$double.eps,
        1,
        0.5 * pchisq(lr_X2, df = 1, lower.tail = FALSE)
      )
      lrts <- cbind(lr_converge, lr_X2, lr_p)
      dimnames(lrts) <- list("h2_g = 0", c("converge", "X2", "Pr(>X2)"))
      res <- list(
        f_sizes = mod_data$f_sizes,
        converge = optres$convergence,
        message = optres$message,
        max_grad = max(report$gradient.fixed),
        deviance = 2 * optres$objective,
        parm_ests = summary(report, p.value = TRUE),
        vc_pl_cis = vc_pl_cis,
        lrts = lrts
      )
      res
    },
    # Plotting methods ========================================================
    plot_pedigree = function(famid) {
      if (!is.null(private$ped_list)) plot(private$ped_list[famid])
    }
  ),
  # Private members ===========================================================
  private = list(
    # Data members ============================================================
    data = NULL,
    data_name = NULL,
    ped_list = NULL,
    phi = NULL,
    var_map = NULL,
    # Utility methods =========================================================
    make_model_mats = function(formula) {
      if (private$data[pr == 1, .N, keyby = fmid][, max(N)] > 1) {
        stop(
          "Cannot have more than one proband per family under single ",
          "ascertainment"
        )
      }
      non_na <- which(apply(
        model.frame(formula, private$data, na.action = na.pass),
        1,
        function(x) all(!is.na(x))
      ))
      fam_w_proband <- private$data[non_na][pr == 1, fmid]
      fam_w_rels <- private$data[
        non_na,
        .(n_rel = sum(pr == 0)),
        keyby = fmid
      ][n_rel > 0, fmid]
      excl_fam <- setdiff(
        unique(private$data$fmid),
        intersect(fam_w_proband, fam_w_rels)
      )
      incl <- intersect(
        non_na,
        which(private$data$fmid %in% intersect(fam_w_proband, fam_w_rels))
      )
      mf <- model.frame(
        formula,
        private$data[incl],
        na.action = na.fail,
        drop.unused.levels = TRUE,
        proband = pr
      )
      probands <- which(mf$`(proband)` == 1)
      y <- as.numeric(model.response(mf))
      X <- model.matrix(formula, mf)
      f_sum <- private$data[incl, .N, keyby = fmid]
      f_sizes <- f_sum$N
      names(f_sizes) <- f_sum$fmid
      list(
        y = y,
        X = X,
        y_pr = y[probands],
        X_pr = X[probands, ],
        phi = private$phi[incl, incl],
        f_sizes = f_sizes,
        incl = incl
      )
    }
  )
)
