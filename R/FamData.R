# FamData.R - Definition for FamData R6 class and methods

#' Family data for modeling
#'
#' @description A \code{R6} class that stores family data and a kinship
#'   matrix. Methods of this class can then be used to fit regression models
#'   such as polygenic mixed models. The class currently supports only family
#'   data with one observation per subject and will check that this requirement
#'   is met.
#'
#' @export
FamData <- R6Class(
  "FamData",
  parent_env = getNamespace("FamModel"),
  lock_class = TRUE,

  # Public members ============================================================

  public = list(

    # Constructor =============================================================

    #' @description Constructs a new instance of this class. All constructor
    #'   arguments other than \code{data} \emph{must be named}. There are two
    #'   possible constructor signatures:
    #' \preformatted{
    #' FamData$new(data, family_id, indiv_id, proband, sex,
    #'             maternal_id, paternal_id, mzgrp, dzgrp)
    #' }
    #' \preformatted{FamData$new(data, family_id, indiv_id, proband, phi)}
    #' @param data A \code{data.frame} or \code{data.table} containing all
    #'   phenotypes and genetic variables of interest.
    #' @param ... The arguments required for a particular method signature;
    #'   see below for detailed descriptions.
    #' @param family_id Character string containing the name of a numeric or
    #'   character column in \code{data} containing a family identifier.
    #' @param indiv_id Character string containing the name of a numeric or
    #'   character column in \code{data} containing a individual identifier.
    #' @param proband Character string containing the name of a numeric column
    #'   in \code{data} containing 1 for probands and 0 otherwise. \code{NA} is
    #'   not permitted.
    #' @param sex Character string containing the name of a numeric column in
    #'   \code{data} containing the individual's sex (1 = "male", 2 = "female")
    #'   or \code{NA} if unknown.
    #' @param maternal_id Character string containing the name of a column in
    #'   \code{data} of the same type as \code{indiv_id} containing the
    #'   \code{indiv_id} of the individual's mother or 0 (numeric), ""
    #'   (character), or \code{NA} if a founder. See documentation for
    #'   \code{\link[kinship2]{pedigree}}.
    #' @param paternal_id Character string containing the name of a column in
    #'   \code{data} of the same type as \code{indiv_id} containing the
    #'   \code{indiv_id} of the individual's father or or 0 (numeric), ""
    #'   (character), or \code{NA} if a founder. See documentation for
    #'   \code{\link[kinship2]{pedigree}}.
    #' @param mzgrp Character string containing the name of a numeric,
    #'   character, or factor column in \code{data} containing the same value
    #'   for all members of a monozygotic twin group. \code{NA} should be used
    #'   for all individuals who are not monozygotic twins.
    #' @param dzgrp Character string containing the name of a numeric,
    #'   character, or factor column in \code{data} containing the same value
    #'   for all members of a dizygotic twin group. \code{NA} should be used for
    #'   all individuals who are not dizygotic twins.
    #' @param phi A matrix coercible to a
    #'   \code{\link[Matrix:dsCMatrix-class]{dsCMatrix}} containing the pairwise
    #'   kinship coefficients for individuals. Note that \code{phi[i, j]} must
    #'   contain the kinship coefficient between the individuals in
    #'   \code{data[i, ]} and \code{data[j, ]}.
    initialize = function(data, ...) {
      private$data_name <- deparse(substitute(data))
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
          fmid = args[["family_id"]],
          id = args[["indiv_id"]],
          mid = args[["maternal_id"]],
          pid = args[["paternal_id"]],
          sex = args[["sex"]],
          pr = args[["proband"]],
          mzid = args[["mzgrp"]],
          dzid = args[["dzgrp"]]
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
        # Make sure that there is only one row per individual within each fmid
        if (anyDuplicated(private$data[, .(fmid, id)])) {
          stop(
            "Multiple rows found for the same individual within a fmid"
          )
        }
        # Check variables
        if (
          !is.numeric(private$data[, fmid]) &
          !is.character(private$data[, fmid])
        ) {
          stop(
            "family_id column contains improperly formatted data. ",
            "See ?FamData."
          )
        }
        if (
          !is.numeric(private$data[, id]) &
          !is.character(private$data[, id])
        ) {
          stop(
            "indiv_id column contains improperly formatted data. See ?FamData."
          )
        }
        if (class(private$data[, mid]) != class(private$data[, id])) {
          stop(
            "maternal_id column contains improperly formatted data. ",
            "See ?FamData."
          )
        }
        if (class(private$data[, pid]) != class(private$data[, id])) {
          stop(
            "paternal_id column contains improperly formatted data. ",
            "See ?FamData."
          )
        }
        if (
          !is.numeric(private$data[, pr]) |
          !setequal(unique(private$data[, pr]), c(0, 1))
        ) {
          stop(
            "proband column contains improperly formatted data. See ?FamData."
          )
        }
        if (
          !is.numeric(private$data[, sex]) |
          length(setdiff(unique(private$data[, sex]), c(1, 2, NA))) > 0
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
          ped_args[["relation"]] <- rel
        }
        private$ped_list <-
          do.call(kinship2::pedigree, ped_args)
        # Calculate autosomal kinship matrix
        private$phi <-
          kinship2::kinship(private$ped_list, chrtype = "autosome")
        # Data set and kinship matrix constructor -----------------------------
      } else if (
        setequal(
          names(args),
          c("family_id", "indiv_id", "proband", "phi")
        )
      ) {
        # Create variable name map
        private$var_map <- c(
          fmid = args[["family_id"]],
          id = args[["indiv_id"]],
          pr = args[["proband"]]
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
          !is.numeric(private$data[, fmid]) &
          !is.character(private$data[, fmid])
        ) {
          stop(
            "family_id column contains improperly formatted data. ",
            "See ?FamData."
          )
        }
        if (
          !is.numeric(private$data[, id]) &
          !is.character(private$data[, id])
        ) {
          stop(
            "indiv_id column contains improperly formatted data. See ?FamData."
          )
        }
        if (
          !is.numeric(private$data[, pr]) |
          !setequal(unique(private$data[, pr]), c(0, 1))
        ) {
          stop(
            "proband column contains improperly formatted data. See ?FamData."
          )
        }
        # Make sure that there is only one row per individual within each fmid
        if (anyDuplicated(private$data[, .(fmid, id)])) {
          stop(
            "Multiple rows found for the same individual within a fmid"
          )
        }
        # Cast kinship matrix argument to appropriate type
        phi <- as(args[["phi"]], "symmetricMatrix")
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

    #' @description Returns data member.
    get_data = function() private$data,

    #' @description Returns name of input \code{data.frame} or
    #'   \code{data.table}.
    get_data_name = function() private$data_name,

    #' @description Returns kinship matrix member.
    get_phi = function() private$phi,

    #' @description Returns list mapping variable names in input data set in
    #'   \code{data} argument to those in data member.
    get_var_map = function() private$var_map,

    # Model fitting methods ===================================================

    #' @description Fits a linear mixed model to ascertained families
    #'
    #' @details Fits a linear mixed model to family data assuming that the
    #'    family has been ascertained through a single proband meeting certain
    #'    criteria. To do this, it uses the appropriate multivariate normal
    #'    likelihood conditional on the observed data in the proband (see the
    #'    Beaty et al. reference). Families without probands or more than one
    #'    proband are excluded. The model assumes an additive polygenic effect,
    #'    parameterized in terms of heritability, and no shared environmental
    #'    effect. The heritability and total variance can be allowed to vary
    #'    across groups of families (see \code{formula} argument), but, within
    #'    each group, all founders in the family are assumed to be drawn from
    #'    the same randomly mating population. See the Lange and Boerwinkle et
    #'    al. references for additional details on this class of mixed models.
    #'
    #' @param formula A \code{\link[Formula]{Formula}} object describing the
    #'   model for a phenotype. The formula is of the form
    #'   \code{y ~ mean | group}, where \code{y} (required) is the outcome,
    #'   \code{mean} (optional) specifies the mean model that includes an
    #'   intercept by default, and \code{| group} (optional) specifies a factor
    #'   term that can be used assign families to homogeneous groups with
    #'   different variance components. Multiple variables can be combined into
    #'   a single factor using \code{:}.
    #' @param ... Additional parameters to pass to the optimization function
    #'   (currently \code{\link{nlminb}}).
    #'
    #' @return Returns a \code{\link{FamLMMFit}} object.
    #'
    #' @references
    #' Beaty TH, Liang KY, Rao DC. Robust inference for variance components
    #'   models in families ascertained through probands: I. Conditioning on
    #'   proband's phenotype. \emph{Genet Epidemiol}. 1987;4(3):203-10.
    #'
    #' Boerwinkle E, Chakraborty R, Sing CF. The use of measured genotype
    #'   information in the analysis of quantitative phenotypes in man. I.
    #'   Models and analytical methods. \emph{Ann Hum Genet}.
    #'   1986;50(Pt 2):181-94.
    #'
    #' Lange K. \emph{Mathematical and statistical methods for genetic
    #'   analysis}. 2nd ed. New York: Springer; 2002.
    lmm = function(formula, ...) {
      if (!Formula::is.Formula(formula)) formula <- Formula::Formula(formula)
      formula_dims <- length(formula)
      if (formula_dims[1] != 1) {
        stop(
          "A single response variable must be specified in the supplied ",
          "formula"
        )
      } else if (formula_dims[2] > 2) {
        stop(
          "More than 2 RHS components cannot be specified in formula"
        )
      } else if (formula_dims[2] == 2) {
        formula <- update(formula, . ~ . | 0 + .)
      } else if (formula_dims[2] == 1) {
        formula <- update(formula, . ~ . | 1)
      }
      mod_data <- private$make_model_mats_lmm(formula)
      init_fits <- lm(
        formula(formula, rhs = 1),
        data = private$data[mod_data[["incl"]]]
      )
      parameters <- list(
        betas = coef(init_fits),
        h2_a = rep(0, ncol(mod_data[["f_pops"]])),
        sigma2 = rep(
          summary(init_fits)[["sigma"]] ^ 2,
          ncol(mod_data[["f_pops"]])
        )
      )
      if (ncol(mod_data[["f_pops"]]) > 1) {
        names(parameters[["h2_a"]]) <- names(parameters[["sigma2"]]) <-
          colnames(mod_data[["f_pops"]])
      }
      objfun <- TMB::MakeADFun(
        mod_data,
        parameters,
        DLL = "sing_asc_lmm",
        method = NULL,
        silent = TRUE
      )
      attr(objfun, "type") <- "TMB"
      parm_lower <- rep(-Inf, length(objfun[["par"]]))
      parm_upper <- rep(Inf, length(objfun[["par"]]))
      names(parm_lower) <- names(parm_upper) <- names(objfun[["par"]])
      parm_lower[names(parm_lower) == "h2_a"] <- 0
      parm_upper[names(parm_upper) == "h2_a"] <- 1
      parm_lower[names(parm_lower) == "sigma2"] <-
        sqrt(.Machine[["double.eps"]])
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
      attr(optres, "optfun") <- "nlminb"
      # Return new FamLMMFit object with results
      FamLMMFit$new(self, formula, objfun, optres)
    },

    # Plotting methods ========================================================

    #' @description For objects initialized with pedigree data, allows plotting
    #'   of individual pedigrees using \pkg{kinship2}.
    #' @param famid Argument matching the value of a single family in the
    #'   \code{family_id} column used in the constructor.
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

    make_model_mats_lmm = function(formula) {
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
        unique(private$data[, fmid]),
        intersect(fam_w_proband, fam_w_rels)
      )
      incl <- sort(intersect(
        non_na,
        private$data[, which(fmid %in% intersect(fam_w_proband, fam_w_rels))]
      ))
      mf <- model.frame(
        formula,
        private$data[incl],
        na.action = na.fail,
        drop.unused.levels = TRUE,
        fmid = fmid,
        id = id,
        proband = pr
      )
      # Check that selection of observations has maintained original sort order
      # of data member, which is ascending fmid, then ascending id
      if (
        !isTRUE(all.equal(
          order(mf[["(fmid)"]], mf[["(id)"]]),
          seq(1, nrow(mf))
        ))
      ) {
        stop("Model frame not ordered by ascending fmid, then ascending id")
      }
      y <- as.numeric(model.response(mf))
      X <- model.matrix(formula, mf, rhs = 1)
      pops <- model.matrix(formula, mf, rhs = 2)
      if (any(apply(pops, 1, sum) != 1)) {
        stop(
          "Grouping formula for variance parameters does not result in a ",
          "mutually exclusive and exhaustive partition into groups."
        )
      }
      f_sum <- data.table(
        fmid = mf[["(fmid)"]],
        proband = mf[["(proband)"]],
        pops
      )[,
        .(pr_idx = .I[proband == 1], .N),
        keyby = c("fmid", colnames(pops))
      ]
      if (anyDuplicated(f_sum[, fmid])) {
        stop(
          "The following families have heterogeneous population IDs: ",
          paste0(f_sum[duplicated(f_sum[, fmid]), fmid], collapse = ", "),
          ".\nAll members of a family included in the model must come from ",
          "the same population."
        )
      }
      f_sizes <- f_sum[, N]
      f_pr_idxs <- f_sum[, pr_idx]
      names(f_pr_idxs) <- names(f_sizes) <- f_sum[, fmid]
      f_pops <- as.matrix(
        f_sum[, colnames(pops), with = FALSE],
        nrow = nrow(f_sum),
        ncol = ncol(pops),
        rownames = f_sum[, fmid],
        colnames = colnames(pops)
      )
      list(
        y = y,
        X = X,
        phi = private$phi[incl, incl],
        f_sizes = f_sizes,
        f_pr_idxs = f_pr_idxs,
        f_pops = f_pops,
        incl = incl
      )
    }
  )
)
