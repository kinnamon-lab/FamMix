# FamData.R - Definition for FamData R6 class and methods

#' Family data for modeling
#'
#' @description An R6 class that stores family data and a kinship matrix.
#'
#' @details Methods of this class can be used to fit regression models such as
#'   polygenic mixed models. The class currently supports only family data with
#'   one observation per subject and will check that this requirement is met.
#'
#' Note that the private `data` member of the instance (returned by the
#'   [`FamData$get_data()`](#method-get_data) method) will always be sorted by
#'   ascending `family_id` and ascending `indiv_id` within `family_id`
#'   regardless of the sort order of the input data set. Also, pedigree columns
#'   in the input data set will be renamed in the `data` member; this mapping
#'   can be returned using the [`FamData$get_var_map()`](#method-get_var_map)
#'   method.
#'
#' @export
FamData <- R6Class(
  "FamData",
  parent_env = getNamespace("FamModel"),
  lock_class = TRUE,

  # Public members ============================================================

  public = list(

    # Constructor =============================================================

    #' @description Constructs a new instance of this class.
    #'
    #' @details All constructor arguments other than `data`
    #'   *must be named*. There are two possible constructor signatures:
    #'   ```
    #'   FamData$new(data, family_id, indiv_id, proband, sex,
    #'               maternal_id, paternal_id, mzgrp, dzgrp)
    #'
    #'   FamData$new(data, family_id, indiv_id, proband, phi)
    #'   ```
    #' @param data A `data.frame` or [`data.table`] containing all
    #'   phenotypes and genetic variables of interest.
    #' @param ... The arguments required for a particular method signature;
    #'   see below for detailed descriptions.
    #' @param family_id Character string containing the name of a numeric or
    #'   character column in `data` containing a family identifier.
    #' @param indiv_id Character string containing the name of a numeric or
    #'   character column in `data` containing a individual identifier.
    #' @param proband Character string containing the name of a numeric column
    #'   in `data` containing 1 for probands and 0 otherwise. `NA` is  not
    #'   permitted.
    #' @param sex Character string containing the name of a numeric column in
    #'   `data` containing the individual's sex (1 = "male", 2 = "female")
    #'   or `NA` if unknown.
    #' @param maternal_id Character string containing the name of a column in
    #'   `data` of the same type as `indiv_id` containing the `indiv_id` of the
    #'   individual's mother or 0 (numeric), "" (character), or `NA` if a
    #'   founder. See documentation for [`kinship2::pedigree`].
    #' @param paternal_id Character string containing the name of a column in
    #'   `data` of the same type as `indiv_id` containing the `indiv_id` of the
    #'   individual's father or or 0 (numeric), "" (character), or `NA` if a
    #'   founder. See documentation for [`kinship2::pedigree`].
    #' @param mzgrp Character string containing the name of a numeric,
    #'   character, or factor column in `data` containing the same value for all
    #'   members of a monozygotic twin group. `NA` should be used for all
    #'   individuals who are not monozygotic twins.
    #' @param dzgrp Character string containing the name of a numeric,
    #'   character, or factor column in `data` containing the same value for all
    #'   members of a dizygotic twin group. `NA` should be used for all
    #'   individuals who are not dizygotic twins.
    #' @param phi A matrix coercible to a [`Matrix::dsCMatrix-class`] containing
    #'   the pairwise kinship coefficients for individuals. Note that
    #'   `phi[i, j]` must contain the kinship coefficient between the
    #'   individuals in `data[i, ]` and `data[j, ]`.
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
        private$data <- if (is.data.table(data)) {
          copy(data)
        } else {
          data.table(data)
        }
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
        private$unique_id <- !anyDuplicated(private$data[, id])
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
        bad_peds <- data.table(
          with(
            private$data,
            kinship2::familycheck(
              famid = fmid,
              id = if (private$unique_id) id else paste(fmid, id, sep = "/"),
              father.id = if (private$unique_id) {
                pid
              } else {
                paste(fmid, pid, sep = "/")
              },
              mother.id = if (private$unique_id) {
                mid
              } else {
                paste(fmid, mid, sep = "/")
              }
            )
          )
        )[unrelated > 0 | split != 1 | join != 0]
        if (nrow(bad_peds) > 0) {
          cat("Incorrectly specified pedigrees:\n")
          print(bad_peds)
          stop("Incorrectly specified pedigrees found. See above output.")
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
        # Check that order of individuals in kinship matrix is the same as in
        # data member. Note that kinship2::kinship uses "id" unless it is not
        # unique across families, in which case it uses "fmid/id"
        if (
          !(
            isTRUE(all.equal(
              private$data[, as.character(id)], rownames(private$phi)
            )) &&
            isTRUE(all.equal(
              private$data[, as.character(id)], colnames(private$phi)
            ))
          ) &&
          !(
            isTRUE(all.equal(
              private$data[,  paste(fmid, id, sep = "/")], rownames(private$phi)
            )) &&
            isTRUE(all.equal(
              private$data[,  paste(fmid, id, sep = "/")], colnames(private$phi)
            ))
          )
        ) {
          stop("Subject order in phi matrix does not match data member")
        }
        # Store IDs of families with consanguinity, determined as those
        # containing individuals with phi_ii > 0.5. This follows because phi_ii
        # = 0.5(1 + f_i) and f_i > 0 iff the individual is inbred (i.e., parents
        # are related)
        private$consang <- private$data[
          Matrix::diag(private$phi) > 0.5,
          unique(fmid)
        ]
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
        private$unique_id <- !anyDuplicated(private$data[, id])
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
        # Get "id" or "fmid/id" in original sort order from data
        # member. Mimicking behavior of kinship2::kinship, which uses "id"
        # unless it is not unique across families
        ids <- if (private$unique_id) {
          private$data[, as.character(id)]
        } else {
          private$data[, paste(fmid, id, sep = "/")]
        }
        # Make "id" or "fmid/id" the row and column names of kinship matrix
        dimnames(phi) <- list(ids, ids)
        # Sort data member member by fmid, then id
        setkey(private$data, fmid, id)
        # Get "id" or "fmid/id" from data member in new sort order
        sorted_ids <- if (private$unique_id) {
          private$data[, as.character(id)]
        } else {
          private$data[, paste(fmid, id, sep = "/")]
        }
        # Permute rows and columns of kinship matrix to be in this order and
        # assign to kinship matrix member
        private$phi <- phi[sorted_ids, sorted_ids]
        # Store IDs of families with consanguinity, determined as those
        # containing individuals with phi_ii > 0.5. This follows because phi_ii
        # = 0.5(1 + f_i) and f_i > 0 iff the individual is inbred (i.e., parents
        # are related)
        private$consang <- private$data[
          Matrix::diag(private$phi) > 0.5,
          unique(fmid)
        ]
        # Uncrecognized constructor -------------------------------------------
      } else {
        stop("Unrecognized constructor signature. Please see documentation.")
      }
    },

    # Accessors ===============================================================

    #' @description Returns a copy of data member (to prevent accidental
    #'   modification of data member [`data.table`] by reference).
    get_data = function() copy(private$data),

    #' @description Returns name of input `data.frame` or [`data.table`].
    get_data_name = function() private$data_name,

    #' @description Returns vector containing family IDs in which consanguinity
    #'   was found (based on individual self-kinship coefficient greater than
    #'   0.5).
    get_consang = function() private$consang,

    #' @description Returns kinship matrix member.
    get_phi = function() private$phi,

    #' @description Returns list mapping variable names in input data set in
    #'   `data` argument to those in data member.
    get_var_map = function() private$var_map,

    # Print method ============================================================

    #' @description Prints information about contents of `FamData` objects
    #'   with nice formatting.
    #'
    print = function() {
      cat("\n==========================\n")
      cat("FamData OBJECT INFORMATION\n")
      cat("==========================\n")
      cat("Input data set: ", private$data_name, "\n", sep = "")
      cat("Families: ", private$data[, uniqueN(fmid)], "\n", sep = "")
      cat("Subjects: ", uniqueN(private$data[, .(fmid, id)]), "\n", sep = "")
      cat(
        "Kinship matrix: ",
        if ("mid" %in% names(private$var_map)) {
          "Calculated from data set"
        } else {
          "Supplied to constructor"
        },
        "\n", sep = ""
      )
      cat("Consanguinous families: ")
      if (length(private$consang) == 0) {
        cat("None\n")
      } else {
        cat("\n")
        cat(private$consang, sep = "\n")
      }
      cat("Pedigree variable mapping:\n")
      cat(
        sprintf("  %-6s : %s", names(private$var_map), private$var_map),
        sep = "\n"
      )
      cat(
        "Unique subject ID: ", if (private$unique_id) "id" else "fmid + id",
        "\n", sep = ""
      )
      cat("Available variables:\n")
      cat(
        strwrap(paste0(
          setdiff(names(private$data), names(private$var_map)), collapse = ", "
        )),
        sep = "\n"
      )
    },

    # Model fitting methods ===================================================

    #' @description Fits a linear mixed model to ascertained families
    #'
    #' @details See `vignette("linear_mixed_models")` for background,
    #'   implementation details, and references.
    #'
    #' @param formula A [Formula::Formula()] object describing the model for a
    #'   phenotype. The formula is of the form `y ~ mean | group`, where `y`
    #'   (required) is the outcome, `mean` (optional) specifies the mean model
    #'   that includes an intercept by default, and `| group` (optional)
    #'   specifies a factor term that can be used assign families to homogeneous
    #'   groups with different variance components. Multiple variables can be
    #'   combined into a single factor using `:`.
    #' @param ... Additional parameters to pass to the `control` list for
    #'   [optim()] with `method = "L-BGFS-B"`. Note that `parscale` and
    #'   `fnscale` cannot be modified.
    #'
    #' @return Returns a [`FamLMMFit`] object.
    #'
    lmm = function(formula, ...) {
      if (!Formula::is.Formula(formula)) formula <- Formula::Formula(formula)
      # Set Formula environment to base namespace to prevent using variables not
      # explicitly included in data argument of model.frame
      environment(formula) <- .BaseNamespaceEnv
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
      mod_data[["model"]] <- "sing_asc_lmm"
      init_fits <- lm(
        formula(formula, rhs = 1),
        data = private$data[mod_data[["incl_ids"]], on = .(fmid, id)][pr == 0]
      )
      # Make sure name patterns reserved for variance components are not used in
      # mean model
      if (any(grepl("^(h2_a|sigma)(\\..+)?$", names(coef(init_fits))))) {
        stop(
          "Variables in mean model cannot have names matching the pattern ",
          "'^(h2_a|sigma)(\\..+)?$'"
        )
      }
      parameters <- list(
        betas = coef(init_fits),
        h2_a = rep(0, ncol(mod_data[["f_pops"]])),
        sigma = rep(
          summary(init_fits)[["sigma"]],
          ncol(mod_data[["f_pops"]])
        )
      )
      if (ncol(mod_data[["f_pops"]]) > 1) {
        names(parameters[["h2_a"]]) <- names(parameters[["sigma"]]) <-
          colnames(mod_data[["f_pops"]])
      }
      objfun <- TMB::MakeADFun(
        mod_data,
        parameters,
        DLL = "FamModel",
        method = NULL,
        silent = TRUE
      )
      attr(objfun, "type") <- "LMM-SA"
      optres <- lmm_optim(objfun, ...)
      # Return new FamLMMFit object with results
      FamLMMFit$new(self, formula, objfun, optres)
    },

    # Plotting methods ========================================================

    #' @description For objects initialized with pedigree data, allows plotting
    #'   of individual pedigrees using [`kinship2::plot.pedigree`].
    #' @param famid Argument matching the value of a single family in the
    #'   `family_id` column used in the constructor.
    plot_pedigree = function(famid) {
      if (!is.null(private$ped_list)) plot(private$ped_list[famid])
    }
  ),

  # Private members ===========================================================

  private = list(

    # Data members ============================================================

    consang = NULL,
    data = NULL,
    data_name = NULL,
    ped_list = NULL,
    phi = NULL,
    unique_id = NULL,
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
      incl_ids <- mf[, c("(fmid)", "(id)")]
      names(incl_ids) <- c("fmid", "id")
      incl_phi_ids <- if (private$unique_id) {
        with(mf, as.character(`(id)`))
      } else {
        with(mf, paste(`(fmid)`, `(id)`, sep = "/"))
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
        phi = private$phi[incl_phi_ids, incl_phi_ids],
        f_sizes = f_sizes,
        f_pr_idxs = f_pr_idxs,
        f_pops = f_pops,
        incl_ids = as.list(incl_ids)
      )
    }
  )
)
