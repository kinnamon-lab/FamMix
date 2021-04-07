# Tests for linear mixed model functionality
library(future.apply, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(vdiffr, quietly = TRUE)

# Note that Mendel 16.0 uses N rather than N-1 for SD divisor in standardization
lmna_data <- copy(lmna_nonseg)[,
  `:=`(
    female = as.integer(sex == 2),
    age_echo_std = (age_echo_yrs - mean(age_echo_yrs, na.rm = TRUE)) /
      (
        sd(age_echo_yrs, na.rm = TRUE) *
          sqrt((sum(!is.na(age_echo_yrs)) - 1) / sum(!is.na(age_echo_yrs)))
      )
  )
]
lmna_fd <- FamData$new(
  lmna_data,
  family_id = "family_ID",
  indiv_id = "individual_ID",
  proband = "proband",
  sex = "sex",
  maternal_id = "maternal_ID",
  paternal_id = "paternal_ID",
  mzgrp = "mzpair",
  dzgrp = "dzpair"
)
lmna_lvef_model <- lmna_fd$lmm(
  lvef ~ female + age_echo_std + I(n_lmna_vars > 0) + I(n_oth_vars > 0)
)
lmna_lvedd_z_model <- lmna_fd$lmm(
  lvedd_z ~ female + age_echo_std + I(n_lmna_vars > 0) + I(n_oth_vars > 0)
)

test_that("Can reproduce univariate model results from Mendel 16.0", {
  # Always check serialized objects underlying output previously checked against
  # Mendel 16.0 results with appropriate numerical tolerance
  check_model <- function(model) {
    # Optimization results
    expect_snapshot_value(
      model$get_optres(), style = "serialize", cran = TRUE,
      tolerance = testthat_tolerance()
    )
    # Parameter estimates
    expect_snapshot_value(
      model$get_theta_hat(), style = "serialize", cran = TRUE,
      tolerance = testthat_tolerance()
    )
    # Covariance matrix
    expect_snapshot_value(
      model$get_V_theta_hat(), style = "serialize", cran = TRUE,
      tolerance = testthat_tolerance()
    )
    # Likelihood ratio tests, skipping items relating to gradient criteria
    # that can have relative differences larger than tolerance even when
    # all other model results are nearly identical
    lrts <- model$get_h2_a_lrts(print = FALSE)
    expect_snapshot_value(
      lrts[,
        grep("_grad", names(lrts), invert = TRUE, value = TRUE),
        with = FALSE
      ],
      style = "serialize", cran = TRUE, tolerance = testthat_tolerance()
    )
    res <- model$get_model_res()[
      pr == 0, .(fmid, id, r_star_hat, c_star_hat, c_star_hat_df, p_c_star_hat)
    ]
    # Goodness-of-fit statsitics
    expect_snapshot_value(
      res[, .(fmid, id, r_star_hat_2 = r_star_hat ^ 2)],
      style = "serialize",
      cran = TRUE,
      tolerance = testthat_tolerance()
    )
    expect_snapshot_value(
      unique(res[, .(fmid, c_star_hat, c_star_hat_df, p_c_star_hat)]),
      style = "serialize",
      cran = TRUE,
      tolerance = testthat_tolerance()
    )
  }

  check_model(lmna_lvef_model)
  check_model(lmna_lvedd_z_model)
})

test_that("Can reproduce univariate model output from Mendel 16.0", {
  # Check against snapshot output and data.tables of r_star_hat ^ 2 and
  # c_star_hat previously checked against Mendel 16.0 results in development
  # environment. Only snapshots of serialized objects underlying these accepted
  # snapshot outputs are compared on CI and CRAN to avoid failures due to minor
  # numerical differences across platforms
  skip_on_ci()
  skip_on_cran()

  check_model_output <- function(model) {
    expect_snapshot_output(model$print())
    res <- model$get_model_res()[
      pr == 0, .(fmid, id, r_star_hat, c_star_hat, c_star_hat_df, p_c_star_hat)
    ]
    expect_snapshot_output(
      print(res[, .(fmid, id, r_star_hat_2 = r_star_hat ^ 2)])
    )
    expect_snapshot_output(
      print(unique(res[, .(fmid, c_star_hat, c_star_hat_df, p_c_star_hat)]))
    )
  }

  check_model_output(lmna_lvef_model)
  check_model_output(lmna_lvedd_z_model)
})

# Skip on CRAN or if FAMMODEL_NO_SIMS environment variable is set to
# a value recognized as TRUE
skip_on_cran()
skip_if(
  isTRUE(as.logical(Sys.getenv("FAMMODEL_NO_SIMS"))),
  "FAMMODEL_NO_SIMS environment variable is TRUE; skipping simulation tests"
)

# Running testthat.R as part of R CMD CHECK sets the R_CMD_CHECK environment
# variable to "true" to signify that a temporary install of the package in a
# temporary library is not necessary for future.apply. For interactive use,
# this temporary install is necessary
if (!isTRUE(as.logical(Sys.getenv("R_CMD_CHECK")))) {
  withr::local_temp_libpaths()
  devtools::install(
    reload = FALSE, quick = TRUE, build = TRUE, quiet = TRUE,
    upgrade = "never"
  )
}

# Make sure R CMD CHECK core limit is not enforced by future
Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "false")

# Only perform the first 2 simulation tests when on CI
n_reps <- 1000
n_tests <- if (isTRUE(as.logical(Sys.getenv("CI")))) 2 else 5
RNGkind("L'Ecuyer-CMRG")
set.seed(271828182)
plan(multisession)

# nolint start
# Simulation uses a population model with exactly 2 populations. In a random
# individual selected from a particular population:
# # very rare variant alleles: n_rv ~ Pois(lambda_rv)
# Random covariate (iid across family members): x ~ N(0, 1)
# Phenotype: y | n_rv, x ~ N(beta_rv * n_rv + beta_x * x, sigma ^ 2)
# => y | n_rv ~ N(beta_rv * n_rv, sigma ^ 2 + beta_x ^ 2)
# Ascertainment threshold: The bottom p_thresh*100% of each population in
# terms of the marginal distribution of y qualifies as probands
# Residual narrow-sense heritability: h2_a
# nolint end
pop_parms <- data.table(
  pop = c(1, 2),
  lambda_rv = c(0.35, 0.45),
  beta_rv = c(-0.4, -0.3),
  beta_x = c(0.2, 0.1),
  sigma = c(0.8, 1.2),
  p_thresh = c(0.05, 0.10),
  h2_a = c(0, 0.10)
)

# p_thresh can be used to solve for the actual thresholds, tau, in y for
# each population
pop_parms[,
  tau := {
    pop_tau <- function(tau) {
      n_rv <- 0:100
      mu_n_rv <- beta_rv * n_rv
      sum(
        exp(
          pnorm(tau, mu_n_rv, sqrt(sigma ^ 2 + beta_x ^ 2), log.p = TRUE) +
            dpois(n_rv, lambda_rv, log = TRUE)
        )
      ) - p_thresh
    }
    uniroot(
      pop_tau, interval = c(qnorm(1e-4), qnorm(1 - 1e-4)),
      tol = sqrt(.Machine[["double.eps"]])
    )[["root"]]
  },
  by = pop
]

# Prepare test data comprising pedigrees with FDRs and proband spouses only
data(minnbreast, package = "kinship2", envir = environment())
mb <- data.table(minnbreast)[,
  .(famid, id, fatherid, motherid, proband, sex, mzgrp = NA, dzgrp = NA)
]
probands <- mb[proband == 1, id]
mothers_pr <- mb[proband == 1, motherid]
fathers_pr <- mb[proband == 1, fatherid]
full_sibs_pr <- mb[
  mb[proband == 1, .(famid, motherid, fatherid)],
  .(id, proband),
  on = .(famid, motherid, fatherid)
][proband == 0, id]
children_pr <- mb[
  motherid %in% probands | fatherid %in% probands,
  id
]
spouses_pr <- mb[
  id %in% children_pr,
  unique(ifelse(motherid %in% probands, fatherid, motherid))
]
mb[,
  `:=`(
    mother_pr = as.numeric(id %in% mothers_pr),
    father_pr = as.numeric(id %in% fathers_pr),
    full_sib_pr = as.numeric(id %in% full_sibs_pr),
    child_pr = as.numeric(id %in% children_pr),
    spouse_pr = as.numeric(id %in% spouses_pr)
  )
]
if (
  mb[,
    !all(
    (proband + mother_pr + father_pr + full_sib_pr + child_pr + spouse_pr)
    %in% c(0, 1)
    )
  ]
) {
  stop("Problem with proband relative classification variables")
}
mb_fdrs_sp_pr <- mb[
  proband + mother_pr + father_pr + full_sib_pr + child_pr + spouse_pr == 1
]

# Create one copy of mb_fdrs_pr_sp for each population in pop_parms for test
# data set
test_data <- rbindlist(
  lapply(seq(1, nrow(pop_parms)), function(x) {
    copy(mb_fdrs_sp_pr)[,
      `:=`(
        famid = sprintf("%1d_%03d", x, famid),
        id = sprintf("%1d_%06d", x, id),
        fatherid = fifelse(
          mother_pr + father_pr + spouse_pr == 1,
          "", sprintf("%1d_%06d", x, fatherid)
        ),
        motherid = fifelse(
          mother_pr + father_pr + spouse_pr == 1,
          "", sprintf("%1d_%06d", x, motherid)
        ),
        pop = x
      )
    ]
  })
)

# Merge in population parameters
test_data <- test_data[pop_parms, on = .(pop)]
setkey(test_data, famid, id)

# Make kinship matrix and verify
test_ped <- with(
  test_data,
  kinship2::pedigree(id, fatherid, motherid, sex, famid = famid)
)
test_phi <- kinship2::kinship(test_ped)
if (
  any(test_data[, id] != rownames(test_phi)) ||
    any(test_data[, id] != colnames(test_phi))
) {
  stop("Test data set and kinship matrix not in same order")
}
if (any(Matrix::diag(test_phi) > 0.5)) {
  stop("There should be no consanguinity in the test data")
}
test_fam_ids <- test_data[famid %in% c("1_004", "1_005"), id]
test_fam_phi <- test_phi[test_fam_ids, test_fam_ids]

# Make population covariance matrix
h2_a <- Matrix::Diagonal(x = test_data[, h2_a])
sigma2 <- Matrix::Diagonal(x =  test_data[, sigma ^ 2])
sigma_mat <- as(
  as(
    2 * (sigma2 * h2_a) %*% test_phi + sigma2 - sigma2 * h2_a,
    "symmetricMatrix"
  ),
  "dsCMatrix"
)

# Obtain components used to calculate distribution conditional on probands'
# realized phenotypes
pr_ids <- test_data[proband == 1, id]
non_pr_ids <- test_data[proband == 0, id]
sigma_mat_pr_pr <- sigma_mat[pr_ids, pr_ids]
sigma_mat_non_pr_pr <- sigma_mat[non_pr_ids, pr_ids]
# Let L_pr_pr refer to the lower Cholesky factor of sigma_mat_pr_pr. Use
# of Matrix::crossprod with L_pr_pr^{-1} sigma_mat_non_pr_pr' ensures
# result of
# sigma_mat_non_pr_pr L_pr_pr^{-1}' L_pr_pr^{-1} sigma_mat_non_pr_pr'
# = sigma_mat_non_pr_pr (L_pr_pr L_pr_pr')^{-1} sigma_mat_non_pr_pr'
# = sigma_mat_non_pr_pr sigma_mat_pr_pr^{-1} sigma_mat_non_pr_pr'
# and omega_mat are recognized as sparse symmetric
L_pr_pr_inv_non_pr_pr_t <- Matrix::solve(
  Matrix::t(Matrix::chol(sigma_mat_pr_pr)),
  Matrix::t(sigma_mat_non_pr_pr)
)
omega_mat <- sigma_mat[non_pr_ids, non_pr_ids] -
  Matrix::crossprod(L_pr_pr_inv_non_pr_pr_t)
L_omega_mat <- Matrix::t(Matrix::chol(omega_mat))

sim_rep <- function(r) {

  setDTthreads(1)

  # --- Generate replicate data ---
  rep_data <- copy(test_data)[,
    `:=`(
      pop = factor(pop),
      sex = sapply(sex, FUN = switch, M = 1, F = 2, NA)
    )
  ][,
    c("n_rv", "x", "mu", "y") := {
      # Sample proband data from population distributions given above until
      # quantitative trait below threshold is obtained
      repeat{
        n_rv_pr <- rpois(1, lambda_rv[proband == 1])
        x_pr <- rnorm(1)
        mu_pr <- beta_rv[proband == 1] * n_rv_pr +
          beta_x[proband == 1] * x_pr
        y_pr <- rnorm(1, mu_pr, sd = sigma)
        if (y_pr <= tau[proband == 1]) break
      }
      # nolint start
      # A proband with n_rv_pr = c is assumed to be a REF/ALT heterozygote at
      # c unlinked loci each with a very rare ALT allele, in which case:
      # 1) Each ALT allele independently has a 50% chance of coming from the
      #    mother, so the number inherited from the mother is
      #    n_rv_m = Bin(c, 0.5). The rarity assumption implies that the mother
      #    would be REF/ALT and the father would be REF/REF at each of the
      #    transmitted loci.
      # 2) The remaining n_rv_f = c - n_rv_m ALT alleles must have come from
      #    the father. Under the same rarity assumption, the father would be
      #    REF/ALT and the mother REF/REF at each of these loci.
      # 3) Full sibs independently receive from mother and father
      #    Bin(n_rv_m, 0.5) + Bin(n_rv_f, 0.5) = Bin(c, 0.5)
      #    ALT alleles.
      # 4) Children of the proband receive Bin(c, 0.5) ALT alleles
      #    because the spouse is REF/REF under the rarity assumption.
      # nolint end
      n_rv <- rbinom(.N, n_rv_pr, 0.5)
      n_rv[proband == 1] <- n_rv_pr
      if (any(mother_pr == 1) && any(father_pr == 1)) {
        n_rv[father_pr == 1] <- n_rv_pr - n_rv[mother_pr == 1]
      }
      # Spouse data are set to NA to exclude from analysis without affecting
      # correspondence between data.table rows and kinship matrix rows/columns
      n_rv[spouse_pr == 1] <- NA
      # Covariate that is iid N(0, 1) across family members
      x <- rnorm(.N)
      x[proband == 1] <- x_pr
      x[spouse_pr == 1] <- NA
      mu <- beta_rv * n_rv + beta_x * x
      y <- rep(NA, .N)
      y[proband == 1] <- y_pr
      list(n_rv, x, mu, y)
    },
    by = famid
  ]
  # Note that rep_data retains famid, id key of test_data at this point, but
  # we still check ordering to be sure that correspondence with ordering of
  # sigma_mat partitions is exact
  if (
    any(pr_ids != rep_data[proband == 1, id]) ||
      any(non_pr_ids != rep_data[proband == 0, id])
  ) {
    stop("Ordering of rep_data and sigma_mat partitions differs")
  }
  rep_data[
    proband == 0,
    c("eta", "y") := {
      # NA mu for spouses ensures that their quantitative traits are
      # eventually NA
      eta <- mu + as.numeric(
        sigma_mat_non_pr_pr %*%
          Matrix::solve(sigma_mat_pr_pr, rep_data[proband == 1, y - mu])
      )
      y <- eta + as.numeric(L_omega_mat %*% rnorm(.N))
      list(eta, y)
    }
  ]
  if (rep_data[proband == 1, any(y >= tau)]) {
    stop("Proband trait exceeded threshold")
  }
  if (rep_data[spouse_pr + proband == 0, any(is.na(y))]) {
    stop("FDR missing trait value")
  }
  if (rep_data[spouse_pr == 1, any(!is.na(y))]) {
    stop("Trait data provided for proband's spouse(s)")
  }

  # --- Fit true model ---
  fd <- FamData$new(
    rep_data, family_id = "famid", indiv_id = "id", proband = "proband",
    sex = "sex", maternal_id = "motherid", paternal_id = "fatherid",
    mzgrp = "mzgrp", dzgrp = "dzgrp"
  )
  lmm <- fd$lmm(y ~ 0 + pop + n_rv:pop + x:pop | pop)

  # --- Get test quantities ---

  # nolint start
  # The p-value for the Wald chi-square contrast comparing all parameters to
  # their true values should be approximately U(0, 1) in the populations where
  # no parameters are on the boundary of the parmeter space (population 2) if
  # everything is working correctly. Note that this result holds with the
  # orthogonal parameterization used but not necessarily with
  # parameterizations where parameters are shared between populations. See
  # vignette("linear_mixed_models")
  # nolint end
  id_mat <- diag(length(lmm$get_theta_hat()))
  rownames(id_mat) <- colnames(id_mat) <- names(lmm$get_theta_hat())
  nb_parm_cx <- lmm$contrast(
    id_mat[c("pop2", "pop2:n_rv", "pop2:x", "h2_a.pop2", "sigma.pop2"), ],
    pop_parms[pop == 2, c(0, beta_rv, beta_x, h2_a, sigma)]
  )

  # It can be shown that the LRT p-value for pop1 with h2_a = 0 has a CDF of
  # the form F(z) = z for z in [0, 0.5), 0.5 for z in [0.5, 1), and 1 for z =
  # 1. Thus, the EDF should behave similarly. We can also check that the
  # p-values for h2_a.pop2 are stochastically less than U(0, 1)
  h2_a_lrts <- lmm$get_h2_a_lrts(print = FALSE)

  # The EDFs for various types of residuals and goodness-of-fit diagnostics
  # should be asymptotically unbiased estimates of the theoretical marginal
  # CDFs at any t, so we can verify that they do not appear to deviate
  # systematically (i.e., on average) from the theoretical marginal CDFs at
  # any t. Note that EDF-based tests and confidence intervals are *not* valid
  # for this purpose because of variance shrinkage arising from substitution
  # of estimated parameters. See vignette("linear_mixed_models")
  mod_res <- lmm$get_model_res()

  # Calculate tests described above and package results for return
  res <- list(
    rep = r,
    pvals = data.table(
      p_nb_parms = nb_parm_cx$get_p_X2(),
      p_h2_a_0_pop1 = h2_a_lrts[h_0 == "h2_a.pop1 = 0", p_X2],
      p_h2_a_0_pop2 = h2_a_lrts[h_0 == "h2_a.pop2 = 0", p_X2]
    ),
    resids = mod_res[pr == 0, .(r_c_hat, r_star_hat)],
    p_c_star_hat = mod_res[
      pr == 0,
      .(p_c_star_hat = unique(p_c_star_hat)),
      by = fmid
    ]
  )
  if (nrow(res[["p_c_star_hat"]]) != mod_res[, uniqueN(fmid)]) {
    stop("Multiple p_c_star_hat values found in a single fmid")
  }

  res
}

# Run simulations in parallel and return data.table of results
sim_res <- future_lapply(
  seq(1, n_tests * n_reps),
  sim_rep,
  future.seed = TRUE,
  future.packages = "FamModel"
)

p_val_edf_cdf_check <- function(colname, data) {
  edf <- ecdf(data[[colname]])
  # Reference CDF should be U(0, 1) except for p_h2_a_0_pop1, which has the
  # special CDF described above
  cdf <- if (colname == "p_h2_a_0_pop1") {
    function(x) pmin(x, 0.50) + 0.5 * as.numeric(x == 1)
  } else {
    function(x) x
  }
  t <- c(0, unique(data[[colname]]), 1)
  edf_col <- edf(t)
  # Minimal abscissae necessary to draw both CDFs
  t_cdf <- c(0, 0.5, 1 - .Machine[["double.eps"]], 1)
  cdf_col <- cdf(t_cdf)
  # 95% DKW confidence band from Wasserman (2006), All of Nonparametric
  # Statistics, p. 15
  eps_n <- sqrt(log(2 / 0.05) / (2 * length(data[[colname]])))
  u_edf_col <- pmin(edf_col + eps_n, 1)
  l_edf_col <- pmax(edf_col - eps_n, 0)
  plot_data <- data.table(
    x = c(rep(t, 3), t_cdf),
    y = c(edf_col, u_edf_col, l_edf_col, cdf_col),
    type = rep(
      c("EDF", "EDF 95% UCL", "EDF 95% LCL", "CDF"),
      times = c(rep(length(t), 3), length(t_cdf))
    )
  )

  plot_data
}

resids_dx_edfs_cdf_check <- function(colname, data) {
  # Reference CDF should be N(0, 1) except for p_c_star_hat, which should be
  # U(0, 1). EDF and its difference from the CDF are evaluated at points
  # ranging from -4 to 4 by 0.05, inclusive, for N(0,1) or from 0 to 1 by
  # 0.01, incusive, for U(0, 1) in each replicate. Sampling at a smaller
  # number of points is necessary to keep vector snapshot file sizes down. The
  # mean of the EDF across replicates is taken at each of these points
  if (colname == "p_c_star_hat") {
    cdf <- punif
    t <- seq(0, 1, by = 0.01)
  } else {
    cdf <- pnorm
    t <- seq(-4, 4, by = 0.05)
  }
  eval(bquote(
    plot_data <- data[,
      list(t, edf = ecdf(.(as.name(colname)))(t), cdf = cdf(t)),
      by = rep
    ]
  ))
  mean_data <- plot_data[, .(mean_edf = mean(edf)), by = t][, cdf := cdf(t)]
  if (!identical(mean_data[["t"]], t)) {
    stop("t values were not identical for all replicates")
  }

  list(plot_data = plot_data, mean_data = mean_data)
}

pvals_edf_cdf_data <- list()
resids_dx_edfs_cdf_data <- list()
for (test in seq(1, n_tests)) {
  test_res <- sim_res[seq((test - 1) * n_reps + 1, test * n_reps)]

  pvals <- rbindlist(
    lapply(test_res, function(x) data.table(rep = x[["rep"]], x[["pvals"]]))
  )
  pvals_edf_cdf_data[[test]] <- sapply(
    grep("^rep$", names(pvals), invert = TRUE, value = TRUE),
    p_val_edf_cdf_check,
    data = pvals,
    simplify = FALSE
  )
  resids <- rbindlist(
    lapply(test_res, function(x) data.table(rep = x[["rep"]], x[["resids"]]))
  )
  resids_dx_edfs_cdf_data[[test]] <- sapply(
    grep("^rep$", names(resids), invert = TRUE, value = TRUE),
    resids_dx_edfs_cdf_check,
    data = resids,
    simplify = FALSE
  )
  p_c_star_hat <- rbindlist(
    lapply(
      test_res,
      function(x) data.table(rep = x[["rep"]], x[["p_c_star_hat"]])
    )
  )
  resids_dx_edfs_cdf_data[[test]][["p_c_star_hat"]] <-
    resids_dx_edfs_cdf_check("p_c_star_hat", p_c_star_hat)
  rm(test_res, pvals, resids, p_c_star_hat)
}

test_that("Simulation results match those that produced approved graphs", {
  # Snapshot test kinship matrix verified against plotted pedigrees for
  # families 4 and 5 in population 1
  expect_snapshot_value(
    test_fam_phi, style = "serialize", cran = TRUE,
    tolerance = testthat_tolerance()
  )

  # Always check serialized objects underlying snapshot graphs approved in
  # development environment with appropriate numerical tolerance. Need to loop
  # through by test to allow comparison of only first two tests on CI
  for (test in seq(1, n_tests)) {
    expect_snapshot_value(
      pvals_edf_cdf_data[[test]], style = "serialize", cran = TRUE,
      tolerance = testthat_tolerance()
    )
    expect_snapshot_value(
      resids_dx_edfs_cdf_data[[test]], style = "serialize", cran = TRUE,
      tolerance = testthat_tolerance()
    )
  }
})

test_that("Simulation results graphs match expectations", {
  # Perform output/graphical comparison only in development environment. Only
  # snapshots of serialized objects underlying these accepted snapshot outputs
  # are compared on CI and CRAN to avoid failures due to minor numerical
  # differences across platforms
  skip_on_ci()
  skip_on_cran()

  # Snapshot test kinship matrix verified against plotted pedigrees for
  # families 4 and 5 in population 1
  expect_snapshot_output(Matrix::print(test_fam_phi, col.names = TRUE))

  for (test in seq(1, n_tests)) {

    # Produce snapshot graph comparing p-value EDF with 95% DKW confidence bands
    # to appropriate CDF
    for (colname in names(pvals_edf_cdf_data[[test]])) {
      plot_data <- pvals_edf_cdf_data[[test]][[colname]]
      edf_cdf_plot <- ggplot() +
        geom_step(
          aes(x = x, y = y, lty = type, color = type),
          data = plot_data[type != "CDF"]
        ) +
        geom_line(
          aes(x = x, y = y, lty = type, color = type),
          data = plot_data[type == "CDF"]
        ) +
        scale_color_manual(
          name = NULL,
          values = c(
            CDF = "black", EDF = "red", `EDF 95% LCL` = "red",
            `EDF 95% UCL` = "red"
          )
        ) +
        scale_linetype_manual(
          name = NULL,
          values = c(CDF = 1, EDF = 1, `EDF 95% LCL` = 2, `EDF 95% UCL` = 2)
        ) + theme_bw() + ylab("F(t)") + xlab("t")
      expect_doppelganger(
        paste0("Test ", test, ": ", colname, " EDF vs. CDF"),
        edf_cdf_plot,
        path = ""
      )
      rm(plot_data, edf_cdf_plot)
    }

    # Produce snapshot graphs comparing EDFs of residuals or diagnostic in each
    # replicate sample and their mean across replicate to appropriate CDF
    for (colname in names(resids_dx_edfs_cdf_data[[test]])) {
      plot_data <- resids_dx_edfs_cdf_data[[test]][[colname]][["plot_data"]]
      mean_data <- resids_dx_edfs_cdf_data[[test]][[colname]][["mean_data"]]
      edfs_cdf_plot <- ggplot() +
        geom_line(
          aes(
            x = t, y = edf - cdf, group = rep, color = "Individual\nReplicate"
          ),
          data = plot_data,
          alpha = 0.10
        ) +
        geom_abline(slope = 0, intercept = 0) +
        geom_line(
          aes(x = t, y = mean_edf - cdf, color = "Mean"),
          data = mean_data,
          alpha = 0.8,
          size = 2
        ) + theme_bw() +
        scale_color_manual(
          name = NULL,
          values = c(`Individual\nReplicate` = "red", Mean = "purple")
        ) +
        ylab(paste(
          "EDF - Standard",
          if (colname == "p_c_star_hat") "Uniform" else "Normal",
          "CDF"
        ))
      expect_doppelganger(
        paste0("Test ", test, ": ", colname, " EDFs vs. CDF"),
        edfs_cdf_plot,
        path = ""
      )
      rm(plot_data, mean_data, edfs_cdf_plot)
    }
  }

})
