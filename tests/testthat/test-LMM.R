# Tests for linear mixed model functionality

test_that("Can reproduce univariate model results from Mendel 16.0", {
  # Note that mendel uses N rather than N-1 for SD divisor in standardization
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

  # Check against snapshot output and data.tables of r_star_hat ^ 2 and
  # c_star_hat previously checked against Mendel 16.0 results
  check_model <- function(model) {
    expect_snapshot_output(model$print(), cran = TRUE)
    res <- model$get_model_res()[
      pr == 0, .(fmid, id, r_star_hat, c_star_hat, c_star_hat_df, p_c_star_hat)
    ]
    expect_snapshot_output(
      print(res[, .(fmid, id, r_star_hat_2 = r_star_hat ^ 2)]),
      cran = TRUE
    )
    expect_snapshot_output(
      print(unique(res[, .(fmid, c_star_hat, c_star_hat_df, p_c_star_hat)])),
      cran = TRUE
    )
  }

  check_model(lmna_lvef_model)
  check_model(lmna_lvedd_z_model)
})
