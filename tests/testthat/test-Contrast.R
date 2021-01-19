# Tests of Contrast class functionality

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

test_that("Constructor catches incorrect inputs", {

  # Error on model_fit argument with wrong object type
  expect_error(
    Contrast$new(lmna_fd, diag(7)),
    regexp = "Argument model_fit must be an object inheriting from FamModelFit"
  )

  # Error on bad L_mat
  expect_error(
    Contrast$new(lmna_lvef_model, Matrix::Diagonal(7)),
    regexp = "Argument L_mat must be a numeric vector or matrix"
  )
  expect_error(
    Contrast$new(lmna_lvef_model, diag(6)),
    regexp = paste0(
      "Argument L_mat must have exactly as many elements \\(vector\\) or ",
      "columns \\(matrix\\) as model parameters"
    )
  )
  expect_error(
    Contrast$new(lmna_lvef_model, rep(1, 6)),
    regexp = paste0(
      "Argument L_mat must have exactly as many elements \\(vector\\) or ",
      "columns \\(matrix\\) as model parameters"
    )
  )
  expect_error(
    Contrast$new(lmna_lvef_model, rbind(diag(7), diag(7))),
    regexp = "Argument L_mat must be of full row rank"
  )

  # Error on bad m
  expect_error(
    Contrast$new(lmna_lvef_model, c(1, rep(0, 6)), rep(0, 2)),
    regexp = paste0(
      "Argument m must be a numeric vector with as many elements as there are ",
      "rows of L_mat"
    )
  )
  expect_error(
    Contrast$new(lmna_lvef_model, diag(7), rep(0, 6)),
    regexp = paste0(
      "Argument m must be a numeric vector with as many elements as there are ",
      "rows of L_mat"
    )
  )
})

test_that("Identity contrasts reproduce inputs", {

  # Identity contrast with m = 0
  id_contrast <- Contrast$new(lmna_lvef_model, diag(7))
  expect_equal(
    id_contrast$get_L_theta_hat_m(), lmna_lvef_model$get_theta_hat(),
    ignore_attr = "names"
  )
  expect_equal(
    id_contrast$get_V_L_theta_hat(), lmna_lvef_model$get_V_theta_hat(),
    ignore_attr = "dimnames"
  )
  expect_equal(id_contrast$get_m(), rep(0, 7))
  expect_equal(id_contrast$get_df_X2(), 7)

  # Identity contrast with m = theta_hat
  id_zero_contrast <- Contrast$new(
    lmna_lvef_model, diag(7), lmna_lvef_model$get_theta_hat()
  )
  expect_equal(
    id_zero_contrast$get_L_theta_hat_m(), rep(0, 7),
    ignore_attr = "names"
  )
  expect_equal(
    id_contrast$get_V_L_theta_hat(), id_zero_contrast$get_V_L_theta_hat()
  )
  expect_equal(id_zero_contrast$get_X2(), 0)
  expect_equal(id_zero_contrast$get_p_X2(), 1)

  # Single parameter identity contrast
  id_1p_contrast <- Contrast$new(lmna_lvef_model, c(0, 0, 1, rep(0, 4)))
  expect_equal(
    id_1p_contrast$get_L_theta_hat_m(), lmna_lvef_model$get_theta_hat()[3],
    ignore_attr = "names"
  )
  expect_equal(
    id_1p_contrast$get_V_L_theta_hat(),
    lmna_lvef_model$get_V_theta_hat()[3, 3, drop = FALSE],
    ignore_attr = "dimnames"
  )
  expect_equal(
    id_1p_contrast$get_X2(),
    (lmna_lvef_model$get_theta_hat()[3] ^ 2) /
      lmna_lvef_model$get_V_theta_hat()[3, 3],
    ignore_attr = "names"
  )
  expect_equal(id_1p_contrast$get_df_X2(), 1)
  expect_equal(
    id_1p_contrast$get_p_X2(),
    pchisq(id_1p_contrast$get_X2(), 1, lower.tail = FALSE)
  )
})
