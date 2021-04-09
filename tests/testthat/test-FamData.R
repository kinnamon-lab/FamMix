# Tests for FamData class functionality

test_that("Constructor works properly", {

  # Set up testbed with randomly ordered input data set
  lmna_nonseg_perm <- lmna_nonseg[sample(.N)]
  phi_ids_perm <- lmna_nonseg_perm[, paste(family_ID, individual_ID, sep = "/")]
  fd_sig1 <- FamData$new(
    lmna_nonseg_perm, family_id = "family_ID", indiv_id = "individual_ID",
    proband = "proband", sex = "sex", maternal_id = "maternal_ID",
    paternal_id = "paternal_ID", mzgrp = "mzpair", dzgrp = "dzpair"
  )
  fd_sig1_phi_perm <- fd_sig1$get_phi()[phi_ids_perm, phi_ids_perm]
  # This constructor signature requires that rows/columns of phi matrix are in
  # the same order as subjects in the data set
  fd_sig2 <- FamData$new(
    lmna_nonseg_perm, family_id = "family_ID", indiv_id = "individual_ID",
    proband = "proband", phi = fd_sig1_phi_perm
  )

  # Should throw error with unrecognized constructor signature
  expect_error(
    FamData$new(lmna_nonseg, family_id = "family_ID"),
    regexp = "Unrecognized constructor signature"
  )

  # Both recognized constructor signatures should work
  expect_s3_class(fd_sig1, c("FamData", "R6"), exact = TRUE)
  expect_s3_class(fd_sig2, c("FamData", "R6"), exact = TRUE)

  # Data set name should be captured
  expect_identical(fd_sig1$get_data_name(), "lmna_nonseg_perm")
  expect_identical(fd_sig2$get_data_name(), "lmna_nonseg_perm")

  # Data member should be keyed by fmid, then id for both signatures
  expect_identical(key(fd_sig1$get_data()), c("fmid", "id"))
  expect_identical(key(fd_sig2$get_data()), c("fmid", "id"))

  # Data members should be identical to original data set ordered by family_ID,
  # then individual_ID after using the object's var_map member to restore
  # original column names and converting all integer columns in original data
  # set to numeric
  lmna_nonseg_orig <- copy(lmna_nonseg)
  is_int_col <- lmna_nonseg_orig[, sapply(.SD, is.integer)]
  int_cols <- names(is_int_col)[is_int_col]
  lmna_nonseg_orig[,
    (int_cols) := lapply(.SD, as.numeric),
    .SDcols = int_cols
  ]
  setkey(lmna_nonseg_orig, family_ID, individual_ID)
  fd_sig1_data <- fd_sig1$get_data()
  fd_sig1_vm <- fd_sig1$get_var_map()
  setnames(fd_sig1_data, names(fd_sig1_vm), fd_sig1_vm)
  expect_identical(fd_sig1_data, lmna_nonseg_orig)
  fd_sig2_data <- fd_sig2$get_data()
  fd_sig2_vm <- fd_sig2$get_var_map()
  setnames(fd_sig2_data, names(fd_sig2_vm), fd_sig2_vm)
  expect_identical(fd_sig2_data, lmna_nonseg_orig)

  # Phi members should be equal for both signatures
  expect_identical(fd_sig1$get_phi(), fd_sig2$get_phi())

  # Multiple rows for a fmid + id combo should trigger an error
  expect_error(
    FamData$new(
      rbind(lmna_nonseg, lmna_nonseg[1]), family_id = "family_ID",
      indiv_id = "individual_ID", proband = "proband", sex = "sex",
      maternal_id = "maternal_ID", paternal_id = "paternal_ID",
      mzgrp = "mzpair", dzgrp = "dzpair"
    ),
    regexp = "Multiple rows found for the same individual"
  )
  expect_error(
    FamData$new(
      rbind(lmna_nonseg, lmna_nonseg[1]), family_id = "family_ID",
      indiv_id = "individual_ID", proband = "proband",
      phi = fd_sig1_phi_perm
    ),
    regexp = "Multiple rows found for the same individual"
  )

  # A disjoint pedigree should trigger error with pedfile signature
  disjoint_peds <- copy(lmna_nonseg)[family_ID == "L", family_ID := "N"]
  expect_error(capture_output(
    FamData$new(
      disjoint_peds, family_id = "family_ID", indiv_id = "individual_ID",
      proband = "proband", sex = "sex", maternal_id = "maternal_ID",
      paternal_id = "paternal_ID", mzgrp = "mzpair", dzgrp = "dzpair"
    ),
    regexp = "Incorrectly specified pedigrees"
  ))
})

test_that("Twins are handled properly", {

  # Family 1 has 2 DZ and 2 MZ twins; family 2 has MZ triplets; family 3 has
  # no twins
  twin_test <- data.table(
    famid = rep(c(1, 2, 3), times = c(6, 5, 4)),
    id = as.numeric(c(seq(1, 6), seq(1, 5), seq(1, 4))),
    pr = c(rep(0, 5), 1, rep(0, 4), 1, rep(0, 3), 1),
    sex = c(1, 2, 1, 2, rep(1, 2), 1, 2, rep(2, 3), 1, 2, rep(2, 2)),
    mid = c(0, 0, rep(2, 4), 0, 0, rep(2, 3), 0, 0, rep(2, 2)),
    fid = c(0, 0, rep(1, 4), 0, 0, rep(1, 3), 0, 0, rep(1, 2)),
    mzgrp = c(rep(NA, 4), rep("A", 2), NA, NA, rep("A", 3), rep(NA, 4)),
    dzgrp = c(NA, NA, rep("A", 2), rep(NA, 11))
  )

  twin_test_fd <- FamData$new(
    twin_test, family_id = "famid", indiv_id = "id", proband = "pr",
    sex = "sex", maternal_id = "mid", paternal_id = "fid", mzgrp = "mzgrp",
    dzgrp = "dzgrp"
  )

  expect_snapshot_output(
    Matrix::print(twin_test_fd$get_phi(), col.names = TRUE)
  )

  # Should also work correctly with numeric or factor twin group identifiers
  twin_test_fd2 <- FamData$new(
    copy(twin_test)[,
      `:=`(
        mzgrp = as.numeric(mzgrp == "A"),
        dzgrp = factor(dzgrp)
      )
    ],
    family_id = "famid", indiv_id = "id", proband = "pr", sex = "sex",
    maternal_id = "mid", paternal_id = "fid", mzgrp = "mzgrp", dzgrp = "dzgrp"
  )

  expect_equal(twin_test_fd$get_phi(), twin_test_fd2$get_phi())

})

test_that("Consanguinity is detected", {

  # Family 1 is consanguinous; family 2 is not
  consang_test <- data.table(
    famid = rep(c(1, 2), each = 6),
    id = as.numeric(rep(seq(1, 6), times = 2)),
    pr = rep(c(rep(0, 5), 1), times = 2),
    sex = rep(c(1, 2), times = 6),
    mid = c(0, 0, 2, 2, 4, 4, 0, 0, 2, 2, 2, 2),
    fid = c(0, 0, 1, 1, 3, 3, 0, 0, 1, 1, 1, 1),
    mzgrp = NA,
    dzgrp = NA
  )

  consang_test_fd <- FamData$new(
    consang_test, family_id = "famid", indiv_id = "id", proband = "pr",
    sex = "sex", maternal_id = "mid", paternal_id = "fid", mzgrp = "mzgrp",
    dzgrp = "dzgrp"
  )

  expect_equal(consang_test_fd$get_consang(), 1)

  # Should get same result with apporpirate type with character family_id column
  consang_test_fd2 <- FamData$new(
    consang_test[, famid := as.character(famid)], family_id = "famid",
    indiv_id = "id", proband = "pr", sex = "sex", maternal_id = "mid",
    paternal_id = "fid", mzgrp = "mzgrp", dzgrp = "dzgrp"
  )

  expect_equal(consang_test_fd2$get_consang(), "1")

})
