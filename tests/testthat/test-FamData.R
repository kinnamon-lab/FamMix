context("FamData class")

test_that("Constructor works properly", {

  # Set up testbed with randomly ordered input data set
  lmna_nonseg_orig <- copy(lmna_nonseg)
  setkey(lmna_nonseg_orig, family_ID, individual_ID)
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
    "Unrecognized constructor signature"
  )

  # Both recognized constructor signatures should work
  expect_is(fd_sig1, c("FamData", "R6"))
  expect_is(fd_sig2, c("FamData", "R6"))

  # Data set name should be captured
  expect_identical(fd_sig1$get_data_name(), "lmna_nonseg_perm")
  expect_identical(fd_sig2$get_data_name(), "lmna_nonseg_perm")

  # Data member should be keyed by fmid, then id for both signatures
  expect_identical(key(fd_sig1$get_data()), c("fmid", "id"))
  expect_identical(key(fd_sig2$get_data()), c("fmid", "id"))

  # Data members should be identical to original data set ordered by family_ID,
  # then individual_ID after using the object's var_map member to restore
  # original column names
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
    "Multiple rows found for the same individual"
  )

  # A disjoint pedigree should trigger error with pedfile signature
  disjoint_peds <- copy(lmna_nonseg)[family_ID == "L", family_ID := "N"]
  expect_error(capture_output(
    FamData$new(
      disjoint_peds, family_id = "family_ID", indiv_id = "individual_ID",
      proband = "proband", sex = "sex", maternal_id = "maternal_ID",
      paternal_id = "paternal_ID", mzgrp = "mzpair", dzgrp = "dzpair"
    ),
    "Incorrectly specified pedigrees"
  ))
})

# Test kinship matrices in a couple of known cases
#test_that("Kinship matrices are correct", {
#  expect_snapshot_value(fd_sig1$get_phi())
#})
