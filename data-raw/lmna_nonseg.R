# Fetch and prepare lmna_nonseg example data
library(data.table)
lmna_nonseg <-
  fread(
    paste0(
      "https://raw.githubusercontent.com/kinnamon-lab/lmna_nonseg/",
      "469141d54f0276722b1c81a93018ef46a9f8e304/data/lmna_nonseg_data.csv"
    ),
    na.strings = "n/a",
    data.table = TRUE
  )
# Fix reversed paternal and maternal ID
lmna_nonseg[family_ID == "O" &
  individual_ID == 10, `:=`(maternal_ID = 7, paternal_ID = 8)]
lmna_nonseg[, `:=`(mzpair = NA, dzpair = NA)]
usethis::use_data(lmna_nonseg, overwrite = TRUE)
