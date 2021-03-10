library(testthat, quietly = TRUE)
library(FamModel, quietly = TRUE)

withr::with_envvar(new = c("R_CMD_CHECK" = "true"), test_check("FamModel"))
