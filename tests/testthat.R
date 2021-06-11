Sys.setenv("R_TESTS" = "") 
library(testthat)
library(ENMeval)

testthat::test_check("ENMeval")
