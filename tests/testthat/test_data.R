library(tmod)
context("Testing data structures")

test_that("tmod is sane", {
  data(tmod)
  expect_is(tmod, "tmod")
  expect_output(print(tmod), 'An object of class "tmod"' )
  expect_output(print(tmod), '606 modules, 12712 genes')
  expect_equal(nrow(tmod$MODULES), 606)
  expect_identical(names(tmod), c( "MODULES", "GENES", "MODULES2GENES", "GENES2MODULES" ))
})
