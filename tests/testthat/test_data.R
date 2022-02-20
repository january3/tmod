library(tmod)
context("Testing data structures")

test_that("tmod is sane", {
  data(tmod)
  expect_is(tmod, "tmodGS")
  expect_output(print(tmod), 'An object of class "tmodGS"' )
  expect_output(print(tmod), '606 gene sets, 12712 genes')
  expect_equal(nrow(tmod$gs), 606)
  expect_identical(names(tmod), c( "gs", "gs2gv", "gv", "info", "weights" ))
})
