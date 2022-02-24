library(tmod)
context("Testing data structures")

test_that("tmod is sane", {
  data(tmod)
  expect_is(tmod, "tmodGS")
  expect_output(print(tmod), 'An object of class "tmodGS"' )
  expect_output(print(tmod), '606 gene sets, 12712 genes')
  expect_equal(nrow(tmod$gs), 606)
  expect_identical(names(tmod), c( "gs", "gs2gv", "gv", "info", "weights" ))
  expect_equal(nrow(tmod$gs), length(tmod$gs2gv))

  mset <- tmod[c("LI.M37.0", "LI.M75", "LI.M3")]
  expect_s3_class(mset, "tmodGS")
  expect_output(print(mset), '3 gene sets, 393 genes')
  expect_identical(names(mset), c( "gs", "gs2gv", "gv", "info", "weights" ))
  expect_equal(nrow(mset$gs), length(mset$gs2gv))

})
