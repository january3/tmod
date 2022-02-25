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

test_that("makeTmodGS works", {

  gs <- data.frame(ID=letters[1:3], Title=LETTERS[1:3])
  gs2gv <- list(a=c("g1", "g2"), b=c("g3", "g4"), c=c("g1", "g2", "g4"))
  mset <- makeTmodGS(gs2gene=gs2gv, gs=gs)
  expect_s3_class(mset, "tmodGS")
  expect_setequal(names(mset), c( "gs", "gs2gv", "gv", "info", "weights" ))
  expect_equal(length(mset$gv), 4)

  mods <- tmod_ids(mset)
  expect_setequal(mods, letters[1:3])

  expect_setequal(mset$gv[ mset$gs2gv[[ 1 ]] ], c("g1", "g2"))
  expect_equal(length(mset$gs2gv), nrow(mset$gs))

})
