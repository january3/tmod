library(tmod)
context("Testing importing and exporting functions")

data(tmod)
mset <- tmod[c("LI.M37.0", "LI.M75", "LI.M3")]

test_that("MsigDB import", {

  fn <- "msigdb_test.xml"
  mset <- tmodImportMSigDB(fn)
  expect_is(mset, "tmodGS")
  expect_output(print(mset), '5 gene sets, 609 genes')
  expect_setequal(names(mset), c( "gs", "gs2gv", "gv"))
  expect_equal(nrow(mset$gs), length(mset$gs2gv))
})

test_that("working makeTmodFromDataFrame", {

  df <- data.frame(gene_id=LETTERS[1:10], geneset_id=rep(letters[1:2], each=5))
  mset <- makeTmodFromDataFrame(df)
  expect_is(mset, "tmodGS")
  expect_output(print(mset), '2 gene sets, 10 genes')
  expect_equal(nrow(mset$gs), length(mset$gs2gv))
  expect_setequal(names(mset), c( "gs", "gs2gv", "gv"))
})

test_that("tmod2DataFrame works", {
  df1 <- tmod2DataFrame(mset)
  expect_is(df1, "tbl_df")
  expect_equal(nrow(df1), nrow(mset$gs))

  df2 <- tmod2DataFrame(mset, rows="features")
  expect_equal(nrow(df2), length(mset$gv))

})

