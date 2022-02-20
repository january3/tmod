library(tmod)
context("Statistical tests")

tt <- read.csv( "sample_data.csv", stringsAsFactors=FALSE)[,1]

testRes <- function(res, res_e, test_auc=TRUE) {
  expect_equal(nrow(res), nrow(res_e))
  expect_equal(ncol(res), ncol(res_e))

  expect_identical(res$ID, res_e$ID)
  expect_identical(res$Title, res_e$Title)
  if(test_auc) {
    expect_equal(res$AUC, res_e$AUC)
  }
}

test_that("tmodCERNO", {
  res   <- tmodCERNOtest(tt)
  res_e <- read.csv("res_CERNO.csv", stringsAsFactors=FALSE, row.names=1)
  testRes(res, res_e)
})

test_that("tmodUtest", {
  res   <- tmodUtest(tt)
  res_e <- read.csv("res_U.csv", stringsAsFactors=FALSE, row.names=1)
  testRes(res, res_e)
})

test_that("tmodHGtest", {
  fg    <- head(tt, 100)
  res   <- tmodHGtest(fg, tt)
  res_e <- read.csv("res_HG.csv", stringsAsFactors=FALSE, row.names=1)
  testRes(res, res_e)
})

test_that("tmodLEA", {
  res_e <- read.csv("res_LEA.csv", stringsAsFactors=FALSE, row.names=1)
  res  <- tmodLEASummary(tmodLEA(tt, c("LI.M37.0", "LI.M75", "LI.M3")))
  testRes(res, res_e, test_auc=FALSE)
  expect_equal(res$N, res_e$N)
  expect_equal(res$fraction, res_e$fraction)
  expect_equal(res$LEA, res_e$LEA)
})

