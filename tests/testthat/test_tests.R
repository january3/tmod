library(tmod)
context("Statistical tests")

tt <- read.csv( "sample_data.csv", stringsAsFactors=FALSE)[,1]

testRes <- function(res, res_e) {
  expect_equal(nrow(res), nrow(res_e))
  expect_equal(ncol(res), ncol(res_e))

  expect_identical(res$Title, res_e$Title)
  expect_equal(res$AUC, res_e$AUC)
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



