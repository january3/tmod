library(tmod)
data(Egambia)
data(tmod)
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
  res   <- tmodCERNOtest(tt, mset="LI")
  res_e <- read.csv("res_CERNO.csv", stringsAsFactors=FALSE, row.names=1)
  testRes(res, res_e)
})

test_that("tmodUtest", {
  res   <- tmodUtest(tt, mset="LI")
  res_e <- read.csv("res_U.csv", stringsAsFactors=FALSE, row.names=1)
  testRes(res, res_e)
})

test_that("tmodHGtest", {
  fg    <- head(tt, 100)
  res   <- tmodHGtest(fg, tt, mset="LI")
  res_e <- read.csv("res_HG.csv", stringsAsFactors=FALSE, row.names=1)
  testRes(res, res_e)
})

test_that("tmodLEA", {
  res_e <- read.csv("res_LEA.csv", stringsAsFactors=FALSE, row.names=1)
  res  <- tmodLEASummary(tmodLEA(tt, c("LI.M37.0", "LI.M75", "LI.M3"), mset="LI"))
  testRes(res, res_e, test_auc=FALSE)
  expect_equal(res$N, res_e$N)
  expect_equal(res$fraction, res_e$fraction)
  expect_equal(res$LEA, res_e$LEA)
})

test_that("tmodPLAGEtest", {
  res_e <- read.csv("res_PLAGE.csv", stringsAsFactors=FALSE, row.names=1)
  x <- Egambia[ , -1:-3 ]
  l <- Egambia$GENE_SYMBOL
  group <- gsub("\\..*", "", colnames(x))
  res <- tmodPLAGEtest(l, x, group, mset="LI")
  testRes(res, res_e, test_auc=FALSE)
})

test_that("tmodZtest", {
  res_e <- read.csv("res_Z.csv", stringsAsFactors=FALSE, row.names=1)
  res <- tmodZtest(tt, mset="LI")
  testRes(res, res_e, test_auc=FALSE)
})

test_that("tmodGeneSetTest", {
  res_e <- read.csv("res_GeneSetTest.csv", stringsAsFactors=FALSE, row.names=1)
  x <- Egambia[ , -1:-3 ]
  l <- Egambia$GENE_SYMBOL
  group <- gsub("\\..*", "", colnames(x))

  dd <- rowMeans(x[ , group == "NID" ]) - rowMeans(x[ , group == "TB" ])
  dd <- dd / apply(x, 1, sd)

  res <- tmodGeneSetTest(l, dd, Nsim=100, mset="LI")

  expect_equal(ncol(res), ncol(res_e))

  #rr <- merge(res, res_e, by=c("ID", "Title"))
  #expect_gte(cor(rr$D.x, rr$D.y), .8)
  #expect_gte(cor(rr$P.Value.x, rr$P.Value.y, method="s"), .8)

})


test_that("tmodAUC", {

  res_e <- read.csv("res_AUC.csv", stringsAsFactors=FALSE, row.names=1)
  res <- data.frame(tmodAUC(tt, 1:length(tt), mset="LI"))

  expect_equal(nrow(res), sum(tmod$gs$SourceID == "LI"))
  expect_equal(ncol(res), 1)
  expect_equal(res[,1], res_e[,1])


})

test_that("tmodDecideTests", {
  res_tt <- read.csv("sample_results.csv")
  res_e <- read.csv("res_tmodDecideTests.csv", row.names=1)
  res   <- tmodDecideTests(res_tt$GENE_SYMBOL, lfc=res_tt$logFC, pval=res_tt$adj.P.Val)
  expect_equal(length(res), 1)
  expect_is(res, "list")
  res   <- as.data.frame(res[[1]])
  expect_equal(res_e, res)



})
