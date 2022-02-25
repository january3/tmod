library(tmod)
context("Testing eigengenes and gene set overlaps")

data(Egambia)
data(tmod)

x <- Egambia[ , -1:-3 ]
g <- Egambia$GENE_SYMBOL 
mset <- tmod[c("LI.M37.0", "LI.M75", "LI.M3")]

test_that("Checking eigengene calculations", {

  eig <- eigengene(x, g, mset=mset)

  expect_equal(nrow(eig), 3)
  expect_equal(ncol(eig), ncol(x))
  expect_equal(eig[2,3], -2.4923429)

  mcors <- modcors(x, g, mset)
  expect_equal(mcors[2, 1], 0.5338431)

})


test_that("Checking modjacard", {
  mjac <- modjaccard(mset)
  expect_equal(mjac[1, 1], 1)
  expect_equal(mjac[2, 1], 0.008196721)

  g <- tmod$gv[ tmod$gs2gv[[ match("LI.M37.0", tmod$gs$ID) ]] ]
  mjac <- modjaccard(mset, g)

  expect_equal(mjac[2, 1], 0.008645533)
})

test_that("Checking modOverlaps", {
  mo1 <- modOverlaps(mset=mset)
  mo2 <- modOverlaps(mset$gs$ID, mset=mset)
  tmp <- mset$gs2gv
  names(tmp) <- mset$gs$ID
  mo3 <- modOverlaps(tmp)

  expect_identical(mo1, mo2)
  expect_identical(mo1, mo3)
  expect_equal(mo1[1, 1], 1)
  expect_equal(mo1[2, 1], 0.008196721)

  mymods <- list(A=LETTERS[1:10], B=LETTERS[5:15], c=LETTERS[5:25])
  mo4 <- modOverlaps(mymods, stat="number")
  expect_equal(mo4[2,1], 6)
  expect_equal(mo4[2,3], 11)


})

test_that("modGroups – finding groups of gene sets", {

  res_e <- read.csv("res_CERNO.csv", stringsAsFactors=FALSE, row.names=1)
  mg <- modGroups(res_e$ID)
  expect_equal(length(mg), 2)
  mg <- modGroups(res_e$ID, min.overlap=5)
  expect_equal(length(mg), 6)
  mg <- modGroups(res_e$ID, min.overlap=0.05, stat="jac")
  expect_equal(length(mg), 7)
  mg <- modGroups(res_e$ID, min.overlap=0.1, stat="soerensen")
  expect_equal(length(mg), 7)
})

