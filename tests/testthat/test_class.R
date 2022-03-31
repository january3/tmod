library(tmod)

context("Testing tmodGS class")

check_obj <- function(x, gs, gs2gv) {
  expect_s3_class(x, "tmodGS")
  expect_type(x, "list")
  expect_equal(length(x), nrow(gs))

  expect_setequal(tmod_ids(x), gs$ID)

  for(id in gs$ID) {
    mm <- getModuleMembers(id, mset=x)
    expect_setequal(gs2gv[[id]], mm[[1]])
  }
}

test_that("makeTmod works", {

  df <- data.frame(ID=LETTERS[1:4], Title=letters[1:4])
  m2g <- list(A = 1:10, B=11:20, C=21:30, D=31:40)
  tgs <- makeTmodGS(m2g)

  check_obj(tgs, df, m2g)

  m2g2 <- m2g[c(3,1,2,4)]
  tgs <- makeTmodGS(m2g2)

  check_obj(tgs, df, m2g)

  tgs <- makeTmodGS(m2g, gs=df)
  check_obj(tgs, df, m2g)

  df2 <- data.frame(ID=letters[1:4], Title=letters[1:4])
  expect_error(makeTmodGS(m2g, gs=df2))

})


test_that("Class operators work", {

  df <- data.frame(ID=LETTERS[1:4], Title=letters[1:4])
  m2g <- list(A = 1:10, B=11:20, C=21:30, D=31:40)
  tgs <- makeTmodGS(m2g, gs=df)


  tgs2 <- tgs[c("A", "B")]
  check_obj(tgs2, df[1:2, ], m2g[1:2])

  tmod_ids(tgs) <- tolower(df$ID)
  df2 <- data.frame(ID=letters[1:4], Title=letters[1:4])
  m2g2 <- m2g
  names(m2g2) <- tolower(names(m2g2))
  check_obj(tgs, df2, m2g2)

  titles2 <- LETTERS[5:8]
  tmod_titles(tgs) <- titles2
  expect_setequal(tmod_titles(tgs), titles2)

  

})
