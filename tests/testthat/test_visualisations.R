library(tmod)
data(Egambia)
data(tmod)
context("Visualisations")


test_that("ggEvidencePlot works", {
  skip_if_not_installed("vdiffr")

  genesets <- data.frame(
    gene_symbol = paste0("gene", 1:10),
    geneset = c(rep("diseaseA", 3), rep("diseaseB", 3), rep("diseaseC", 4))
  )
  ranked_genes <- paste0("gene", 1:10)

  tmod_db <- makeTmodFromDataFrame(genesets)

  p <- ggEvidencePlot(
    ranked_genes,
    m = c("diseaseA", "diseaseB"),
    mset = tmod_db,
    gene.labels = FALSE
  )

  vdiffr::expect_doppelganger("ggEvidencePlot", p)
})
