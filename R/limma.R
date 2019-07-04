#' Run tmod enrichment tests directly on a limma object
#'
#' Order the genes according to each of the coefficient
#' found in a limma object and run an enrichment test on the ordered list.
#'
#' For each coefficient in the fit returned by the eBayes / lmFit functions
#' from the limma package, tmodLimmaTest will order the genes run an enrichment test and return
#' the results.
#'
#' The ordering of the genes according to a certain metric is the fundament for
#' gene enrichment analysis. tmodLimmaTest allows three orderings: p-values, "MSD" and log fold changes.
#' The default MSD ("minimal significant difference") is the lower boundary of the 95%
#' confidence interval for positive log fold changes, and 0 minus the upper boundary
#' of the 95% confidence interval for negative log fold changes. MSD performs
#' better than ordering by p-value or by log fold change. See discussion in
#' the package vignette.
#'
#' @return A list with length equal to the number of coeffients. Each
#' element is the value returned by tmod test function. The list can be
#' directly passed to the functions tmodSummary and tmodPanelPlot.
#' @param f result of linear model fit produced by limma functions lmFit and eBayes
#' @param genes Either the name of the column in f$genes which
#' contains the gene symbols corresponding to the gene set collection used, or
#' a character vector with gene symbols
#' @param sort.by How the gene names should be ordered: "msd" (default), "pval" or "lfc"
#' @param tmodFunc The function to run the enrichment tests. Either tmodCERNOtest or tmodUtest
#' @param coef If not NULL, only run tmod on these coefficients
#' @param ... Further parameters passed to the tmod test function
#' @examples
#' data(Egambia)
#' design <- cbind(Intercept=rep(1, 30), TB=rep(c(0,1), each= 15))
#' if(require(limma)) {
#'   fit <- eBayes( lmFit(Egambia[,-c(1:3)], design))
#'   ret <- tmodLimmaTest(fit, genes=Egambia$GENE_SYMBOL)
#'   tmodSummary(ret)
#'   tmodPanelPlot(ret)
#' }
#' @seealso tmodCERNOtest, tmodUtest, tmodPlotPanel, tmodSummary
#' @export
tmodLimmaTest <- function(f,
  genes,
  sort.by="msd",
  tmodFunc=tmodCERNOtest,
  coef=NULL,
  ...
) {

  # sanity checks
  .check.limma(f)

  sort.by <- match.arg(sort.by, c("msd", "pval", "lfc"))

  N  <- ncol(f$coefficients)
  cn <- coef

  if(is.null(cn)) {
    cn <- colnames(f$coefficients)
  }

  if(is.null(cn)) {
    cn <- 1:ncol(f$coefficients)
  }
  
  gl <- genes

  # gene list in the f$genes
  if(length(gl) == 1) {
    if(!gl %in% colnames(f$genes)) 
      stop(sprintf("column %s not in f$genes", gl))
    gl <- f$genes[,gl]
  }

  # calculate MSD
  if(sort.by == "msd") {
    if(!all(c("s2.post", "stdev.unscaled", "df.total") %in% names(f)))
      stop("Incorrect f object, did you use eBayes?")

    # this is from limma's toptable() function
    msd <- .limma.msd(f)
  }

  # define the test function
  tmodCoefTest <- function(cn) {
    ord <- switch( sort.by, 
      msd=order(msd[,cn], decreasing=TRUE),
      lfc=order(abs(f$coefficients[,cn]), decreasing=TRUE),
      pval=order(f$p.value[,cn]))
    tmodFunc( gl[ord], ... )
  }

  # run the tests
  res <- lapply(cn, tmodCoefTest)
  names(res) <- cn
  res
}


#' tmod's replacement for the limma topTable function
#'
#' Produce a data frame for all or for selected coefficients of a linear
#' fit object, including log fold changes, q-values, confidence intervals and
#' MSD.
#'
#' Produce a data frame for all or for selected coefficients of a linear
#' fit object, including log fold changes, q-values, confidence intervals and
#' MSD. For each coefficient, these four columns will be created in the
#' output file, with the name consisting of a prefix indicating the type of
#' the column ("msd", "logFC", "qval", "SE", "ci.L", "ci.R") and the name of the
#' coefficient.
#' @return A data frame with all genes.
#' @param f result of linear model fit produced by limma functions lmFit and eBayes
#' @param genelist A data frame or a character vector with additional
#' information on the genes / probes
#' @param coef Which coefficients to extract
#' @param adjust.method Which method of p-value adjustment; see "p.adjust()"
#' @param confint Confidence interval to be calculated
#' @seealso tmodLimmaTest
#' @export
tmodLimmaTopTable <- function(f, genelist=NULL, coef=NULL, 
  adjust.method="BH", confint=0.95) {

  # sanity checks
  .check.limma(f)

  N  <- ncol(f$coefficients)
  cn <- colnames(f$coefficients)
  Nr <- nrow(f$coefficients)

  if(is.null(cn)) 
    cn <- 1:ncol(f$coefficients)

  if(!is.null(coef)) {
    if(!all(coef %in% cn)) stop("Unknown coefficients")
    cn <- coef
  }

  pval <- f$p.value
  pval <- apply(pval, 2, function(x) p.adjust(x, method=adjust.method))

  msdci <- .limma.ci.msd(f, confint)

  if(is.null(genelist)) genelist <- f$genes

  ret <- data.frame(placeholder=1:Nr, stringsAsFactors=FALSE)

  if(!is.null(genelist)) 
    ret <- data.frame(ret, genelist, stringsAsFactors=FALSE)

  for(c in cn) {
    ret[,paste0("logFC.", c)] <- f$coefficients[,c]
    ret[,paste0("t.", c)]     <- f$t[,c]
    ret[,paste0("msd.", c)]   <- msdci$msd[,c]
    ret[,paste0("SE.", c)]    <- msdci$se[,c]
    ret[,paste0("d.", c)]     <- msdci$d[,c]
    ret[,paste0("ciL.", c)]   <- msdci$ci.l[,c]
    ret[,paste0("ciR.", c)]   <- msdci$ci.r[,c]
    ret[,paste0("qval.", c)]  <- pval[,c]
  }

  rownames(ret) <- rownames(f)
  ret$placeholder <- NULL
  ret
}

#' Up- and down-regulated genes in modules based on limma object
#' 
#' For each module in mset and for each coefficient in f$coefficients, this
#' function calculates the numbers of significantly up- and down-regulated genes.
#' 
#' For an f object returned by eBayes(), tmodLimmaDecideTests considers every
#' coefficient in this model (every column of f$coefficients). For each such
#' coefficient, tmodLimmaDecideTests calculates, for each module, the number
#' of genes which are up- or down-regulated.
#'
#' In short, tmodLimmaDecideTests is the equivalent of tmodDecideTests, but
#' for limma objects returned by eBayes().
#' 
#' @param f result of linear model fit produced by limma functions lmFit and eBayes
#' @param genes Either the name of the column in f$genes which
#' contains the gene symbols corresponding to the gene set collection used, or
#' a character vector with gene symbols
#' @param lfc.thr log fold change threshold
#' @param pval.thr p-value threshold
#' @param adjust.method method used to adjust the p-values for multiple testing. See p.adjust(). Default:BH.
#' @param mset Which module set to use (see tmodUtest for details)
#' @param filter.unknown If TRUE, modules with no annotation will be omitted
#' @return A list with as many elements as there were coefficients in f. Each element of the list is a data 
#' frame with the columns "Down", "Zero" and "Up" giving the number of the
#' down-, not- and up-regulated genes respectively. Rows of the data frame
#' correspond to module IDs. The object can directly be used in
#' tmodPanelPlot as the pie parameter.
#' @examples
#' data(Egambia)
#' design <- cbind(Intercept=rep(1, 30), TB=rep(c(0,1), each= 15))
#' if(require(limma)) {
#'   fit <- eBayes( lmFit(Egambia[,-c(1:3)], design))
#'   ret <- tmodLimmaTest(fit, Egambia$GENE_SYMBOL)
#'   pie <- tmodLimmaDecideTests(fit, Egambia$GENE_SYMBOL)
#'   tmodPanelPlot(ret, pie=pie)
#' }
#' @seealso tmodDecideTests, tmodLimmaTest, tmodPanelPlot
#' @export
tmodLimmaDecideTests <- function(f, genes, lfc.thr=0.5, pval.thr=0.05, 
  filter.unknown=FALSE,
  adjust.method="BH",
  mset="all") {

  # sanity checks
  .check.limma(f)

  cn <- colnames(f$coefficients)
  if(is.null(cn)) cn <- 1:ncol(f$coefficients)

  if(length(genes) == 1) {
    if(!genes %in% colnames(f$genes)) 
      stop(sprintf("column %s not in f$genes", genes))
    genes <- f$genes[,genes]
  }

  lfc  <- f$coefficients
  pval <- f$p.value
  pval <- apply(pval, 2, function(x) p.adjust(x, method=adjust.method))

  tmodDecideTests(genes, lfc=lfc, pval=pval, 
    lfc.thr=lfc.thr, 
    pval.thr=pval.thr,
    labels=cn,
    filter.unknown=filter.unknown,
    mset=mset)

}


.limma.ci.msd  <- function(f, confint=0.95) {

  qq <- (1+confint)/2

  # Bayes posterior standard deviation
  s <- sqrt(f$s2.post)

  # this is from limma's toptable() function
  se  <- s * f$stdev.unscaled

  # Cohen's d 
  d   <- s * f$coefficients

  margin.error <- se * qt(qq, df=f$df.total)

  # "minimum significant difference"
  msd <- abs(f$coefficients) - margin.error
  
  # confidence intervals
  ci.l <- f$coefficients - margin.error
  ci.r <- f$coefficients + margin.error

  list(msd=msd, se=se, d=d, ci.l=ci.l, ci.r=ci.r)
}

.limma.msd <- function(f) {

  # this is from limma's toptable() function
  margin.error <- sqrt(f$s2.post) * f$stdev.unscaled * qt(0.975, df=f$df.total)
  msd <- abs(f$coefficients) - margin.error
  msd
}

.check.limma <- function(f) {

  if(!is(f, "MArrayLM")) stop("f must be the result of the limma functions lmFit/eBayes")
  if(is.null(f$p.value)) stop("use eBayes from limma to calculate the p-values first")
  if(!all(c("p.value", "s2.post", "stdev.unscaled", "df.total") %in% names(f)))
      stop("Incorrect f object, did you use eBayes?")

}

