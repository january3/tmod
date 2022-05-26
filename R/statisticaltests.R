## this is the function that does the actual testing and concocting the
## results
.tmodTest <- function(mod.test, post.test=NULL, qval= 0.05, order.by= "pval", mset=NULL, 
                      cols=c("ID", "Description", "Title")) {
  if(is.null(mset)) stop("Something went wrong in .tmodTest")
  gs_ids <- mset$gs$ID

  order.by <- match.arg(order.by, c("pval", "none", "auc"))

  # --------------------------------------------------------------
  #                  The actual test
  ret <- lapply(1:length(gs_ids), mod.test)
  ret <- tibble(as.data.frame(Reduce(rbind, ret)))
  if(!is.null(post.test)) ret <- post.test(ret)
  # --------------------------------------------------------------


  cols <- intersect(c("ID", cols), colnames(mset$gs))
  ret2 <- tibble(mset$gs[ ret$n_id, cols, drop=FALSE ])

  ret <- cbind(ret2, ret)
  ret$n_id <- NULL

  ret$adj.P.Val <- p.adjust(ret$P.Value, method= "fdr")
  rownames(ret) <- ret$ID
  ret <- ret[ ret$adj.P.Val < qval,, drop= FALSE ]

  if(order.by == "pval") ret <- ret[ order(ret$P.Value), ]
  if(order.by == "auc")  ret <- ret[ order(ret$AUC), ]
  class(ret) <- c("tmodReport", "colorDF", class(ret))

  # set colorDF column type
  col_type(ret, "P.Value")   <- "pval"
  col_type(ret, "adj.P.Val") <- "pval"

  ret

}


#' Perform a statistical test of module expression
#'
#' Perform a statistical test of module expression
#'
#' Performs a test on either on an ordered list of genes (tmodUtest,
#' tmodCERNOtest, tmodZtest) or on two groups of genes (tmodHGtest).
#' tmodUtest is a U test on ranks of genes that are contained in a module.
#'
#' tmodCERNOtest is also a nonparametric test working on gene ranks, but it
#' originates from Fisher's combined probability test. This test weights
#' genes with lower ranks more, the resulting p-values better correspond to
#' the observed effect size. In effect, modules with small effect but many
#' genes get higher p-values than in case of the U-test.
#'
#' tmodPLAGEtest is based on the PLAGE, "Pathway level analysis of gene
#' expression" published by Tomfohr, Lu and Kepler (2005), doi 10.1186/1471-2105-6-225.
#' In essence it is just a t-test run on module eigengenes, but it
#' performs really well. This approach can be used with any complex linear
#' model; for this, use the function eigengene(). See users guide for details.
#'
#' tmodZtest works very much like tmodCERNOtest, but instead of combining
#' the rank-derived p-values using Fisher's method, it uses the Stouffer
#' method (known also as the Z-transform test). 
#'
#' tmodGeneSetTest is an implementation of the function geneSetTest from
#' the limma package (note that tmodUtest is equivalent to the limma's
#' wilcoxGST function).
#'
#' For a discussion of the above three methods, read M. C. Whitlock,
#' "Combining probability from independent tests: the weighted Z-method is
#' superior to Fisher's approach", J. Evol. Biol. 2005 (doi:
#' 10.1111/j.1420-9101.2005.00917.x) for further details.
#'
#' tmodHGtest is simply a hypergeometric test.
#'
#' In tmod, two module sets can be used, "LI" (from Li et al. 2013), or
#' "DC" (from Chaussabel et al. 2008). Using the parameter "mset", the
#' module set can be selected, or, if mset is "all", both of sets are used.
#'
#' @section Custom module definitions:
#'
#' Custom and arbitrary module, gene set or pathway definitions can be also provided through the mset
#' option, if the parameter is a list rather than a character vector. The
#' list parameter to mset must contain the following members: "MODULES",
#' "MODULES2GENES" and "GENES". 
#' 
#' "MODULES" and "GENES" are data frames. It is required that MODULES
#' contains the following columns: "ID", specifying a unique identifier of a
#' module, and "Title", containing the description of the module. The
#' data frame "GENES" must contain the column "ID".
#'
#' The list MODULES2GENES is a mapping between modules and genes. The names
#' of the list must correspond to the ID column of the MODULES data frame. The members of the
#' list are character vectors, and the values of these vectors must correspond
#' to the ID column of the GENES data frame.
#'
#' @return The statistical tests return a data frame with module names, additional statistic (e.g.
#' enrichment or AUC, depending on the test), P value and FDR q-value (P value corrected for multiple
#' testing using the p.adjust function and Benjamini-Hochberg correction.
#' The data frame has class 'colorDF' (see package colorDF for details),
#' but except for printing using colors on the terminal behaves just like an
#' ordinary data.frame. To strip the coloring, use [colorDF::uncolor()].
#' @param l sorted list of HGNC gene identifiers
#' @param fg foreground gene set for the HG test
#' @param bg background gene set for the HG test
#' @param x Expression matrix for the tmodPLAGEtest; a vector for tmodGeneSetTest
#' @param group group assignments for the tmodPLAGEtest
#' @param modules optional list of modules for which to make the test
#' @param qval Threshold FDR value to report
#' @param nodups Remove duplicate gene names in l and corresponding
#' rows from ranks
#' @param order.by Order by P value ("pval") or none ("none")
#' @param filter Remove gene names which have no module assignments
#' @param mset Which module set to use. Either a character vector ("LI", "DC" or "all", default: LI) or an object of class tmod (see "Custom module definitions" below)
#' @param cols Which columns from the MODULES data frame should be included in resulsts
#' @param useR use the R \code{wilcox.test} function; slow, but with exact p-values for small samples
#' @param Nsim for tmodGeneSetTest, number of replicates for the randomization test
#' @seealso tmod-package
#' @import colorDF
#' @examples 
#' data(tmod)
#' fg <- tmod$MODULES2GENES[["LI.M127"]]
#' bg <- tmod$GENES$ID
#' result <- tmodHGtest( fg, bg )
#' 
#' ## A more sophisticated example
#' ## Gene set enrichment in TB patients compared to 
#' ## healthy controls (Egambia data set)
#' \dontrun{
#' data(Egambia)
#' library(limma)
#' design <- cbind(Intercept=rep(1, 30), TB=rep(c(0,1), each= 15))
#' fit <- eBayes( lmFit(Egambia[,-c(1:3)], design))
#' tt <- topTable(fit, coef=2, number=Inf, genelist=Egambia[,1:3] )
#' tmodUtest(tt$GENE_SYMBOL)
#' tmodCERNOtest(tt$GENE_SYMBOL)
#' }
#' @import stats
#' @export
tmodUtest <- function( l, modules=NULL, qval= 0.05, 
  order.by= "pval", filter= FALSE, mset="LI", 
  cols="Title",
  useR=FALSE, nodups=TRUE) {

  # process mset parameter
  mset <- .getmodules_gs(modules, mset)

  # prepare the variables specific to that test
  l <- .prep_list(l, mset, filter=filter, nodups=nodups)

  if(is.null(l)) {
    warning( "No genes in l match genes in GENES" )
    return(NULL)
  }

  N <- length(l)

  # set up the test function
mod.test <- function(m) {
    x <- l %in% mset$gs2gv[[m]]
    pos <- which(x)
    N1 <- length(pos)
    R1 <- sum(pos)
    if (!useR)
        return(c(n_id = m, N1 = N1, R1 = R1))

    # if useR:
    neg <- which(!x)
    AUC <- sum(1 - pos/N)/N1

    if (N1 == 0) {
        ret <- c(n_id = m, U = NA, N1 = NA, AUC = NA, P.Value = 1)
    } else {
        tt <- wilcox.test(pos, neg, alternative = "less")
        ret <- c(n_id = m, U = unname(tt$statistic), N1 = N1, AUC = AUC, P.Value = tt$p.value)
    }
    ret
}

  # homebrew U test; post.test function
  post.test <- function(ret) {
    if(useR) return(ret)

    # if useR is FALSE
    N1 <- ret$N1
    N2 <- N - ret$N1
    R1 <- ret$R1
    U  <- N1*N2+N1*(N1+1)/2-R1
    Mu <- N1*N2/2
    Su <- sqrt( N1*N2*(N+1)/12 )
    z  <- (-abs(U-Mu)+0.5)/Su
    P  <- pnorm(z)
    AUC <- U/(N1*N2)
    return(data.frame(n_id=ret$n_id, U=U, N1=N1, AUC=AUC, P.Value=P, stringsAsFactors=FALSE))
  }

  ret <- .tmodTest(mod.test, post.test, qval=qval, order.by=order.by, mset=mset, cols=cols)
  attr(ret, "effect_size") <- "AUC"
  attr(ret, "pvalue")      <- "adj.P.Val"
  ret
}

#' @name tmodUtest
#' @export
tmodGeneSetTest <- function(l, x, modules=NULL, qval= 0.05, order.by= "pval", filter= FALSE, mset="LI", 
  cols="Title", Nsim=1000, nodups=TRUE) {

  # process mset parameter
  mset <- .getmodules_gs(modules, mset)

  # prepare the variables specific to that test
  tmp <- .prep_list(l, mset, x=x, filter=filter, nodups=nodups)
  if(is.null(tmp)) {
    warning( "No genes in l match genes in GENES" )
    return(NULL)
  }

  l <- tmp$l
  x <- tmp$x

  N <- length(l)

  # generate N samples
  x.rand <- replicate(Nsim, sample(x))

  # set up the test function
  mod.test <- function(m) {
    xi <- l %in% mset$gs2gv[[m]]
    N1 <- sum(xi)
    mm <- mean(x[xi])
    mm.rand <- apply(x.rand[xi, , drop=FALSE], 2, mean)
    p <- sum(mm.rand > mm)/Nsim
    d <- (mm - mean(mm.rand))/sd(mm.rand)
    ranks <- c(1:N)[xi]
    
    ret <- c(n_id=m, D=d, M=mm, R1=sum(ranks), N1=N1, P.Value=p)
    ret
  }

  post.test <- function(ret) {
    ret <- data.frame(ret, stringsAsFactors=FALSE)
    N1 <- ret$N1
    N2 <- N - N1
    R1 <- ret$R1
    U  <- N1*N2+N1*(N1+1)/2-R1

    ret$AUC <- U/(N1*N2)
    ret <- ret[ , c("n_id", "D", "M", "N1", "AUC", "P.Value" ) ]
    ret$P.Value[ is.na(ret$P.Value) ] <- 1
    ret
  }

  ret <- .tmodTest( mod.test, post.test, qval=qval, order.by=order.by, mset=mset, cols=cols )
  attr(ret, "effect_size") <- "D"
  attr(ret, "pvalue")      <- "adj.P.Val"
  ret
}


## prepare the list of the genes for the analysis
## returns NULL if no genes match the genes in mset
## @param l list of genes (character vector)
## @param mset tmodGS object
## @param nodups whether duplicate genes should be removed
## @param filter whether only genes present in mset should be kept
.prep_list <- function(l, mset, x=NULL, nodups=TRUE, filter=FALSE) {

  # prepare the variables specific to that test
  if(nodups) {
    sel <- !duplicated(l)
    l <- l[ sel ]
    if(!is.null(x)) {
      x <- subset(x, sel)
    }
  }

  l_ret <- match(l, mset$gv)

  if(all(is.na(l_ret))) {
    return(NULL)
  }

  .max_n <- length(mset$gv) + 1
  nas <- is.na(l_ret)
  l_ret[ nas ] <- seq(.max_n, length.out=sum(nas))

  if(filter) {
    sel <- l_ret < .max_n
    l_ret <- l_ret[ sel ]
    if(!is.null(x)) {
      x <- subset(x, sel)
    }
  }

  if(!is.null(x)) {
    return(list(l=l_ret, x=x))
  }

  l_ret
}


#' @name tmodUtest
#' @export
tmodCERNOtest <- function(l, modules=NULL, qval= 0.05, order.by= "pval", filter= FALSE, mset="LI", 
  cols="Title", nodups=TRUE) {

  # process mset parameter
  mset <- .getmodules_gs(modules, mset)

  l <- .prep_list(l, mset, nodups=nodups, filter=filter)

  if(is.null(l)) {
    warning( "No genes in l match genes in GENES" )
  }

  N <- length(l)

  # set up the test function
  mod.test <- function(m) {
    x <- l %in% mset$gs2gv[[m]]
    N1 <- sum(x)

    ranks <- c(1:N)[x]
    cerno <- -2 * sum( log(ranks/N) )
    cES <- cerno/(2*N1)
    ret <- c(n_id=m, cerno=cerno, N1=N1, R1=sum(ranks), cES=cES, P.Value=1)
    ret
  }

  post.test <- function(ret) {
    ret <- data.frame(ret, stringsAsFactors=FALSE)
    N1 <- ret$N1
    N2 <- N - N1
    R1 <- ret$R1
    U  <- N1*N2+N1*(N1+1)/2-R1
    ret$AUC <- U/(N1*N2)
    ret <- ret[ , c("n_id", "cerno", "N1", "AUC", "cES" ) ]
    ret$P.Value= pchisq(ret$cerno, 2*ret$N1, lower.tail=FALSE)
    ret
  }

  ret <- .tmodTest(mod.test, post.test, qval=qval, order.by=order.by, mset=mset, cols=cols)
  attr(ret, "effect_size") <- "AUC"
  attr(ret, "pvalue")      <- "adj.P.Val"
  ret
}


.cohensd <- function(x, g) {

  M <- tapply(x, g, mean)
  VAR <- tapply(x, g, sd)^2
  N <- table(g)

  SDpool <- sqrt( (VAR[1]*(N[1]-1) + VAR[2]*(N[2]-1)) / (sum(N) -2 ) )

  (M[1] - M[2])/SDpool
}


#' @name tmodUtest
#' @export
tmodPLAGEtest <- function(l, x, group, modules=NULL, qval=0.05, order.by="pval", 
                          mset="LI", cols="Title", filter=FALSE, nodups=TRUE) {

  # process mset parameter
  mset <- .getmodules_gs(modules, mset)

  tmp <- .prep_list(l, mset, x=x, nodups=nodups, filter=filter)
  if(is.null(tmp)) {
    warning( "No genes in l match genes in GENES" )
    return(NULL)
  }

  x <- tmp$x
  l <- tmp$l

  if(!is.factor(group)) {
    message("Converting group to factor")
    group <- factor(group)
  }

  if(length(levels(group)) != 2) 
    stop("group must be a vector with exactly two different values")


  if(length(l) != nrow(x))
    stop("length of l must be equal to number of rows in x")
  if(length(group) != ncol(x))
    stop("length of group must be equal to number of columns in x")

  N <- length(l)

  message("Calculating eigengenes...")
  eig <- eigengene(x, mset$gv[l], mset=mset, k=1)

  mset <- mset[ rownames(eig) ]

  mod.test <- function(m) {
    x <- l %in% mset$gs2gv[[m]]
    if(sum(x) < 3) return(c(n_id=m, t=NA, D=NA, AbsD=NA, P.Value=1))

    tt <- t.test(eig[m,] ~ group)
    d <- .cohensd(eig[m,], group)
    ret <- c(n_id=m, t=tt$statistic, D=d, AbsD=abs(d), P.Value=tt$p.value)
    names(ret) <- c("n_id", "t", "D", "AbsD", "P.Value")
    ret
  }

  ret <- .tmodTest(mod.test, qval=qval, order.by=order.by, mset=mset, cols=cols)
  attr(ret, "effect_size") <- "D"
  attr(ret, "pvalue")      <- "adj.P.Val"
  ret
}





#' @name tmodUtest
#' @export
tmodZtest <- function(l, modules=NULL, qval= 0.05, order.by= "pval", 
                      filter= FALSE, mset="LI", cols="Title", nodups=TRUE) {

  # process mset parameter
  mset <- .getmodules_gs(modules, mset)

  # prepare the variables specific to that test
  l <- .prep_list(l, mset, filter=filter, nodups=nodups)
  if(is.null(l)) {
    warning( "No genes in l match genes in GENES" )
    return(NULL)
  }

  N <- length(l)

  # set up the test function
  mod.test <- function(m) {
    x <- l %in% mset$gs2gv[[m]]
    N1 <- sum(x)

    ranks <- c(1:N)[x]
    pvals <- ranks/N
    ret <- c(n_id=m, z=sum(qnorm(pvals)), N1=N1, R1=sum(ranks), P.Value=1)
    ret
  }

  post.test <- function(ret) {
    ret <- data.frame(ret, stringsAsFactors=FALSE)
    N1 <- ret$N1
    N2 <- N - N1
    R1 <- ret$R1
    U  <- N1*N2+N1*(N1+1)/2-R1
    ret$AUC <- U/(N1*N2)
    ret$z <- ret$z / sqrt(ret$N1)
    ret <- ret[ , c( "n_id", "z", "N1", "AUC" ) ]
    ret$P.Value=pnorm(ret$z)
    ret$P.Value[is.nan(ret$P.Value)] <- 1
    ret
  }


  ret <- .tmodTest( mod.test, post.test, qval=qval, order.by=order.by, mset=mset, cols=cols )
  attr(ret, "effect_size") <- "AUC"
  attr(ret, "pvalue")      <- "adj.P.Val"
  ret
}



## tmodWZtest is the weighted version of the tmodZtest. The weights can
## be provided as a named numeric vector (the "weights" parameter). In this case, the weights of any
## genes in the parameter l which are not in the names of "weights" are set to
## 0. These weights will be constant for a given gene in all modules.
## Alternatively, weights associated with the "mset" parameter (an object
## of class tmod) in the "WEIGHTS" member of the object (if present).
##
## @name tmodUtest
### tmodWZtest <- function( l, modules=NULL, weights=NULL, qval= 0.05, order.by= "pval", filter= FALSE, mset="LI", 
###   cols="Title") {
### 
###   # process mset parameter
###   mset <- .getmodules2(modules, mset)
### 
###   if(is.null(weights) && is.null(mset$WEIGHTS)) 
###     stop("If no weights are found in mset, the weights parameter must not be NULL")
### 
###   # prepare the variables specific to that test
###   l <- as.character( unique( l ) )
### 
###   if( sum( l %in% mset$GENES$ID ) == 0 ) {
###     warning( "No genes in l match genes in GENES" )
###     return(NULL)
###   }
### 
###   if( filter ) l <- l[ l %in% mset$GENES$ID ]
###   N <- length( l )
### 
###   if(!is.null(weights)) {
###     ww <- rep(0, N)
###     names(ww) <- l
###     weights <- weights[names(weights) %in% l]
###     ww[names(weights)] <- weights
###     weights <- ww
###   }
### 
### 
###   # set up the test function
###   mod.test <- function( m ) {
###     x <- l %in% mset$MODULES2GENES[[m]]
###     N1 <- sum(x)
### 
###     ranks <- c(1:N)[x]
###     pvals <- ranks/N
###     if(is.null(weights)) {
###       w <- mset$WEIGHTS[[m]][ l[x] ]
###     } else {
###       w <- weights[l[x]]
###     }
###     print(w)
###     z <- sum(qnorm(pvals) * w)/ sqrt(sum(w^2))
###     print(z)
###     ret <- c(z=z, N1=N1, R1=sum(ranks), P.Value=1)
###     ret
###   }
### 
###   post.test <- function(ret) {
###     ret <- data.frame(ret, stringsAsFactors=FALSE)
###     N1 <- ret$N1
###     N2 <- N - N1
###     R1 <- ret$R1
###     U  <- N1*N2+N1*(N1+1)/2-R1
###     ret$AUC <- U/(N1*N2)
###     ret <- ret[ , c( "z", "N1", "AUC" ) ]
###     ret$P.Value=pnorm(ret$z)
###     ret$P.Value[is.nan(ret$P.Value)] <- 1
###     ret
###   }
### 
### 
###   ret <- .tmodTest( mod.test, post.test, qval=qval, order.by=order.by, mset=mset, cols=cols )
###   ret
### }
### 

#' Calculate AUC
#'
#' Calculate AUC
#'
#' tmodAUC calculates the AUC and U statistics. The main purpose of this
#' function is the use in randomization tests. While tmodCERNOtest and
#' tmodUtest both calculate, for each module, the enrichment in a single
#' sorted list of genes, tmodAUC takes any number of such sorted lists. Or,
#' actually, sortings -- vectors with ranks of the genes in each replicate.
#' 
#' Note that the input for this function
#' is different from tmodUtest and related functions: the ordering of l
#' and the matrix ranks does not matter, as long as the matrix ranks
#' contains the actual rankings. Each column in the ranks matrix is treated as
#' a separate sample.
#'
#' Also, the `nodups` parameter which is available (and TRUE by default) for 
#' other tests cannot be used here. This means that the AUCs calculated here
#' might be slightly different from the AUCs calculated with default
#' parameters in tests such as the [tmodCERNOtest()]. Use `nodups=FALSE`
#' with [tmodCERNOtest()] to obtain identical results as with `tmodAUC`.
#' @return A matrix with the same number of columns as "ranks" and as many
#' rows as there were modules to be tested.
#' @param l List of gene names corresponding to rows from the ranks matrix
#' @param ranks a matrix with ranks, where columns correspond to samples
#' and rows to genes from the l list
#' @param modules optional list of modules for which to make the test
#' @param filter Remove gene names which have no module assignments
#' @param stat Which statistics to generate. Default: AUC
#' @param recalculate.ranks Filtering and removing duplicates will also
#' remove ranks, so that they should be recalculated. Use FALSE if you don't
#' want this behavior. If unsure, stay with TRUE
#' @param mset Which module set to use. "LI", "DC" or "all" (default: LI)
#' @seealso tmod-package
#' @examples 
#' data(tmod)
#' l <- tmod_ids(tmod)
#' ranks <- 1:length(l)
#' res <- tmodAUC(l, ranks)
#' head(res)
#' @export
tmodAUC <- function(l, ranks, modules=NULL, 
  stat="AUC",
  recalculate.ranks=TRUE, 
  filter=FALSE, mset="LI") {

  stat <- match.arg(stat, c( "AUC", "U"))

  # prepare the variables 
  if(is.vector(ranks)) ranks <- matrix(ranks, ncol=1)
  if(is.null(colnames(ranks))) {
    colnames(ranks) <- paste0("ID", 1:ncol(ranks))
  }

  stopifnot(length(l) == nrow(ranks))

  mset    <- .getmodules_gs(modules, mset)
  modules <- mset$gs$ID

  tmp <- .prep_list(l, mset, x = ranks, filter=filter, nodups=FALSE)
 
  if(is.null(tmp)) {
    warning( "No genes in l match genes in GENES" )
    return(NULL)
  }

  l     <- tmp$l
  ranks <- tmp$x
 
  N <- length(l)
 
  if(recalculate.ranks) { ranks <- apply(ranks, 2, rank) }

  ret <- sapply(1:length(modules), function(m) {
    x <- l %in% mset$gs2gv[[m]]
    N1 <- sum(x)
    R1 <- apply(ranks[x, , drop=F], 2, sum )
    c(N1, R1)
  })

  ret <- t(ret)
  #colnames(ret) <- c("n_id", "N1", "R1")
  rownames(ret) <- mset$gs$ID
  N1  <- ret[, 1 ]
  ret <- ret[, -1, drop=F]

  N2 <- N - N1
  if(stat == "AUC" ) {
    ret <- apply(ret, 2, function(r1) 1+(N1*(N1+1)/2-r1)/(N1*N2))
  } else if(stat == "U" ) {
    ret <- apply(ret, 2, function(r1) 1+(N1*(N1+1)/2-r1))
  }

  rownames(ret) <- modules

  ret
}


#' @name tmodUtest
#' @export
tmodHGtest <- function(fg, bg, modules=NULL, qval= 0.05, order.by= "pval", filter= FALSE, mset="LI",
  cols="Title", nodups=TRUE) {

  mset <- .getmodules_gs(modules, mset)

  fg <- .prep_list(fg, mset, filter=filter, nodups=nodups)
  bg <- .prep_list(bg, mset, filter=filter, nodups=nodups)

  if(is.null(bg)) {
    warning( "No genes in bg match any of the genes in the GENES" )
  }

  if(is.null(fg)) {
    warning( "No genes in fg match any of the genes in the GENES" )
    return( NULL )
  }

  bg <- setdiff(bg, fg)
  if(length(bg) == 0) stop( "All features from bg in fg." )

  tot <- unique(c( fg, bg ))
  n   <- length(tot)
  k   <- length(fg)

  mod.test <- function(n_id) {
    mg <- mset$gs2gv[[n_id]]
    q <- sum( fg %in% mg ) 
    m <- sum( tot %in% mg )

    if( m == 0 ) { E <- NA } 
    else { E <- (q/k)/(m/n) }

    if( q == 0 || m == 0 ) return( c(n_id=n_id, b=q, B=m, n=k, N=n, E=E, P.Value=1 ) )

    pv <- phyper( q-1, m, n-m, k, lower.tail= FALSE )
    c(n_id=n_id, b=q, B=m, n=k, N=n, E=E, P.Value=pv)
  }


  ret <- .tmodTest( mod.test, NULL, qval=qval, order.by=order.by, mset=mset, cols=cols )
  attr(ret, "effect_size") <- "E"
  attr(ret, "pvalue")      <- "adj.P.Val"
  ret
}



#' Count the Up- or Down-regulated genes per module
#' 
#' For each module in a set, calculate the number of genes which are in
#' that module and which are significantly up- or down-regulated.
#'
#' This function can be used to decide whether a module, as a whole, is up- or
#' down regulated. For each module, it calculates the number of genes which
#' are up-, down- or not regulated. 
#' A gene is considered to be up- regulated
#' if the associated p-value is smaller than pval.thr and the associated log
#' fold change is greater than lfc.thr.
#' A gene is considered to be down- regulated
#' if the associated p-value is smaller than pval.thr and the associated log
#' fold change is smaller than lfc.thr.
#'
#' Note that unlike decideTests from limma, tmodDecideTests does not correct
#' the p-values for multiple testing -- therefore, the p-values should already
#' be corrected.
#'
#' @param g a character vector with gene symbols
#' @param lfc a numeric vector or a matrix with log fold changes
#' @param pval a numeric vector or a matrix with p-values. Must have the same dimensions as lfc
#' @param lfc.thr log fold change threshold
#' @param pval.thr p-value threshold
#' @param labels Names of the comparisons. Either NULL or a character
#' vector of length equal to the number of columns in lfc and pval.
#' @param filter.unknown If TRUE, modules with no annotation will be omitted
#' @param mset Which module set to use. Either a character vector ("LI", "DC" or "all", default: LI) or a list (see "Custom module definitions" below)
#' @return A list with as many elements as there were comparisons (columns in lfc and pval). Each element of the list is a data 
#' frame with the columns "Down", "Zero" and "Up" giving the number of the
#' down-, not- and up-regulated genes respectively. Rows of the data frame
#' correspond to module IDs.
#' @seealso tmodSummary, tmodPanelPlot, tmodDecideTestsLimma
#' @export
tmodDecideTests <- function(g, lfc=NULL, pval=NULL, lfc.thr=0.5, pval.thr=0.05, 
  labels=NULL,
  filter.unknown=FALSE, mset="all") {

  mset <- .getmodules_gs(NULL, mset, known.only=filter.unknown)

  g <- .prep_list(g, mset=mset, filter=FALSE, nodups=FALSE)

  # just count the number of genes in each module
  if(is.null(lfc) && is.null(pval)) {
    N <- sapply(1:length(mset$gs$ID), function(m) sum( g %in% mset$gs2gv[[m]] ))
    return(data.frame(ID=mset$gs$ID, N=N))
  }

  # distinguish the one-dimensional case from the multi-dimensional case
  if( (!is.null(lfc)  && is.null(dim(lfc))) ||
      (!is.null(pval) && is.null(dim(pval))) ) {
    one.dim <- TRUE
  } else {
    one.dim <- FALSE
  }

  nr <- NULL; nc <- NULL ;

  # boring sanity checks
  if(!is.null(pval)) {
    pval <- as.matrix(pval)
    if(!is.null(colnames(pval)) && is.null(labels)) 
        labels <- colnames(pval)
    nr <- nrow(pval)
    nc <- ncol(pval)
    if(is.null(lfc)) {
      lfc <- matrix(Inf, ncol=nc, nrow=nr)
    }
  }

  pval[ is.na(pval) ] <- 1

  if(!is.null(lfc)) {
    lfc  <- as.matrix(lfc)
    if(is.null(labels) && !is.null(colnames(lfc))) 
        labels <- colnames(lfc)
    if(is.null(nr)) nr <- nrow(lfc)
    if(is.null(nc)) nc <- ncol(lfc)
    if(is.null(pval)) {
      pval <- matrix(0, ncol=nc, nrow=nr)
    }
  }

  lfc[ is.na(lfc) ] <- 0

  if(length(pval.thr) == 1) pval.thr <- rep(pval.thr, nc)
  if(length(lfc.thr) == 1)  lfc.thr  <- rep(lfc.thr, nc)
  if(length(pval.thr) != nc) stop(sprintf("pval.thr must have length 1 or %d", nc))
  if(length(lfc.thr) != nc)  stop(sprintf("lfc.thr must have length 1 or %d", nc))

  if(is.null(labels)) labels <- paste0( "X.", 1:nc)

  if(!identical(dim(lfc), dim(pval)) || nrow(lfc) != length(g)) {
    stop("Dimensions of lfc and pval must be identical")
  }

  if(length(labels) != ncol(lfc)) {
    stop("Length of labels must be equal to ncol(lfc) and ncol(pval)")
  }

  lfc.a <- abs(lfc)

  signif.up   <- sapply(1:nc, function(i) pval[,i] < pval.thr[i] & lfc[,i] > lfc.thr[i])
  signif.down <- sapply(1:nc, function(i) pval[,i] < pval.thr[i] & lfc[,i] < -lfc.thr[i])

  # set up the function for counting 
  count.m <- function(m) {
    sel  <- which(g %in% mset$gs2gv[[m]])
      up   <- colSums(signif.up[sel,,drop=F  ])
      down <- colSums(signif.down[sel,,drop=F ])

    N <- length(sel)
    return(cbind(down=down, N=N - (up + down), up=up))
  }

  res  <- lapply(1:length(mset$gs$ID), count.m)
  nmod <- length(res)
  res <- simplify2array(res)

  ret <- lapply(1:nc, function(i) { .x <- t(res[i,,]) ; 
                rownames(.x) <- mset$gs$ID ; .x })

  names(ret) <- labels
  return(ret)
}




