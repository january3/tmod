#' Calculate the eigengene of a module from a data set
#'
#' Calculate the eigengene of a module from a data set
#'
#' The eigengene of a module is here defined as the first principal
#' component of a PCA on the gene expression of all genes from a module.
#'
#' @param x data; genes in rows, samples in columns
#' @param g genes -- a vector gene IDs corresponding to annotation in modules
#' @param mset -- a module set; eigengenes will be calculated for each module in the set
#' @param k which component defines the eigengene (default: 1)
#' @return A numeric matrix with rows corresponding to modules. If there
#' was not a sufficient number of genes in a module corresponding to the data
#' set, the row will contain only NA's.
#' @examples
#' data(Egambia)
#' data(tmod)
#' x <- Egambia[ , -c(1:3) ]
#' ifns <- tmod[ grep("[Ii]nterferon", tmod$MODULES$Title) ]
#' eigv <- eigengene(x, Egambia$GENE_SYMBOL, ifns)
#' plot(eigv["LI.M127", ], eigv["DC.M1.2",])
#'
#' # which interferon modules are correlated
#' cor(eigv) 
#' @export
eigengene <- function(x, g, mset=NULL, k=1) {
  mset <- .getmodules2(NULL, mset)

  x <- t(x)

  # filter zero variance
  vars <- apply(x, 2, var)
  sel <- vars > 0
  if(sum(sel) == 0) stop("No variables remain after filtering (variances all zero)")


  x <- x[, sel, drop=F]
  g <- g[sel]
  x <- scale(x)
  n <- nrow(x)  # number of samples

  .eig <- function(y) {
    pca <- prcomp(y, scale.=FALSE, center=FALSE, rank.=k)
    ret <- pca$x[,k]
    if(sum(sign(pca$rotation[,k])) < 0) ret <- -ret
    ret
  }
 
  ret <- sapply(mset$MODULES$ID, function(id) {
    sel <- g %in% mset$MODULES2GENES[[id]]
    if(sum(sel) < k) {
      ret <- NULL
    } else {
      ret <- .eig(x[, sel, drop=FALSE ])
    }
    ret
  }, simplify=FALSE)

  ret <- ret[!sapply(ret, is.null)]
  ret <- t(simplify2array(ret))

  ret
}



#' Module correlation
#'
#' Calculate the correlation between modules
#'
#' The correlation between modules are defined as 
#' correlation coefficient between the modules eigengenes.
#' These are based on a particular gene expression data set.
#' 
#' This function is a simple wrapper combining eigengene() and cor().
#' @return a matrix of module correlation coefficients
#' @param x a data set, with variables (e.g. genes) in rows and samples in columns
#' @param g a vector of variable idenitifiers which correspond to the definition of modules
#' @param mset a module set
#' @param ... any further parameters will be passed to the cor() function
#' @export
modcors <- function(x, g, mset=NULL, ...) {
  mset <- .getmodules2(NULL, mset)

  eigengenes <- eigengene(x, g, mset)

  cor(t(eigengenes), ...)
}



#' Jaccard index for modules
#'
#' Jaccard index for modules
#'
#' For each pair of modules in mset, calculate the Jacard index between
#' these modules.
#' @return matrix with Jaccard index for each pair of modules in mset
#' @param mset set of modules for which to calculate the Jaccard index. If
#' NULL, the default tmod module set will be used.
#' @param g a list of genes. If NULL, the list of genes from the mset will be
#' used.
#' @export
modjaccard <- function(mset=NULL, g=NULL) {
  mset <- .getmodules2(NULL, mset)
  if(length(mset) < 2) stop("mset must contain at least 2 modules")

  if(is.null(g)) g <- mset$GENES$ID

  n <- length(mset$MODULES2GENES)
  mat <- sapply(mset$MODULES2GENES, function(x) g %in% x)
  sums <- apply(mat, 2, sum)
  crmat <- crossprod(mat)

  mm <- matrix(sums, nrow=n, ncol=n)

  crmat/(mm + t(mm) - crmat)
}


#' Calculate overlaps of the modules
#' @param modules either a character vector with module IDs from mset, or a list which
#'        contains the module members
#' @param stat Type of statistics to return. "jaccard": Jaccard index;
#' number: number of overlapping genes;
#'        otherwise, total number of common genes will be returned
#' @inheritParams tmodUtest
#' @export
modOverlaps <- function(modules, mset=NULL, stat="jaccard") {
  if(!is.list(modules)) {
    mset <- .getmodules2(modules, mset)
    modules <- mset$MODULES2GENES[modules]
  }

  stat <- match.arg(stat, c("jaccard", "number"))

  g <- unique(unlist(modules))
  mat <- sapply(modules, function(x) g %in% x)
  crmat <- crossprod(mat)
  sums <- diag(crossprod(mat))
  n <- length(sums)
  mm <- matrix(sums, nrow=n, ncol=n)

  if(stat == "jaccard") crmat <- crmat/(mm + t(mm) - crmat)

  return(crmat)
}




#' Find group of modules 
#'
#' Find group of modules  based on shared genes
#'
#' Split the modules into groups based on the overlapping items.
#'
#' The first argument, modules, is either a character vector of module identifiers from `mset`
#' (NULL mset indicates the default mset of tmod) or a list. If it is a
#' list, then each element is assumed to be a character vector with module
#' IDs.
#' @examples
#' mymods <- list(A=c(1, 2, 3), B=c(2, 3, 4), C=c(5, 6, 7))
#' modGroups(mymods)
#' @param min.overlap Minimum number of overlapping items if stat ==
#'        number, minimum jaccard index if stat == jaccard.
#' @param modules Either a list of modules or character vector. 
#' @param stat statistics to use vor 
#' @inheritParams tmodUtest
#' @export
modGroups <- function(modules, mset=NULL, min.overlap=2, stat="number") {
  stat <- match.arg(stat, c("number", "jaccard"))

  crmat <- modOverlaps(modules, mset, stat=stat)
  crmat[ crmat < min.overlap ] <- 0

  if(!is.list(modules)) {
    mset <- .getmodules2(modules, mset)
    modules <- mset$MODULES2GENES[modules]
  }

  colnames(crmat) <- rownames(crmat) <- names(modules)

  groups <- list()

  .recursive_find <- function(current, m, crmat) {
    current <- c(current, m)

    totest <- colnames(crmat)[ crmat[m, ] > 0 ]
    totest <- setdiff(totest, current)
    while(length(totest) > 0) {
      mm <- totest[1]
      current <- .recursive_find(current, mm, crmat)
      totest <- setdiff(totest, current)
    }

    return(current)
  }

  while(nrow(crmat) > 0) {

    m <- rownames(crmat)[1]
    group <- .recursive_find(c(), m, crmat)
    groups[[m]] <- group

    crmat <- crmat[ setdiff(rownames(crmat), group), , drop=FALSE]
  }

  return(groups)
}



