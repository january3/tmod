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
#' ifns <- tmod[ grep("[Ii]nterferon", tmod$gs$Title) ]
#' eigv <- eigengene(x, Egambia$GENE_SYMBOL, ifns)
#' plot(eigv["LI.M127", ], eigv["DC.M1.2",])
#'
#' # which interferon modules are correlated
#' cor(eigv) 
#' @export
eigengene <- function(x, g, mset=NULL, k=1) {
  mset <- .getmodules_gs(NULL, mset)

  x <- t(x)

  # filter zero variance
  vars <- apply(x, 2, var)
  sel <- vars > 0
  if(sum(sel) == 0) stop("No variables remain after filtering (variances all zero)")


  x <- x[, sel, drop=F]
  g <- g[sel]

  g <- .prep_list(g, mset, filter=FALSE, nodups=FALSE)

  x <- scale(x)
  n <- nrow(x)  # number of samples

  .eig <- function(y) {
    pca <- prcomp(y, scale.=FALSE, center=FALSE, rank.=k)
    ret <- pca$x[,k]
    if(sum(sign(pca$rotation[,k])) < 0) ret <- -ret
    ret
  }
 
  ret <- lapply(1:nrow(mset$gs), function(n_id) {
    sel <- g %in% mset$gs2gv[[n_id]]
    if(sum(sel) < k) {
      ret <- NULL
    } else {
      ret <- .eig(x[, sel, drop=FALSE ])
    }
    ret
  })

  names(ret) <- mset$gs$ID

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
  mset <- .getmodules_gs(NULL, mset)
  if(length(mset) < 2) stop("mset must contain at least 2 modules")

  if(is.null(g)) {
    g <- 1:length(mset$gv)
  } else {
    g <- match(g, mset$gv)
    g <- g[ !is.na(g) ]
    if(!length(g)) {
      stop("No genes in g match genes in mset")
    }
  }

  n <- length(mset$gs2gv)
  mat <- sapply(mset$gs2gv, function(x) g %in% x)
  sums <- apply(mat, 2, sum)
  crmat <- crossprod(mat)

  mm <- matrix(sums, nrow=n, ncol=n)

  crmat/(mm + t(mm) - crmat)
}


#' Calculate overlaps of the modules
#'
#' Calculate overlaps of the modules
#'
#' For a set of modules (aka gene sets) determine the similarity between
#' these. You can run modOverlaps either on a character vector of module
#' IDs or on a list. In the first case, the module elements are taken from
#' `mset`, and if that is NULL, from the default tmod module set. 
#'
#' Alternatively, you can provide a list in which each element is a
#' character vector. In this case, the names of the list are the module IDs,
#' and the character vectors contain the associated elements.
#'
#' The different statistics available are:
#'  * "number": total number of common genes (size of the overlap)
#'  * "jaccard": Jaccard index, i.e. \eqn{\frac{|A \cap B|}{|A \cup B|}}
#'    (number of common elements divided by the total number of unique elements);
#'  * "soerensen": Soerensen-Dice coefficient, defined as \eqn{\frac{2 \cdot |A \cap B|}{|A| + |B|}} – number of common genes in relation to the total number of elements (divided by two, such that the maximum is 1)
#'    (number of common elements divided by the average size of both gene sets)
#'  * "overlap": Szymkiewicz-Simpson coefficient, defined as \eqn{\frac{|A \cap B|}{\min(|A|, |B|)}} – this is the number of common genes scaled by the size of the smaller of the two gene sets
#'    (number of common elements divided by the size of the smaller gene set)
#' 
#' @param modules either a character vector with module IDs from mset, or a list which
#'        contains the module members
#' @param stat Type of statistics to return. 
#'        "jaccard": Jaccard index (default);
#'        "number": number of common genes;
#'        "soerensen": Soerensen-Dice coefficient;
#'        "overlap": Szymkiewicz-Simpson coefficient.
#' @inheritParams tmodUtest
#' @export
modOverlaps <- function(modules=NULL, mset=NULL, stat="jaccard") {

  if(is.null(modules)) {
    if(is.null(mset)) {
      stop("Either mset or modules must be provided")
    }

    modules <- mset$gs2gv
    names(modules) <- mset$gs$ID
  } else if(!is.list(modules)) {
    mset <- .getmodules_gs(modules, mset)
    if(!length(modules)) {
      stop("None of the modules are found in the mset")
    }
    tmp <- mset$gs2gv[ match(modules, mset$gs$ID) ]
    names(tmp) <- modules
    modules <- tmp
  }

  stat <- match.arg(stat, c("jaccard", "number", "soerensen", "overlap"))

  g <- unique(unlist(modules))
  mat <- sapply(modules, function(x) g %in% x)
  crmat <- crossprod(mat)
  sums <- diag(crossprod(mat))
  n <- length(sums)
  mm <- matrix(sums, nrow=n, ncol=n)

  if(stat == "jaccard")   crmat <- crmat/(mm + t(mm) - crmat)
  if(stat == "soerensen") crmat <- 2 * crmat / (mm + t(mm))
  if(stat == "overlap")   {
    minmat <- pmin(mm, t(mm))
    crmat <- crmat / minmat
  }

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
#' mymods <- list(A=c(1, 2, 3), B=c(2, 3, 4, 5), C=c(5, 6, 7))
#' modGroups(mymods)
#' @param min.overlap Minimum number of overlapping items if stat ==
#'        number, minimum jaccard index if stat == jaccard etc.
#' @param modules Either a list of modules or character vector. 
#' @inheritParams tmodUtest
#' @inheritParams modOverlaps
#' @export
modGroups <- function(modules, mset=NULL, min.overlap=2, stat="number") {
  stat <- match.arg(stat, c("number", "jaccard", "soerensen", "overlap"))

  crmat <- modOverlaps(modules, mset, stat=stat)
  crmat[ crmat < min.overlap ] <- 0

  if(!is.list(modules)) {
    mset <- .getmodules_gs(modules, mset)
    tmp <- mset$gs2gv[ match(modules, mset$gs$ID) ]
    names(tmp) <- modules
    modules <- tmp
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


#' Get genes belonging to a gene set
#'
#' Get genes belonging to a gene set
#'
#' Create a data frame mapping each module to a comma separated list of
#' genes. If genelist is provided, then only genes in that list will be
#' shown. An optional column, "fg" informs which genes are in the "foreground"
#' data set.
#' @return data frame containing module to gene mapping, or a list (if
#' as.list == TRUE
#' @param gs gene set IDs; if NULL, returns all genes from all gene sets
#' @param genes character vector with gene IDs. If not NULL, only genes
#'        from this parameter will be considered.
#' @param mset module set to use
#' @param fg genes which are in the foreground set
#' @param as.list should a list of genes rather than a data frame be returned
#' @export
getGenes <- function(gs=NULL, genes=NULL, fg=NULL, mset="LI", as.list=FALSE) {
  mset <- .getmodules_gs(gs, mset)

  g <- .prep_list(gs, mset, FALSE, FALSE)

  if(!is.null(genes)) 
    mset$gs2gv <- lapply(mset$gs2gv, function(x) x[ mset$gv[x] %in% genes ])

  if(as.list) {
    ret <- mset$gs2gv
    names(ret) <- mset$gs$ID
    ret <- lapply(ret, function(x) mset$gv[ x ])
    return(ret)
  }

  ret <- data.frame(ID=mset$gs$ID)
  rownames(ret) <- ret$ID
  ret$N <- sapply(mset$gs2gv, length)
  ret$Genes <- sapply(mset$gs2gv, function(x) paste(mset$gv[ x ], collapse=","))

  if(!is.null(fg)) {
    ret$fg <- sapply(mset$gs2gv, 
      function(x) {
        x <- mset$gv[ x ]
        paste(x[x %in% fg], collapse=",")
    })
  }
  ret
}

#' Filter by genes belonging to a gene set from a data frame
#'
#' Filter a data frame or vector by genes belonging to a gene set
#' 
#' filterGS filters a vector of gene IDs based on whether the IDs belong to
#' a given set of gene sets, returning a logical vector.
#' 
#' The showModule function is deprecated and will be removed in future.
#'
#' @return filterGS returns a logical vector of length equal to genes, with
#' TRUE indicating that the given gene is a member of the gene sets in `gs`.
#' @param x a data frame or a vector
#' @param genes a character vector with gene IDs
#' @param gs a character vector corresponding to the IDs of the gene sets to be shown
#' @param mset Module set to use; see "tmodUtest" for details
#' @param extra no longer used.
#' @examples
#' data(Egambia)
#' ## LI.M127 – type I interferon response
#' sel <- filterGS("LI.M127", Egambia$GENE_SYMBOL)
#' head(Egambia[sel, ])
#' @export
filterGS <- function(genes, gs, mset="all") {
  mset <- .getmodules_gs(gs, mset)

  genes %in% mset$gv
}

#' @rdname filterGS
#' @export
showModule <- function(x, genes, gs, mset="all", extra=NULL) {
  if(is.factor(genes)) genes <- as.character(genes)
  sel <- filterGS(genes, gs, mset)

  subset(x, sel)
}


