#' Leading Edge Analysis
#'
#' For each module, return a list of genes on the leading edge
#'
#' Given a vector of ordered gene identifiers and a vector of module IDs,
#' for each module, return the genes which are on the up-slope of the
#' GSEA-style evidence plot. That is, return the genes that are driving the
#' enrichment.
#' @param l list of genes
#' @param modules character vector of module IDs for which to run the LEA
#' @param mset Which module set to use. Either a character vector ("LI", "DC" or "all", default: LI) or an object of class tmod
#' @inheritParams tmodUtest
#' @export
tmodLEA <- function(l, modules, mset="all", nodups=TRUE, filter=FALSE) {
  # prepare the variables specific to that test
  l <- as.character(l)
  if(nodups) l <- unique(l)

  mset <- .getmodules2(modules, mset)

  if( filter ) l <- l[ l %in% mset$GENES$ID ]

  if( sum( l %in% mset$GENES$ID ) == 0 ) {
    warning( "No genes in l match genes in GENES" )
    return(NULL)
  }

  ret <- lapply(modules, function(m) {
    x <- l %in% mset$MODULES2GENES[[m]] 
    mi <- which.max(cumsum(x)/sum(x) - cumsum(!x)/sum(!x))
    #mi <- which.max(cumsum(x)/1:length(l))
    .ret <- l[1:mi][ x[1:mi] ]
    attr(.ret, "N") <- sum(x)
    attr(.ret, "LEA") <- mi
    attr(.ret, "LEA.frac") <- sum(x[1:mi])/sum(x)
    return(.ret)
  })

  names(ret) <- modules
  class(ret) <- c("tmodLEA", class(ret))

  return(ret)
}

#' Summary stats of a leading edge analysis
#'
#' Summary stats of a leading edge analysis
#' @param lea result of `tmodLEA`
#' @param labels labels to add to the result; if NULL, the labels will be
#'        taken from `mset`
#' @param genes if TRUE, the gene identifiers of the leading edge (joined by commas) will be appended.
#' @inheritParams tmodUtest
#' @return data frame with summary stats
#' @export 
tmodLEASummary <- function(lea, genes=FALSE, labels=NULL, mset=NULL) {
  if(!is(lea, "tmodLEA")) 
    stop("parameter lea must be of class `tmodLEA`")

  ret <- data.frame(ID=names(lea))
  if(is.null(labels)) {
    mset <- .getmodules2(names(lea), mset)
    labels <- mset$MODULES$Title
  }
  labels[is.na(labels)] <- ""
  ret$Title <- labels
  ret$N <- sapply(lea, attr, "N")
  ret$LEA <- sapply(lea, attr, "LEA")
  ret$fraction <- sapply(lea, attr, "LEA.frac")

  if(genes) {
    ret$genes <- sapply(lea, paste, collapse=",")
  }


  return(ret)
}
