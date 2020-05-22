#' Leading Edge Analysis

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
    attr(.ret, "LEA") <- mi
    attr(.ret, "LEA.frac") <- sum(x[1:mi])/sum(x)
    return(.ret)
  })

  names(ret) <- modules

  return(ret)
}
