# helper wrappers around sprintf
.catf   <- function( ... ) cat( sprintf( ... ) )
.printf <- function( ... ) print( sprintf( ... ) )



#' Return the contents of a gene set 
#'
#' Return the contents of a gene set
#'
#' This function returns the selected gene sets from a collection.
#' @examples 
#' # show the interferon related modules
#' getModuleMembers(c("LI.M127", "LI.M158.0", "LI.M158.0"))
#' @return A list of gene sets
#' @param x a character vector of gene set names 
#' @param mset optional, a gene set collection
#' @examples
#' getModuleMembers("LI.M127")
#' @export
getModuleMembers <- function(x, mset="all") {

  mset <- .getmodules_gs(x, mset)

  ret <- lapply(mset$gs2gv, function(x) mset$gv[ x ])
  names(ret) <- tmod_ids(mset)
  ret
}
