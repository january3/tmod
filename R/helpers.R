# helper wrappers around sprintf
.catf   <- function( ... ) cat( sprintf( ... ) )
.printf <- function( ... ) print( sprintf( ... ) )



#' Return the contents of a module 
#'
#' Return the contents of a module
#'
#' This function returns the selected modules from a module set.
#' @examples 
#' # show the interferon related modules
#' getModuleMembers(c("LI.M127", "LI.M158.0", "LI.M158.0"))
#' @return A list of modules
#' @param x a character vector of module names 
#' @param mset optional, a module set
#' @export
getModuleMembers <- function(x, mset="all") {

  mset <- .getmodules2(x, mset)

  mset$MODULES2GENES
}
