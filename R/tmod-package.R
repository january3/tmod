#' Transcriptional Module Analysis
#'
#' Transcriptional Module Analysis
#'
#' The primary role of this package is to provide published module
#' assignments between genes and transcriptional modules, as well as tools
#' to analyse and visualize the modules.
#' @seealso \code{\link{tmodHGtest}}, \code{\link{tmodUtest}}
#' @name tmod-package
NULL


.onAttach <- function(libname, pkgname) {
  packageStartupMessage('For tmod user guide, type `vignette("tmod")`')
}
