#' Cell type signatures
#'
#' Cell type signatures
#'
#' @format An object of class tmodGS
#' @examples
#' ## to use cell signatures, type
#' data(cell_signatures)
#' data(vaccination)
#' gl <- vaccination$GeneName[ order(vaccination$qval.F.D1) ]
#' tmodCERNOtest(gl, mset=cell_signatures)
#' @source CIBERSORT, CellMarkers, PanglaoDB
#' @name cell_signatures
NULL
