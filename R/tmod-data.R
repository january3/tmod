#' Default gene expression module data
#'
#' Gene expression module data from Chaussabel et al. (2008) and Li et al. (2014)
#'
#' The tmod package includes one data set of class tmod which can be loaded with
#' data(tmod). This data set is derived from two studies (see package
#' vignette for details). By default, enrichment analysis with tmod uses
#' this data set; however, it is not loaded into user workspace by default.
#' @references 
#' Chaussabel, Damien, Charles Quinn, Jing Shen, Pinakeen Patel, Casey Glaser, Nicole Baldwin, Dorothee Stichweh, et al. 2008. 
#' "A Modular Analysis Framework for Blood Genomics Studies: Application to Systemic Lupus Erythematosus." Immunity 29(1):150-64.
#'
#' Li, Shuzhao, Nadine Rouphael, Sai Duraisingham, Sandra Romero-Steiner, Scott Presnell, Carl Davis, Daniel
#' S Schmidt, et al. 2014. "Molecular Signatures of Antibody Responses Derived from a Systems Biology Study
#' of Five Human Vaccines." Nature Immunology 15(2):195-204. 
#' @seealso tmod-class, modmetabo
#' @examples
#' # list of first 10 modules
#' data(tmod)
#' tmod
#' tmod$MODULES[1:10, ]
#' tmod[1:10]
#' @name tmod-data
NULL

#' @name tmod
#' @rdname tmod-data
NULL
