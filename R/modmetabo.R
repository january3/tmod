#' Modules for metabolic profiling
#'
#' Feature and data sets for metabolic profiling
#'
#' The module set "modmetabo" can be used with tmod to analyse metabolic profiling
#' data. The clusters defined in this set are based on hierarchical clustering
#' of metabolic compounds from human serum and have been published in a paper on 
#' metabolic profiling in tuberculosis by
#' Weiner et al. (2012).
#' 
#' For an example analysis, "tbmprof" is a data set containing metabolic profiles
#' of serum isolated from tuberculosis (TB) patients and healthy individuals. The tbmprof is
#' a data frame containing observations in rows and metabolite id's
#' (corresponding to the id's in the modmetabo object). See examples below.
#' @references
#' Weiner et al. "Biomarkers of inflammation, immunosuppression
#' and stress with active disease are revealed by metabolomic profiling of
#' tuberculosis patients." PloS one 7.7 (2012): e40221.
#' @seealso tmod-data
#' @examples
#' data(modmetabo)  # module definitions
#' data(tbmprof)    # example data set
#' ids <- rownames(tbmprof)
#' tb  <- factor(gsub("\\..*", "", ids))
#'
#' ## use Wilcoxon test to calculate significant differences
#' wcx <- apply(tbmprof, 2, function(x) wilcox.test(x ~ tb)$p.value)
#' wcx <- sort(wcx)
#' tmodCERNOtest(names(wcx), mset=modmetabo)
#' @name modmetabo
NULL

#' @name tbmprof
#' @rdname modmetabo
NULL
