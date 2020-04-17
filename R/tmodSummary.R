





#' Create a summary of multiple tmod analyses
#'
#' Create a summary of multiple tmod analyses
#'
#' This function is useful if you run an analysis for several conditions or
#' time points and would like to summarize the information in a single data
#' frame. You can use lapply() to generate a list with tmod results and use
#' tmodSummary to convert it to a data frame.
#' @param x list, in which each element has been generated with a tmod test function
#' @param clust whether, in the resulting data frame, the modules should be
#' ordered by clustering them with either q-values ("qval") or the effect size
#' ("effect"). If "sort" or NULL, the modules are sorted alphabetically by their ID.
#' If "keep", then the order of the modules is kept.
#' @param select a character vector of module IDs to show. If clust == "keep", then in that particular
#' order.
#' @param filter.empty If TRUE, all elements (columns) with no significant enrichment will be removed
#' @param filter.unknown If TRUE, modules with no annotation will be omitted
#' @return a data frame with a line for each module encountered anywhere in
#' the list x, two columns describing the module (ID and module title), and
#' two columns(effect size and q value) for each element of list x.
#' @param effect.col The name of the column with the effect size
#' @param pval.col The name of the p-value column
#' @seealso tmodPanelPlot
#' @examples
#' data(Egambia)
#' E <- Egambia[,-c(1:3)]
#' pca <- prcomp(t(E), scale.=TRUE)
#'
#' # Calculate enrichment for each component
#' gs   <- Egambia$GENE_SYMBOL
#' gn.f <- function(r) {
#'     tmodCERNOtest(gs[order(abs(r), 
#'                 decreasing=TRUE)], 
#'                 qval=0.01)
#' }
#' x <- apply(pca$rotation, 2, gn.f)
#' tmodSummary(x)
#' @export
tmodSummary <- function(x, clust=NULL, filter.empty=FALSE, filter.unknown=TRUE,
  select=NULL, 
  effect.col="AUC", pval.col="adj.P.Val") {

  if(!is.null(clust)) 
    clust <- match.arg(clust, c("qval", "effect", "keep", "sort"))
  if(!is(x, "list")) stop( "x must be a list object")

  rid <- names(x)
  if(is.null(rid)) {
    rid <- paste0("X.", 1:length(x))
    names(x) <- rid
  }

  if(!is.null(select)) {
    all.mods <- select 
  } else {
    all.mods <- unique(unlist(lapply(x, `[[`, "ID")))
    if(is.null(clust) || clust == "sort") all.mods <- sort(all.mods)
  }

  ret <- data.frame( ID=all.mods, Title=rep(NA, length(all.mods)), stringsAsFactors=FALSE)
  rownames(ret) <- ret$ID

  # collect the Title, effect and q-value information
  for(n in rid) {
    .x <- x[[n]]

    if(!all(c(effect.col, pval.col) %in% colnames(.x)))
      stop(sprintf("colnames for %s lack either column %s or column %s", n, effect.col, pval.col))

    if(filter.unknown) {
      .x <- .x[ ! .x$Title %in% c(NA, "", "Unknown", "Undetermined", "TBA" ), ] 
    }

    if(filter.empty && nrow(.x) == 0) next ;

    sel <- all.mods %in% .x$ID
    mm  <- match(all.mods[sel], .x$ID)

    ## select the effect size column - depending on the test type
    #effect.col <- which(colnames(.x) %in% c("AUC", "E"))[1]
    #effect.col <- colnames(.x)[effect.col]

    ret[sel, "Title"] <- .x[mm, "Title", drop=TRUE]
    an <- paste0(effect.col, ".", n)
    ret[sel,an] <- .x[mm, effect.col, drop=TRUE]
    qn <- paste0("q.", n)
    ret[sel,qn]   <- .x[mm, pval.col, drop=TRUE]
  }

  # Remove these which are still empty
  ret <- ret[ !is.na(ret$Title), , drop=FALSE]

  # reorder rows if clustering enabled
  if( !is.null(clust) && clust %in% c("effect", "qval")) {
    m <- ret[,-c(1:2)]
    Ncol <- ncol(m)/2
    if(clust == "effect") {
      m <- m[,(1:Ncol)*2-1]
      m[is.na(m)] <- 0.5
    } else {
      m <- m[(1:Ncol)*2]
      m[is.na(m)] <- 1
    }
    o <- hclust(dist(m))$order
    ret <- ret[o,]
  }

#  attr(ret, "tmodSummary") <- TRUE
#  attr(ret, "pval.col") <- pval.col
#  attr(ret, "effect.col") <- effect.col
#  attr(ret, "rid") <- rid
  ret <- new("tmodSummary", ret, pval.col=pval.col, effect.col=effect.col, rid=rid)
  ret
}


