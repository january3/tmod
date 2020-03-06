#' Tag cloud based on tmod results
#'
#' Plot a tag (word) cloud based on results from tmod enrichment. 
#'
#' The tags will be generated based on results from tmod or any other
#' suitable data frame. The data frame must contain two numeric columns,
#' specified with "weights.col" and "pval.col", which will be
#' used to calculate the size and shade of the tags, respectively.
#' Furthermore, it has to contain a column with tags (parameter "tag.col",
#' by default "Title").
#' 
#' Any data frame can be used as long as it contains the specified columns.
#' @param results data frame produced by one of the tmod enrichment tests
#' @param filter Whether redundant and not annotated modules should be removed
#' @param simplify Whether module names should be simplified
#' @param tag.col Which column from results should be used as tags on the plot
#' @param weights.col Which column from results should be used as weights for the tag cloud
#' @param pval.col Which column contains the P values which will be used to shade the tags
#' @param plot Should the tag cloud be plotted or only returned
#' @param min.auc Minimal AUC to show (default: 0.5)
#' @param max.qval Maximal adjusted p value to show (default: 0.05)
#' @param maxn Maximum number of gene set enrichment terms shown on the plot (if NULL – default – all terms will be shown)
#' @param ... Any further parameters are passed to the tagcloud function
#' @return Either NULL or whatever tagcloud returns
#' @examples
#' data(tmod)
#' fg <- tmod$MODULES2GENES[["LI.M127"]]
#' bg <- tmod$GENES$ID
#' result <- tmodHGtest( fg, bg )
#' tmodTagcloud(result)
#' @export
tmodTagcloud <- function( results, filter=TRUE, simplify=TRUE, tag.col="Title", 
  min.auc=.5, max.qval=.05, plot=TRUE,
  weights.col="auto", pval.col="P.Value", maxn=NULL, ... ) {

  res <- results

  if(!is(res, "data.frame")) stop("res must be a data frame")
  if(!is.null(maxn) && nrow(res) > maxn) res <- res[1:maxn, ] 

  res <- res[ res$AUC > min.auc, ]
  res <- res[ res$adj.P.Val < max.qval, ]

  if( weights.col == "auto" ) {
    if( "AUC" %in% colnames(res)) {
      weights.col <- "AUC"
    } else if( "E" %in% colnames(res)) {
      weights.col <- "E"
    } else {
      stop( "Neither E nor AUC columns present in results" )
    }
  }

  if( ! pval.col %in% colnames(res)) {
    stop(sprintf("Column %s not present in results", pval.col ))
  }

  if( ! tag.col %in% colnames(res)) {
    stop(sprintf("Column %s not present in results", tag.col ))
  }

  if( ! weights.col %in% colnames(res)) {
    stop(sprintf("Column %s not present in results", weights.col ))
  }


  ## remove uninteresting results
  if(filter) {
    res <- res[ ! res$Title %in% c( "TBA", "Undetermined", "Not Determined"), ]
  }

  ## simplify the output
  if(simplify) {
    res$Title <- gsub( " \\([IV]*\\)$", "", res$Title )
    res <- res[ ! duplicated( res$Title ), ]
  }

  ret <- NULL
  if(nrow(res) < 2 ) {
    warning("Less than 2 results found, not generating plot")
    if(plot) plot.new()
  } else {
    ret <- tagcloud(strmultline(res$Title), weights=res[,weights.col], col=smoothPalette(-log10(res$P.Value), max=5), ceiling=5, 
      wmin=.5, wmax=1.0, plot=plot,
      ...)
  }

  return(invisible(ret))
}



#' PCA plot annotated with tmod
#'
#' Generate a PCA plot on which each dimension is annotated by a tag cloud
#' based on tmod enrichment test.
#'
#' There are three types of plots that can be generated (parameter "mode"):
#' simple, leftbottom and cross. In the "simple" mode, two enrichments are
#' run, on on each component, sorted by absolute loadings of the PCA
#' components.
#' Both "leftbottom" and "cross" run two enrichment analyses on each component,
#' one on the loadings sorted from lowest to largest, and one on the loadings
#' sorted from largetst to lowest. Thus, two tag clouds are displayed per
#' component. In the "leftbottom" mode, the tag clouds are displayed to the
#' left and below the PCA plot. In the "cross" mode, the tag clouds are
#' displayed on each of the four sides of the plot.
#'
#' By default, the plotting function is pca2d from the pca3d package. Any
#' additional parametrs for pca2d can be passed on using the plot.params
#' parameter. You can define your own function instead of pca2d, however, mind
#' that in any case, there will be two parameters passed to it on the first
#' two positions: pca and components, named "pca" and "components"
#' respectively.
#' @param pca Object returned by prcomp or a matrix of PCA coordinates. In the latter case, a
#' loading matrix must be provided separately.
#' @param loadings A matrix with loadings
#' @param genes A character vector with gene identifiers
#' @param tmodfunc Name of the tmod enrichment test function to use. Either
#' @param plotfunc Function for plotting the PCA plot. See Details
#' @param plot.params A list of parameters to be passed to the plotting function. See
#' Details
#' @param mode Type of the plot to generate; see Details. tmodCERNOtest or tmodUtest (tmodHGtest is not suitable)
#' @param components integer vector of length two: components which components to show on the plot. Must be smaller than the number of columns in pca.
#' @param filter Whether "uninteresting" modules (with no annotation)
#' should be removed from the tag cloud
#' @param simplify Whether the names of the modules should be simplified
#' @param legend whether a legend should be shown
#' @param maxn Maximum number of gene set enrichment terms shown on the plot (if NULL – default – all terms will be shown)
#' @param plot if FALSE, no plot will be shown, but the enrichments will be calculated and returned invisibly
#' @param ... Any further parameters passed to the tmod test function
#' @importFrom tagcloud tagcloud strmultline smoothPalette
#' @return A list containing the calculated enrichments as well as the
#' return value from the plotting function
#' @examples
#' data(Egambia)
#' E <- as.matrix(Egambia[,-c(1:3)])
#' pca <- prcomp(t(E), scale.=TRUE)
#' group <- rep(c("CTRL", "TB"), each=15)
#' tmodPCA(pca, 
#'   genes=Egambia$GENE_SYMBOL, 
#'   components=4:3,
#'   plot.params=list(group=group))
#' @export
tmodPCA <- function(pca, loadings=NULL, genes, 
                     tmodfunc="tmodCERNOtest", 
                     plotfunc=pcaplot,
                     mode="simple", 
                     components=c(1,2), 
                     plot.params=NULL,
                     filter=TRUE, simplify=TRUE, 
                     legend=FALSE, 
                     maxn=NULL,
                     plot=TRUE,
                     ...) {

  mode <- match.arg( mode, c( "simple", "leftbottom", "cross" ) )
  tmodfunc <- match.arg(tmodfunc, c( "tmodCERNOtest", "tmodUtest" ))
  tfunc <- switch(tmodfunc, tmodCERNOtest=tmodCERNOtest, tmodUtest=tmodUtest)

  if(plot) {
    oldpar <- par("mfrow") # try to restore the screen after layout()
    on.exit(par(oldpar))

    if( mode == "simple" ) { 
      layout(matrix(c(2,3,
                      4,1),2,2,byrow=TRUE), 
             widths=c(0.3, 0.7), heights=c(0.7, 0.3))
    } else if( mode == "leftbottom" ) {
      layout(matrix( c(4, 5, 5, 
                       3, 5, 5, 
                       6, 1, 2), 3, 3, byrow=TRUE), 
             widths=rep(1/3, 3), heights=rep(1/3, 3))
    } else if( mode == "cross" ) {
      layout(matrix( c( 6, 4, 0,
                        1, 5, 2,
                        0, 3, 0), 3, 3, byrow=TRUE), 
             widths=c( 1/4, 1/2, 1/4), heights=c( 1/4, 1/2, 1/4 ))
    }

    oldpar <- par( mar=c( 1, 2, 0, 0 ) )
    on.exit( par(oldpar) )
  }

  cc <- components

  ret <- list()
  ret$enrichments <- list()

  if( mode == "simple" ) {
    ## calculate enrichments
    names(res) <- paste0( "Component", cc )
    res <- lapply(1:2, function(i) tfunc( genes[ order( abs( pca$rotation[,cc[i]] ), decreasing=T ) ], ... ))
    names(res) <- ids <- paste0( "Component", cc )
    ret$enrichments <- res
    fverts <- c(0, 0)
    algorithm <- "oval"

    #if(plot) tmodTagcloud( res, filter=filter, simplify=simplify, maxn=maxn )
    #if(plot) tmodTagcloud( res, filter=filter, simplify=simplify, fvert=1, maxn=maxn)
  } else {
    ids <- paste0( "Component", cc )
    ids <- c( paste0(ids[1], c(".left", ".right")), paste0(ids[2], c(".bottom", ".top")))

    fverts <- c( 0, 0, 1, 1 )
    if(mode == "cross") fverts <- c( 0, 0, 0, 0 )
    comp_i <- cc[c(1, 1, 2, 2)]
    decreasing <- c(FALSE, TRUE, FALSE, TRUE)

    ## calculate enrichments, both directions for each component
    res <- lapply(1:4, function(i) 
      tfunc( genes[ order( pca$rotation[,comp_i[i]], decreasing=decreasing[i] ) ], ... )
    )
    names(res) <- ids
    ret$enrichments <- res
    algorithm <- "fill"

    #if(plot) tmodTagcloud( res, filter=filter, simplify=simplify, fvert=fverts[4], algorithm="fill", maxn=maxn )
  }

  if(plot) {

    n <- length(res)
    l_i <- which.max(lapply(res, nrow)) # largest result list

    ## first, calculate a rough scale
    foo <- tmodTagcloud(res[[l_i]], filter=filter, simplify=simplify, fvert=fverts[l_i], algorithm=algorithm, maxn=maxn, plot=FALSE)
    scale <- attr(foo, "scale")
    tagclouds <- lapply(1:n, function(i) {
      tmodTagcloud(res[[i]], filter=filter, simplify=simplify, fvert=fverts[i], algorithm=algorithm, maxn=maxn, scale=scale)
    })
    names(tagclouds) <- ids
    ret$tagclouds <- tagclouds




  }

  plot.params <- c(list(pca=pca, components=components), plot.params)

  # ------------ main plot --------------------
  if(plot) ret$plot.return <- do.call(plotfunc, plot.params )
  # ------------ main plot --------------------

  if(plot && legend) {
    `%.n%` <- function(x, y) if(is.na(x) || is.null(x)) y else x

    par(mar=c(0,0,0,0), usr=c(0,1,0,1))
    plot.new()
    legend("topleft",
      as.character(ret$plot.return$groups),
      col=ret$plot.return$colors %.n% ret$plot.return$col %.n% "black",
      pch=ret$plot.return$shapes %.n% ret$plot.return$pch %.n% 19,
      bty="n",
      cex=1.5)
  }
  #plot.new()

  return(invisible(ret))
}
