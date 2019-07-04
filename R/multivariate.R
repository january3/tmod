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
#' @param ... Any further parameters are passed to the tagcloud function
#' @return Either NULL or whatever tagcloud returns
#' @examples
#' data(tmod)
#' fg <- tmod$MODULES2GENES[["LI.M127"]]
#' bg <- tmod$GENES$ID
#' result <- tmodHGtest( fg, bg )
#' tmodTagcloud(result)
#' @export
tmodTagcloud <- function( results, filter=TRUE, simplify=TRUE, tag.col="Title", weights.col="auto", pval.col="P.Value", ... ) {

  res <- results

  if(!is(res, "data.frame")) stop("res must be a data frame")

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

  if(nrow(res) < 2 ) {
    plot.new()
  } else {
    tagcloud(strmultline(res$Title), weights=res[,weights.col], col=smoothPalette(-log10(res$P.Value), max=5), ceiling=5, ...)
  }

  return(NULL)
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
                     legend=FALSE, ...) {

  mode <- match.arg( mode, c( "simple", "leftbottom", "cross" ) )
  tmodfunc <- match.arg(tmodfunc, c( "tmodCERNOtest", "tmodUtest" ))
  tfunc <- switch(tmodfunc, tmodCERNOtest=tmodCERNOtest, tmodUtest=tmodUtest)
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

  cc <- components

  ret <- list()
  ret$enrichments <- list()

  if( mode == "simple" ) {
    ids <- paste0( "Component", cc )
    res <- tfunc( genes[ order( abs( pca$rotation[,cc[1]] ), decreasing=T ) ], ... )
    tmodTagcloud( res, filter=filter, simplify=simplify )
    ret$enrichments[[ids[1]]] <- res

    res <- tfunc( genes[ order( abs( pca$rotation[,cc[2]] ), decreasing=T ) ], ... )
    tmodTagcloud( res, filter=filter, simplify=simplify, fvert=1 )
    ret$enrichments[[ids[2]]] <- res
  } else {
    ids <- paste0( "Component", cc )
    ids <- c( paste0(ids[1], c(".left", ".right")), paste0(ids[2], c(".bottom", ".top")))

    fverts <- c( 0, 0, 1, 1 )
    if(mode == "cross") fverts <- c( 0, 0, 0, 0 )
    
    res <- tfunc( genes[ order( pca$rotation[,cc[1]], decreasing=F ) ], ... )
    tmodTagcloud( res, filter=filter, simplify=simplify, fvert=fverts[1], algorithm="fill" )
    ret$enrichments[[ids[1]]] <- res
    res <- tfunc( genes[ order( pca$rotation[,cc[1]], decreasing=T ) ], ... )
    tmodTagcloud( res, filter=filter, simplify=simplify, fvert=fverts[2], algorithm="fill" )
    ret$enrichments[[ids[2]]] <- res
    res <- tfunc( genes[ order( pca$rotation[,cc[2]], decreasing=F ) ], ... )
    tmodTagcloud( res, filter=filter, simplify=simplify, fvert=fverts[3], algorithm="fill" )
    ret$enrichments[[ids[3]]] <- res
    res <- tfunc( genes[ order( pca$rotation[,cc[2]], decreasing=T ) ], ... )
    tmodTagcloud( res, filter=filter, simplify=simplify, fvert=fverts[4], algorithm="fill" )
    ret$enrichments[[ids[4]]] <- res
  }

  plot.params <- c(list(pca=pca, components=components), plot.params)

  # ------------ main plot --------------------
  ret$plot.return <- do.call(plotfunc, plot.params )
  # ------------ main plot --------------------

  if(legend) {
    par(mar=c(0,0,0,0), usr=c(0,1,0,1))
    plot.new()
    legend("topleft",
      ret$plot.return$groups,
      col=ret$plot.return$colors,
      pch=ret$plot.return$shapes,
      bty="n",
      cex=1.5)
  }
  #plot.new()

  return(invisible(ret))
}
