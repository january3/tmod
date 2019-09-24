## adds a grid behind the picture
.plotGrid <- function(grid, x.vec, y.vec, row.h, col.w, grid.col="#33333333", scale.min=NULL, scale.max=NULL, text.cex=1) {

  if(grid == "none") return ;

  Nr <- length(y.vec)
  Nc <- length(x.vec)

  if(grid == "at") {
    segments(x.vec, y.vec[1] - row.h/2, 
             x.vec, y.vec[Nr] + row.h/2, col=grid.col )
    segments(x.vec[1] - col.w/2, y.vec,
             x.vec[Nc] + col.w/2, y.vec, col=grid.col )
  }

  if(grid == "between" || grid == "scale") {
    x.vec <- x.vec - col.w / 2
    x.vec2 <- c(x.vec, x.vec[Nc] + col.w)

    y.vec <- y.vec - row.h/2
    y.vec2 <- c(y.vec, y.vec[Nr] + row.h)

    segments(x.vec2, y.vec2[1], 
             x.vec2, y.vec2[Nr+1], col=grid.col )
    segments(x.vec2[1], y.vec2,
             x.vec2[Nc+1], y.vec2, col=grid.col )

  }


  if(grid == "scale" && !(is.null(scale.min) || is.null(scale.max))) {
    ticks <- axisTicks(c(scale.min, scale.max), log=F)
    tick.height <- row.h / 5
    n.t <- length(ticks)

    x.vec3 <- rep(x.vec, each=n.t) + rep(0:(n.t-1), length(x.vec)) * col.w / (n.t - 1)
    segments(x.vec3, y.vec2[1], x.vec3, y.vec2[Nr + 1], col="#33333333" )

    segments(x.vec2, y.vec2[1], x.vec2, y.vec2[Nr+1], col="#333333cc" )
    segments(x.vec2[1], y.vec2, x.vec2[Nc+1], y.vec2, col="#333333cc" )
    
    labs <- as.character(c(scale.min, scale.max))
    sw <- strwidth(labs, cex=text.cex)
    for(i in 1:Nc) {
      text(c(x.vec[i] + sw[1]/2, x.vec[i] + col.w - sw[2]/2), y.vec[1] - tick.height, labs, pos=1, cex=text.cex, col="#333333cc")
    }
    
    segments(x.vec3, y.vec[1], x.vec3, y.vec[1] - tick.height, col="#333333cc")
  }
}

## return p value formatted for legend
.pformat <- function(p) {
  pp <- 10^-p
  if(pp < 0.001) {
    return(bquote(10^-.(p)))
  }

  if(pp < 0.01) {
    return(sprintf("%.3f", pp))
  }

  return(sprintf("%.2f", pp))
}

## legend for the pvalEffectPlot
.pvalEffectPlotLegend <- function(
  .draw.test, palfunc,
  line.h, col.w, row.h, 
  min.e, max.e,
  min.p, max.p, 
  text.cex,
  legend.style="broad"
 ) {

  N <- 5

  # broad: pvalue and effect size on one row
  # tall: pvalue and effect size one above another
  x.off <- switch(legend.style, broad=0.5, tall=0)
  y.off <- switch(legend.style, broad=0,   tall=line.h * 5 + row.h)
  width <- switch(legend.style, broad=0.5, tall=1)

  # labels
  text(0,             line.h * 2.5 + row.h, "Effect size:", adj=c(0, 0), cex=text.cex[3])
  text(x.off, y.off + line.h * 2.5 + row.h, "P value:", adj=c(0, 0), cex=text.cex[3])

  xc <- seq(col.w, width - col.w, length.out=N)
  yc <- rep(line.h * 2 + row.h/2, N)

  ev <- seq(0, 1, length.out=N)
  ex <- rep(max.p, N)

  # draw the figures for the effect size range (no color)
  .draw.test(NULL, NULL, xc, yc, ev, mq=ex, color=NULL)

  lab.y.pos <- line.h * 0.5

  # smallest and highest effect size
  text(xc[1], lab.y.pos, signif(min.e,2), adj=c(0.5, 0), cex=text.cex[3])
  text(xc[5], lab.y.pos, signif(max.e,2), adj=c(0.5, 0), cex=text.cex[3])

  # a range of significance thresholds
  ex  <- signif(seq(min.p, max.p, length.out=N),2)
  col <- palfunc(1.1*ex)

  # draw the figures for the p values (same size)
  .draw.test(NULL, NULL, x.off + xc, y.off + yc, ev=rep(0.4, N), mq=1.1*ex, color=col)

  for(i in 1:5) {
    text(x.off + xc[i], y.off + lab.y.pos, .pformat(ex[i]), adj=c(0.5, 0), cex=text.cex[3])
  }

}

## returns a function for mapping a value onto a color palette
.getColFunc <- function(min, max, end.col, null.col="#FFFFFF00", start.col="auto", alpha="99", pal.l=12) {

  if(start.col == "auto") {
    ppp <- colorRampPalette(c("white", end.col))(10)
    start.col <- ppp[3]
  }

  pal <- c( null.col, paste0(colorRampPalette(c(start.col, end.col))(pal.l), alpha))
  pal <- c(pal, pal[length(pal)])

  iv <- seq(min, max, length.out=pal.l+1)

  ret <- function(x) {
    ind <- findInterval(x, iv, rightmost.closed=TRUE) + 1
    pal[ind]
  }

  ret
}





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
  select=NULL) {

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
    all.mods <- unique(unlist(lapply(x, function(.) .$ID)))
    if(is.null(clust) || clust == "sort") all.mods <- sort(all.mods)
  }

  ret <- data.frame( ID=all.mods, Title=NA, stringsAsFactors=FALSE)
  rownames(ret) <- ret$ID


  # collect the Title, effect and q-value information
  for(n in rid) {
    if(filter.unknown) {
      x[[n]] <- x[[n]][ ! x[[n]]$Title %in% c( "Undetermined", "TBA" ), ] 
    }

    if(filter.empty && nrow(x[[n]]) == 0) next ;
    sel <- all.mods %in% x[[n]]$ID
    mm  <- match(all.mods[sel], x[[n]]$ID)

    ## select the effect size column - depending on the test type
    effect.col <- which(colnames(x[[n]]) %in% c("AUC", "E"))[1]
    effect.col <- colnames(x[[n]])[effect.col]

    ret[sel,"Title"] <- x[[n]][mm,]$Title
    an <- paste0(effect.col, ".", n)
    ret[sel,an] <- x[[n]][mm,effect.col]
    qn <- paste0("q.", n)
    ret[sel,qn]   <- x[[n]][mm,]$adj.P.Val
  }

  # Remove these which are still empty
  ret <- ret[ !is.na(ret$Title), ]

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

  ret
}




#' Create an effect size / p-value plot
#'
#' Create a heatmap-like plot showing information about both effect size
#' and p-values.
#'
#' pvalEffectPlot shows a heatmap-like plot. Each row corresponds to one
#' series of tests (e.g. one module), and each column corresponds to the time points or conditions for which
#' a given analysis was run. Each significant result is shown as
#' a red dot. Size of the dot corresponds to the effect size (or any
#' arbitrary value), and intensity of the color corresponds to the log10 of
#' p-value.
#'
#' Just like a heatmap corresponds to a single numeric matrix, the pvalue /
#' effect plot corresponds to two matrices: one with the effect size, and
#' another one with the p-values. Each cell in the matrix corresponds to the
#' results of a single statistical test.
#'
#' For example, a number of genes or transcriptional modules might be tested
#' for differential expression or enrichment, respectively, in several conditions. 
#'
#' By default, each test outcome is represented by a dot of varying size
#' and color. Alternatively, a function may be specified with the parameter
#' 'plot.func'. It will be called for each test result to be drawn. The plot.func function must take the following arguments:
#' \itemize{
#' \item{row, col}{either row / column number or the id of the row / column to plot; NULL if drawing legend}
#' \item{x, y}{user coordinates of the result to visualize}
#' \item{w, h}{width and height of the item to plot}
#' \item{e}{Enrichment -- a relative value between 0 and 1, where 0 is the minimum and 1 is the maximum enrichment found}
#' \item{p}{P-value -- an absolute value between 0 and 1}
#' }
#' For the purposes of drawing a legend, the function must accept NULL
#' p-value or a NULL enrichment parameter.
#' @return Invisibly returns a NULL value.
#' @param e matrix with effect sizes
#' @param p matrix with probabilities
#' @param col.labels Labels for the columns. If NULL, names of the elements
#' of the list x will be used.
#' @param pval.thr The p-value must be this or lower in order for a test result to be visualized
#' @param pval.cutoff On visual scale, all p-values below pval.cutoff will be replaced by pval.cutoff
#' @param row.labels Labels for the modules. This must be a named vector, with module IDs as vector names. If NULL, module titles from
#' the analyses results will be used.
#' @param grid Style of a light-grey grid to be plotted; can be "none", "at" and "between"
#' @param grid.color Color of the grid to be plotted (default: light grey)
#' @param col.labels.style Style of column names: "top" (default), "bottom", "both", "none"
#' @param plot.cex a numerical value giving the amount by which the plot
#' symbols will be maginfied 
#' @param text.cex a numerical value giving the amount by which the plot
#' text will be magnified, or a vector containing three cex values for row labels, column labels and legend, respectively
#' @param plot.func Optionally, a function to be used to draw the dots. See "Details"
#' @param legend.style Style of the legend: "auto" -- automatic; "broad":
#' pval legend side by side with effect size legend; "tall": effect size
#' legend above pval legend; "none" -- no legend.
#' @param min.e,max.e scale limits for the effect size
#' @import grDevices
#' @import graphics
#' @import plotwidgets
#' @export
pvalEffectPlot <- function(e, p, 
  pval.thr=0.01, pval.cutoff=1e-6,
  row.labels=NULL, col.labels=NULL, 
  plot.func=NULL,
  grid="at", grid.color="#33333333",
  plot.cex=1, text.cex=1, 
  col.labels.style="top",
  legend.style="auto",
  min.e=NULL,max.e=NULL)  {

   
  # ---------------------------------------
  # initializations
  grid <- match.arg(grid, c("none", "at", "between", "scale"))

  me <- e ; mq <- -log10(p)

  Nc <- ncol(me) ; Nr <- nrow(me)
  if(Nc == 0 || Nr == 0) stop("Nothing to plot")


  min.p <- -log10(pval.thr)
  max.p <- -log10(pval.cutoff)

  if(is.null(min.e)) { min.e <- min(me) }
  if(is.null(max.e)) { max.e <- max(me) }

  if(is.null(row.labels) || length(row.labels) != Nr) row.labels <- 1:Nr
  if(is.null(col.labels) || length(col.labels) != Nc) col.labels <- 1:Nc

  col.labels.style <- match.arg( col.labels.style, c( "top", "bottom", "both", "none" ) )
  if(length(text.cex) == 1) text.cex <- rep(text.cex, 3)
  else if(length(text.cex) != 3) stop("Incorrect text.cex parameter: it should be either one or three numbers" )

  # -----------------------------------------
  # Legend style and parameters
  legend.style <- match.arg(legend.style, c("auto", "tall", "broad", "none"))
  if(legend.style == "auto") {
    if(Nc < 5) {
      legend.style <- "tall"
    } else {
      legend.style <- "broad"
    }
  }

  legend.ysize <- switch(legend.style, broad=1, tall=2, none=.2)


  # ---------------------------------------
  # plotting
  plot.new()
  oldpar <- par(mar=rep(1,4), usr=rep(c(-0.04, 1.04), 2), xpd=NA)
  on.exit(par(oldpar))

  # ---------------------------------------
  # calculate plotting parameters
  #
  # maximum row text height
  line.h <- max(strheight(row.labels, cex=text.cex[1]))
  legend.l.h <- strheight("PpJk", cex=text.cex[3])

  pp <- par("pin")

  # maximum height of the text in the top row
  # with columns, we need to recalculate due to rotation
  colnht <- max(strwidth(col.labels, cex=text.cex[2]))  * pp[1] / pp[2] 
  colnwd <- max(strheight(col.labels, cex=text.cex[2])) * pp[2] / pp[1]

  #print(sprintf("rownwd=%.2f colnwd=%.2f\n", rownwd, colnwd))

  # maximum row text width ("row name width = rownwd")
  rownwd <- max(strwidth(row.labels, cex=text.cex[1])) * 1.1
  if(rownwd > 1) 
    warning("Figure too narrow, text will not fit in; use smaller text.cex")

  col.w <- 1/(Nc) *(1-rownwd-colnwd)
  if(col.w < colnwd) warning( "Figure too narrow, the labels will overlap.\nConsider using smaller text.cex" )
  message(sprintf("rownwd=%.2f, col.w=%.2f, colnwd=%.2f\n", rownwd, col.w, colnwd))

  n.cn <- switch(col.labels.style, top=1, bottom=1, both=2, none=0)
  col.lab.space <- 1 - n.cn * (line.h + colnht)

  row.h <- 1/(Nr + legend.ysize) * (col.lab.space - 5 * legend.ysize * legend.l.h)

  if(row.h < 0 ) warning("Figure too short, text will not fit in; use smaller text.cex")
  if(row.h < line.h) warning( "Figure too short, the labels will overlap.\nConsider using smaller text.cex" )

  # columns tend to get too wide if there are too few
  # if(col.w > 3 * row.h * pp[2] / pp[1]) {
  #   col.w <- 3 * row.h * pp[2] / pp[1]
  # }

  # reserve space at the bottom
  bottom.h <- (legend.l.h * 5 + row.h) * legend.ysize

  if(col.labels.style %in% c("bottom", "both"))  {
    bottom.h <- bottom.h + line.h + colnht
  }

  x.vec <- rownwd + colnwd + ((1:Nc) - 0.5) * col.w
  y.vec <- bottom.h + ((1:Nr) - 0.5) * row.h 
  #abline(h=y.vec[1])
  #abline(h=bottom.h, col="red")

  message(sprintf("bottom.h=%.3f, y.vec[1]=%.3f", bottom.h, y.vec[1]))


  # setup the internal drawing function
  .draw.test <- function(row, col, xc, yc, ev, mq, color) {
    # if a plotting function is provided, we will use it
    if(!is.null(plot.func)) {
      dev.hold()
      on.exit(dev.flush())
      for(i in 1:length(xc)) { plot.func(row=row[i], col=col[i], x=xc[i], y=yc[i], col.w, row.h, e=ev[i], p=mq[i]) }
    } else {
      if(is.null(color)) color <- "#33333333" ;
      cex <- (1 + 2*ev) * plot.cex
      points(xc, yc, col=color, pch=19, cex=cex)
    }
  }

  # plot the labels
  text(rownwd + colnwd * 0.5, rev(y.vec), row.labels, pos=2, cex=text.cex[1] )

  # column labels
  if(col.labels.style %in% c("both", "top")) {
    coln.y.pos <- rep( bottom.h + Nr * row.h + line.h, Nc)
    text( x.vec, coln.y.pos, col.labels, srt=90, adj=c(0, 0.5), cex=text.cex[2])
  }

  if(col.labels.style %in% c("both", "bottom")) {
    coln.y.pos <- rep( bottom.h - line.h, Nc)
    text( x.vec, coln.y.pos, col.labels, srt=90, adj=c(1, 0.5), cex=text.cex[2])
  }


  # ---------------------------------------
  # add light grey grid
  .plotGrid(grid, x.vec, y.vec, row.h, col.w, grid.col=grid.color, scale.max=max.e, scale.min=min.e, text.cex=text.cex[1])

  # ---------------------------------------
  # plot the test results
  mq[mq > max.p] <- max.p
  signif <- mq >= min.p       # not significant points
  if(all(!signif)) { # nothing significant at all
    stop(sprintf("Nothing to do: no p-values below %.2e", 10^-min.p))
  }

  pmin2 <- min(mq[signif])  # lowest p but not above pval

  palfunc <- .getColFunc(pmin2, max(mq), "red")
  col <- palfunc(mq)

  dme <- max.e - min.e
  if(dme == 0) dme <- 1
  ev  <- (me - min.e)/dme

  xc <- x.vec[col(me)]
  yc <- rev(y.vec)[row(me)]
  
  .draw.test(row(me)[signif], col(me)[signif], xc[signif], yc[signif], ev[signif], mq[signif], col[signif])

  lines(c(x.vec[1], x.vec[2]), rep(y.vec[1], 2), col="red", lwd=2)

  # ---------------------------------------
  # add a legend

  if(legend.style != "none") 
    .pvalEffectPlotLegend(.draw.test, palfunc,
      legend.l.h, col.w, row.h,
      min.e, max.e,
      min.p, max.p, 
      text.cex,
      legend.style=legend.style

    ) 

  return(invisible(NULL))
}


#' Plot a summary of multiple tmod analyses
#'
#' Plot a summary of multiple tmod analyses
#'
#' This function is useful if you run an analysis for several conditions or
#' time points and would like to summarize the information on a plot.
#' You can use lapply() to generate a list with tmod results and use
#' tmodPanelPlot to visualize it.
#' 
#' tmodPanelPlot shows a heatmap-like plot. Each row corresponds to one
#' module, and columns correspond to the time points or conditions for which
#' the tmod analyses were run. Each significantly enriched module is shown as
#' a red dot. Size of the dot corresponds to the effect size (for example, AUC
#' in the CERNO test), and intensity of the color corresponds to the q-value.
#'
#' By default, tmodPanelPlot visualizes each the results of a single
#' statistical test by a red dot. However, it is often interesting to know how
#' many of the genes in a module are significantly up- or down regulated.
#' tmodPanelPlot can draw a pie chart based on the optional argument "pie".
#' The argument must be a list of length equal to the length of x.
#' Note also that the names of the pie list must be equal to the names of x.
#' Objects returned by the function tmodDecideTests can be directly used
#' here. The rownames of either the data frame or the array must be the
#' module IDs.
#' 
#' @param x list, in which each element has been generated with a tmod test function
#' @param pie a list of data frames with information for drawing a pie chart
#' @param clust whether, in the resulting data frame, the modules should be
#' ordered by clustering them with either q-values ("qval") or the effect size
#' ("effect"). If "sort" or NULL, the modules are sorted alphabetically by their ID.
#' If "keep", then the order of the modules is kept.
#' @param select a character vector of module IDs to show. If clust == "keep", then in that particular
#' order.
#' @param filter.empty.cols If TRUE, all elements (columns) with no enrichment below pval.thr in any row will be removed
#' @param filter.empty.rows If TRUE, all modules (rows) with no enrichment below pval.thr in any column will be removed
#' @param filter.unknown If TRUE, modules with no annotation will be omitted
#' @param filter.rows.pval Rows in which no p value is below this threshold will be omitted
#' @param filter.rows.auc Rows in which no AUC value is above this threshold will be omitted
#' @param filter.by.id if provided, show only modules with IDs in this character vector
#' @param pval.thr Results with p-value above pval.thr will not be shown
#' @param pval.thr.lower Results with p-value below pval.thr.lower will look identical on the plot
#' @param col.labels Labels for the columns. If NULL, names of the elements
#' of the list x will be used.
#' @param row.labels Labels for the modules. This must be a named vector, with module IDs as vector names. If NULL, module titles from
#' the analyses results will be used.
#' @param row.labels.auto Automatic generation of row labels from module
#' data: "both" (default, ID and title), "id" (only ID), "title" (only title),
#' "none" (no row label)
#' @param pie.colors character vector of length equal to the cardinality of the third dimension of the pie argument. By default: blue, grey and red.
#' @param grid Style of a light-grey grid to be plotted; can be "none", "at" and "between"
#' @param plot.cex a numerical value giving the amount by which the plot
#' symbols will be maginfied 
#' @param text.cex a numerical value giving the amount by which the plot
#' text will be magnified, or a vector containing three cex values for row labels, column labels and legend, respectively
#' @param plot.func Optionally, a function to be used to draw the dots. See "pvalEffectPlot"
#' @param pie.style Can be "pie", "boxpie" or "rug"
#' @param col.labels.style Style of column names: "top" (default), "bottom", "both", "none"
#' @param legend.style Style of the legend: "auto" -- automatic; "broad":
#' pval legend side by side with effect size legend; "tall": effect size
#' legend above pval legend
#' @param min.e,max.e scale limits for the effect size (default: 0.5 and 1.0)
#' @param ... Any further arguments will be passed to the pvalEffectPlot function (for example, grid.color)
#' @return a data frame with a line for each module encountered anywhere in
#' the list x, two columns describing the module (ID and module title), and
#' two columns(effect size and q value) for each element of list x.
#' @seealso tmodDecideTests, tmodSummary, pvalEffectPlot, simplePie
#' @examples
#' data(Egambia)
#' E <- Egambia[,-c(1:3)]
#' pca <- prcomp(t(E), scale.=TRUE)
#'
#' # Calculate enrichment for first 5 PCs
#' gs   <- Egambia$GENE_SYMBOL
#' gn.f <- function(r) {
#'     o <- order(abs(r), decreasing=TRUE)
#'     tmodCERNOtest(gs[o], 
#'                 qval=0.01)
#' }
#' x <- apply(pca$rotation[,3:4], 2, gn.f)
#' tmodPanelPlot(x, text.cex=0.7)
#' @export
tmodPanelPlot <- function(x, pie=NULL, clust="qval", select=NULL,
  filter.empty.cols=FALSE, filter.empty.rows=TRUE, filter.unknown=TRUE,
  filter.rows.pval=0.05,
  filter.rows.auc=0.5,
  filter.by.id=NULL,
  col.labels=NULL, 
  col.labels.style="top",
  row.labels=NULL, 
  row.labels.auto="both",
  pval.thr=10^-2,
  pval.thr.lower=10^-6,
  plot.func=NULL,
  grid="at", 
  pie.colors=c("#0000FF", "#cccccc", "#FF0000" ),
  plot.cex=1, text.cex=1, 
  pie.style="rug", 
  min.e=.5, max.e=1,
  legend.style="auto", ... ) {

  x <- .xcheck(x)

  # check whether the pie object is correct
  if(!is.null(pie)) pie <- .piecheck(pie, x)

  # create a summary
  df <- tmodSummary(x, clust=clust, filter.empty=FALSE, filter.unknown=filter.unknown, select=select)

  if(!is.null(filter.by.id)) 
    df <- df[ df$ID %in% filter.by.id, , drop=FALSE ]


  # check the column labels
  if(is.null(col.labels)) col.labels <- names(x)
  if(is.null(col.labels)) col.labels <- paste0("X.", 1:length(x))

  # handle row labels
  row.labels.auto <- match.arg(row.labels.auto, c("title", "id", "both", "none"))
  if(is.null(row.labels)) {
    row.labels <- switch(row.labels.auto,
      title=df$Title,
      id=df$ID,
      both=sprintf("%s (%s)", df$Title, df$ID),
      none=rep("", nrow(df)))
  } else {
    if(!all(df$ID %in% names(row.labels)))
      stop("row.labels must be a named vector with all module IDs")
    row.labels <- row.labels[ df$ID ]
  }

  min.p <- -log10(pval.thr)
  max.p <- -log10(pval.thr.lower)

  # calculate q and e matrices
  m <- as.matrix(df[,-c(1:2),drop=FALSE])
  Nc <- ncol(m)/2
  Nr <- nrow(m)

  # split the matrix in two
  me <- m[,(1:Nc)*2-1,drop=FALSE]
  me[is.na(me)] <- 0.5

  mq <- m[,(1:Nc)*2,drop=FALSE]
  mq[is.na(mq)] <- 1
  mq <- -log10(mq)

  # ----- FILTERING ------------
  row.ids <- df$ID

  # remove rows with no at least one pval below filter.rows.pval
  minps <- apply(mq, 1, max)
  sel   <- minps > -log10(filter.rows.pval)

  # remove rows with at least one AUC above filter.rows.auc
  maxaucs <- apply(me, 1, max)
  sel     <- sel & maxaucs >= filter.rows.auc

  if(sum(!sel) == nrow(mq)) stop("No rows remain after filtering")

  mq    <- mq[sel,,drop=F]
  me    <- me[sel,,drop=F]
  df    <- df[sel,,drop=F]
  row.labels <- row.labels[sel]
  row.ids    <- row.ids[sel]


  # remove rows w/o a significant p-value
  if(filter.empty.rows) {
    minps <- apply(mq, 1, max)
    sel   <- minps > min.p
    mq <- mq[sel,,drop=F]
    me <- me[sel,,drop=F]
    row.labels <- row.labels[sel]
    row.ids    <- row.ids[sel]
  }

  # remove cols w/o a significant p-value
  if(filter.empty.cols) {
    minps <- apply(mq, 2, max)
    sel   <- minps > min.p
    mq <- mq[,sel,drop=F]
    me <- me[,sel,drop=F]
    col.labels <- col.labels[sel]
  }

  # ------------ PLOTTING ----------------

  # prepare the pie plotting function
  if(!is.null(pie)) {
    plot.func <- .preparePlotFunc(row.ids, 
      pie, min.p, max.p, pie.colors, plot.cex, style=pie.style)
  }

  # make the actual plot -- pass it on to pvalEffectPlot
  pvalEffectPlot(me, 10^-mq, 
    row.labels=row.labels, col.labels=col.labels, pval.thr=pval.thr,
    grid=grid, plot.cex=plot.cex, text.cex=text.cex, 
    plot.func=plot.func, legend.style=legend.style, col.labels.style=col.labels.style,
    min.e=min.e, max.e=max.e,
    ...)
  return(invisible(df))
}


## select a function which is used to recalculate the bounding box and its
## position, depending on the pie style
.selCalcXYWH <- function(style="pie") {

  ret <- NULL

  if(style %in% c( "pie", "boxpie")) {
    # for pie and boxpie, the size of the plotting widget is proportional to
    # the effect size, scaled in vertical as well as horizontal, centered
    ret <- function(x, y, w, h, e, plot.cex) {
      w <- w * plot.cex * (0.5 + e/2)
      h <- h * plot.cex * (0.5 + e/2)
      pp <- par("pin")
      w2 <- w * pp[1] / pp[2]
      if(w2 < h) {
        h <- w2
      } else {
        w <- h * pp[2] / pp[1]
      }
      return(list(x=x, y=y, w=w, h=h))
    }
  } else if(style == "rug") {
    # rug takes up all the available vertical space, but is scaled
    # horizontally to reflect the effect size, left-adjusted
    ret <- function(x, y, w, h, e, plot.cex) {
      x <- x - w /2
      w <- w * (0.5 + e )/1.5
      x <- x + w / 2
      h <- h * 0.9
      return(list(x=x, y=y, w=w, h=h))
    }
  }

  return(ret) 
}


## choose the right function for plotting
.selPiePlotFunc <- function(style) {

  switch(style,
    pie=wgPie,
    boxpie=wgBoxpie,
    rug=wgRug)

}


## return a plotting function for use as parameter in pvalEffectPlot
.preparePlotFunc  <- function(row.ids, pie, min.p, max.p, pie.colors, plot.cex, style="pie") {

  style <- match.arg(style, c("pie", "rug", "boxpie"))

  # prepare the pie coloring function
  col.funcs <- lapply(pie.colors,
    function(xx) .getColFunc( min.p, max.p, xx, alpha="FF" )
    )

  # different pie styles call for different bounding boxes, positions and
  # plotting functions. Here is the right place to decide which to take.
  calcXYWH.func <- .selCalcXYWH(style)
  piePlot.func  <- .selPiePlotFunc(style)

  # here we generate the plotting function that will be called both for
  # drawing the legend and drawing the pie charts on the panel plot
  plot.func <- function(row, col, x, y, w, h, e, p) {

    # this is where we take the additional arguments from the "pie" object
    # and turn it into arguments for the pie widget plotting function
    if(is.null(row)) {
      v <- c(10, 10, 10) # row is null when called from legend!
    } else {
      id <- row.ids[row]
      v <- pie[[col]][id,]
    }

    if(sum(v) == 0) {
      warning(sprintf("For module %s, condition %d, sum of pie data is zero", row.ids[row], col))
      return()
    }

    mycols <- sapply(col.funcs, function(x) x(p))

    xywh <- calcXYWH.func(x, y, w, h, e, plot.cex)

    with(xywh, piePlot.func(x, y, w, h, v=v, col=mycols))

  }

  return(plot.func)
}


## ---------------- argument checking for tmodPanelPlot ---------------------------


## check whether the "pie" argument is correct
.piecheck <- function(pie, x) {
  if(!is(pie, "list")) 
    stop( "pie must be a list. Make sure that is.list(pie) == TRUE")

  if(is.null(names(pie))) {
    warning("names(pie) is NULL. Generating default names")
    names(pie) <- paste0("X.", 1:length(pie))
  }

  if(! all(names(x) %in% names(pie))) 
    stop("All named elements of x must be found in pie. Please make sure that all(names(x) %in% names(pie))")

  pie <- pie[names(x)]

  pie <- sapply(pie, as.matrix, simplify=FALSE)
  pie
}


## check whether the x argument is correct
.xcheck <- function(x) {
  if(!is(x, "list")) stop( "x must be a list object. Make sure that is.list(x) == TRUE")

  z <- sapply(x, is.data.frame)
  if(any(!z)) {
    warning(sprintf("Some elements of x are not data frames, removing %d elements", sum(!z)))
    x <- x[ z ]
  }

  if(length(x) < 1)
    stop("No usable elements of x")

  if(is.null(names(x))) {
    warning("names(x) is NULL. Generating default names")
    names(x) <- paste0("X.", 1:length(x))
  }

  x
}
