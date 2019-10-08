
## adds an alpha chanel to a color
.alphacol <- function(x, suff="99") {
  x <- col2rgb(x)
  x <- apply(x, 2, function(xx) paste0("#", paste(as.hexmode(xx), collapse= ""), suff) )
  x
}


## This function implements a friendly palette to be used in graphics
.mypalette <- function (n = NULL, transparent = "99") {


    pal="E69F00 56B4E9 009E73 F0E442 0072B2 D55E00 CC79A7 669900 660099 996600 990066 666633 666600 aa3366 5B4E85 FF6A5C ADAEA3 A0A376 FF8040 A2D6DA DA9CA5"

    pal <- unlist(strsplit(pal, " "))
    pal <- paste("#", pal, transparent, sep = "")
    if (!is.null(n)) {
        if (n > length(pal)) {
            pal <- rep(pal, ceiling(n/length(pal)))
        } else {
            pal <- pal[1:n]
        }
    }
    return(pal)
}

## return a gradient palette
.gradientpal <- function(n=NULL, set="bluewhitered", transparent="99") {
  if(is.null(n)) n <- 3

  cols <- switch(set,
    bwr=c("blue", "white", "red"),
    rwb=c("red", "white", "blue"),
    ckp=c("cyan", "black", "purple"),
    rwb=c("purple", "black", "cyan")
    )

  pal <- colorRampPalette(cols)(n)
  paste0(pal, transparent)
}


#' A selection of color palettes
#'
#' Return a preset selection of colors, adjusted by alpha
#'
#' A few palettes have been predefined in tmod, and this function can be
#' used to extract them. The following palettes have been defined:
#' * friendly -- a set of distinct, colorblind-friendly colors
#' * bwr, rwb, ckp, pkc -- gradients (b-blue, r-red, w-white, c-cyan, k-blacK, p-purple)
#' By default, either all colors are returned, or, if it is a gradient
#' palette, only three.
#'
#' @param n Number of colors to return (default: all for "friendly", 3 for everything else)
#' @param set Which palette set (see Details).
#' @param alpha 0 for maximum transparency, 1 for no transparency.
#' @param func if TRUE, the returned object will be a function rather than
#'        a character vector
#' @return Either a character vector, or, when the func parameter is TRUE,
#'         a function that takes only one argument (a single number)
#' @export
tmodPal <- function(n=NULL, set="friendly", alpha=0.7, func=FALSE) {

  if(alpha > 1) alpha <- 1
  if(alpha < 0) alpha <- 0

  set <- match.arg(set,
    c("friendly", "bwr", "rwb", "ckp", "pkc"))

  transparent <- sprintf("%X", as.integer(alpha * 255))

  if(func) {
    if(set == "friendly") pal <- function(n) .mypalette(n, transparent)
    else pal <- function(n) .gradientpal(n, set, transparent)
  } else {
    if(set == "friendly") pal <- .mypalette(n, transparent)
    else pal <- .gradientpal(n, set, transparent)
  }

  return(pal)
}


#' Select genes belonging to a module from a data frame
#'
#' Select genes belonging to a module from a data frame or vector
#' 
#' showModule filters a data frame or a vector such that only genes from a module are
#' shown. Use it, for example, to show a subset of topTable result from
#' limma in order to see which genes from a module are significantly
#' regulated. In essence, this is just a wrapper around "subset()".
#' @return Either a filtered vector or (if extra==TRUE) a data frame.
#' @param x a data frame or a vector
#' @param genes a character vector with gene IDs
#' @param module a single character value, ID of the module to be shown
#' @param mset Module set to use; see "tmodUtest" for details
#' @param extra If TRUE, additional information about the features will be shown
#' @export
showModule <- function(x, genes, module, mset="all", extra=TRUE) {
  mset <- .getmodules2(NULL, mset)

  if(!is(x, "data.frame") && !is(x, "vector") && !is(x, "factor")) 
    stop( "x must be a data frame, a factor or a vector" )

  is.df <- is(x, "data.frame")

  if(!module %in% mset$MODULES$ID) stop("no such module")
  genes <- as.character(genes)

  if(is.df && length(genes) != nrow(x)) 
    stop("genes must be of length equal to the number of rows of x")

  if(!is.df && length(genes) != length(x))
    stop("genes must be of length equal to the length of x")

  sel <- genes %in% mset$MODULES2GENES[[module]]

  if(sum(sel) == 0) warning("No genes belonging to that module found")

  x   <- subset(x, sel)
  genes <- genes[sel]

  if(extra) x <- cbind(x, mset$GENES[genes,])

  x
}


#' A combined beeswarm / boxplot
#'
#' A combined beeswarm / boxplot
#'
#' This is just a simple wrapper around the beeswarm() and boxplot()
#' commands.
#' @param data a vector of numeric values to be plotted
#' @param group factor describing the groups
#' @param main title of the plot
#' @param pch character to plot the points
#' @param xlab,ylab x and y axis labels
#' @param las see par()
#' @param pwcol colors of the points (see beeswarm)
#' @param ... any additional parameters to be passed to the beeswarm command
#' @examples
#' data(Egambia)
#' E <- as.matrix(Egambia[,-c(1:3)])
#' showGene(E["20799",], rep(c("CTRL", "TB"), each=15))
#' @importFrom beeswarm beeswarm
#' @export
showGene <- function( data, group, main= "", pch= 19, 
                         xlab= "", ylab= "log2 expression", las= 2, pwcol= NULL, ... ) {
  group  <- factor( group )
  pal    <- .mypalette( n= length( unique( group ) ) )

  if(! is.null(pwcol) & length( pwcol ) == 1 ) pwcol <- rep( pwcol, length( group ) )
  if(  is.null(pwcol)) pwcol= pal[ group ]

  
  beeswarm( data ~ group, 
    pch= pch, xlab= xlab, ylab= ylab, main= main, las= las, 
    pwcol= pwcol, bty="n",
    ... )

  boxplot( data ~ group, col= "#ffffff00", add= T, yaxt= "n", xaxt= "n", main= "", outline= FALSE, frame= FALSE )

  return( invisible( list( groups= levels( group ), col= pal ) ) )
}


#' Get genes belonging to a module
#'
#' Get genes belonging to a module
#'
#' Create a data frame mapping each module to a comma separated list of
#' genes. If genelist is provided, then only genes in that list will be
#' shown. An optional column, "fg" informs which genes are in the "foreground"
#' data set.
#' @return data frame containing module to gene mapping, or a list (if
#' as.list == TRUE
#' @param modules module IDs 
#' @param genelist list of genes 
#' @param mset module set to use
#' @param fg genes which are in the foreground set
#' @param as.list should a list of genes rather than a data frame be returned
#' @export
getGenes <- function(modules, genelist=NULL, fg=NULL, mset="LI", as.list=FALSE) {
  mset <- .getmodules2(modules, mset)

  if(!is.null(genelist)) 
    mset$MODULES2GENES <- lapply(mset$MODULES2GENES, function(x) x[ x %in% genelist ])

  if(as.list) return(mset$MODULES2GENES)

  ret <- data.frame(ID=mset$MODULES$ID)
  rownames(ret) <- ret$ID
  ret$N <- sapply(mset$MODULES2GENES, length)
  ret$Genes <- sapply(mset$MODULES2GENES, function(x) paste(x, collapse=","))

  if(!is.null(fg)) {
    ret$fg <- sapply(mset$MODULES2GENES, function(x) paste(x[x %in% fg], collapse=","))
  }
  ret
}




#' Create a visualisation of enrichment
#'
#' Create a visualisation of enrichment
#' 
#' This functions creates a barplot visualizing the enrichment of a
#' module in the foreground (fg) set as compared to the background (bg) set.
#' It is the counterpart 
#' @param fg the foreground set of genes
#' @param bg the background set of genes (gene universe)
#' @param m module for which the plot should be created
#' @param mset Which module set to use (see tmodUtest for details)
#' @param ... additional parameters to be passed to the plotting function
#' @seealso \code{\link{tmod-package}}, \code{\link{evidencePlot}}
#' @examples 
#' set.seed(123)
#' data(tmod)
#' bg <- tmod$GENES$ID
#' fg <- sample(c(tmod$MODULES2GENES[["LI.M127"]], bg[1:1000]))
#' hgEnrichmentPlot(fg, bg, "LI.M127")
#' @export
hgEnrichmentPlot <- function(fg, bg, m, mset="all", ...) {

  mset <- .getmodules2(NULL, mset)

  if(is.null(mset$MODULES2GENES[[m]])) stop("No such module")

  fg <- unique(fg)
  bg <- unique(c(fg, bg))
  b <- sum(fg %in% mset$MODULES2GENES[[m]])
  n <- length(fg)
  B <- sum(bg %in% mset$MODULES2GENES[[m]])
  N <- length(bg)

  mm <- matrix(c(b/n, 1-b/n, B/N, 1-B/N), byrow= F, nrow= 2)

  barplot(mm,
  xlim=c(0, 5), legend.text=c("In module", "Other"), bty="n", names.arg= c(
  "Foreground", "Background"), col= c("#E69F0099",  "#56B4E999"), ...)

    
}



#' Create an evidence plot for a module
#'
#' Create an evidence plot for a module
#'
#' This function creates an evidence plot for a module, based on an
#' ordered list of genes. The plot shows the receiving operator
#' characteristic (ROC) curve and a rug below, which indicates the distribution of the
#' module genes in the sorted list.
#' @param l sorted list of HGNC gene identifiers
#' @param m character vector of modules for which the plot should be created
#' @param mset Which module set to use (see tmodUtest for details)
#' @param scaled if TRUE, the cumulative sums will be divided by the total sum (default)
#' @param rug if TRUE, draw a rug-plot beneath the ROC curve
#' @param roc if TRUE, draw a ROC curve above the rug-plot
#' @param filter if TRUE, genes not defined in the module set will be removed
#' @param unique if TRUE, duplicates will be removed
#' @param add if TRUE, the plot will be added to the existing plot
#' @param gene.labels if TRUE, gene names are shown; alternatively, a named character vector with gene labels to be shown, or NULL (default) for no labels (option evaluated only if rug is plotted)
#' @param gene.colors NULL (default) or a character vectors indicating the color for each gene. Either a named vector or a vector with the same order of genes as `l`.
#' @param gene.lines a number or a vector of numbers; line width for marking the genes on the rug (default=1). If the vector is named, the names should be gene ids.
#' @param gl.cex Text cex (magnification) for gene labels
#' @param style "roc" for receiver-operator characteristic curve (default), and "gsea" for GSEA-style (Kaplan-Meier like plot)
#' @param col a character vector color to be used
#' @param col.rug a character value specifying the color of the rug
#' @param lwd line width (see par())
#' @param lty line type (see par())
#' @param rug.size fraction of the plot that should show the rug. If rug.size is 0, rug is not drawn. If rug.size is 1, ROC curve is not drawn.
#' @param legend position of the legend. If NULL, no legend will be drawn
#' @param ... Further parameters passed to the plotting function
#' @seealso \code{\link{tmod-package}}, \code{\link{hgEnrichmentPlot}}
#' @examples 
#' # artificially enriched list of genes
#' set.seed(123)
#' data(tmod)
#' bg <- tmod$GENES$ID
#' fg <- sample(c(tmod$MODULES2GENES[["LI.M127"]], bg[1:1000]))
#' l <- unique(c(fg, bg))
#' evidencePlot(l, "LI.M127")
#' evidencePlot(l, filter=TRUE, "LI.M127")
#' @export
evidencePlot <- function(l, m, mset="all", scaled= TRUE, rug=TRUE, roc=TRUE,
  filter= FALSE, unique=TRUE, add= FALSE, col="black", 
  col.rug="#eeeeee",
  gene.labels=NULL, 
  gene.colors=NULL,
  gene.lines=1,
  gl.cex=1,
  style="roc",
  lwd=1, lty=1, rug.size=0.2,
  legend=NULL, ...) {

  if(rug.size == 0) rug <- FALSE
  if(rug.size == 1) roc <- FALSE
  if(!rug && !roc) stop("Both rug and roc are FALSE, nothing to plot")
  if(!roc) rug.size <- 1
  if(!rug) rug.size <- 0

  style <- match.arg(style, c("roc", "gsea"))

  # standardize the modules
  m <- as.character(m)
  mset <- .getmodules2(NULL, mset)

  if(! all(m %in% names(mset$MODULES2GENES))) stop("No such module")

  #if(filter) l <- l[ l %in% mset$GENES$ID ]
  keep <- rep(TRUE, length(l))
  if(filter) keep[ ! l %in% mset$GENES$ID ] <- FALSE
  #if(unique) l <- unique(l)
  if(unique) keep[ duplicated(l) ] <- FALSE
  l <- l[ keep ]

  # gene colors
  if(!is.null(gene.colors)) {
    if(is.null(names(gene.colors))) {
      gene.colors <- gene.colors[ keep ]
      stopifnot(length(gene.colors) == length(l))
      names(gene.colors) <- l
    }
  }

  # gene line widths

  if(is.null(names(gene.lines))) {
    if(length(gene.lines) == 1) {
      gene.lines <- rep(gene.lines, length(l))
      names(gene.lines) <- l
    } else if(length(gene.lines) != length(keep)) {
      stop("if gene.lines is shorter than l, then it must be named")
    }
    gene.lines <- gene.lines[keep]
    names(gene.lines) <- l
  }


  n <- length(l)
  Nm <- length(m)

  # check graphical parameters and fill up if necessary
  col <- col[(((0:(Nm-1)) ) %% length(col))+1]
  lty <- lty[(((0:(Nm-1)) ) %% length(lty))+1]
  lwd <- lwd[(((0:(Nm-1)) ) %% length(lwd))+1]

  # check gene labels
  if(!is.null(gene.labels)) {

    if(is.logical(gene.labels) && gene.labels) {
      sel <- unlist(unique(lapply(m, function(mm) l[ l %in% mset$MODULES2GENES[[mm]] ])))
      gene.labels <- sel
      names(gene.labels) <- sel
    }

    gene.labels <- unique(as.character(gene.labels))
    
    if(is.null(names(gene.labels))) {
      names(gene.labels) <- gene.labels
    }

  }

  gene.labels <- gene.labels[ names(gene.labels) %in% l ]
  gene.labels <- gene.labels[ order(match(names(gene.labels), l)) ]

  # for each module to draw, find out which genes from l are found in that
  # module
  x <- sapply(m, function(mm) l %in% mset$MODULES2GENES[[mm]])

  # cumulative sum or scaled cumulative sum
  if(roc) {
    if(scaled) {
      cfunc <- function(xx) cumsum(xx) / sum(xx)
      r <- c(0, 1)
      ylab="Fraction of genes in module"
    } else {
      cfunc <- function(xx) cumsum(xx)
      r <- range(xcs)
      ylab <- "Number of genes in module"
    }

    cfunc2 <- switch(style,
      roc=cfunc,
      gsea=function(xx) { cfunc(xx) - cfunc(!xx) }
      )

    xcs <- apply(x, 2, cfunc2)

    # additional space for the rug
    r[1] <- - rug.size * (r[2]-r[1]) / (1-rug.size)
  } else {
    ylab <- ""
    r <- c(-1, 0)
  }

  if(!add) {
    plot(NULL, xlim= c(1, n), ylim= r, bty="n", xlab="List of genes", ylab=ylab, yaxt="n", ...)
    if(roc) axis(side=2, at= axisTicks(range(xcs), log=FALSE))
  }

  if(roc) {
    # plot the actual ROC curves
    for(i in 1:Nm) { lines(1:n, xcs[,i], col= col[i], lty=lty[i], lwd=lwd[i]) }

    # draw the diagonal
    if(style == "roc") segments(0, 0, n, r[2], col= "grey")
  }

  # plot the rug(s)
  if(rug) { 
    step <- r[1] / Nm

    # draw the rug background
    rect(0, (1:Nm) * step - step * 0.8, n, (1:Nm) * step, col= "#eeeeee", border= NA) 

    # draw each individual rug
    for(i in 1:Nm) {
      w <- which(x[,i]) # which genes to show for module m[i]
      gi <- l[w] # gene ids to show
      print(gi)

      if(is.null(gene.colors)) {
        cols <- col[i]
      } else {
        cols <- ifelse(is.na(gene.colors[gi]), col[i], gene.colors[gi])
      }

      lwds <- ifelse(is.na(gene.lines[gi]), lwd[i], gene.lines[gi])

      print(cols)
      segments(w, step * i, w, step * (i-0.8), col= cols, lty=lty[i], lwd=lwds)
    }

    # add gene labels
    if(!is.null(gene.labels)) {
      lw <- max(strheight(gene.labels, units="i", cex=gl.cex))
      lw <- grconvertX( c(0, lw), from="i", to="u")
      lw <- lw[2] - lw[1]
      lpos <- lw * 1.5 * (1:length(gene.labels))
      gi <- names(gene.labels)
      rpos <- match(gi, l)

      if(is.null(gene.colors)) {
        cols <- "black"
      } else {
        cols <- ifelse(is.na(gene.colors[gi]), "black", gene.colors[gi])
      }

      lwds <- ifelse(is.na(gene.lines[gi]), lwd[1], gene.lines[gi])


      for(i in 1:length(lpos)) {
        if(rpos[i] + lw > lpos[i]) {
          sel <- i:length(lpos)
          lpos[sel] <- lpos[sel] + (rpos[i] - lpos[i] + lw)
        }
      }

      text(lpos, 0.1, gene.labels, srt=90, pos=4, col=cols)
      segments(rpos, step * 0.2, lpos + lw/2, 0.08, col=cols, lwd=lwds)
    }
  }

  # add a legend
  if(!is.null(legend)) {
    legend(legend, m, col=col, lty=lty, lwd=lwd, bty= "n")
  }
}



#' Plot a PCA object returned by prcomp
#'
#' Plot a PCA object returned by prcomp
#' @param pca PCA object returned by prcomp
#' @param components a vector of length two indicating the components to plot
#' @param pch Type of character to plot (default: 19)
#' @param col Color for plotting (default: grey)
#' @param group a factor determining shapes of the points to show (unless
#'        overriden by pch=...)
#' @param ... any further parameters will be passed to the plot() function
#'        (e.g. col, cex, ...)
#' @return If group is NULL, then NULL; else a data frame containing
#'         colors and shapes matching each group
#' @export
pcaplot <- function(pca, components=1:2, group=NULL, col="black", pch=19, ...) {

  args <- list(...)

  x <- pca$x[, components[1]]
  y <- pca$x[, components[2]]

  if(!is.null(group)) {
    group <- factor(group)
    gr <- as.numeric(factor(group))
    pch <- c(15:19)[ gr %% 5L ]
    if(length(col) == 1L) col <- rep(col, length(x))
  }

  default.args <- list(pch=pch, xlim=range(x), ylim=range(y), bty="n", pch=pch, 
    col=col,
    xlab=paste("PC", components[1]), ylab=paste("PC", components[2]))
  args <- c(args, default.args)
  args <- c(list(x, y), args[unique(names(args))])
  do.call(plot, args)

  abline(h=0, col="grey")
  abline(v=0, col="grey")

  if(is.null(group)) return(NULL)

  sel <- !duplicated(group)

  return(list(groups=group[sel], pch=args$pch[sel], colors=args$col[sel]))
}
