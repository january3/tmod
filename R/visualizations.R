
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





#' Create a visualisation of enrichment
#'
#' Create a visualisation of enrichment
#' 
#' This functions creates a barplot visualizing the enrichment of a
#' module in the foreground (fg) set as compared to the background (bg) set.
#' It is the counterpart 
#' @param fg the foreground set of genes
#' @param bg the background set of genes (gene universe)
#' @param m gene set for which the plot should be created
#' @param mset Which module set to use (see tmodUtest for details)
#' @param ... additional parameters to be passed to the plotting function
#' @seealso \code{\link{tmod-package}}, \code{\link{evidencePlot}}
#' @examples 
#' set.seed(123)
#' data(tmod)
#' bg <- tmod$gv
#' fg <- getGenes("LI.M127", as.list=TRUE)[[1]]
#' fg <- sample(c(fg, bg[1:100]))
#' hgEnrichmentPlot(fg, bg, "LI.M127")
#' @export
hgEnrichmentPlot <- function(fg, bg, m, mset="all", ...) {

  mset <- .getmodules_gs(NULL, mset)

  if(!m %in% mset$gs$ID) {
    stop("No such gene set")
  }

  m <- match(m, mset$gs$ID)

  fg <- unique(fg)
  bg <- unique(c(fg, bg))

  fg <- .prep_list(fg, mset=mset, filter=FALSE, nodups=FALSE)
  bg <- .prep_list(bg, mset=mset, filter=FALSE, nodups=FALSE)


  b <- sum(fg %in% mset$gs2gv[[m]])
  n <- length(fg)
  B <- sum(bg %in% mset$gs2gv[[m]])
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
#' ordered list of genes. By default, the plot shows the receiving operator
#' characteristic (ROC) curve and a rug below, which indicates the distribution of the
#' module genes in the sorted list.
#'
#' Several styles of the evidence plot are possible:
#'  * roc (default): a receiver-operator characteristic like curve; the
#'    area under the curve corresponds to the effect size (AUC)
#'  * roc_absolute: same as above, but the values are not scaled by the
#'    total number of genes in a module
#'  * gsea
#'  * enrichment: the curve shows relative enrichment at the given position
#'
#' 
#' @param l sorted list of HGNC gene identifiers
#' @param m character vector of modules for which the plot should be created
#' @param mset Which module set to use (see tmodUtest for details)
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
#' bg <- sample(tmod$gv)
#' fg <- getGenes("LI.M127", as.list=TRUE)[[1]]
#' fg <- sample(c(fg, bg[1:1000]))
#' l <- unique(c(fg, bg))
#' evidencePlot(l, "LI.M127")
#' evidencePlot(l, filter=TRUE, "LI.M127")
#' @export
evidencePlot <- function(l, m, mset="all", rug=TRUE, roc=TRUE,
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

  style <- match.arg(style, c("roc", "roc_absolute", "gsea", "enrichment"))

  # standardize the modules
  m <- as.character(m)
  mset <- .getmodules_gs(NULL, mset)

  if(! all(m %in% tmod_ids(mset))) stop("No such gene set")
  m_orig <- m
  m <- match(m, tmod_ids(mset))

  # #if(filter) l <- l[ l %in% mset$GENES$ID ]
  # keep <- rep(TRUE, length(l))
  # if(filter) keep[ ! l %in% mset$GENES$ID ] <- FALSE
  # #if(unique) l <- unique(l)
  # if(unique) keep[ duplicated(l) ] <- FALSE
  # l <- l[ keep ]
  l_orig <- l
  tmp    <- .prep_list(l, x=l_orig, mset=mset, filter=filter, nodups=unique)
  keep <- l_orig %in% tmp$x
  l_orig <- tmp$x
  l      <- tmp$l

  # gene colors
  if(!is.null(gene.colors)) {
    if(is.null(names(gene.colors)) && length(gene.colors) == 1) {
      gene.colors <- rep(gene.colors, length(l_orig))
      names(gene.colors) <- l_orig
    } else {
      if(is.null(names(gene.colors))) {
        gene.colors <- gene.colors[ keep ]
        stopifnot(length(gene.colors) == length(l))
        names(gene.colors) <- l_orig
      }
    }
  }

  # gene line widths

  if(is.null(names(gene.lines))) {
    if(length(gene.lines) == 1) {
      gene.lines <- rep(gene.lines, length(l))
      names(gene.lines) <- l_orig
    } else if(length(gene.lines) != length(keep)) {
      stop("if gene.lines is shorter than l, then it must be named")
    }
    gene.lines <- gene.lines[keep]
    names(gene.lines) <- l_orig
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
      sel <- unlist(unique(lapply(m, function(mm) l[ l %in% mset$gs2gv[[mm]] ])))
      gene.labels <- mset$gv[ sel ]
      names(gene.labels) <- gene.labels
    }

    #gene.labels <- unique(as.character(gene.labels))
    
    if(is.null(names(gene.labels))) {
      names(gene.labels) <- gene.labels
    }

  }

  gene.labels <- gene.labels[ names(gene.labels) %in% l_orig ]
  gene.labels <- gene.labels[ !duplicated(names(gene.labels)) ]
  gene.labels <- gene.labels[ order(match(names(gene.labels), l_orig)) ]

  # for each module to draw, find out which genes from l are found in that
  # module
  x <- sapply(m, function(mm) l %in% mset$gs2gv[[mm]])

  # cumulative sum or scaled cumulative sum
  if(roc) {
    cfunc2 <- switch(style,
      roc=function(xx) cumsum(xx)/sum(xx),
      gsea=function(xx) { cumsum(xx)/sum(xx) - cumsum(!xx)/sum(!xx) },
      roc_absolute=cumsum,
      enrichment=function(xx) (cumsum(xx)/(1:length(xx)))/(sum(xx)/length(xx))
      )

    xcs <- apply(x, 2, cfunc2)

    ylab <- switch(style,
      roc="Fraction of genes in module",
      roc_absolute="Number of genes in module",
      gsea="Relative enrichment",
      enrichment="Enrichment")

    r <- switch(style,
      roc=c(0, 1),
      roc_absolute=c(0, max(apply(x, 2, sum))),
      gsea=c(0, max(xcs)),
      enrichment=c(0, max(xcs)))

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
      gi <- l_orig[w] # gene ids to show

      if(is.null(gene.colors)) {
        cols <- col[i]
      } else {
        cols <- ifelse(is.na(gene.colors[gi]), col[i], gene.colors[gi])
      }

      lwds <- ifelse(is.na(gene.lines[gi]), lwd[i], gene.lines[gi])

      segments(w, step * i, w, step * (i-0.8), col= cols, lty=lty[i], lwd=lwds)
    }

    # add gene labels
    if(!is.null(gene.labels)) {
      lw <- max(strheight(gene.labels, units="i", cex=gl.cex))
      lw <- grconvertX( c(0, lw), from="i", to="u")
      lw <- lw[2] - lw[1]
      lpos <- lw * 1.5 * (1:length(gene.labels))
      gi <- names(gene.labels)
      rpos <- match(gi, l_orig)

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
    legend(legend, m_orig, col=col, lty=lty, lwd=lwd, bty= "n")
  }
}



#' Plot a PCA object returned by prcomp
#'
#' Plot a PCA object returned by prcomp
#'
#' This is a simplistic function. A much better way is to use the pca2d
#' function from the pca3d package.
#' @param pca PCA object returned by prcomp
#' @param components a vector of length two indicating the components to plot
#' @param pch Type of character to plot (default: 19)
#' @param col Color for plotting (default: internal palette)
#' @param group a factor determining shapes of the points to show (unless
#'        overriden by pch=...)
#' @param legend draw a legend? If legend is a position (eg. "topright"), then a legend
#'        is drawn. If NULL or if the group parameter is NULL, then not.
#' @param cex size of the symbols used for plotting
#' @param ... any further parameters will be passed to the plot() function
#'        (e.g. col, cex, ...)
#' @return If group is NULL, then NULL; else a data frame containing
#'         colors and shapes matching each group
#' @export
pcaplot <- function(pca, components=1:2, group=NULL, col=NULL, pch=19, cex=2, legend=NULL, ...) {

  args <- list(...)
  

  x <- pca$x[, components[1]]
  y <- pca$x[, components[2]]

  if(!is.null(group)) {
    group <- factor(group)
    gr <- as.numeric(factor(group))
    pch <- c(15:19)[ gr %% 5L ]
    if(is.null(col)) col <- .mypalette(n=length(levels(group)))
    col <- col[ gr ]
    if(length(col) == 1L) col <- rep(col, length(x))
  }

  if(is.null(col)) col <- "grey"

  default.args <- list(pch=pch, xlim=range(x), ylim=range(y), bty="n", pch=pch, 
    cex=cex,
    col=col,
    xlab=paste("PC", components[1]), ylab=paste("PC", components[2]))
  args <- c(args, default.args)
  args <- c(list(x, y), args[unique(names(args))])
  do.call(plot, args)

  abline(h=0, col="grey")
  abline(v=0, col="grey")

  if(!is.null(legend) && !is.null(group)) {
      n <- length(levels(group))
      mm <- match(levels(group), group)
      legend(legend, levels(group), pch=pch[mm], col=col[mm], bty="n")
  }
  if(is.null(group)) return(NULL)

  sel <- !duplicated(group)

  return(list(groups=group[sel], pch=args$pch[sel], colors=args$col[sel]))
}


.auto_labels <- function(labels, modules, mset) {

  if(is.null(labels)) {
    mset <- .getmodules_gs(NULL, mset)
    
    mn     <- tmod_ids(mset)
    titles <- tmod_titles(mset)
    titles[is.na(titles)] <- mn[is.na(titles)]

    labels <- titles
    ne <- titles != mn
    labels[ ne ] <- sprintf("%s (%s)", labels[ne], mn[ne])
  }

  if(is.null(labels)) labels <- names(modules)
  names(labels) <- names(modules)

  return(labels)
}

#' Plot a correlation heatmap for modules
#'
#' Plot a correlation heatmap for modules
#' @param upper.cutoff Specify upper cutoff for the color palette
#' @param labels Labels for the modules (if NULL, labels will be retrieved from `mset`)
#' @param heatmap_func function drawing the heatmap
#' @param ... Any further parameters are passed to the heatmap function (by
#'        default, [pheatmap()].
#' @return Returns the return value of heatmap_func (by default, a pheatmap object).
#' @inheritParams modOverlaps
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @export
modCorPlot <- function(modules, mset=NULL, heatmap_func=pheatmap, labels=NULL, stat="jaccard", upper.cutoff=NULL, ...) {
  n <- 13

  mod_orig <- modules
  
  if(!is.list(modules)) {
    mset <- .getmodules_gs(modules, mset)
    m_label <- modules
    modules <- mset$gs2gv[modules]
    labels <- .auto_labels(labels, modules, mset)
  } else {
    mset <- .getmodules_gs(NULL, mset)
    modules <- lapply(modules, function(x) match(x, mset$gv))
    if(is.null(labels)) {
      labels <- names(modules)
    }
  }


  crmat <- modOverlaps(mod_orig, mset, stat=stat)
  if(is.null(upper.cutoff)) {
    upper.cutoff <- max(crmat)
  }

  heatmap.args <- list(...)
  if(identical(pheatmap, heatmap_func)) {
  default.args <- list(
    trace="n", 
    color=colorRampPalette(c("black", "cyan"))(n),
    breaks=seq(0, upper.cutoff, length.out=n + 1),
    scale="none",
    treeheight_col=0,
    labels_row=labels,
    labels_col=rep("", length(labels)),
    dendrogram="r")
  } else {
    default.args <- list()
  }

  heatmap.args <- c(heatmap.args, default.args)
  heatmap.args <- heatmap.args[ !duplicated(names(heatmap.args)) ]

  do.call(heatmap_func, c(list(crmat), heatmap.args))
}



## return all combinations of modules sharing a certain number of elements
## min.size: minimum number of modules in a combination
## min.overlap: minimum number of genes to share between modules
.combn_upset <- function(set, min.size=1, min.overlap=1, max.comb=NULL) {

  n <- length(set)
  setN <- names(set)
  if(min.size > n) min.size <- n
  if(!is.null(max.comb) && n > max.comb) n <- max.comb

  ret <- lapply(min.size:n, function(i) {
    combs <- combn(setN, i)

    .ret <- lapply(1:ncol(combs), function(j) {
      x <- combs[,j]
      common <- Reduce(intersect, set[x])
      full <- Reduce(union, set[x])
      total <- sum(sapply(set[x], length))
      attr(common, "modules") <- x
      attr(common, "jaccard") <- length(common)/length(full)
      attr(common, "soerensen") <- length(set[x])*length(common)/total
      attr(common, "number") <- length(common)
      attr(common, "overlap") <- length(common)/min(sapply(set[x], length))
      return(common)
    })
    names(.ret) <- apply(combs, 2, paste, collapse="_:_")
    return(.ret)
  })

  ret <- unlist(ret, recursive=FALSE)
  ret.l <- sapply(ret, length)
  ret <- ret[ ret.l > min.overlap - 1 ]

  return(ret)
}

## generate combinations for upset
.upset_generate_combinations <- function(modules, modgroups, min.size, min.overlap, max.comb, value) {

  message("generating combinations")
  ups <- lapply(modgroups, function(g) {
    message("generating combinations for group: ")
    print(g)
    .ups <- .combn_upset(modules[g], min.size=min.size, min.overlap=min.overlap, max.comb=max.comb)
    vals <- sapply(.ups, attr, value)
    .ups <- .ups[ order(-vals) ]
    attr(.ups, "maxval") <- max(vals)
    .ups
  })
  message("done")

  return(ups)
}

## find module groups for upset
.upset_find_groups <- function(modules, mset, group.cutoff, min.group, group.stat="number") {
  message(sprintf("finding groups stat=%s cutoff=%.2f", group.stat, group.cutoff))
  modgroups <- modGroups(modules, mset, min.overlap=group.cutoff, stat=group.stat)
  min.group.orig <- min.group

  while(sum(sapply(modgroups, length) >= min.group) < 1) {
    min.group <- min.group - 1
  }

  if(min.group < min.group.orig) {
    warning(sprintf("No groups found for group size of %d;\nchanging parameter to %d",
      min.group.orig, min.group))
  }

  modgroups <- modgroups[ sapply(modgroups, length) >= min.group ]
  modgroups <- lapply(modgroups, function(x) match(x, tmod_ids(mset)))
  message("done")

  return(modgroups)
}


## automagically adapt font size to the given vertical and horizontal space
.adapt_cex <- function(labels, lab.cex, horiz=1, vert=1, vert.expand=1.3, vert.mar=1, step=.95) {

  message("adapting fonts")
  while((maxwidth <- max(sapply(labels, strwidth, cex=lab.cex))) < horiz ||
        (length(labels) + vert.mar) * vert.expand * max(sapply(labels, strheight, cex=lab.cex)) < vert)  {
    lab.cex <- lab.cex / step
  }
  message(sprintf("lab.cex=%.2f", lab.cex))
  while((maxwidth <- max(sapply(labels, strwidth, cex=lab.cex))) > horiz ||
        (length(labels) + vert.mar) * vert.expand * max(sapply(labels, strheight, cex=lab.cex)) > vert)  {
    lab.cex <- lab.cex * step
  }
  message(sprintf("lab.cex=%.2f", lab.cex))
  return(lab.cex)
}


#' Upset plot
#'
#' Upset plots help to interpret the results of gene set enrichment.
#' 
#' The plot consists of three parts. The main part shows the overlaps
#' between the different modules (module can be a gene set, for example).
#' Each row corresponds to one module. Each column corresponds to an
#' intersection of one or more gene sets. Dots show which gene sets are in
#' that combination. Which combinations are shown depends on the parameters
#' `min.overlap` (which is the cutoff for the similarity measure specified by
#' the `value` parameter), the parameter `min.group` which specifies the
#' minimum number of modules in a group and the parameter `max.comb` which
#' specifies the maximum number of combinations tested (too many combinations
#' are messing the plot).
#'
#' Above the intersections, you see a plot showing a similarity measure of
#' the intersected gene sets. By default it is the number of module members
#' (genes in case of a gene set), but several
#' other measures (e.g. the Jaccard index) are also implemented.
#' 
#' To the left are the module descriptions (parameter `label`; if label is
#' empty, the labels are taken from the mset object provided or, if that is
#' NULL, from the default tmod module set). The function attempts to scale
#' the text in such a way that all labels are visible. 
#'
#' By default, upset attempts to group the modules. This is done by
#' defining a similarity measure (by default the Jaccard index, parameter
#' `group.stat`) and a cutoff threshold (parameter `group.cutoff`).
#' @return upset returns invisibly the identified module groups: a list of
#'         character vectors.
#' @param min.size minimal number of modules in a comparison to show
#' @param min.overlap smallest overlap (number of elements) between two modules to plot
#' @param min.group Minimum number of modules in a group. Group with a
#'        smaller number of members will be ignored. Change this value to 
#'        1 to see also modules which could not be grouped.
#' @param group Should the modules be grouped by the overlap?
#' @param max.comb Maximum number of combinations to show (i.e., number of
#'        dots on every vertical segment in the upset plot)
#' @param value what to show on the plot: "number" (number of common
#'        elements; default), "soerensen" (Sørensen–Dice coefficient), 
#'        "overlap" (Szymkiewicz–Simpson coefficient) or "jaccard" (Jaccard index)
#' @param cutoff Combinations with the `value` below cutoff will not be shown.
#' @param labels Labels for the modules. Character vector with the same
#'        length as `modules`
#' @param group.cutoff cutoff for group statistics 
#' @param group.stat Statistics for finding groups 
#'        (can be "number", "overlap", "soerensen" or "jaccard"; see function modOverlaps)
#' @param lab.cex Initial cex (font size) for labels
#' @param pal Color palette to show the groups. 
#' @seealso modGroups, modOverlaps
#' @inheritParams tmodUtest
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils combn
#' @examples
#' \dontrun{
#' data(Egambia)
#' design <- cbind(Intercept=rep(1, 30), TB=rep(c(0,1), each= 15))
#' library(limma)
#' fit <- eBayes( lmFit(Egambia[,-c(1:3)], design))
#' tt <- topTable(fit, coef=2, number=Inf, genelist=Egambia[,1:3] )
#' res <- tmodCERNOtest(tt$GENE_SYMBOL)
#'
#' upset(res$ID, group.cutoff=.1, value="jaccard")
#' }
#' @export
upset <- function(modules, mset=NULL, min.size=2, min.overlap=2, max.comb=4, 
  min.group=2, value="number", cutoff=NULL, labels=NULL, group.stat="jaccard",
  group.cutoff=.1,
  group=TRUE,
  pal=brewer.pal(8, "Dark2"),
  lab.cex=1) {



  value <- match.arg(value, c("number", "jaccard", "soerensen", "overlap"))
  group.stat <- match.arg(group.stat, c("number", "jaccard", "soerensen", "overlap"))

  if(!is.list(modules)) {
    mset <- .getmodules_gs(modules, mset)
    labels <- .auto_labels(labels, modules, mset)
    names(labels) <- modules

    modules_n <- modules
    modules <- mset$gs2gv
    names(modules) <- modules_n
  } else {
    mset <- .getmodules_gs(NULL, mset)
    labels <- modules
    names(labels) <- modules
  }


  if(group) {
    modgroups <- .upset_find_groups(modules_n, mset, group.cutoff, min.group, group.stat)
    group.n <- length(modgroups)
    pal <- rep(pal, ceiling(group.n/length(pal)))
    pal <- pal[1:group.n]
    names(pal) <- names(modgroups)
  } else {
    modgroups <- list(all=modules_n)
    pal <- c(all="#333333")
    pal.bar <- c(all="#333333")
  }

  message(sprintf("Number of groups: %d", length(modgroups)))
  ## translate mod IDs to numeric values

  ups <- .upset_generate_combinations(modules, modgroups, min.size, min.overlap, max.comb, value)
  ord <- order(-sapply(ups, attr, "maxval"))
  ups <- ups[ ord ]
  modgroups <- modgroups[ ord ] # reorder by best value from ups
  ups.col <- unlist(lapply(names(modgroups), function(mgn) { rep(pal[mgn], length(ups[[mgn]])) }))
  modules <- modules[ unlist(modgroups) ] # reorder by group

  ups <- unlist(ups, recursive=FALSE)
  ups.l <- sapply(ups, attr, value)
  if(!is.null(cutoff)) {
    sel <- ups.l > cutoff
    ups <- ups[ sel ]
    ups.l <- ups.l[ sel ]
    ups.col <- ups.col[ sel ]
  }
  ups.n <- length(ups)

  modlist <- unlist(lapply(ups, attr, "modules"))
  modules <- modules[ names(modules) %in% modlist ]

  # ----------- plotting ----------------------------------

  div.right <- ups.n + 1
  div.left   <- .3 * 2 * div.right
 #if((.mwidth <- max(sapply(labels, strwidth, cex=lab.cex))) < div.left) {
 #  div.left <- 1.1 * .mwidth
 #}
  div.top <- 1.1 * max(ups.l)
  div.bottom <- 1.4 * div.top

  xlim <- c(-div.left, ups.n + 1)
  ylim <- c(-div.bottom, 1.1 * max(ups.l))

  dev.hold()
  plot(NULL, xlim=xlim, ylim=ylim,
    xaxt="n", yaxt="n",
    bty="n", xlab="", ylab="")

  #rect(xlim[1], ylim[1], xlim[2], ylim[2])

  sw <- strwidth("XX")
  ylab <- switch(value, number="Number of members", jaccard="Jaccard Index", 
    soerensen="Soerensen-Dice coefficient", overlap="Overlap coefficient")
  text(-sw * 3, max(ups.l)/2, ylab, srt=90)

  ticks <- axisTicks(c(0, max(ups.l)), log=FALSE)
  axis(2, at=ticks, pos=0)

  ## ------ upper barplot
  rect(
    1:ups.n - .4,
    0, 
    1:ups.n + .4,
    ups.l,
    col=ups.col, border=NA)

  segments(0, ticks, ups.n + 1, ticks, col="white")

  ## ------ left bottom labels

  labels <- labels[ names(modules) ]
  lab.cex <- .adapt_cex(labels, lab.cex, div.left, div.bottom)
  modules.n <- length(modules)
  step <- div.bottom/(modules.n + 1)

  vpos <- (1:modules.n) * -step - step/2
  names(vpos) <- names(modules)
  text(0, vpos, labels, pos=2, cex=lab.cex)
 
  ## ------ right bottom upset plot

  sel <- rep(c(TRUE, FALSE), floor(modules.n/2))
  if(length(sel) > 0) {
    rect(.5, vpos[sel] - step/2, ups.n + .5, vpos[sel] + step/2, col="#33333333", border=NA)
  }
  segments(.5 + 1:ups.n, max(vpos) + step/2, .5 + 1:ups.n, min(vpos) - step/2, col="white")
 
  for(i in 1:ups.n) {
    mm <- attr(ups[[i]], "modules")
    mm.v <- vpos[mm]
    segments(i, min(mm.v), i, max(mm.v), lwd=4, col=ups.col[i])
    points(rep(i, length(mm.v)), mm.v, pch=19, cex=1.5, col=ups.col[i])
  }
  dev.flush()
  return(invisible(modgroups))
}


upsetPlot <- function(items, labels, overlaps, add=FALSE, coords=NULL) {
  if(is.null(names(labels))) {
    names(labels) <- items
  }

  n <- length(items)

  sel <- rep(c(TRUE, FALSE), floor(n/2))




}
