
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

## cluster modules by overlaps
.cluster_mods <- function(ids, mset=NULL) {

  mset <- .getmodules_gs(ids, mset)
  if(!all(ids %in% tmod_ids(mset))) {
    warning("gene set IDs not in mset, not clustering")
    return(rev(sort(ids)))
  }

  ovr  <- modOverlaps(tmod_ids(mset), mset)
  hc <- hclust(as.dist(1 - ovr))

  ids[ hc$order ]    
}

## decide what to show on the ggPanelplot as gene set names
.get_modnames <- function(resS, add_ids) {

  if(any(resS$Title != resS$ID)) {
    modnames <- resS$Title
    if(add_ids) {
      modnames <- paste(resS$ID, modnames)
    }
  } else {
    modnames <- resS$ID
  }

  modnames
}

## sort the gene sets; cluster if necessary
.get_id_ord <- function(resS, id_order, cluster, mset) {

  if(!is.null(id_order)) {
    id_ord <- id_order[ id_order %in% resS$ID ]
    id_ord <- c(id_ord, setdiff(resS$ID, id_ord))
    id_ord <- rev(id_ord)
  } else {
    if(cluster) {
      id_ord <- .cluster_mods(resS$ID, mset)
    } else {
      id_ord <- sort(resS$ID)
    }
  }

  id_ord
}

## determine which of the columns of the results is the effect size
.get_effect_size <- function(res, effect_size) {

  if(!is.null(effect_size) && length(effect_size) == 1L && effect_size != "auto") {
    es <- effect_size[ which(effect_size %in% colnames(res))[1] ]
    if(is.na(es)) {
      stop(paste0("Effect size `", effect_size, "` not found in columns of results"))
    }
  } else {
    es <- attr(res, "effect_size")
    if(is.null(es)) {
      effect_size <- c("AUC", "D", "E")
      es <- effect_size[ which(effect_size %in% colnames(res))[1] ]
      #warning(paste0("Automatically selecting column `", es, "` as effect size"))
    }
  }

  es
}


#' Create a tmod panel plot using ggplot
#'
#' Create a tmod panel plot using ggplot
#'
#' Panel plot is a kind of a heatmap. This is the most compact way of
#' representing the results of several gene set enrichment analyses. Each row of a
#' panel plot shows the result for one gene set, and each column shows
#' corresponds to one analysis. For example, if one tests gene set enrichment
#' for a number of different contrasts, then each contrast will be represented
#' in a separate column.
#' 
#' Each cell of a panel plot shows both the effect size and the p-value.
#' The p-value is encoded as transparency: the enrichments with a lower
#' p-value have stronger colors. The size of the bar corresponds to the effect
#' size, however it is defined. For example, in case of the tmodCERNOtest,
#' tmodZtest or tmodUtest it is the area under curve, AUC.
#'
#' In addition, the bars may also encode information about the number of
#' up- or down-regulated genes. For this, an object must be created using the
#' function tmodDecideTests. This object provides  information about which
#' genes in a particular gene set are regulated in which direction.
#'
#' The order of the gene sets displayed is, by default, determined by
#' clustering the gene sets based on their overlaps. For this to work,
#' ggPanelplot must know about what genes are contained in which gene sets.
#' This is provided by the parameter `mset`. By default (when mset is NULL)
#' this is the built-in list of gene sets. If, however, the results of gene
#' set enrichment come from a different set of gene sets, you need to specify
#' it with the mset parameter. See Examples.
#' @param res a list with tmod results (each element of the list is a data
#' frame returned by a tmod test function)
#' @param sgenes a list with summaries of significantly DE genes by gene set.
#' Each a element of the list is a matrix returned by tmodDecideTests. If
#' NULL, the bars on the plot will be monochromatic.
#' @seealso [tmodDecideTests()], [tmodPanelPlot()]
#' @param colors character vector with at least 1 (when sgenes is NULL) or 3
#'        (when sgenes is not NULL) elements
#' @param effect_size name of the column which contains the effect sizes;
#'        by default, the name of this column is taken from the "effect_size"
#'        attribute of the first result table.
#' @param auc_thr gene sets enrichments with AUC (or other effect size) below `auc_thr` will not be shown
#' @param q_thr gene sets enrichments with q (adjusted P value) above `q_thr` will not be shown
#' @param filter_row_q Gene sets will only be shown if at least in one
#' contrast q is smaller than `filter_row_q`
#' @param filter_row_auc Gene sets will only be shown if at least in one
#' contrast AUC (or other effect size if specified) is larger than `filter_row_auc`
#' @param q_cutoff Any q value below `q_cutoff` will be considered equal to
#' `q_cutoff`
#' @param label_angle The angle at which column labels are shown
#' @param add_ids add IDs of gene sets to their titles on the plot
#' @param cluster whether to cluster the IDs by similarity
#' @param id_order character vector specifying the order of IDs to be
#' shown. This needs not to contain all IDs shown, but whatever IDs are in this
#' vector, they will be shown on top in the given order.
#' @param mset an object of the type 'tmodGS'. If the option `cluster` is
#' TRUE, the mset object is used to cluster the gene sets. By default, the
#' built-in transcription modules are used. See details.
#' @importFrom tidyr pivot_longer pivot_wider 
#' @importFrom purrr imap map
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes aes_string geom_bar geom_boxplot geom_point guides guide_legend lims ylab 
#' @importFrom ggplot2 ylab xlab element_text element_blank coord_flip
#' @importFrom ggplot2 scale_x_discrete scale_y_continuous scale_fill_manual facet_wrap theme
#' @importFrom tidyselect all_of starts_with
#' @return The object returned is a ggplot2 object which can be further
#'         modified the usual way.
#' @examples
#' ## prepare a set of results
#' data(Egambia)
#' genes <- Egambia$GENE_SYMBOL
#' exprs <- Egambia[ , -1:-4 ]
#' group <- gsub("\\..*", "", colnames(exprs))
#' ## test differential expression using limma
#' design <- cbind(Intercept=rep(1, 30), TB=rep(c(0,1), each= 15))
#' \dontrun{
#' library(limma)
#' fit <- eBayes( lmFit(Egambia[,-c(1:3)], design))
#' tt <- topTable(fit, coef=2, number=Inf, genelist=Egambia[,1:3] )
#' res <- tmodCERNOtest(tt$GENE_SYMBOL)
#' ## show the results using a panel plot
#' ggPanelplot(list(limma=res))
#' ## add information about the significant genes
#' sgenes <- tmodDecideTests(tt$GENE_SYMBOL, lfc=tt$logFC, pval=tt$adj.P.Val)
#' names(sgenes) <- "limma"
#' ggPanelplot(list(limma=res), sgenes=sgenes)
#' ## we will now compare the results of enrichments for different types of
#' ## differential expression tests on the data
#' res_utest <- apply(exprs, 1, function(x) wilcox.test(x ~ group)$p.value)
#' res_ttest <- apply(exprs, 1, function(x) t.test(x ~ group)$p.value)
#' ## Calculate the gene set enrichment analysis results for each of the
#' ## different types of tests
#' res_tmod <- list()
#' res_tmod$limma <- res
#' res_tmod$utest <- tmodCERNOtest(genes[order(res_utest)])
#' res_tmod$ttest <- tmodCERNOtest(genes[order(res_ttest)])
#' ggPanelplot(res_tmod)
#' ## Using the `mset` parameter
#' ## First, we generate results using a different set of gene sets
#' data(cell_signatures)
#' res_cs <- tmodCERNOtest(tt$GENE_SYMBOL, mset=cell_signatures)
#' ## the following will triger a warning that no clustering is possible
#' ## because ggPanelplot doesn't have the information about the gene set
#' ## contents
#' ggPanelplot(list(res=res_cs))
#' ## if we use the mset parameter, clustering is available
#' ggPanelplot(list(res=res_cs), mset=cell_signatures)
#' }
#' @export
ggPanelplot <- function(res, sgenes=NULL, auc_thr=.5, q_thr=.05,
                         filter_row_q=.01,
                         filter_row_auc=.65,
                         q_cutoff=1e-12,
                         cluster=TRUE,
                         id_order=NULL,
                         effect_size="auto",
                         colors=c("red", "grey", "blue"),
                         label_angle=45, add_ids=TRUE,
                         mset=NULL) {

  if(is(res, "tmodReport")) {
    res <- list(default=res)
  }

  label_angle=as.numeric(label_angle)
  resS <- tmodSummary(res) 

  modnames <- .get_modnames(resS, add_ids)
  names(modnames) <- resS$ID

  id_ord <- .get_id_ord(resS, id_order, cluster, mset)

  ## choose the column for effect size
  effect_size <- .get_effect_size(res[[1]], effect_size)

  if(is.na(effect_size)) {
    stop(sprintf("the effect size column(s) '%s' is not found in the results",
    paste(effect_size, collapse=", ")))
  }
    

  resS_l <- pivot_longer(resS, starts_with(c(effect_size, "q")), 
                 names_to=c("Param", "Contrast"), 
                 names_sep="\\.", 
                 values_to="Value") 
  resS_l <- pivot_wider(resS_l, all_of(c("ID", "Title", "Contrast")), 
                names_from="Param", 
                values_from="Value") 

  resS_l[["q"]] <- ifelse(is.na(resS_l[["q"]]), 1, resS_l[["q"]])
  resS_l[[effect_size]] <- ifelse(is.na(resS_l[[effect_size]]), .5, resS_l[[effect_size]])
  resS_l <- resS_l[ resS_l[[effect_size]] > auc_thr & resS_l[["q"]] < q_thr, ]

  if(!is.na(q_cutoff)) {
    resS_l[["q"]] <- ifelse(resS_l[["q"]] < q_cutoff, q_cutoff, resS_l[["q"]])
  }

  resS_l[["alpha"]] <- -log10(resS_l[["q"]]) / -log10(q_cutoff)
  selMod <- resS$ID

  if(!is.na(filter_row_auc)) {
    .s <- resS_l$ID[ resS_l[[effect_size]] > filter_row_auc ]
    selMod <- intersect(selMod, .s)
  }

  if(!is.na(filter_row_q)) {
    .s <- resS_l$ID[ resS_l[["q"]] < filter_row_q ]
    selMod <- intersect(selMod, .s)
  }

  resS_l <- resS_l[ resS_l[["ID"]] %in% selMod, ]

  # in order to keep also the empty contrasts and have them shown on the plot
  # in precisely the same order as in the input list
  resS_l$Contrast <- factor(resS_l$Contrast, levels=names(res))

  minq <- max(-log10(resS_l[["q"]]))

  if(is.null(sgenes)) {
    resS_l$ID <- factor(resS_l$ID, levels=id_ord)

  ret <- ggplot(resS_l, aes(x=.data[["ID"]], y=.data[[effect_size]], 
                 contrast=.data[["Contrast"]],
                 id=.data[["ID"]],
                 alpha=-log10(.data[["q"]]))) + 
    facet_wrap(~ .data[["Contrast"]], nrow=1, drop=FALSE) + 
    geom_bar(stat="identity", fill=colors[1]) + 
    coord_flip() +
    scale_x_discrete(breaks=names(modnames), labels=modnames) +
    theme(strip.text.x = element_text(angle=label_angle), strip.background=element_blank(),
          axis.text.x=element_text(angle=90), axis.title.y=element_blank()) +
    scale_y_continuous(breaks=c(0, .5, 1), labels=c("0", ".5", "1"), limits=c(0, 1)) +
    guides(alpha = guide_legend(override.aes = list(fill = colors[1]))) +
    lims(alpha = c(0, minq)) +
    ylab("Effect size")
 

    #stop("Don't know how to do it yet :-(")
    return(ret)
  }

  pieS <- imap(sgenes, ~ { 
               colnames(.x) <- paste0(.y, '.', colnames(.x))
               .x <- as.data.frame(.x)
               .x$ID <- rownames(.x)
               .x
                })

  pieS <- Reduce(function(x, y) merge(x, y, all=TRUE, by="ID"), pieS) 

  pieS <-  pivot_longer(pieS, -which(colnames(pieS) == "ID"), names_to=c("Contrast", "Direction"), 
                 names_sep="\\.", 
                 values_to="Number")

  ## add pie information about significant genes
  df <- merge(resS_l, pieS, by=c("ID", "Contrast"), all.x=TRUE) 
  df$id.cntr <- paste(df$ID, df$Contrast)

  tots <- tapply(df$Number, df$id.cntr, sum)
  df$Tot <- tots[ df$id.cntr ]
  df$Number <- df$Number * df[[effect_size]] / df$Tot
  df$Direction <- factor(df$Direction, levels=c("up", "N", "down"))

 #df <- df %>% group_by(paste0(.data[["Contrast"]], .data[["Direction"]])) %>%
 #  slice(match(resS$ID, .data[["ID"]])) %>%
 #  ungroup()

  df$Contrast <- factor(df$Contrast, levels=names(res))

  names(colors) <- levels(df$Direction)

  df$ID <- factor(df$ID, levels=rev(id_ord))

  ggplot(df, aes(x=.data[["ID"]], y=.data[["Number"]], 
                 fill=.data[["Direction"]],
                 contrast=.data[["Contrast"]],
                 id=.data[["ID"]],
                 alpha=-log10(.data[["q"]]))) + 
    facet_wrap(~ .data[["Contrast"]], nrow=1, drop=FALSE) + 
    geom_bar(stat="identity") + 
    coord_flip() +
    scale_fill_manual(values=colors) +
    scale_x_discrete(breaks=names(modnames), labels=modnames) +
    theme(strip.text.x = element_text(angle=label_angle), strip.background=element_blank(),
          axis.text.x=element_text(angle=90), axis.title.y=element_blank()) +
    scale_y_continuous(breaks=c(0, .5, 1), labels=c("0", ".5", "1"), limits=c(0, 1)) +
    guides(alpha = guide_legend(override.aes = list(fill = "grey"))) +
    lims(alpha = c(0, minq)) +
    ylab("Effect size")
    
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
#' @seealso \code{\link{tmod-package}}, [evidencePlot()]
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


#' Create an evidence plot for a module (ggplot2 version)
#' 
#' @inheritParams evidencePlot
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 geom_line geom_polygon geom_segment
#' @importFrom ggplot2 scale_color_manual geom_abline
#' @importFrom purrr imap_dfr map_dfr 
#' @export
ggEvidencePlot <- function(l, m, mset=NULL, filter=FALSE, unique=TRUE,
  gene.labels=NULL,
  gene.colors=NULL
  ) {

  mset <- .getmodules_gs(m, mset)

  l_orig <- l
  tmp    <- .prep_list(l, x=l_orig, mset=mset, filter=filter, nodups=unique)
  #keep <- match(tmp$x, l_orig) #l_orig %in% tmp$x
  l_orig <- tmp$x
  l      <- tmp$l

  mods <- tmod_ids(mset)
  names(mods) <- mods

  mm <- getModuleMembers(mods, mset)  
  mm <- map(mm, ~ .x[ .x %in% l_orig ])

  coords <- imap_dfr(mm, ~ {
    pos <- sort(match(.x, l_orig))
    pos <- rep(pos, each=2)
    pos <- c(1, pos, length(l_orig))
    y <- rep(c(0:length(.x)), each=2) / length(.x)
    data.frame(x=pos, y=y, "mod"=.y)
  })

  coords_segm <- imap_dfr(mm, ~ {
    .match <- match(.x, l_orig)
    .ord   <- order(.match)
    x <- sort(.match)
    ret <- data.frame(x=x, y=-.1, xend=x, yend=0, "mod"=.y, label=.x[.ord], gene=.x[.ord])
    if(!is.null(gene.labels) && !gene.labels) {
      ret$label <- ""
    }
    ret
  })




  if(length(mm) > 1) {
    pal <- tmodPal(length(mm), alpha=1)
  } else {
    pal <- "black"
  }
  names(pal) <- names(mm)

  if(!is.null(gene.colors)) {
    pal <- c(pal, gene.colors)
    names(pal) <- c(names(mm), names(gene.colors))
  }


  ret <- ggplot(coords, aes_string(x="x", y="y", color="mod")) 

  ret <- ret + 
    geom_segment(aes(x=0, y=0, xend=length(l_orig), yend=1), col="darkgrey") +
    geom_line() + 
    geom_polygon(data=data.frame(x=c(1, 1, length(l_orig), length(l_orig)), 
                                   y=c(-.1, 0, 0, -.1)), aes_string(x="x", y="y"), 
                                   color="grey",
                                   fill="LightGrey") +
    scale_y_continuous(breaks=seq(0, 1, by=.25)) +
    xlab("Gene list") +
    ylab("Fraction of gene set")

  if(!is.null(gene.colors)) {
    ret <- ret +
    geom_segment(data=coords_segm, aes_string(x="x", y="y", xend="xend", yend="yend", color="gene")) +
    geom_text_repel(data=coords_segm, aes_string(x="x", y="yend", label="label", color="gene"), 
      max.overlaps=999, direction="x", angle=90,
      hjust="left", min.segment.length=.01, xlim=c(1, length(l_orig)), ylim=c(.05, .5))
  } else {
    ret <- ret +
    geom_segment(data=coords_segm, aes_string(x="x", y="y", xend="xend", yend="yend", color="mod")) +
    geom_text_repel(data=coords_segm, aes_string(x="x", y="yend", label="label", color="mod"), 
      max.overlaps=999, direction="x", angle=90,
      hjust="left", min.segment.length=.01, xlim=c(1, length(l_orig)), ylim=c(.05, .5))
  }

  ret <- ret + scale_color_manual(values=pal, breaks=names(pal[1:length(mm)])) +
    guides(color=guide_legend(title="Gene set")) 

  if(length(mm) < 2) {
    ret <- ret + theme(legend.position = "none")
  }

  ret

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
#' @seealso [tmod-package()], [hgEnrichmentPlot()]
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
  keep <- match(tmp$x, l_orig) #l_orig %in% tmp$x
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
      # calculate the actual size of the labels in user units
      lw <- max(strheight(gene.labels, units="i", cex=gl.cex))
      lw <- grconvertX( c(0, lw), from="i", to="u")
      lw <- lw[2] - lw[1]

      # initial positions of the labels, stacked one on the other
      lpos <- lw * 1.2 * (1:length(gene.labels))

      # remove labels which are too far to the right
      maxsel <- lpos < n
      lpos <- lpos[maxsel]

      # find the position of the gene in the original list
      # rpos is the x coordinate of the line from the rug to the label
      gi <- names(gene.labels)[maxsel]
      rpos <- match(gi, l_orig)

      if(is.null(gene.colors)) {
        cols <- "black"
      } else {
        cols <- ifelse(is.na(gene.colors[gi]), "black", gene.colors[gi])
        cols <- cols[maxsel]
      }

      lwds <- ifelse(is.na(gene.lines[gi]), lwd[1], gene.lines[gi])
      # this needs to be tested more
      if(length(lwds) > 1) lwds <- lwds[maxsel]

      # shift the position of the labels to the left if there is still some space
      for(i in 1:length(lpos)) {
        if(rpos[i] + lw > lpos[i]) {
          sel <- i:length(lpos)
          lpos[sel] <- lpos[sel] + (rpos[i] - lpos[i] + lw)
        }
      }

      text(lpos, 0.1, gene.labels[maxsel], srt=90, pos=4, col=cols)
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
#' @seealso [modGroups()], [modOverlaps()]
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
