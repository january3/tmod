## we need to use the tmod data set internally
.myDataEnv <- new.env(parent=emptyenv())

if(!exists("tmod", .myDataEnv)) {
  data("tmod", package="tmod", envir=.myDataEnv)
}

.gettmod <- function() {
  .myDataEnv[["tmod"]]
}


## check user provided mset for sanity
.mset_sanity_check <- function(mset, modules=NULL) {

  # sanity checks
  if(!all( c( "MODULES", "MODULES2GENES", "GENES") %in% names(mset)))
    stop("Required members missing from the list mset parameter")

  for( i in c( "MODULES", "GENES" )) {
    if(is.null(mset[[i]])) 
      stop(sprintf("Member %s of mset is NULL", i))
    if(!is(mset[[i]], "data.frame")) 
      stop(sprintf("Member %s of mset is not a data frame", i))
  }

  if(!is(mset[["MODULES2GENES"]], "list"))
    stop("Member MODULES2GENES of mset is not a list")

  if(!all(c("ID", "Title") %in% colnames(mset[["MODULES"]])))
    stop("Required columns missing from member MODULES")

  if(any(duplicated(mset[["MODULES"]]$ID))) 
    stop("Module IDs must not be duplicated")

  mset[["MODULES"]]$ID <- as.character(mset[["MODULES"]]$ID)
  rownames(mset[["MODULES"]]) <- mset[["MODULES"]]$ID

  if(!"ID" %in% colnames(mset[["GENES"]]))
    stop("Required column ID missing from member GENES")
  mset[["GENES"]]$ID <- as.character(mset[["GENES"]]$ID)

  missing <- !mset[["MODULES"]]$ID %in% names(mset[["MODULES2GENES"]])
  if(any(missing)) {
    stop(sprintf("Modules from MODULES member are missing in MODULES2GENES, first is %s", mset[["MODULES"]]$ID[missing][1] ))
  }

  if(!is.null(modules)) {

    if(!all(modules %in% mset[["MODULES"]]$ID )) {
      stop("Modules specified with the modules parameter are missing from definitions in the mset parameter")
    }
  }

  if(!is(mset, "tmod")) mset <- new("tmod", mset)
  mset
}


## prepare the modules set
.getmodules2 <- function(modules=NULL, mset="all", known.only=FALSE, skipcheck=FALSE) {

  # user provided mset
  if(is(mset, "list")) {
    mset <- .mset_sanity_check(mset, modules)
    if(is.null(modules)) modules <- mset[["MODULES"]]$ID
  } else {

    tmod <- .gettmod()

    mset <- match.arg( mset, c( "all", unique( tmod$MODULES$SourceID )) )
    if( mset != "all" ) { 
      if( is.null( modules ) ) modules <- tmod$MODULES$ID
      modules <- modules[ tmod$MODULES[modules,]$SourceID == mset ]
    }

    mset <- tmod
  }

  # filter the modules if hand-picked
  if(!is.null(modules)) {
    mset <- mset[modules,]
  }

  if(known.only) {
    mset <- mset[ ! is.na(mset$MODULES$Title) & ! mset$MODULES$Title %in% c( "TBA", "Undetermined"), ]
  }

  mset
}


#' Convert between matrix representation of modules and tmod objects 
#'
#' Converts a matrix where columns correspond to modules and rows to
#' individual variables into an tmod object ("mset") represenation and back
#'
#' A matrix mapping genes onto modules and vice versa can be converted with
#' mtx2mset into a tmod object. The numeric matrix either only contains '0'
#' and '1' values indicating presence or absence of a variable (gene) in a module
#' (gene set), or any other numbers which are treated as weights. In the latter case, '0' and only '0' is
#' interpreted as absence of a variable in a module; any other value is
#' interpreted as a weight and stored in the WEIGHTS slot of the tmod object.
#'
#' The mset2mtx function does the reverse.
#' @seealso tmod-class
#' @param x for mtx2mset, a numeric matrix with named rows and columns; for
#'        mset2mtx, an object of the class tmod
#' @export
mtx2mset <- function(x) {

  x <- x[ apply(x, 1, function(xx) any(xx > 0)),, drop=F ]

  modules <- colnames(x)
  genes   <- rownames(x)

  if(is.null(modules))     stop("x must have row names")
  if(is.null(genes))       stop("x must have column names")

  genes2modules <- sapply(genes, function(g) modules[ x[g,] > 0 ], simplify=FALSE)
  modules2genes <- sapply(modules, function(m) genes[ x[,m] > 0 ], simplify=FALSE)
  weights <- NULL

  if(!all(x %in% c(0,1))) {
    weights <- sapply(modules, function(m) { 
      print(m)
      sel <- x[,m] > 0
      ret <- x[sel,m]
  }, simplify=FALSE)
  }

  modules <- data.frame(ID=modules, Title=modules, row.names=modules)
  genes   <- data.frame(ID=genes, row.names=genes)

  new( "tmod", list( 
    MODULES=modules,
    GENES=genes,
    MODULES2GENES=modules2genes,
    GENES2MODULES=genes2modules,
    WEIGHTS=weights ))
}


#' @rdname mtx2mset
#' @export
mset2mtx <- function(x) {

  x <- .mset_sanity_check(x)

  ng <- nrow(x$GENES)
  nm <- nrow(x$MODULES)

  ret <- matrix(0, nrow=ng, ncol=nm)
  rownames(ret) <- x$GENES$ID
  colnames(ret) <- x$MODULES$ID

  if(!is.null(x$WEIGHTS)) {
    for(m in x$MODULES$ID) {
      sel <- x$MODULES2GENES[[m]]
      ret[sel,m] <- x$WEIGHTS[[m]][sel]
    }
  } else {
    for(m in x$MODULES$ID) {
      ret[x$MODULES2GENES[[m]],m] <- 1
    }
  }

  ret
}


## checking a new tmod object
check_tmod <- function(object) {

  if(is.null(object$MODULES))       return("Required list member MODULES missing")
  if(is.null(object$MODULES2GENES)) return("Required list member MODULES2GENES missing")

  if(!is(object$MODULES, "data.frame")) return("MODULES must be a data frame")
  if(!is(object$MODULES2GENES, "list")) return("MODULES2GENES must be a list")

  if(is.null(object$MODULES$ID))    return("MODULES must have columns ID and Title")
  if(is.null(object$MODULES$Title)) return("MODULES must have columns ID and Title")

  if(!all(object$MODULES$ID %in% names(object$MODULES2GENES)))
    return("All MODULES$ID must be found in names of MODULES2GENES")


  if(!is.null(object$WEIGHTS)) {
    if(!all(object$MODULES$ID %in% names(object$WEIGHTS)))
      return("All MODULES$ID must be found in names of WEIGHTS")
  }

  TRUE

}

.invert_hash <- function(l) {
  if(is.null(names(l))) names(l) <- seq_along(l)
  ret <- split(rep(names(l), lengths(l)), unlist(l))
  ret
}


#' Rename module IDs
#'
#' Change the IDs of modules
#'
#' Rename the modules by replacing the selected IDs with new IDs.
#' @param mset a module set (object of class tmod)
#' @param newids a named vector of unique new IDs
#' @return object of class tmod with renamed IDs
#' @export
renameMods <- function(mset, newids) {
  mset <- .getmodules2(NULL, mset)

  if(any(duplicated(newids))) 
    stop("new IDs must be unique")

  from <- names(newids)
  if(is.null(from)) stop("newids must be a named vector")

  all <- mset$MODULES$ID

  from <- from[ from %in% all ]
  if(!length(from)) stop("no names(newids) match the IDs from mset")

  newids <- newids[from]

  remain <- setdiff(all, from)
  newids <- c(newids, setNames(remain, remain))

  ## this should not happen
  if(!setequal(names(newids), all)) 
    stop("ooops, something went wrong")

  mset$MODULES$ID <- rownames(mset$MODULES) <- newids[mset$MODULES$ID]
  names(mset$MODULES2GENES) <- newids[names(mset$MODULES2GENES)]

  mset$GENES2MODULES <- .invert_hash(mset$MODULES2GENES)

  mset
}

#' @param modules A data frame with at least columns ID and Title
#' @param modules2genes A list with module IDs as names. Each member of the list is a character vector with IDs of genes contained in that module
#' @param genes2modules A list with gene IDs as names. Each member of the list is a character vector with IDs of modules in 
#'        which that gene is contained. This object will be created automatically if the provided parameter is NULL
#' @param genes A data frame with meta-information on genes. Must contain the column ID. If NULL, then a data frame with only one column (ID) will be created automatically.
#' @rdname tmod-class
#' @export
makeTmod <- function(modules, modules2genes, genes2modules=NULL, genes=NULL) {


  if(!is(modules, "data.frame")) stop("MODULES must be a data frame")
  if(!is(modules2genes, "list")) stop("MODULES2GENES must be a list")

  if(is.null(modules$ID))    stop("MODULES must have column ID")
  if(is.null(modules$Title)) stop("MODULES must have column Title")

  if(is.null(genes)) {
    genes <- unique(unlist(modules2genes))
    genes <- data.frame(ID=genes, stringsAsFactors=FALSE)
  } else {
    if(is.null(genes$ID)) {
      stop("genes must have an ID column")
    }
  }

  genes$ID     <- as.character(genes$ID)
  modules$ID   <- as.character(modules$ID)
  rownames(genes)   <- genes$ID
  rownames(modules) <- modules$ID

  if(!all(modules$ID %in% names(modules2genes)))
    stop("All MODULES$ID must be found in names of MODULES2GENES")


  if(is.null(genes2modules)) {
    genes2modules <- .invert_hash(modules2genes)
  } else {
    if(!all(genes$ID %in% names(genes2modules)))
      stop("All genes$ID must be found in names of genes2modules")
  }


  mset <- new("tmod", list(MODULES=modules, MODULES2GENES=modules2genes, GENES2MODULES=genes2modules, GENES=genes))
  mset
}

## bind two modules together
.module_bind <- function(x, y) {
  x <- .getmodules2(NULL, x)
  y <- .getmodules2(NULL, y)
  m.com.cols <- intersect(colnames(x$MODULES), colnames(y$MODULES))
  g.com.cols <- intersect(colnames(x$GENES), colnames(y$GENES))

  if(length(intersect(x$MODULES$ID, y$MODULES$ID)) > 0) 
    stop("module IDs of x and y cannot overlap")

  newM <- rbind(x$MODULES[,m.com.cols, drop=FALSE],
                y$MODULES[,m.com.cols, drop=FALSE])
  newM2G <- c(x$MODULES2GENES, y$MODULES2GENES)

  newG <- rbind(x$GENES[,g.com.cols, drop=FALSE],
                y$GENES[,g.com.cols, drop=FALSE])

  newG <- newG[ !duplicated(newG$ID), , drop=FALSE]

  newG2M <- .invert_hash(newM2G)

  mset <- new("tmod", list(MODULES=newM, MODULES2GENES=newM2G, GENES2MODULES=newG2M, GENES=newG))
  mset
}


#' Bind two or more transcriptional module sets
#' 
#' Bind two or more transcriptional module sets
#'
#' @param x First module to bind
#' @param ... further modules will be concatenated with x
#' @return an object of class tmod which is the result of binding together
#'         all modules
#' @export
mbind <- function(x, ...) {
  args <- list(...)

  args <- c(list(x), args)
  Reduce(.module_bind, args)
}



#' S4 class for tmod
#'
#' S4 class for tmod
#'
#' An object of class tmod contains the module annotations (tmod$MODULES), gene to
#' module (tmod$GENES2MODULES) and module to gene (tmod$MODULES2GENES) mappings
#' and a gene list (tmod$GENES).
#'
#' tmod$MODULES and tmod$GENES are data frames, while tmod$MODULES2GENES and
#' tmod$GENES2MODULES are lists with, respectively, module and gene
#' identifiers as names. The data frames MODULES and GENES must contain the
#' column "ID", and the data frame MODULES must contain additionally the
#' column "Title".
#' 
#' Objects of class tmod should be constructed 
#' by calling the function makeTmod(). This function check the validity and
#' consistency of the provided objects, and, if necessary, automatically
#' creates the members GENES and GENES2MODULES.
#'
#' See the package vignette for more on constructing custom module sets.
#'
#' @rdname tmod-class
#' @importFrom methods setClass setMethod loadMethod is new representation signature
#' @seealso tmod-data
#' @examples
#' # A minimal example
#' m <- data.frame(ID=letters[1:3], Title=LETTERS[1:3])
#' m2g <- list(a=c("g1", "g2"), b=c("g3", "g4"), c=c("g1", "g2", "g4"))
#' mymset <- makeTmod(modules=m, modules2genes=m2g)
#' @exportClass tmod
setClass( "tmod",
  representation("list"),
  validity=check_tmod
)

## how to display a tmod object

#' Shows the tmod object
#'
#' @name show
#' @param object a tmod object
#' @aliases show,tmod-method
#' @rdname extract-methods
#' @docType methods
setMethod( "show", "tmod",
  function(object) {
    .catf( "An object of class \"%s\"\n", class(object) )
    .catf( "\t%d modules, %d genes\n", 
      nrow(object$MODULES),
      nrow(object$GENES) )
  })



#' Shows the length of a tmod object
#'
#' @name length
#' @param x a tmod object
#' @aliases length,tmod-method
#' @rdname length
#' @docType methods
setMethod("length", "tmod",
  function(x) {
    nrow(x$MODULES)
  })

## allow easy subsetting of tmod objects

#' Extracts parts of a tmod object
#'
#' @param x a tmod object
#' @param i indices specifying elements to extract or replace
#' @aliases [,tmod-method
#' @rdname extract-methods
setMethod("[", c(x="tmod", i="ANY"),

  function(x, i) {

    object  <- x
    modules <- i
    modules <- object$MODULES[ modules, "ID" ]
    #modules <- modules[ modules %in% object$MODULES$ID ]

    MODULES       <- object$MODULES[modules,,drop=FALSE]
    MODULES2GENES <- object$MODULES2GENES[modules]

    glist <- unique(unlist(MODULES2GENES))

    GENES <- object$GENES[ object$GENES$ID %in% glist,,drop=FALSE]

    GENES2MODULES <- NULL
    if(!is.null(object$GENES2MODULES)) {
      GENES2MODULES <- object$GENES2MODULES[object$GENES$ID]
    }

    if(!is.null(object$WEIGHTS)) {
      WEIGHTS <- object$WEIGHTS[modules]
    } else {
      WEIGHTS <- NULL
    }

    new( "tmod", list( 
      MODULES=MODULES,
      GENES=GENES,
      MODULES2GENES=MODULES2GENES,
      GENES2MODULES=GENES2MODULES,
      WEIGHTS=WEIGHTS ))
  }
)





#' S4 class for tmodSummary
#'
#' S4 class for tmodSummary
#'
#' An object of class tmodSummary contains results from multiple
#' enrichment tests for the same set of gene sets, for example from running
#' tmodCERNOtest on different contrasts.
#'
#' @rdname tmodSummary-class
#' @importFrom methods setClass setMethod loadMethod is new representation signature
#' @seealso tmodSummary-data
#' @exportClass tmodSummary
setClass(Class="tmodSummary",
  slots=c(pval.col="character", effect.col="character", rid="character"),
  contains="data.frame")



#' Shows a tmodSummary object
#'
#' @name show
#' @param object a tmod object
#' @aliases show,tmodSummary-method
#' @rdname extract-methods
#' @docType methods
setMethod( "show", "tmodSummary",
  function(object) {
    show(as(object, "data.frame"))
  })



