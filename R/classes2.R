## we need to use the tmod data set internally
.myDataEnv <- new.env(parent=emptyenv())

if(!exists("tmod", .myDataEnv)) {
  data("tmod", package="tmod", envir=.myDataEnv)
}

.gettmod <- function() {
  .myDataEnv[["tmod"]]
}

## check user provided mset for sanity
.mset_sanity_check_gs <- function(mset, gs=NULL) {

  required_members <-  c( "gs", "gv", "gs2gv")
  # sanity checks
  if(!all(required_members %in% names(mset))) {
    stop("Required members missing from the list mset parameter")
  }

  if(any(unlist(lapply(required_members, is.null)))) {
    stop("Required members missing from the list mset parameter")
  }

  stopifnot(is(mset$gs, "data.frame"))
  stopifnot("ID" %in% colnames(mset$gs))
  stopifnot(is(mset$gs2gv, "list"))
    
  if(any(duplicated(mset[["gs"]]$ID))) 
    stop("Gene set IDs must not be duplicated")

  stopifnot(length(mset$gs2gv) == nrow(mset$gs)) 

  if(!is.null(gs)) {

    if(!all(gs %in% mset[["gs"]]$ID )) {
      stop("Gene sets specified with the gs parameter are missing from definitions in the mset parameter")
    }
  }

  mset
}



## prepare the modules set
.getmodules_gs <- function(modules=NULL, mset="all", known.only=FALSE, skipcheck=FALSE) {

  # user provided mset
  if(is(mset, "list")) {
    mset <- .mset_sanity_check_gs(mset, modules)
    if(!is.null(modules)) { mset <- mset[ modules ] }
  } else {

    tmod <- .gettmod()

    mset <- match.arg(mset, c( "all", unique(tmod$gs$SourceID)))
    if(mset != "all") { 
      if(is.null(modules)) modules <- tmod$gs$ID
      sel <- tmod$gs$SourceID[ match(modules, tmod$gs$ID) ] == mset
      modules <- modules[ sel ]
    }

    if(!is.null(modules)) {
      modules <- modules[ modules %in% tmod$gs$ID ]
      mset <- tmod[ modules ]
    } else {
      mset <- tmod
    }
  }


  if(known.only && "Title" %in% colnames(mset$gs)) {
    mset <- mset[ ! is.na(mset$gs$Title) & ! mset$gs$Title %in% c( "TBA", "Undetermined", ""), ]
  }

  mset
}



## this is the new and more efficient tmod class

#' Check an object of class tmodGS
#'
#' Check an object of class tmodGS
#' @export
check_tmod_gs <- function(object) {

  if(is.null(object$gs))  return("Required list member `gs` missing")
  if(is.null(object$gs2gv)) return("Required list member `gs2gv` missing")

  if(!is(object$gs, "data.frame")) return("List member `gs` must be a data frame")
  if(!is(object$gs2gv, "list")) return("List member `gs2gv` must be a list")

  if(is.null(object$gs$ID))    return("List member `gs` must have column `ID`")

  if(nrow(object$gs) != length(object$gs2gv)) 
    return("Number of gene sets in `gs2gv` different from number of rows of `gs`")

  if(!is.null(object$weights)) {
    if(nrow(object$gs) != length(object$weights)) {
      return("Number of gene sets in `weights` different from number of rows of `gs`")
    }
  }

  TRUE

}



#' S4 class for tmod gene set collections
#'
#' S4 class for tmod gene set collections
#'
#' An object of class tmod contains the gene set annotations (`tmod$gs`), 
#' a character vector of gene identifiers (`tmod$gv`)
#' and a mapping between gene sets and gene identifiers (`tmod$gs2gv`).
#' Optionally, a vector of numeric weights of the same length as `gs2gv` may
#' be provided (`tmod$weights`).
#'
#' Furthermore, it may contain additional information about the gene set
#' (`tmod$info`).
#'
#' `tmod$gs` is a data frame which must contain the column "ID". Additional
#' optional columns `Title` and `Description` are recognized by some
#' functions. Any further columns may contain additional information on the
#' gene sets. The number of rows of that data frame is equal to the number
#' of gene sets in a gene set collection.
#'
#' Each element of the tmod$g2m list corresponds to the respective row of
#' the `tmod$gs` data frame. Each element is an integer vector containing
#' the positions of the gene identifiers in the `tmod$gv` character vector.
#' 
#' Objects of class tmodGS should be constructed 
#' by calling the function makeTmodGS(). This function check the validity and
#' consistency of the provided objects.
#'
#' See the package vignette for more on constructing custom module sets.
#'
#' @rdname tmodGS-class
#' @importFrom methods setClass setMethod loadMethod is new representation signature
#' @seealso tmod-data
#' @examples
#' # A minimal example
#' gs <- data.frame(ID=letters[1:3], Title=LETTERS[1:3])
#' gs2gv <- list(a=c("g1", "g2"), b=c("g3", "g4"), c=c("g1", "g2", "g4"))
#' mymset <- makeTmodGS(gs=gs, gs2gv=gs2gv)
#' @exportClass tmod
setClass( "tmodGS",
  representation("list"),
  validity=check_tmod_gs
)

## how to display a tmod object

#' Shows the tmod object
#'
#' @name show
#' @param object a tmod object
#' @aliases show,tmod-method
#' @rdname extract-methods
#' @docType methods
setMethod( "show", "tmodGS",
  function(object) {
    .catf( "An object of class \"%s\"\n", class(object) )
    .catf( "\t%d gene sets, %d genes\n", 
      nrow(object$gs),
      length(object$gv) )
  })



#' Shows the length of a tmod object
#'
#' @name length
#' @param x a tmod object
#' @aliases length,tmod-method
#' @rdname length
#' @docType methods
setMethod("length", "tmodGS",
  function(x) {
    nrow(x$gs)
  })

## allow easy subsetting of tmod objects

#' Extracts parts of a tmod object
#'
#' @param x a tmod object
#' @param i indices specifying elements to extract or replace
#' @aliases [,tmod-method
#' @rdname extract-methods
setMethod("[", c(x="tmodGS", i="ANY"),

  function(x, i) {

    object  <- x
    modules <- i

    if(is.numeric(modules)) {
      modules <- object$gs[ modules, "ID" ]
    } 

    sel <- match(modules, object$gs$ID)
    #modules <- modules[ modules %in% object$MODULES$ID ]

    gs  <- object$gs[sel, , drop=FALSE]
    gs2gv <- object$gs2gv[sel]

    ## gene list: complicated

    gv <- object$gv[ unique(unlist(gs2gv)) ]

    gs2gv <- lapply(gs2gv, function(x) {
                    orig_gv <- object$gv[ x ]
                    new_gv  <- match(orig_gv, gv)
    })


    if(!is.null(object$weights)) {
      weights <- object$weights[sel]
    } else {
      weights <- NULL
    }

    new( "tmodGS", list( 
      gs=gs,
      gs2gv=gs2gv,
      gv=gv,
      info=object$info,
      weights=weights))
  }
)

#' Convert the old tmod objects to the tmodGS objects
#'
#' Convert the old tmod objects to the tmodGS objects
#'
#' The old tmod representation was very inefficient. This function converts
#' the objects to a new representation which is both allowing faster
#' computations and more memory efficient.
#' @param x an object of class tmod
#' @return Returns an object of class tmodGS.
#' @importFrom tibble tibble
#' @export
tmod2tmodGS <- function(x) {
  stopifnot(check_tmod(x))

  gs <- tibble(x$MODULES)
  gs2gv <- x$MODULES2GENES[ gs$ID ]
  gv <- sort(unique(unlist(gs2gv)))

  gs2gv <- lapply(gs2gv, function(.) {
                  match(., gv)
  })
  names(gs2gv) <- NULL

  if(!is.null(x$WEIGHTS)) {
    weights <- x$WEIGHTS[ gs$ID ]
  } else {
    weights <- NULL
  }

  new( "tmodGS", list( 
    gs=gs,
    gs2gv=gs2gv,
    gv=gv,
    info=x$info,
    weights=weights))
}


#' @param gs2gene A list with module IDs as names. Each member of the list is a character vector with IDs of genes contained in that module
#' @param gs [Optional] A data frame with at least columns ID and Title
#' @param weights [Optional] a named numeric vector of weights for each gene set
#' @param info [Optional] a list containing meta-information about the gene set collection
#' @rdname tmodGS-class
#' @export
makeTmodGS <- function(gs2gene, gs=NULL, weights=NULL, info=NULL) {

  if(!is(gs2gene, "list")) stop("`gs2gene` must be a list")

  if(!is.null(gs)) {
    if(!is(gs, "data.frame")) {
      stop("`gs` must be a data frame")
    } else {
      gs <- tibble(gs)
    }
    if(is.null(gs$ID))    stop("`gs` must have column ID")
    stopifnot(all(gs$ID %in% names(gs2gene)))
  } 

  if(is.null(names(gs2gene))) {
    names(gs2gene) <- paste0("GS_", sprintf("%05d", 1:length(gs2gene)))
  }


  if(is.null(gs)) {
    gs <- tibble(ID=names(gs2gene))
  } else {
    gs <- gs[ match(gs$ID, names(gs2gene)) ]
  }


  gv <- sort(unique(unlist(gs2gene)))
  gs2gv <- lapply(gs2gene, function(.) { match(., gv) })
  names(gs2gv) <- NULL

  new( "tmodGS", list( 
    gs=gs,
    gs2gv=gs2gv,
    gv=gv,
    info=info,
    weights=weights))

}


