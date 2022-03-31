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

  if(!"tmodGS" %in% class(mset)) {
    class(mset) <- c(class(mset), "tmodGS")
  }

  mset
}



## prepare the modules set
.getmodules_gs <- function(modules=NULL, mset="all", known.only=FALSE, skipcheck=FALSE) {

  if(is(mset, "tmod")) {
    warning(
"You are loading an obsolete version of tmod R object.
The class `tmod` has been retired.
The data will still work, but it will incur a penalty
on the computational time. Please use the `tmod2tmodGS`
function to convert your object to the new tmodGS class.")
    mset <- tmod2tmodGS(mset)
  }


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
#' @param object an object of class tmodGS
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



#' S3 class for tmod gene set collections
#'
#' S3 class for tmod gene set collections
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
#' The makeTmod function remains for compatibility with previous versions
#' of the package. It produces the objects of the new class tmodGS,
#' however.
#'
#' See the package vignette for more on constructing custom module sets.
#'
#' @rdname tmodGS-class
#' @importFrom methods setClass setMethod loadMethod is new representation signature
#' @seealso tmod-data
#' @param ... further arguments passed to `print()`
#' @param gs2gene,modules2genes A list with module IDs as names. Each member of the list is a character vector with IDs of genes contained in that module
#' @param gs,modules [Optional] A data frame with at least columns ID and Title
#' @param genes2modules,genes Ignored
#' @param weights [Optional] a named numeric vector of weights for each gene set
#' @param info [Optional] a list containing meta-information about the gene set collection
#' @examples
#' # A minimal example
#' gs <- data.frame(ID=letters[1:3], Title=LETTERS[1:3])
#' gs2gv <- list(a=c("g1", "g2"), b=c("g3", "g4"), c=c("g1", "g2", "g4"))
#' mymset <- makeTmodGS(gs2gene=gs2gv, gs=gs)
#' str(mymset)
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
    stopifnot(all(names(gs2gene) %in% gs$ID))
    gs <- gs[ match(names(gs2gene), gs$ID), ]
  }


  gv <- sort(unique(unlist(gs2gene)))
  gs2gv <- lapply(gs2gene, function(.) { match(., gv) })
  names(gs2gv) <- NULL

  ret <- list( 
    gs=gs,
    gs2gv=gs2gv,
    gv=gv,
    info=info,
    weights=weights)
  class(ret) <- c(class(ret), "tmodGS")

  ret
}

#' @rdname tmodGS-class
#' @export
makeTmod <- function(modules, modules2genes, genes2modules=NULL, genes=NULL) {
  makeTmodGS(gs2gene=modules2genes, gs=modules, weights=NULL, info=NULL)
}


#' @rdname tmodGS-class
#' @param check_sanity whether the tmodGS object should be tested for correctness
#' @export
as_tmodGS <- function(x, check_sanity=FALSE) {
  stopifnot(is.list(x)) 
  
  if(!"tmodGS" %in% class(x)) {
    class(x) <- c("tmodGS", class(x))
  }

  if(check_sanity) {
    x <- .mset_sanity_check_gs(x)
  }

  x
}


## how to display a tmod object

#' @param x a tmod object
#' @rdname tmodGS-class
#' @export
print.tmodGS <- function(x, ...) {
    .catf("An object of class \"tmodGS\"\n" )
    .catf("\t%d gene sets, %d genes\n", 
      nrow(x$gs),
      length(x$gv))
}

#' Query and set IDs of gene sets in a tmodGS object
#' 
#' Query and set IDs (tmod_id) or Titles (tmod_title) of gene sets in a tmodGS object
#' @param x an object of class tmodGS
#' @param value a character vector of unique IDs
#' @return Returns character vector corresponding to x$gs$ID
#' @examples
#' data(tmod)
#' mset <- tmod[ c("LI.M37.0", "LI.M75", "LI.M3") ]
#' tmod_ids(mset)
#' tmod_ids(mset) <- c("em", "pstrem", "bzdrem")
#' tmod_titles(mset) <- c("foo", "bar", "baz")
#' mset$gs
#' @export
tmod_ids <- function(x) {
  x$gs$ID
}

#' @rdname tmod_ids
#' @export
`tmod_ids<-` <- function(x, value) {
  stopifnot(length(value) == nrow(x$gs))
  stopifnot(all(!duplicated(value)))
  x$gs$ID <- value
  x
}

#' @rdname tmod_ids
#' @export
tmod_titles <- function(x) {
  if("Title" %in% names(x$gs)) {
    return(x$gs$Title)
  } else {
    return(NULL)
  }
}

#' @rdname tmod_ids
#' @export
`tmod_titles<-` <- function(x, value) {
  stopifnot(length(value) == nrow(x$gs))
  stopifnot(all(!duplicated(value)))
  x$gs$Title <- value
  x
}


#' @param x a tmod object
#' @rdname tmodGS-class
#' @export
length.tmodGS <- function(x) {
    nrow(x$gs)
}  

#' @param x a tmod object
#' @param i indices specifying elements to extract or replace
#' @rdname tmodGS-class
#' @export
`[.tmodGS` <- function(x, i) {

  object  <- x

  if(is.logical(i)) {
    stopifnot(length(i) == nrow(object$gs))
    i <- which(i)
  }


  if(!is.numeric(i)) {
    stopifnot(all(i %in% object$gs$ID))
    sel <- match(i, object$gs$ID)
  } else {
    sel <- i
  }

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

  as_tmodGS(list( 
    gs=gs,
    gs2gv=gs2gv,
    gv=gv,
    info=object$info,
    weights=weights), check_sanity=FALSE)
}

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

  as_tmodGS(list( 
    gs=gs,
    gs2gv=gs2gv,
    gv=gv,
    info=x$info,
    weights=weights))
}


.invert_hash <- function(l) {
  if(is.null(names(l))) names(l) <- seq_along(l)
  ret <- split(rep(names(l), lengths(l)), unlist(l))
  ret
}


## checking anold style tmod object
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



