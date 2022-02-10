## this is the new and more efficient tmod class

#' Check an object of class tmodGS
#'
#' Check an object of class tmodGS
#' @export
check_tmod_gs <- function(object) {

  if(is.null(object$gs))  return("Required list member `gs` missing")
  if(is.null(object$m2g)) return("Required list member `m2g` missing")

  if(!is(object$gs, "data.frame")) return("List member `gs` must be a data frame")
  if(!is(object$m2g, "list")) return("List member `m2g` must be a list")

  if(is.null(object$gs$ID))    return("List member `gs` must have column `ID`")

  if(nrow(object$gs) != length(object$m2g)) 
    return("Number of gene sets in `m2g` different from number of rows of `gs`")

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
#' and a mapping between gene sets and gene identifiers (`tmod$m2g`).
#' Optionally, a vector of numeric weights of the same length as `m2g` may
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
#' @rdname tmod-class
#' @importFrom methods setClass setMethod loadMethod is new representation signature
#' @seealso tmod-data
#' @examples
#' # A minimal example
#' m <- data.frame(ID=letters[1:3], Title=LETTERS[1:3])
#' m2g <- list(a=c("g1", "g2"), b=c("g3", "g4"), c=c("g1", "g2", "g4"))
#' mymset <- makeTmod(modules=m, modules2genes=m2g)
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
    m2g <- object$m2g[sel]

    ## gene list: complicated

    gv <- object$gv[ unique(unlist(m2g)) ]

    m2g <- lapply(m2g, function(x) {
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
      m2g=m2g,
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
#' @export
tmod2tmodGS <- function(x) {
  stopifnot(check_tmod(x))

  gs <- x$MODULES
  m2g <- x$MODULES2GENES[ gs$ID ]
  gv <- sort(unique(unlist(m2g)))

  m2g <- lapply(m2g, function(.) {
                  match(., gv)
  })

  if(!is.null(x$WEIGHTS)) {
    weights <- x$WEIGHTS[ gs$ID ]
  } else {
    weights <- NULL
  }

  new( "tmodGS", list( 
    gs=gs,
    m2g=m2g,
    gv=gv,
    info=object$info,
    weights=weights))
}



