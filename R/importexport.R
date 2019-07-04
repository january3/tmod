## imports the XML format of MSigDB
.importMsigDBXML <- function(file, fields, organism) {
  msig <- list()

  fields <- unique(c("SYSTEMATIC_NAME", "STANDARD_NAME", "CATEGORY_CODE", "SUB_CATEGORY_CODE", fields))
  field.names <- c( "ID", "Title", "Category", "Subcategory", fields[-c(1:4)] )

  foo <- xmlParse(file)
  foo <- xmlToList(foo)
  
  if(organism != "all") {
    orgs <- sapply(foo, function(x) x["ORGANISM"])
    foo <- foo[ orgs == organism ]
  }

  # remove NULLs
  foo <- foo[ ! sapply(foo, is.null) ]

  msig$MODULES <- t(sapply(foo,
    function(x) x[ fields ]))
  colnames(msig$MODULES) <- field.names
  msig$MODULES <- data.frame(msig$MODULES, stringsAsFactors=FALSE, row.names=NULL)

  if(any(duplicated(msig$MODULES$ID))) {
    warning("Duplicated IDs found; automatic IDs will be generated")
    msig$MODULES$oldID <- msig$MODULES$ID
    msig$MODULES$ID    <- make.unique(as.character(msig$MODULES$ID))
  }

  rownames(msig$MODULES) <- msig$MODULES[,"ID"]

  msig$MODULES2GENES <- lapply(foo, function(x) strsplit( x["MEMBERS_SYMBOLIZED"], "," )[[1]])

  names(msig$MODULES2GENES) <- msig$MODULES$ID
  msig$GENES <- data.frame( ID=unique(unlist(msig$MODULES2GENES)))

  msig <- new("tmod", msig)
  msig
}


## imports the GMT format of MSigDB
.importMsigDBGMT <- function(file) {
  msig <- list()

  con <- file(file, open="r")
  lines <- readLines(con)
  close(con)

  ids   <- gsub( "\t.*", "", lines)
  desc  <- gsub( "^[^\t]*\t([^\t]*)\t.*", "\\1", lines )
  genes <- gsub( "^[^\t]*\t[^\t]*\t(.*)", "\\1", lines )

  msig$MODULES <- data.frame(
    ID=ids, Title=desc, stringsAsFactors=FALSE)
  if(any(duplicated(msig$MODULES$ID))) {
    warning("Duplicated IDs found; automatic IDs will be generated")
    msig$MODULES$oldID <- msig$MODULES$ID
    msig$MODULES$ID    <- make.unique(as.character(msig$MODULES$ID))
  }

  rownames(msig$MODULES) <- msig$MODULES[,"ID"]

  msig$MODULES2GENES <- strsplit(genes, "\t")
  names(msig$MODULES2GENES) <- ids

  msig$GENES <- data.frame( ID=unique(unlist(msig$MODULES2GENES)))
  msig <- new("tmod", msig)
  msig
}



#' Import data from MSigDB
#'
#' Import data from an MSigDB file in either XML or GMT format
#'
#' This command parses a file from MSigDB. Both XML and the MSigDB-specific
#' "GMT" format are supported (however, the latter is discouraged, as it
#' contains less information).
#' @param file The name of the file to parse
#' @param format Format (either "xml" or "gmt")
#' @param organism Select the organism to use. Use "all" for all organisms in the file (only for "xml" format; default: "Homo sapiens")
#' @param fields Which fields to import to the MODULES data frame (only for "xml" format)
#' @importFrom XML xmlParse xmlToList
#' @examples
#' \dontrun{
#' ## First, download the file "msigdb_v5.0.xml" from http://www.broadinstitute.org/gsea/downloads.jsp
#' msig <- tmodImportMSigDB( "msigdb_v5.0.xml" )
#' }
#' @export

tmodImportMSigDB <- function( file=NULL, format="xml", organism="Homo sapiens",
  fields=c( "STANDARD_NAME", "CATEGORY_CODE", "SUB_CATEGORY_CODE", "EXTERNAL_DETAILS_URL") ) {

  if(length(file) != 1) stop("Incorrect file parameter")
  if(!file.exists(file)) stop( sprintf("File %s does not exist", file))

  format <- match.arg(format, c( "xml", "gmt"))
  msig <- switch(format,
    xml=.importMsigDBXML(file, fields, organism),
    gmt=.importMsigDBGMT(file))

  s <- msig$MODULES$Title
  msig$MODULES$Title <- paste0(toupper(substring(s, 1,1)), tolower(substring(s, 2)) )
  msig$MODULES$Title <- gsub( "^Gse([0-9])", "GSE\\1", msig$MODULES$Title )
  msig$MODULES$Title <- gsub( "_", " ", msig$MODULES$Title )

  msig$MODULES$B <- sapply(msig$MODULES2GENES, length)
  msig
}


