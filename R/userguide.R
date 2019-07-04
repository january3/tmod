#' Display the tmod user guide
#'
#' Display the tmod user guide
#' 
#' @return invisibly return the location of the PDF file
#' @export
tmodUserGuide <- function() {

  ## code inspired by tools:::print.vignette
  pdfviewer <- getOption("pdfviewer")

  f <- system.file("doc", "tmod_user_manual.pdf", package = "tmod")
  if(identical(pdfviewer, "false")) 
    stop(sprintf("Cannot display the file %s", f))

  if (.Platform$OS.type == "windows" && 
      identical(pdfviewer, file.path(R.home("bin"), "open.exe"))) {
    shell.exec(f)
  }  else {
    system2(pdfviewer, shQuote(f), wait = FALSE)
  }
    
  return(invisible(f))
}
