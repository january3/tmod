#' Simple Pie Chart
#' 
#' The simplePie function draws a simple pie chart at specified coordinates with
#' specified width, height and color. The simpleRug function draws a
#' corresponding rug plot, while simpleBoxpie creates a "rectangular pie
#' chart" that is considered to be better legible than the regular pie.
#' 
#' simplePie() draws a pie chart with width w and height h at coordinates
#' (x,y). The size of the slices is taken from the numeric vector v, and
#' their color from the character vector col.
#' @param border color of the border. Use NA (default) or NULL for no border
#' @param x,y coordinates at which to draw the plot
#' @param w,h width and height of the plot
#' @param v sizes of the slices
#' @param col colors of the slices
#' @param res resolution (number of polygon edges in a full circle)
#' @param grid boxpie only: the grid over which the areas are distributed.
#'             Should be roughly equal to the number of areas shown.
#' @examples
#' # demonstration of the three widgets
#' plot.new()
#' par(usr=c(0,3,0,3))
#' x <- c(7, 5, 11)
#' col <- tmodPal()
#' b <- "black"
#' simpleRug(0.5, 1.5, 0.8, 0.8, v=x, col=col, border=b)
#' simplePie(1.5, 1.5, 0.8, 0.8, v=x, col=col, border=b)
#' simpleBoxpie(2.5, 1.5, 0.8, 0.8, v=x, col=col, border=b)
#'
#' # using pie as plotting symbol
#' plot(NULL, xlim=1:2, ylim=1:2, xlab="", ylab="")
#' col <- c("#cc000099", "#0000cc99")
#' for(i in 1:125) { 
#'   x <- runif(1) + 1 
#'   y <- runif(1) + 1
#'   simplePie( x, y, 0.05, 0.05, c(x,y), col)
#' }
#'
#' # square filled with box pies
#' n <- 10 
#' w <- h <- 1/(n+1)
#' plot.new()
#' for(i in 1:n) for(j in 1:n) 
#'  simpleBoxpie(1/n*(i-1/2), 1/n*(j-1/2), w, h, 
#'  v=runif(3), col=tmodPal())
#' @export
simplePie <- function(x, y, w, h, v, col, res=100, border=NA) {
  dev.hold()
  on.exit(dev.flush())

  v <- c(0, cumsum(v)/sum(v))
  dv <- diff(v)
  nv <- length(dv)

  nn <- res
  w <- w / 2
  h <- h / 2

  a2xy <- function(a) {
    t <- pi * (0.5  - 2 * a)
    list( x=x + w * cos(t), y=y + h * sin(t) )
  }

  for(i in 1:nv) {
    ni <- max(2, floor(nn * dv[i])) 
    po <- a2xy(seq.int(v[i], v[i+1], length.out=ni))
    polygon(c(x, po$x), c(y, po$y), col=col[i], border=NA)
  }

  if(!(is.null(border) || is.na(border))) {
    po <- a2xy(seq.int(0, 1, length.out=ni))
    lines(c(po$x, po$x[1]), c(po$y, po$y[1]), col=border)
  }

  invisible(NULL)
}


#' @rdname simplePie
#' @export
simpleRug <- function(x, y, w, h, v, col, border=NULL) {
  dev.hold()
  on.exit(dev.flush())

  x <- x - w/2
  y <- y - h/2

  v <- x + c(0, cumsum(v)/sum(v)) * w
  nv <- length(v)

  rect( v[-nv], rep(y, nv-1), 
        v[-1],  rep(y+h, nv-1), col=col, border=NA)

  if(!(is.null(border) || is.na(border))) {
    rect(x, y, x + w, y + h, col=NA, border=border)
  }
}



#' @rdname simplePie
#' @export
simpleBoxpie <- function(x, y, w, h, v, col, border=NA, grid=3) {
  dev.hold()
  on.exit(dev.flush())

  x <- x - w/2
  y <- y - h/2

  v <- grid^2 * v / sum(v)
  cv <- c(0, cumsum(v))

  left <- TRUE

  for(i in 1:length(v)) {
    ret <- .boxpie.calculate.polygon(cv[i], cv[i+1], grid=grid, left=left)
    xp <- ret$vec[,1] * w / grid + x
    yp <- ret$vec[,2] * h / grid + y

    polygon(xp, yp, col=col[i], border=NA)
    left <- ret$left
  }

  if(!(is.null(border) || is.na(border))) 
    rect(x, y, x + w, y + h)

  invisible(NULL)
}


## joins a rectangle (defined by x1, x2, y1, y2) with a polygon (defined by vec)
## vec is a two-column matrix with coordinates of the polygon
.boxpie.addlines <- function(vec, x, y) {
  t <- cbind(rep(x, each=2), c(y, rev(y)))
  return(rbind(t[1:2,], vec, t[3:4,]))
}

## calculates the coordinates of a polygon corresponding to the area of
## (end-start) on a surface grid x grid. The drawing happens in three
## parts:
##                                +------------+                                      
##              starting block -> |    ########|                
##                     /->        |############|                
##         full rows -+---->      |############|                
##                     \->        |############|                
##              ending block ->   |        ####|              
##                                |            |                
##                                +------------+                
##                                                                       
## the "left" parameter decides whether the polygon should be sticking to
## the left, or to the right side of the grid.
## The function returns a list with two elements:
##
## vec  -- a matrix with two columns encoding the coordinates of the
##         resulting polygon
## left -- a T/F value indicating whether the *next* polygon in row should
##         stick to the left, or to the right side of the grid
##
.boxpie.calculate.polygon <- function(start, end, left=TRUE, grid=3) {

  dx <- end - start

  # start drawing from top to bottom (first row to draw is the top-most
  # row)
  row0     <- grid - start %/% grid 
  leftover <- start %% grid 

  # vec will hold the polygon coordinates
  vec <- NULL

  ret.left <- !left

  # starting block
  if(leftover > 0) {
    startbloclen <- min(grid - leftover, dx)

    xx <- c(0, startbloclen) + leftover
    if(left) xx <- grid - rev(xx) 

    vec <- .boxpie.addlines(vec, xx, c(row0 - 1, row0))

    # if the whole polygon fits in that remainder, and if there is still
    # some place left, the empty place will be on the same side as before,
    # so the next polygon should stick to the same side as that one
    if(grid - leftover > dx) ret.left <- left

    row0 <- row0 - 1
    dx <- dx - startbloclen
  } 

  # rows completly filled
  fullrows <- dx %/% grid

  if(fullrows > 0) {
    vec <- .boxpie.addlines(vec, c(0, grid), c(row0 - fullrows, row0))
    row0 <- row0 - fullrows
  }

  # remaining bit to be filled out
  endrest  <- dx %% grid

  if(endrest > 0) {
    xx <- c(0, endrest)
    if(!left) xx <- grid - rev(xx)
    vec <- .boxpie.addlines(vec, xx, c(row0 - 1, row0))
  }

  return(list(vec=vec, left=ret.left))
}



