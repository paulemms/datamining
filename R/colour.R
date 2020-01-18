# Routines for manipulating colors and making color plots
# Tom Minka

# intended for white background
gray.colors <- function(n) gray(rev(0:(n-1))/n)
default.colors <- function(n) ((0:(n-1))%%6+1)
# good for objects which are being outlined (e.g. pies)
default.colors.w <- function(n) c(8,seq(len=n-1))
# these are from
# http://www.personal.psu.edu/faculty/c/a/cab38/ColorBrewerBeta2.html
OrRd.colors <- function(n) {
  if(n == 7) c("#fef0d9","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#99000d")
  else if(n == 6) c("#fef0d9","#fdd49e","#fdbb84","#fc8d59","#e34a33","#b2000d")
  else if(n == 5) c("#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b2000d")
  else if(n == 4) c("#fef0d9","#fdcc8a","#fc8d59","#d7301f")
  else if(n == 3) c("#fee8c8","#fdbb84","#e34a33")
  else if(n == 2) c("#fdbb84","#e34a33")
  else stop("can't make that many levels")
}
YlGnBu.colors <- function(n) {
  if(n == 8) c("#ffffd9","#edf8b0","#c7e9b4","#7fcdbb","#41b6c3","#1d91c0","#225ea8","#0c2c83")
  else if(n == 7) c("#ffffcc","#c7e9b4","#7fcdbb","#41b6c3","#1d91c0","#386cb0","#0c2c83")
  else if(n == 6) c("#ffffcc","#c7e9b4","#7fcdbb","#41b6c3","#2c7fb8","#253494")
  else if(n == 5) c("#ffffcc","#a1dab4","#41b6c3","#2c7fb8","#253494")
  else if(n == 4) c("#ffffcc","#a1dab4","#41b6c3","#386cb0")
  else if(n == 3) c("#edf8b0","#7fcdbb","#2c7fb8")
  else if(n == 2) c("yellow","blue")
  else stop("can't make that many levels")
}
# diverging
BrBg.colors <- function(n) {
  if(n == 7) c("#8c510a","#d8b365","#f6e8c3","#f5f5f5","#c7eae5","#5ab4ac","#01665e")
  else if(n == 6) c("#8c510a","#d8b365","#f6e8c3","#c7eae5","#5ab4ac","#01665e")
  else if(n == 5) c("#a6611a","#dfc27d","#f5f5f5","#80cdc1","#018571")
  else if(n == 4) c("#a6611a","#dfc27d","#80cdc1","#018571")
  else if(n == 3) c("#d8b365","#f5f5f5","#5ab4ac")
  else if(n == 2) c("#d8b365","#5ab4ac")
  else stop("can't make that many levels")
}
RYB.colors <- function(n) {
  if(n == 7) c("#d73d29","#fc8d59","#fee090","#ffffbf","#e0f3f8","#91bfdb","#4575b4")
  else if(n == 6) c("#d73d29","#fc8d59","#fee090","#e0f3f8","#91bfdb","#4575b4")
  else if(n == 5) c("#ca373b","#fdae61","#ffffbf","#abd9e9","#2c7bb6")
  else if(n == 4) c("#ca373b","#fdae61","#abd9e9","#2c7bb6")
  else if(n == 3) c("#fc8d59","#ffffbf","#91bfdb")
  else if(n == 2) c("red","blue")
  else stop("can't make that many levels")
}

# sequential WYRK
YR.colors <- function(k,lt=0.97,dk=0.03) {
  r <- c(lt, rep(lt,k),
         rep(lt,k),
         seq(lt,dk,len=k+1)[-1])
  g <- c(lt, rep(lt,k),
         seq(lt,dk,len=k+1)[-1],
         rep(dk,k))
  b <- c(lt, seq(lt,dk,len=k+1)[-1],
         rep(dk,k),
         rep(dk,k))
  i <- seq(1,length(r),len=k)
  rgb(r[i],g[i],b[i])
}

# diverging
diverging.colors <- function(n,hmax=0,hmin=(0.5+hmax)%%1,v=1) {
  if(n %% 2 == 1) {
    # odd case
    k <- n%/%2
    x <- seq(1,0,len=k+1)
    c(hsv(h = hmin, s = x, v),
      hsv(h = hmax, s = rev(x[1:k]), v))
  } else {
    # even case
    k <- n%/%2
    x <- seq(1,-1,len=n)[1:k]
    c(hsv(h = hmin, s = x, v),
      hsv(h = hmax, s = rev(x), v))
  }
}
RC.colors <- diverging.colors
GM.colors <- function(n) diverging.colors(n,2/6)

# returns a new vector of darker colors
darker <- function(col,scale=0.8) {
  crgb <- col2rgb(col)/255*scale
  if(T) {
    # handle black specially
    i <- which(apply(crgb,2,sum)==0)
    if(length(i)>0) {
      crgb[,i] <- rep.col(c(1,1,1)/10,length(i))
    }
  }
  rgb(crgb[1,],crgb[2,],crgb[3,])
}

# note: pch 21 looks like pch 1, but is actually a filled circle
#plot(y,pch=21,bg=hsv(6/12,1,1))




#' Add a thin key to the top of a plot
#'
#' Similar to \code{\link{legend}} but takes up less space.
#'
#'
#' @param col a vector of colors, one for each group.  Assumed black if none
#' given.
#' @param pch a vector of symbols, one for each group.  If less than two
#' distinct symbols are given, no symbols are put in the key.
#' @param labels a character vector giving the label of each group.
#' @param breaks a numeric vector of boundaries between groups.
#' @param digits number of digits to use for \code{breaks}.
#' @param cex character expansion factor.
#' @return A thin bar is placed at the top of the plot, just inside the
#' plotting area.  The bar is divided into equal-sized lengths, colored
#' according to \code{col}.  If \code{pch} is given, that symbol is plotted in
#' the center of each segment.  If \code{labels} is given, each label is placed
#' at the center of the corresponding segment, just above the plotting area
#' (using \code{\link{mtext}}).  If \code{breaks} is given, the boundary
#' between segments is also labeled with a number.
#' @author Tom Minka
#' @seealso \code{\link{color.plot}}
#' @examples
#'
#' data(iris)
#' y = as.numeric(iris$Species)
#' plot(Sepal.Width ~ Sepal.Length, iris,col=y,pch=y)
#' color.key(1:3,1:3,levels(iris$Species))
#'
color.key <- function(col=NULL,pch=NULL,labels=NULL,breaks=NULL,digits=2,cex=0.75*par("cex")) {
  # don't draw symbols in the key if there is only one type
  if(length(unique(pch)) == 1) pch <- NULL
  if(is.null(col)) col <- rep(1,length(pch))
  rx <- range(par("usr")[1:2])
  ry <- range(par("usr")[3:4])
  if(T && !is.null(pch)) {
    # key is par("csi")*cex inches
    ry[1] <- ry[2] - diff(ry)/par("pin")[2]*par("csi")*cex
  } else {
    # key is 0.06 inches high
    ry[1] <- ry[2] - diff(ry)/par("pin")[2]*0.06
  }
  if(!is.null(labels)) n <- length(labels)
  else if(!is.null(breaks)) n <- length(breaks)-1
  else n <- length(col)
  i <- seq(0,1,len=n+1)
  x1 <- rx[1]+i[-length(i)]*diff(rx)
  x2 <- rx[1]+i[-1]*diff(rx)
  if(n == 0) stop("n == 0")
  col <- rep(col,ceiling(n/length(col)))[1:n]
  rect(x1, ry[1], x2, ry[2], col=col, border=NA)
  if(!is.null(pch))
    points((x1+x2)/2, rep(mean(ry),n), col=8, pch=pch, cex=0.75)
  if(!is.null(labels))
    mtext(labels,at=(x1+x2)/2,cex=cex)
  if(!is.null(breaks)) {
    x <- rx[1] + i*diff(rx)
    # same as format.pretty
    labels <- prettyNum(formatC(breaks,format="fg",digits=digits),big.mark=",")
    #labels = format(breaks,digits=digits,trim=T)
    mtext(labels[1],at=x[1],cex=cex,adj=0)
    mtext(labels[n+1],at=x[n+1],cex=cex,adj=1)
    mtext(labels[2:n],at=x[2:n],cex=cex,adj=0.5)
  }
}

# extends r=c(min,max) by 4%, to mimic plotting limits
extend.pct <- function(r,pct=0.04) {
  m <- diff(r)*pct
  if(m == 0) m <- 0.432*abs(r[1])
  if(m == 0) m <- 1.08
  c(r[1]-m,r[2]+m)
}
extend.inches <- function(r,upper=0,lower=0) {
  # extends the given range by a specified number of inches of each side.
  # upper and lower are actually inches/plot.size
  # constraints:
  # (r.new[2] - r[2])/diff(r.new) = upper
  # (r[1] - r.new[1])/diff(r.new) = lower
  r.new <- r
  a = lower/(1-upper)
  r.new[1] <- (a*r[2] - r[1])/(a-1)
  r.new[2] <- (r[2] - upper*r.new[1])/(1-upper)
  r.new
}

# returns the minimum user limits for a plotting region of size pin
# which would enclose the given labels at the given positions
# with the justification adj at rotation srt.
range.text <- function(x,y,labels,cex=NULL,adj=NULL,srt=0,
                       pin=par("pin")) {
  bad = (is.na(x) | is.na(y) | is.na(labels))
  x = x[!bad]; y = y[!bad]; labels = labels[!bad]
  lab.w <- strwidth(labels,"inches",cex=cex)
  lab.h <- strheight(labels,"inches",cex=cex)

  # bounding box of text, columns are (x1,y1)(x2,y1)(x2,y2)(x1,y2)
  n <- length(labels)
  ix <- seq(1,8,by=2)
  iy <- seq(2,8,by=2)
  hull <- array(0,c(n,8))
  hull[,c(1,7,2,4)] <- rep(0,n)
  hull[,c(3,5,6,8)] <- rep(1,n)
  # put adj point at origin and scale
  if(is.null(adj)) adj <- c(0.5,0.5)
  if(length(adj) == 1) adj <- c(adj,0.5)
  if(is.null(nrow(adj))) {
    hull[,ix] <- (hull[,ix] - adj[1])*rep.col(lab.w,4)
    hull[,iy] <- (hull[,iy] - adj[2])*rep.col(lab.h,4)
  } else {
    # adj is an n by 2 matrix
    adj = adj[!bad,]
    if(nrow(adj) != n) stop("nrow(adj) != length(x)")
    hull[,ix] <- (hull[,ix] - rep.col(adj[,1],4))*rep.col(lab.w,4)
    hull[,iy] <- (hull[,iy] - rep.col(adj[,2],4))*rep.col(lab.h,4)
  }
  # rotate
  srt <- srt/180*pi
  cr <- cos(srt)
  sr <- sin(srt)
  R <- array(c(cr,-sr,sr,cr),c(2,2))
  R <- diag(4) %x% R
  hull <- hull %*% R
  bbox <- list()
  bbox$left <- apply(hull[,ix],1,min)/pin[1]
  bbox$right <- apply(hull[,ix],1,max)/pin[1]
  bbox$bottom <- apply(hull[,iy],1,min)/pin[2]
  bbox$top <- apply(hull[,iy],1,max)/pin[2]

  # solve constraints
  xmid <- mean(range(x))
  ymid <- mean(range(y))
  xlim <- c()
  # these come from the constraints
  # (x-xmin)/(xmax-xmin) + left > 0
  # (x-xmin)/(xmax-xmin) + right < 1
  # where xmax is fixed at 2*xmid - xmin
  min1 <- min((x + 2*bbox$left*xmid)/(1+2*bbox$left))
  min2 <- min((2*(1-bbox$right)*xmid - x)/(1-2*bbox$right))
  xlim[1] <- min(min1,min2)
  xlim[2] <- 2*xmid - xlim[1]
  ylim <- c()
  min1 <- min((y + 2*bbox$bottom*ymid)/(1+2*bbox$bottom))
  min2 <- min((2*(1-bbox$top)*ymid - y)/(1-2*bbox$top))
  ylim[1] <- min(min1,min2)
  ylim[2] <- 2*ymid - ylim[1]
  c(xlim,ylim)
}

range.aspect <- function(asp,lim=par("usr"),pin=par("pin")) {
  # returns c(xlim,ylim), strictly bigger than lim, which satisfy asp
  # simulates the effect of plot.window(asp)
  xlim <- lim[1:2]
  ylim <- lim[3:4]
  dx <- diff(xlim)/pin[1]
  dy <- diff(ylim)/pin[2]
  if(dy > asp*dx) {
    dx <- dy/asp*pin[1]
    xlim <- mean(xlim) + c(-dx/2,dx/2)
  } else {
    dy <- asp*dx*pin[2]
    ylim <- mean(ylim) + c(-dy/2,dy/2)
  }
  c(xlim,ylim)
}

default.pch <- c(1,41,4,20,35,38,2,17)
# use for the screen
sequential.pch <- c(20,41,4,1,38,35,11,19)
# use for printing
sequential.pch.paper <- c(45,41,4,1,38,35,11,19)



#' Plot subgroups as colored points
#'
#' Like \code{\link{plot}} and \code{\link{text.plot}} but colors according to
#' a third variable.
#'
#' Each (x,y) point is plotted with a color determined by \code{z}.  If
#' \code{z} is a factor, each factor level is in a different color.  If
#' \code{z} is numeric, then it is color-coded by assigning each quantile a
#' different color.
#'
#' Despite the name, this function can also make plots using different symbols
#' instead of colors.  For example, if \code{color.palette=1} then all points
#' will be black but use different symbols.
#'
#' When \code{labels != NULL}, the result is equivalent to
#' \code{\link{text.plot}} with colors.
#'
#' If \code{key=TRUE} and there is more than one color or symbol on the plot, a
#' key is displayed at the top of the figure.
#'
#' @aliases color.plot color.plot.default
#' @param x,y numeric vectors of the same length.
#' @param z a numeric vector or factor, same length as \code{x} and \code{y}.
#' @param labels a character vector of labels, same length as \code{z}.  If
#' NULL, cases are plotted as points (default).
#' @param axes If FALSE, no axes are plotted (but there may still be a key).
#' @param key If FALSE, no key is plotted.
#' @param add If FALSE, a new plot is made.  If TRUE, points or labels are
#' added to the existing plot.
#' @param nlevels an integer.  If \code{z} is numeric, it is color-coded using
#' this many levels.  (If \code{z} is a factor, color-coding follows the factor
#' levels.)
#' @param color.palette a vector of colors, arbitrary length, or a function
#' with integer argument which generates a vector of colors (e.g.
#' \code{\link{YlGnBu.colors}}).  Used if \code{col} is not specified.  If
#' shorter than the number of levels, colors will be recycled and the plotting
#' symbol will change.
#' @param col a vector of colors, as in a call to \code{\link{plot}}.  Used to
#' explicitly set the color of each point.
#' @param pch.palette a vector of plotting symbols, arbitrary length.  If
#' \code{labels=NULL}, the plot symbol will rotate through these when there
#' aren't enough colors.
#' @param pch a vector of plotting symbols, as in a call to \code{\link{plot}}.
#' Used to explicitly set the symbol of each point.
#' @param digits the number of digits to use in the color key when \code{z} is
#' numeric.
#' @param bg the background color
#' @param mar figure margins (see \code{par}).  If not specified, appropriate
#' margins will be chosen automatically, which will not necessarily match the
#' current value of \code{par("mar")}.
#' @return A plot is produced.
#' @note This function sets the figure margins permanently, so that you can
#' draw on the color plot.  Unfortunately, this also means future plots will
#' use the same margins, until you change them with \code{par("mar")}.
#' @author Tom Minka
#' @seealso
#' \code{\link{color.plot.data.frame}},\code{\link{color.plot.loess}},\code{\link{color.plot.glm}},\code{\link{color.plot.knn}},\code{\link{color.plot.tree}},\code{\link{YlGnBu.colors}}
#' @examples
#'
#' # See the examples for color.plot.data.frame
#'
color.plot <- function(object, ...) UseMethod("color.plot")
color.plot.formula <- function(formula,data=parent.frame(),...) {
  x <- model.frame.default(formula,data,na.action=na.pass)
  # put response at end (my preferred ordering)
  x = cbind(x[-1],x[1])
  color.plot.data.frame(x,...)
}


#' Plot cases as colored points
#'
#' Like \code{\link{plot}} and \code{\link{text.plot}} but colors according to
#' the response variable.
#'
#' Calls \code{\link{color.plot.default}} with \code{x} as the first predictor
#' in \code{data}, \code{y} as the second predictor, and \code{z} as the
#' response.  To get a different predictor/response breakdown than the default,
#' use \code{color.plot(formula, x, ...)}, which is shorthand for
#' \code{color.plot(model.frame(formula, x), ...)}.
#'
#' Each case is plotted with a color determined by the response.  If the
#' response is a factor, each factor level is in a different color.  If the
#' response is numeric, then it is color-coded by assigning each quantile a
#' different color.
#'
#' @aliases color.plot.data.frame color.plot.formula
#' @param data a data frame.
#' @param formula a formula specifying a response and two predictors from
#' \code{data}
#' @param labels If NULL, cases are plotted as points.  If T, cases are plotted
#' as labels, according to \code{rownames}.
#' @param ... Extra arguments passed to \code{\link{color.plot.default}}.
#' @author Tom Minka
#' @seealso \code{\link{color.plot.default}}
#' @examples
#'
#' data(iris)
#' color.plot(iris)
#' color.plot(Species ~ Petal.Length + Petal.Width, iris)
#' color.plot(Species ~ Petal.Length, iris)
#' color.plot(Species ~ Petal.Length, iris,jitter=T)
#' color.plot(iris, col=1)
#' color.plot(iris, col=c(1,2))
#'
#' data(state)
#' x <- data.frame(state.x77)
#' color.plot(Murder ~ Frost + Illiteracy, x, labels=T, cex=0.5)
#'
color.plot.data.frame <- function(x,z,zlab=NULL,xlab=NULL,ylab=NULL,labels=F,...) {
  if(missing(z)) {
    pred <- predictor.terms(x)
    resp <- response.var(x)
    z <- x[[resp]]
    if(is.null(zlab)) zlab <- resp
  } else {
    pred <- colnames(x)
    if(is.null(zlab)) zlab <- deparse(substitute(z))
  }
  if(length(pred) < 2) stop("must have two predictors to plot")
  if(is.logical(labels)) {
    labels <- if(labels) rownames(x) else NULL
  }
  if(is.null(xlab)) xlab <- pred[1]
  if(is.null(ylab)) ylab <- pred[2]
  color.plot.default(x[[pred[1]]],x[[pred[2]]],z,labels=labels,
                     xlab=xlab,ylab=ylab,zlab=zlab,...)
}
color.plot.default <- function(x,y,z,labels=NULL,data=parent.frame(),
                               xlab=NULL,ylab=NULL,zlab=NULL,main="",
                               xlim=NULL,ylim=NULL,
                               axes=T,key=T,add=F,nlevels=4,pretty=T,
                               color.palette=NULL,col,
                               pch.palette=NULL,pch,
                               jitter=F,digits=2,mar,bg=NULL,
                               cex=par("cex"),yaxt="s",...) {
  # continous responses are automatically quantized into nlevels levels, with
  # colors given by col.
  # jitter is only a suggestion.  jittering is only done when necessary.
  if(is.null(xlab)) xlab <- deparse(substitute(x))
  if(is.null(ylab)) ylab <- deparse(substitute(y))
  if(is.null(zlab)) zlab <- deparse(substitute(z))
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  z <- eval(substitute(z),data)
  # treat character as factor
  if(mode(x) == "character") x <- factor(x,levels=unique(x))
  if(mode(y) == "character") y <- factor(y,levels=unique(x))
  if(is.numeric(z)) {
    b <- break.quantile(z,nlevels,pretty=pretty)
    f <- rcut(z,b)
    if(!any(is.na(z)) && any(is.na(f)))
      stop("internal error: cut introduced NAs")
    if(is.null(color.palette)) {
      color.palette = YlGnBu.colors
      #color.palette = YR.colors
      # pick bg color using color.cone(c(YlGnBu.colors(8),grey(0.65)))
      # don't use grey(0.6) or grey(0.7)
      if(is.null(bg)) bg = grey(0.65)
    }
    if(is.null(pch.palette)) {
      pch.palette = if(mode(color.palette) == "function") 16
      else sequential.pch
    }
  } else {
    f <- factor(z)
    if(is.null(color.palette)) color.palette = default.colors
    if(is.null(pch.palette)) {
      pch.palette = if(mode(color.palette) == "function") 16
      else default.pch
    }
  }
  nlevels <- length(levels(f))
  if(nlevels == 0) stop("0 levels in factor")
  # lev.col is the color of each level
  # lev.pch is the symbol of each level
  # if there are too few colors, recycle them and change pch
  lev <- 1:nlevels
  # col is the color of the individual points
  q <- as.numeric(f)
  if(missing(col)) {
    if(mode(color.palette) == "function")
      color.palette <- color.palette(nlevels)
    lev.col <- color.palette[((lev-1) %% length(color.palette))+1]
    col <- lev.col[q]
  } else {
    lev.col <- NULL
  }
  if(missing(pch)) {
    lev.pch <- pch.palette[(((lev-1) %/% length(color.palette)) %% length(pch.palette))+1]
    pch <- lev.pch[q]
  } else {
    lev.pch <- NULL
  }

  # setup plot window
  if(!add) {
    if(missing(mar)) mar = auto.mar(main=zlab,axes=axes,yaxt=yaxt,key=key)
    # change mar permanently so we can add things on top
    par(mar=mar)
    if(!is.null(bg) && bg != par("bg")) {
      opar <- par(bg=bg)
      on.exit(par(opar))
    }
  }
  # xlim,ylim needed for jittering and label placement
  if(add) {
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]
  } else {
    if(is.null(labels)) {
      if(is.null(xlim)) xlim <- range(as.numeric(x),na.rm=T)
      if(is.null(ylim)) ylim <- range(as.numeric(y),na.rm=T)
      if(key) {
        # make room at the top of the plot for the key
        if(length(unique(lev.pch)) == 1) {
          a = 0.06/par("pin")[2]
        } else {
          a = par("csi")*cex/par("pin")[2]
        }
        ylim = extend.inches(ylim,a)
      }
    } else {
      plot.new()
      lab.lim <- range.text(x,y,labels,cex)
      # make room for labels
      if(is.null(xlim)) xlim <- lab.lim[1:2]
      if(is.null(ylim)) ylim <- lab.lim[3:4]
    }
  }
  if(jitter) {
    # this test needs to be less strict.  need to consider actual figure size.
    if(length(levels(factor(x))) < length(x)) {
      cat("jittering",xlab,"\n")
      jitter1 <- (runif(length(x))-0.5)*diff(extend.pct(xlim))/100
      x <- x+jitter1
    }
    if(length(levels(factor(y))) < length(y)) {
      cat("jittering",ylab,"\n")
      jitter2 <- (runif(length(y))-0.5)*diff(extend.pct(ylim))/100
      y <- y+jitter2
    }
  }
  if(add) {
    if(!is.null(labels)) {
      return(text(x, y, labels=labels, col=col, cex=cex, ...))
    } else {
      return(points(x, y, pch=pch, col=col, cex=cex, ...))
    }
  }
  # from here on, add=F
  if(!is.null(labels)) {
    plot.default(x, y, xlim=xlim, ylim=ylim, type="n",
                 xlab=xlab,ylab=ylab,main="",axes=axes,yaxt=yaxt,...)
    text(x, y, labels=labels, col=col, cex=cex, ...)
  } else {
    plot.default(x, y, xlim=xlim, ylim=ylim, xlab=xlab,ylab=ylab,main="",
                 pch=pch,col=col,axes=axes,yaxt=yaxt,cex=cex,...)
  }
  if(key) {
    if(is.numeric(z))
      color.key(lev.col,lev.pch,breaks=b,digits=digits)
    else
      color.key(lev.col,lev.pch,levels(f))
    if(axes) title(zlab,line=1)
  }
  else if(axes) title(zlab)
}
# tests:
# color.plot(runif(10),runif(10),c(rep("a",5),rep("b",5)))
# color.plot(runif(10),runif(10),c(rep("a",5),rep("b",5)),rep("foo",10))

# this part same as color.plot
text.plot <- function(object,...) UseMethod("text.plot")
text.plot.formula <- function(formula,data=parent.frame(),...) {
  x <- model.frame.default(formula,data,na.action=na.pass)
  text.plot.data.frame(x,...)
}
text.plot.data.frame <- function(x,labels=TRUE,xlab,ylab,...) {
  pred <- predictor.vars(x)[1]
  resp <- response.var(x)
  x1 <- x[[pred]]
  x2 <- x[[resp]]
  if(identical(labels,TRUE)) labels <- rownames(x)
  if(missing(xlab)) xlab = pred
  if(missing(ylab)) ylab = resp
  text.plot.default(x1,x2,labels,xlab=xlab,ylab=ylab,...)
}
text.plot.default <- function(x,y,labels=NULL,xlab,ylab,xlim=NULL,ylim=NULL,
                              main="",asp=NULL,move=F,
                              cex=par("cex"),pch=par("pch"),
                              adj=NULL,srt=0,axes=T,...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(is.null(labels)) labels = names(x)
  if(is.null(labels)) labels = names(y)
  if(length(labels) != length(x)) stop("labels is wrong length")
  plot.new()
  lim <- range.text(x,y,labels,cex,adj,srt)
  # make room for labels
  if(is.null(xlim)) xlim <- lim[1:2]
  if(is.null(ylim)) ylim <- lim[3:4]
  if(F) {
    plot(x,y,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,axes=axes,
         main=main,...,type="n")
  } else {
    plot.window(xlim,ylim,asp=asp)
    if(axes) { axis(1,...); axis(2,...); box() }
    title(main,xlab=xlab,ylab=ylab)
  }
  if(move) {
    points(x,y,cex=cex,pch=pch,...)
    x = inchx(x)
    y = inchy(y)
    w = strwidth(labels,units="inch",cex=cex)
    h = strheight(labels,units="inch",cex=cex)
    n = length(x)
    s = inch.symbol(pch=pch,cex=cex) + 0.1
    w2 = c(w,rep(s,n))
    h2 = c(h,rep(s,n))
    x2 = c(x,x)
    y2 = c(y,y)
    #rect(xinch(x2-w2/2),yinch(y2-h2/2),xinch(x2+w2/2),yinch(y2+h2/2))
    cost = c(rep(1,n),rep(0,n))
    r = move.collisions2(x2,y2,w2,h2,cost)
    #rect(xinch(r[,1]-w2/2),yinch(r[,2]-h2/2),xinch(r[,1]+w2/2),yinch(r[,2]+h2/2))
    r = r[1:n,]
    x = xinch(r[,1])
    y = yinch(r[,2])
  }
  if(!is.null(nrow(adj))) {
    if(nrow(adj) != length(x)) stop("nrow(adj) != length(x)")
    # adj is an n by 2 matrix
    for(i in 1:nrow(adj)) {
      text(x[i],y[i],labels=labels[i],cex=cex,adj=drop(adj[i,]),srt=srt,...)
    }
  } else {
    text(x,y,labels=labels,cex=cex,adj=adj,srt=srt,...)
  }
}
# test: text.plot(runif(10),runif(10),rep("a",10))



#' Contour plot of a regression surface
#'
#'
#' The regression surface is evaluated at all points on a grid, clipped values
#' are set to \code{NA}, and \code{contour.plot} is used to plot the contours.
#'
#' If \code{add=FALSE}, the data is plotted on top using
#' \code{color.plot.data.frame}.
#'
#' @param object a \code{\link{loess}} object.
#' @param data data to use instead of \code{model.frame(object)}.
#' @param res resolution of the sampling grid in each direction.
#' @param fill passed to \code{\link{contour.plot}}.
#' @param add If \code{TRUE}, the contours are added to an existing plot.
#' Otherwise a new plot is created.
#' @param clip a polygon over which the surface is to be defined.  Possible
#' values are \code{FALSE} (no clipping), \code{TRUE} (clip to the convex hull
#' of the data), or a matrix with two columns specifying (x,y) coordinates.
#' @param ... extra arguments to \code{\link{color.plot.data.frame}}
#' @return A plot is produced.
#' @author Tom Minka
#' @seealso \code{\link{contour.plot}}
#' @examples
#'
#' data(Housing)
#' fit = loess(Price ~ Rooms + Low.Status, Housing)
#' color.plot(fit)
#'
color.plot.loess <- function(object,x,res=50,fill=F,add=fill,
                             clip=T,xlab=NULL,ylab=NULL,zlab=NULL, ...) {
  if(missing(x)) x <- model.frame(object)
  else x <- model.frame(terms(object),x)
  resp <- response.var(object)
  pred <- predictor.vars(object)
  if(length(pred) < 2) stop("must have 2 predictors")
  #return(plot.loess(object,x,add=add,lwd=lwd,...))

  if(is.null(xlab)) xlab = pred[1]
  if(is.null(ylab)) ylab = pred[2]
  if(is.null(zlab)) zlab = resp
  if(!add && !fill) color.plot.data.frame(x,xlab=xlab,ylab=ylab,zlab=zlab,...)
  # use the par options set by color.plot.data.frame

  # x is only used to get plotting range
  # this is ugly, but otherwise recursively calling color.plot.data.frame is hard
  if(is.null(...$xlim)) {
    if(add || !fill) xlim <- range(par("usr")[1:2])
    else xlim <- range(x[[pred[1]]])
  } else xlim <- ...$xlim
  x1 <- seq(xlim[1],xlim[2],length=res)
  if(is.null(...$ylim)) {
    if(add || !fill) ylim <- range(par("usr")[3:4])
    else ylim <- range(x[[pred[2]]])
  } else ylim <- ...$ylim
  x2 <- seq(ylim[1],ylim[2],length=res)
  xt <- expand.grid(x1,x2)
  names(xt) <- pred[1:2]
  z <- predict(object,xt)
  if(clip != F && length(pred) >= 2) {
    if(clip == T) {
      i <- chull(x[pred])
      clip <- x[pred][i,]
    }
    z[!in.polygon(clip,xt)] <- NA
  }
  dim(z) <- c(length(x1),length(x2))
  contour.plot(x1,x2,z,fill=fill,level.data=x[[resp]],
               xlab=xlab,ylab=ylab,zlab=zlab,add=(add || !fill),...)
  if(!add && fill) color.plot.data.frame(x,col=1,add=T,...)
}



#' Display contours
#'
#' Create a filled or unfilled contour plot.
#'
#' This is a wrapper function for \code{\link{contour}} and
#' \code{\link{filled.contour}} whose main purpose is to provide a uniform
#' interface and provide a decent automatic choice of levels.
#'
#' If \code{equal=FALSE}, the levels are chosen according to the equal-count
#' algorithm of \code{\link{break.quantile}}.
#'
#' @param x,y locations at which the values in \code{z} are measured.
#' @param z a matrix containing the values to be plotted (NAs are allowed).
#' @param fill If \code{TRUE}, makes a filled contour plot.
#' @param levels numeric vector of levels at which to draw contour lines.  The
#' next few arguments only apply to the case when \code{levels==NULL}.
#' @param nlevels The number of contour lines to draw.
#' @param level.data numeric vector to use for computing \code{levels}.
#' @param zlim The range that the levels should cover.  Defaults to the range
#' of \code{level.data}.
#' @param equal If \code{TRUE}, the levels are equally spaced over \code{zlim}.
#' Otherwise they match the quantiles of \code{level.data}.
#' @param pretty If \code{TRUE}, the levels are chosen to be round numbers.
#' @param lwd width of the contour lines.
#' @param drawlabels If \code{TRUE}, the contour lines are labeled with their
#' level.
#' @param color.palette a vector of colors, or a function which takes a number
#' and returns a vector of that length.
#' @param bg the background color.
#' @param key If \code{TRUE}, a color key is drawn at the top of the plot.  If
#' \code{key==2}, the key of \code{filled.contour} is used.
#' @param main the title of the plot.
#' @param ... extra arguments passed to \code{contour} or
#' \code{filled.contour}.
#' @author Tom Minka
#' @seealso \code{\link{contour}},\code{\link{filled.contour}},
#' \code{\link{color.plot.loess}}
#' @examples
#'
#' # see the examples for color.plot.loess
#'
contour.plot <- function(x,y,z,fill=F,
                         levels=NULL,nlevels=if(is.null(levels)) 4 else length(levels),
                         level.data=z,
                         zlim,equal=F,pretty=T,lwd=2,drawlabels=F,
                         color.palette=YlGnBu.colors,bg=grey(0.65),
                         key=T,zlab="",...) {
  if(is.function(color.palette)) {
    color.palette <- color.palette(nlevels)
  }
  if(missing(zlim)) {
    zlim = range(level.data,na.rm=T)
    z[z < zlim[1]] <- zlim[1]
    z[z > zlim[2]] <- zlim[2]
  }
  if(is.null(levels)) {
    if(equal)
      levels <- seq(zlim[1],zlim[2],len=nlevels+1)
    else {
      r <- unique(level.data)
      r <- r[r >= zlim[1] & r <= zlim[2]]
      r <- c(r,zlim)
      levels <- break.quantile(r,nlevels,pretty=pretty)
    }
  }
  opar = par(bg=bg)
  on.exit(par(opar))
  if(fill) filled.contour(x,y,z,levels=levels,col=color.palette,main="",
                          key=if(key==2) T else F,...)
  else {
    color.palette <- c(color.palette,color.palette[length(color.palette)])
    contour(x,y,z,col=color.palette,levels=levels,drawlabels=drawlabels,lwd=lwd,main="",...)
  }
  if(key && (key != 2)) {
    color.key(color.palette,breaks=levels)
    title(zlab,line=1)
  }
  else title(zlab)
}

color.plot.rtable <- function(rt,main=response.var(rt),
                              xlab=predictor.vars(rt)[1],
                              ylab=predictor.vars(rt)[2],...) {
  x = as.numeric(rownames(rt))
  y = as.numeric(colnames(rt))
  contour.plot(x,y,rt,main=main,xlab=xlab,ylab=ylab,...)
}

color.plot.glm <- function(object,data,add=F,se=F,col=3,lwd=2,...) {
  if(missing(data)) data = model.frame(object)
  if(!add) color.plot.data.frame(data,...)
  # use the par options set by color.plot.data.frame

  w <- coef(object)
  z <- log(0.75/0.25)
  if(length(w) == 3) {
    abline(-w[1]/w[3],-w[2]/w[3],col=col,lwd=lwd,...)
    if(se) {
      abline((-z-w[1])/w[3],-w[2]/w[3],col=col,lty=2,lwd=lwd,...)
      abline((z-w[1])/w[3],-w[2]/w[3],col=col,lty=2,lwd=lwd,...)
    }
  } else {
    stop("use color.plot.knn")
    # two predictors with an interaction term
    # boundary is an ellipse
    xlim <- par("usr")[1:2]
    x <- seq(xlim[1],xlim[2],length=50)
    y <- (0-w[1]-w[2]*x)/(w[3]+w[4]*x)
    i <- (w[3]+w[4]*x < 0)
    lines(x[i],y[i],col=col,...)
    lines(x[!i],y[!i],col=col,...)
    if(se) {
      for(q in c(-z,z)) {
        y <- (q-w[1]-w[2]*x)/(w[3]+w[4]*x)
        lines(x[i],y[i],col=col,lty=2,...)
        lines(x[!i],y[!i],col=col,lty=2,...)
      }
    }
  }
}
color.plot.multinom <- color.plot.glm
color.plot.logistic = color.plot.glm

# Plot classification boundaries.
# If pclass is given, the probability of being in that class is plotted.
#
# levels and col are only used if pclass != NULL
#
color.plot.knn <- function(object,x,pclass=NULL,fill=F,res=50,
                           levels=c(0.5),add=F,axes=T,lwd=2,
                           color.palette=default.colors,col=3,drawlabels=F,...) {
  resp <- response.var(object)
  pred <- predictor.vars(object)
  if(missing(x)) x <- model.frame(object)
  else if(length(pred) == 1) {
    # take an unused variable from x
    pred[2] = setdiff(colnames(x),c(resp,pred))[1]
    fmla = formula(paste(resp,"~",paste(pred,collapse="+")))
    x = model.frame(fmla,x)
  } else x <- model.frame(terms(object),x)
  if(!add && !fill) color.plot.data.frame(x,color.palette=color.palette,...)
  # use the par options set by color.plot.data.frame

  # x is only used to get plotting range
  # this is ugly, but otherwise recursively calling color.plot.data.frame is hard
  if(is.null(...$xlim)) {
    if(add || !fill) xlim <- range(par("usr")[1:2])
    else xlim <- range(x[[pred[1]]])
  } else xlim <- ...$xlim
  x1 <- seq(xlim[1],xlim[2],length=res)
  # allow > 2 for derived features
  if(length(pred) >= 2) {
    if(is.null(...$ylim)) {
      if(add || !fill) ylim <- range(par("usr")[3:4])
      else ylim <- range(x[[pred[2]]])
    } else ylim <- ...$ylim
    x2 <- seq(ylim[1],ylim[2],length=res)
    xt <- expand.grid(x1,x2)
    names(xt) <- pred[1:2]
  } else {
    xt <- data.frame(x1)
  }
  prob.plot <- !is.null(pclass)
  if(prob.plot) {
    # z is the probability of a class
    if(inherits(object,"glm") || inherits(object,"gam") ||
       inherits(object,"nnet")) {
      if(inherits(object,"glm") || inherits(object,"gam")) {
        z <- predict(object,xt,type="response")
      } else {
        z = predict(object,xt)
      }
      y <- model.frame(object)[[resp]]
      if(is.character(pclass)) pclass <- pmatch(pclass, levels(y))
      if(pclass == 1) {
        z <- 1-z
        lev <- levels(y)[1]
        main <- paste("Pr(",resp," = ",lev,")",sep="")
      }
      else {
        if(pclass > length(levels(y))) stop("pclass exceeds the number of classes")
        if(length(levels(y)) == 2) {
          lev <- levels(y)[2]
          main <- paste("Pr(",resp," = ",lev,")",sep="")
        }
        else {
          lev <- levels(y)[1]
          main <- paste("Pr(",resp," != ",lev,")",sep="")
        }
      }
    } else {
      z <- predict(object,xt,type="vector")
      if(is.character(pclass)) pclass <- pmatch(pclass, colnames(z))
      lev <- colnames(z)[pclass]
      main <- paste("Pr(",resp," = ",lev,")",sep="")
      z <- z[,pclass]
    }
  }
  else {
    # z is the class
    if(inherits(object,"glm") || inherits(object,"gam")) {
      p <- predict(object,xt,type="response")
      y <- model.frame(object)[[resp]]
      z <- factor(levels(y)[as.numeric(p > 0.5)+1], levels=levels(y))
    }
    else {
      z <- factor(predict(object,xt,type="class"))
    }
  }
  dim(z) <- c(length(x1),length(x2))
  if(fill) {
    if(prob.plot) {
      filled.contour(x1,x2,z,
                     plot.axes={
                       title(main=main,xlab=pred[1],ylab=pred[2]);axis(1);axis(2);
                       if(!add) color.plot.data.frame(x,col=1,add=T,...)
                     },...)
    }
    else {
      lev <- levels(z)
      if(is.function(color.palette))
        col <- color.palette(length(lev))
      else
        col <- color.palette
      actual.lev <- levels(factor(z))
      icol <- col[lev %in% actual.lev]
      z <- as.numeric(z)
      dim(z) <- c(length(x1),length(x2))

      colorbar <- (length(col) > 1)
      if(!add) {
        # change mar permanently so we can add things on top
        if(axes) {
          mar = auto.mar()
          mar[3] = if(colorbar) 2 else 1
          par(mar=mar)
        } else {
          par(mar=c(0.05,0,0,0.1))
        }
        image(x1,x2,z,col=icol,main="",xlab=pred[1],ylab=pred[2])
        # assume a light color scheme (else should be lighter)
        dcol <- darker(col)
        color.plot.data.frame(x,add=T,col=dcol,...)
      }
      else image(x1,x2,z,col=icol,add=T)
      if(colorbar) {
        color.key(col,lev)
        if(!add) title(resp,line=1)
      } else if(!add) title(resp,line=0)
    }
  }
  else {
    # not filled
    if(prob.plot) {
      cat("Plotting contour for",main,"=",levels,"\n")
      contour(x1,x2,z,add=T, levels=levels,col=col,lwd=lwd,
              drawlabels=drawlabels,...)
    }
    else {
      nlevels <- length(levels(z))
      z <- as.numeric(z)
      dim(z) <- c(length(x1),length(x2))
      levels <- 0.5+(1:(nlevels-1))
      contour(x1,x2,z,add=T,levels=levels,col=col,lwd=lwd,
              drawlabels=drawlabels,...)
    }
  }
}
