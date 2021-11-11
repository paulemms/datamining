# different types of plots for data mining

prototypes <- function(x,y=NULL,fun=mean) {
  # basically a cleaned-up version of "aggregate"
  # that is specialized for one "by" factor
  if(!is.null(y) && mode(y) == "function") {
    fun <- y
    y <- NULL
  }
  if(is.null(y)) {
    resp <- response.var(x)
    y <- x[[resp]]
    if(!is.factor(y)) stop("no class variable")
    x <- not(x,resp)
  } else {
    resp <- "Group.1"
  }
  y = as.factor(y)
  if(identical(fun,sum)) {
    # fast implementation of special case
    as.data.frame(indicators.factor(y) %*% as.matrix(x))
  } else if(identical(fun,mean)) {
    # fast implementation of special case
    f <- indicators.factor(y)
    f <- scale.rows(f, 1/rowSums(f))
    as.data.frame(f %*% as.matrix(x))
  } else if(T) {
    # faster than below
    xs <- split(as.data.frame(x),y)
    xp <- sapply(xs,function(xl) apply(xl,2,fun))
    as.data.frame(t(xp))
  } else {
    yl <- list(y)
    names(yl) <- resp
    x <- aggregate(x,yl,fun)
    rownames(x) <- levels(y)
    x[setdiff(names(x),resp)]
  }
}

hclust.prototypes <- function(x,...) {
  resp <- response.var(x)
  px <- not(prototypes(x),resp)
  hclust(dist(px)^2,members=table(x[[resp]]),...)
}

reorder.svd.data.frame <- function(x,dims=c(1,2)) {
  cat("using svd to order\n")
  i <- sapply(x,is.numeric)
  nx <- x[i]
  nx = na.dummy(nx,0)
  s <- svd(nx - rep.row(colMeans(nx),nrow(nx)))
  v1 = s$u[,1]
  s <- svd(nx - rep.col(rowMeans(nx),ncol(nx)))
  v2 = s$v[,1]
  # remove flip ambiguity by alphabetic sort
  if(rownames(x)[which.min(v1)] > rownames(x)[which.max(v1)]) {
    v1 <- -v1
  }
  if(colnames(x)[which.min(v2)] > colnames(x)[which.max(v2)]) {
    v2 <- -v2
  }
  ord = list(order(v1),order(v2))
  for(d in 1:length(dim(x))) {
    if(d %in% dims) {
    } else {
      ord[[d]] = 1:dim(x)[d]
    }
  }
  if(any(!i)) {
    j = 1:length(i)
    j[i] = which(i)[ord[[2]]]
    ord[[2]] = j
  }
  ord
}
reorder.hc.data.frame <- function(x,dims=c(1,2)) {
  cat("using hclust to order\n")
  ord = list()
  if(1 %in% dims) {
    # rows
    hc = hclust(dist(x)^2,method="ward")
    ord[[1]] = hc$order
  } else {
    ord[[1]] = 1:nrow(x)
  }
  if(2 %in% dims) {
    # columns
    hc = hclust(dist(t(x))^2,method="ward")
    ord[[2]] = hc$order
  } else {
    ord[[2]] = 1:ncol(x)
  }
  ord
}


apply_order <- function(x, ord) {
  x <- x[ord[[1]],]
  x <- x[,ord[[2]]]
  x
}



#' Data Image
#'
#' Each value in a data matrix is represented by a colored pixel.
#'
#' @param x a matrix or data frame.
#' @param reorder If \code{TRUE}, the rows and columns are reordered so that
#' the data along each row and each column follows a linear trend.  Can also be
#' a vector of two logicals, indicating separately whether to reorder rows
#' and/or columns.
#' @param scale If \code{TRUE}, each column is scaled to have minimum 0 and
#' maximum 1.
#' @param col a vector of colors to use for representing values.
#' @param ... extra arguments to \code{\link{image.default}}.
#' @author Tom Minka
#' @seealso \code{\link{YR.colors}}
#' @references M. C. Minnotte and R. W. West. "The Data Image: A Tool For
#' Exploring High Dimensional Data Sets." Proceedings of the ASA Section on
#' Statistical Graphics, 1998. \url{http://citeseer.nj.nec.com/72085.html}
#' \url{http://math.usu.edu/~minnotte/research/pubs.html}
#' @examples
#'
#' data(USJudgeRatings)
#' data_image(USJudgeRatings,col=RC.colors(32))
#'
data_image <- function(x,las=1,reorder=T,scale=T,...) {
  #mar=c(2.5,8.5,0,0.1)
  x = as.data.frame(x)
  x <- x[sapply(x,is.numeric)]
  if(scale) x <- apply.df(x,scale.range)
  if(scale %% 2 == 0) {
    x = apply.df(x,rank)
  }
  if(scale > 2) {
    x = scale(x)
    s = linear.profiles(x)$s
    s = abs(s)
    x = scale.cols(x,s)
  }
  if(length(reorder) == 1) reorder = rep(reorder,2)
  if(any(reorder)) {
    ord <- reorder.hc.data.frame(x,which(reorder))
    x <- apply_order(x, ord)
  }
  image.table(as.matrix(x),las=las,...)
  invisible(ord)
}


#' Reordered and scaled star plot
#'
#' Automatically reorders and scales the variables to make a
#' readable star plot.
#'
#' @param x a matrix or data frame.  Each row makes one star.
#' @param proto If \code{TRUE} and the response variable is a factor,
#'     then only the mean value for each response is plotted.
#'     See \code{\link{prototypes}}.
#' @param reorder If \code{TRUE}, the rows and columns are reordered so
#'     that the data along each row and each column follows a linear trend.
#'     Can also be a vector of two logicals, indicating separately
#'     whether to reorder rows and/or columns.
#' @param scale If \code{TRUE}, each column is scaled to have minimum
#'     0 and maximum 1.
#' @param... additional arguments for \code{\link{stars}}.

#' @author{Tom Minka}
#' @seealso \code{\link{stars}}
#' @examples
#' data(iris)
#' star_plot(iris)
#' @export
star_plot <- function(x,draw.segments=T,proto=T,reorder=T,scale=T,...) {
  # reorder = c(T,F) will only reorder cases
  x = as.data.frame(x)
  resp <- response.var(x)
  if(is.factor(x[[resp]]) && proto) {
    cat("showing median of each",resp,"\n")
    f <- x[resp]
    x <- x[sapply(x,is.numeric)]
    x <- aggregate(x,f,median)
    rownames(x) <- levels(f[[1]])
  }
  x <- x[sapply(x,is.numeric)]
  if(scale) {
    x <- apply.df(x,scale.range)
    x = apply.df(x,sqrt)
  }
  if(length(reorder) == 1) reorder = rep(reorder,2)
  if(any(reorder)) {
    ord <- reorder.svd.data.frame(x,which(reorder))
    x <- apply_order(x, ord)
  }
  #opar <- par(mar=c(0.05,0,0,0.1))
  #on.exit(par(opar))
  stars(x, draw.segments=draw.segments, scale=F, ...)
}

parallel.cases <- function(x,yscale=c("linear","range","none"),
                           xscale=c("equal","linear","none"),
                           yaxt=if(yscale=="none") "s" else "n",...) {
  x <- as.data.frame(x)
  y <- x
  # assume no factors
  y <- data.matrix(y)
  x <- 1:nrow(y)
  yscale <- match.arg(yscale)
  xscale <- match.arg(xscale)
  if(yscale=="linear") {
    # do this to subtract means
    # scaling y here has no effect
    y <- scale(y)
    if(xscale != "none") {
      # want yij*sj = aj + bj*xi
      # aj subtracts mean of yij
      # then xi is top eigenvector
      s = svd(y)
      x = s$u[,1]
      s = abs(1/s$v[,1])
      if(xscale == "equal") {
        x <- rank.stable(x)
        # now act as if x was fixed and solve for scaling again
        xscale <- "none"
      }
      xscale <- "none"
    }
    if(xscale=="none") {
      # solve for sj when xj is known
      xt <- x - mean(x)
      if(T) {
        # scale to standardize residuals
        r = colSums(y*y) - (xt %*% y)^2/sum(xt*xt)
        r = r/nrow(y)
        s = 1/sqrt(r)
      } else {
        # scale to match svd solution
        r = (xt %*% y)/sum(xt*xt)
        s = 1/r
      }
    }
    # never any reversals
    y <- y * rep.row(s,nrow(y))
  } else {
    if(yscale == "range") {
      # map to [0,1]
      y <- y - rep.row(apply(y,2,min),nrow(y))
      y <- y / rep.row(apply(y,2,max),nrow(y))
    }
    if(xscale != "none") {
      # want yij = aj + bj*xi
      # subtract mean in each col to remove aj,
      # then xi is top eigenvector
      s <- svd(y)
      x <- s$u[,1]
      if(xscale == "equal") x <- rank.stable(x)
    }
  }
  # remove x flip ambiguity by alphabetic sort
  if(rownames(y)[which.min(x)] > rownames(y)[which.max(x)]) {
    x <- -x
  }
  i <- order(x)
  y <- reorder.rows(y,i)
  x <- x[i]

  labeled_curves(x,y,xlab="",ylab="",yaxt=yaxt,...)
  i
}

#' Parallel-coordinate plot
#'
#' @description Plots a data matrix using parallel coordinates.
#' @param x a data matrix or frame.  Only numeric columns will be
#'     plotted.
#' @param yscale describes how the columns should be scaled.
#' @param xscale describes how the columns should be placed.
#' @param proto If \code{TRUE} and the response variable is a factor,
#'     then only the mean value for each response is plotted.
#'     See \code{\link{prototypes}}.
#' @param ...Extra arguments to \code{\link{labeled.curves}}.
#' @details
#'   Let the rows of the matrix be \emph{cases} and the columns be
#'   \emph{dimensions}.
#'   Each dimension is laid out as a vertical axis, and the value for each
#'   case is plotted as a point on that axis.  Points corresponding to the
#'   same case are connected by lines.  Thus each case is represented by a
#'   curve or \emph{profile}, moving horizontally from dimension to dimension.
#'
#'   To enhance readability, the dimensions should usually be shifted to have
#'   a common center and scaled to have
#'   comparable units.  This is controlled by \code{yscale}.
#'   If \code{yscale="linear"}, the dimensions are automatically shifted
#'   and scaled
#'   by the linear profiles method.  That is, the profiles are
#'   made to be as straight as possible.
#'   This may result in negative scaling for a dimension,
#'   i.e. the axis is reversed so that moving up means a lower value.
#'   When this happens, the dimension name is prepended with a minus sign.
#'   If \code{yscale="range"}, the dimensions are shifted and scaled (by a
#'                                                                    positive number) so that the data ranges from 0 to 1 on each axis.
#'   If \code{yscale="none"}, the dimensions are plotted as-is.
#'
#'   Another important choice is the ordering and spacing of the axes
#'   on the plot.
#'   If \code{xscale="linear"}, the axes are placed according to the
#'   linear profiles method.
#'   If \code{xscale="equal"}, the axes are placed similarly to
#'   \code{"linear"} but with the constraint that they must be equally
#'   spaced.
#'   If \code{xscale="none"}, the axes are placed in the order
#'   that they appear in the matrix.
#'
#'   If the data frame has a categorical response, the profiles are colored
#'   according to the response.  The color scheme can be adjusted using the
#'   arguments to \code{\link{labeled.curves}}.
#' @author Tom Minka
#' @references A. Inselberg and B. Dimsdale.
#'   "Parallel coordinates: A tool for visualizing multidimensional geometry."
#'   Proc. of Visualization 90, p. 361-78, 1990.
#'
#' E. J. Wegman.
#' "Hyperdimensional Data Analysis Using Parallel Coordinates."
#' JASA 85:664-675, 1990.
#' @seealso \code{\link{parallel.cases}},\code{\link{star.plot}}
#' @examples
#' data(iris)
#' parallel_plot(iris)
#' parallel_plot(iris,proto=F,labels=NULL)
#' @export
parallel_plot <- function(x,yscale=c("linear","range","none"),flipy=F,
                          xscale=c("equal","linear","none"),proto=T,flipx=F,
                          yaxt=if(yscale=="none") "s" else "n",type=1,...) {
  x <- as.data.frame(x)
  resp <- response.var(x)
  if(is.factor(x[[resp]])) {
    if(proto) {
      cat("showing median of each",resp,"\n")
      f <- x[resp]
      x <- x[sapply(x,is.numeric)]
      x <- aggregate(x,f,median)
      rownames(x) <- levels(factor(f[[1]]))
      group <- 1:nrow(x)
    } else {
      # plot cases, but color by class
      group <- as.numeric(factor(x[[resp]]))
    }
    x <- not(x,resp)
  } else {
    group <- 1:nrow(x)
  }
  y <- x

  # assume no factors
  y <- data.matrix(y)
  x <- 1:ncol(y)

  yscale <- match.arg(yscale)
  xscale <- match.arg(xscale)
  if(yscale=="linear") {
    # do this to subtract means
    # scaling y here has no effect
    y <- scale(y)
    # bugfix for scale: NaNs are created when a column has zero variance
    y[is.nan(y)] = 0
    if(xscale != "none") {
      r = linear.profiles(y,type=type); x = r$x; s = r$s
      if(xscale == "equal") {
        x <- rank.stable(x)
        # now act as if x was fixed and solve for scaling again
        xscale <- "none"
      }
    }
    if(xscale=="none") {
      s = linear.profiles(y,x=x,type=type)$s
    }
    # minimize number of reversals
    if(sum(s<0) > sum(s>0)) s <- -s
    if(flipy) s = -s
    y <- y * rep.row(s,nrow(y))
    reversed <- (s<0)
  } else {
    if(yscale == "range") {
      # map to [0,1]
      y <- y - rep.row(apply(y,2,min),nrow(y))
      y <- y / rep.row(apply(y,2,max)+eps,nrow(y))
    }
    if(xscale != "none") {
      x = linear.profiles(y,s=1,type=type)$x
      if(any(is.na(x))) stop("internal: x is nan")
      if(xscale == "equal") x <- rank.stable(x)
    }
    reversed = F
  }
  # remove x flip ambiguity by alphabetic sort
  if(colnames(y)[which.min(x)] > colnames(y)[which.max(x)]) {
    x <- -x
  }
  if(flipx) x = -x
  # must do this after the sort
  if(any(reversed)) {
    cat("axis reversed for",colnames(y)[reversed],"\n")
    colnames(y)[reversed] <- sapply(colnames(y)[reversed], function(s) paste("-",s,sep=""))
  }
  i <- order(x)
  y <- reorder.cols(y,i)
  x <- x[i]

  labeled_curves(x,t(y),group=group,xlab="",ylab="",yaxt=yaxt,...)
}

improve.s <- function(A,s,n) {
  d = length(s)
  if(nrow(A) != d || ncol(A) != d) error("A is wrong size")
  for(j in 1:d) {
    b = drop(A[j,-j] %*% s[-j])/A[j,j]
    s1 = 0.5*(-b + sqrt(b^2 + 4*n/A[j,j]))
    s2 = 0.5*(-b - sqrt(b^2 + 4*n/A[j,j]))
    j1 = A[j,j]*(s1 + b)^2 - n*log(s1^2)
    j2 = A[j,j]*(s2 + b)^2 - n*log(s2^2)
    print(c(j1,j2))
    s[j] = if(j1 < j2) s1 else s2
  }
  s
}

linear.profiles2 <- function(y,x=NULL,s=NULL) {
  y = as.matrix(y)
  n = nrow(y)
  d <- ncol(y)
  r = linear.profiles(y,x,s); x = r$x; s = r$s
  s = rnorm(d)
  #print(s)
  #print(x)
  C <- diag(d) - array(1/d,c(d,d))
  for(iter in 1:10) {
    old.x = x
    for(iter2 in 1:100) {
      old.s = s
      xt <- x - mean(x)
      P = C - outer(xt,xt)/sum(xt*xt)
      # Hadamard product
      A = P * (t(y) %*% y)
      # function value
      J = t(s) %*% A %*% s - n*sum(log(s^2))
      print(J[1])
      s = improve.s(A,s,nrow(y))
      activity = max(abs(s - old.s))
      #print(activity)
      if(activity < 1e-2) break
    }
    if(is.null(x)) {
      # improve x
      x = svd(scale.cols(y,s) %*% C)$v[,1]
      #print(x)
    }
    activity = max(abs(x - old.x))
    if(activity < 1e-4) break
  }
  s = s * sqrt(sum(1/s^2))
  cat(iter,"iters\n")
  list(x=x,s=s)
}


linear.profiles <- function(y,x=NULL,s=NULL,type=1) {
  y = as.matrix(y)
  #y = scale(y)
  if(is.null(x) && is.null(s)) {
    # Hartigan's linear profile algorithm
    # want yij*sj = ai + bi*xj
    # yij = ai*(1/sj) + bi*(xj/sj)
    # so 1/sj = v[j,1], xj/sj = v[j,2]
    # the constraint is sum_j 1/sj^2 = 1
    # scaling by singular values is irrelevant
    sv <- svd(y)
    x <- sv$v[,2]/sv$v[,1]
    s <- 1/sv$v[,1]
  } else if(type==1 && is.null(s)) {
    # solve for sj when xj is known
    sv = svd(y)
    s = 1/sv$v[,1]
  } else if(is.null(s)) {
    # solve for sj when xj is known
    # this comes from solving for (ai,bi) analytically and substituting
    # to get eigenvector problem for sj
    # the constraint is sum_j sj^2 = 1
    d <- ncol(y)
    xt <- x - mean(x)
    a <- diag(d) - array(1/d,c(d,d)) - outer(xt,xt)/sum(xt*xt)
    # Hadamard product
    a <- a * (t(y) %*% y)
    s <- svd(a)
    # want the smallest eigenvector of a
    s <- s$u[,d]
  } else if(is.null(x)) {
    # want yij = ai + bi*xj
    # subtract mean in each row to remove ai,
    # then xj is top eigenvector
    sv <- svd(y - rep.col(rowMeans(y),ncol(y)))
    x <- sv$v[,1]
  }
  list(x=x,s=s)
}


# abbreviate strings so that max width is w, and strings remain unique
abbreviate.width <- function(s,w,cex=par("cex")) {
  old.max <- w
  repeat {
    wid <- strwidth(s,cex=cex)
    if(max(wid)<=w) break
    i <- wid > w
    j <- which.max(wid)
    # this may cause a collision
    s[i] <- abbreviate(s[i],nchar(s[j])-1)
  }
  s
}

# minka: use auto.mar and auto.layout
stars <-
  function(x, full = TRUE, scale = TRUE, radius = TRUE,
           labels = dimnames(x)[[1]], locations = NULL,
           layout = NULL, len = 1,
           key.loc = NA, key.labels = dimnames(x)[[2]], key.xpd = TRUE,
           xlim = NULL, ylim = NULL, flip.labels = NULL,
           draw.segments = FALSE, col.segments = 1:n.seg,
           col.stars = NA,
           axes = FALSE, frame.plot = axes,
           main = NULL, sub = NULL, xlab = "", ylab = "",
           cex = 0.8, lwd = 0.25, lty = par("lty"), xpd = FALSE,
           mar = pmin(par("mar"),auto.mar(main,xlab,ylab,axes)),
           ...)
  {
    # technically, strheight doesn't work until we've made axes
    #char.height <- max(strheight(labels,cex=cex))
    char.height <- cex*par("cxy")[2]
    cell.size <- function() {
      if(is.null(flip.labels)) flip.labels = T
      ff <- c(2.1,2.1)
      # minka: adjust loc for labels
      if(!is.null(labels)) {
        ff[2] <- ff[2] + char.height*(1+flip.labels)
      }
      ff
    }
    if (is.data.frame(x))
      x <- data.matrix(x)
    else if (!is.matrix(x))
      stop("x must be a matrix or a data frame")
    if (!is.numeric(x))
      stop("data in x must be numeric")

    n.loc <- nrow(x)
    n.seg <- ncol(x)

    if (is.null(locations)) { ## Default (x,y) locations matrix
      # minka: use auto.layout
      if(is.null(layout)) {
        ff = cell.size()
        asp = ff[2]/ff[1]
        layout = auto.layout(n.loc,asp)
      }
      if(!is.null(labels) && is.null(flip.labels))
        flip.labels <- layout[2] * mean(nchar(labels)) > 30
      ff = cell.size()
      locations <- expand.grid(ff[1] * 1:layout[2],
                               ff[2] * layout[1]:1)[1:n.loc, ]
    }
    else {
      if (is.numeric(locations) && length(locations) == 2) {
        ## all stars around the same origin
        locations <- cbind(rep(locations[1],n.loc),
                           rep(locations[2],n.loc))
        if(!missing(labels) && n.loc > 1)
          warning("labels don't make sense for a single location")
        else labels <- NULL
      }
      else {
        if (is.data.frame(locations))
          locations <- data.matrix(locations)
        if (!is.matrix(locations) || ncol(locations) != 2)
          stop("locations must be a 2-column matrix.")
        if (n.loc != nrow(locations))
          stop("number of rows of locations and x must be equal.")
      }
      if(is.null(flip.labels))
        flip.labels <- FALSE # have no grid
    }
    xloc <- locations[,1]
    yloc <- locations[,2]
    ## Angles start at zero and pace around the circle counter
    ## clock-wise in equal increments.
    angles <-
      if(full)
        seq(0, 2*pi, length=n.seg+1)[-(n.seg+1)]
    else if (draw.segments)
      seq(0, pi, length=n.seg+1)[-(n.seg+1)]
    else
      seq(0, pi, length=n.seg)

    if (length(angles) != n.seg)
      stop("length(angles) must be the same as ncol(x)")

    ## Missing values are treated as 0
    if (scale) {
      x <- apply(x, 2, function(x)
        (x - min(x, na.rm = TRUE))/diff(range(x, na.rm = TRUE)))
      # minka: sqrt because we are depicting area
      x = sqrt(x)
    } else {
      # just scale so max is 1, as assumed by rest of code
      x = x/max(x)
    }
    ## Missing values are treated as 0
    x[is.na(x)] <- 0
    mx <- max(x <- x * len)

    # minka: place key automatically at top
    # do this now for computing ylim
    if(!is.null(key.loc) && is.na(key.loc)) {

      key.loc <- c(mean(xloc), max(yloc) + 2.1 + char.height)
    }
    if(is.null(xlim)) {
      xlim <- range(xloc) + c(-mx, mx)
      # minka: should make room for labels too - use range.text
      if(!is.null(key.loc)) {
        # minka: make room for the key
        w <- max(strwidth(key.labels,cex=cex,units="inches"))/par("pin")[1]
        xlim.key = extend.inches(c(key.loc[1]-len,key.loc[1]+len),w,w)
        xlim[1] <- min(xlim[1],xlim.key[1])
        xlim[2] <- max(xlim[2],xlim.key[2])
      }
    }
    if(is.null(ylim)) {
      ylim <- range(yloc) + c(-mx, mx)
      if(!is.null(labels)) {
        # minka: make room for the bottom row labels
        ylim[1] <- min(yloc - len - char.height)
      }
      if(!is.null(key.loc)) {
        # minka: make room for the key
        ylim[1] <- min(ylim[1],key.loc[2]-len-char.height)
        ylim[2] <- max(ylim[2],key.loc[2]+len+char.height)
      }
    }

    deg <- pi / 180

    ## The asp argument keeps everything (the symbols!) square
    op <- par(mar = mar, xpd = xpd) ; on.exit(par(op))
    plot(0, type="n", ..., xlim=xlim, ylim=ylim,
         main = main, sub = sub, xlab = xlab, ylab=ylab, asp = 1, axes = axes)

    s.x <- xloc + x * rep(cos(angles), rep(n.loc,n.seg))
    s.y <- yloc + x * rep(sin(angles), rep(n.loc,n.seg))

    if ( draw.segments ) {
      aangl <- c(angles, if(full)2*pi else pi)
      for (i in 1:n.loc) { ## for each location, draw a segment diagram
        px <- py <- numeric()
        for (j in 1:n.seg) {
          k <- seq(from = aangl[j], to = aangl[j+1], by = 1*deg)
          px <- c(px, xloc[i], s.x[i,j], x[i,j]*cos(k) + xloc[i], NA)
          py <- c(py, yloc[i], s.y[i,j], x[i,j]*sin(k) + yloc[i], NA)
        }
        polygon(px, py, col = col.segments, lwd=lwd, lty=lty)
      }
    } # Segment diagrams

    else { # Draw stars instead
      for (i in 1:n.loc) {
        polygon(s.x[i,], s.y[i,], lwd=lwd, lty=lty, col = col.stars[i])
        if (radius)
          segments(rep(xloc[i],n.seg),
                   rep(yloc[i],n.seg),
                   s.x[i,], s.y[i,], lwd=lwd, lty=lty)
      }
    }

    if(!is.null(labels)) {
      ## vertical text offset from center
      y.off <- mx * (if(full) 1 else 0.1)
      if(flip.labels)
        # minka: bug fix
        y.off <- y.off + cex*par("cxy")[2] * ((seq(n.loc)-1)%%2)
      ##DBG cat("mx=",format(mx),"y.off:"); str(y.off)
      text(xloc, yloc - y.off, labels, cex=cex, adj=c(0.5, 1))
    }

    if ( !is.null(key.loc) ) { ## Draw unit key
      ## usually allow drawing outside plot region:
      par(xpd = key.xpd) # had `xpd' already above
      key.x <- len * cos(angles) + key.loc[1]
      key.y <- len * sin(angles) + key.loc[2]
      if (draw.segments) {
        px <- py <- numeric()
        for (j in 1:n.seg) {
          k <- seq(from = aangl[j], to = aangl[j+1], by = 1*deg)
          px <- c(px, key.loc[1], key.x[j], len * cos(k) + key.loc[1], NA)
          py <- c(py, key.loc[2], key.y[j], len * sin(k) + key.loc[2], NA)
        }
        polygon(px, py, col = col.segments, lwd=lwd, lty=lty)
      }
      else { # draw unit star
        polygon(key.x, key.y, lwd=lwd, lty=lty)
        if (radius)
          segments(rep(key.loc[1],n.seg), rep(key.loc[2],n.seg),
                   key.x, key.y, lwd=lwd, lty=lty)
      }

      ## Radial Labeling -- should this be a standalone function ?
      lab.angl <- angles +
        if(draw.segments) (angles[2] - angles[1]) / 2 else 0
      label.x <- 1.1 * len * cos(lab.angl) + key.loc[1]
      label.y <- 1.1 * len * sin(lab.angl) + key.loc[2]
      ## Maybe do the following without loop {need not use adj but ..)!
      for (k in 1:n.seg) {
        text.adj <-
          c(## horizontal
            if      (lab.angl[k] < 90*deg || lab.angl[k] > 270*deg) 0
            else if (lab.angl[k] > 90*deg && lab.angl[k] < 270*deg) 1
            else 0.5,
            ## vertical
            if (lab.angl[k] <= 90*deg) (1 - lab.angl[k] / (90*deg)) /2
            else if (lab.angl[k] <= 270*deg)
              (lab.angl[k] - 90*deg) / (180*deg)
            else ## lab.angl[k] > 270*deg
              1 - (lab.angl[k] - 270*deg) / (180*deg)
          )
        text(label.x[k], label.y[k],
             labels= key.labels[k], cex = cex, adj = text.adj)
      }
    } # Unit key is drawn and labelled

    if (frame.plot) box(...)

    invisible(locations)
  }

# modified by Tom Minka
# abbreviation, better spacing, better handling of key
my.stars <-
  function(x, full = TRUE, scale = TRUE, radius = TRUE,
           labels = dimnames(x)[[1]],
           locations = NULL, xlim = NULL, ylim = NULL, len = 1,
           colors = NULL, abbrev=T,
           key.loc = NA, key.labels = NULL, key.xpd = TRUE,
           draw.segments = FALSE, axes = FALSE,
           cex = 0.8, lwd = 0.25, ...)
  {
    if (is.data.frame(x))
      x <- as.matrix(x)
    else if (!is.matrix(x))
      stop("x must be a matrix or a data frame")
    if (!is.numeric(x))
      stop("data in x must be numeric")

    n.loc <- nrow(x)
    n.seg <- ncol(x)
    deg <- pi / 180			# segments only

    seg.colors <- if(!is.null(colors)) colors else 1:n.seg

    if (is.null(locations)) {		# make loc matrix
      md <- ceiling(sqrt(n.loc)) # =>  md^2 >= n.loc
      loc.x <- 2.1*((1:md)-1)
      loc.y <- rev(loc.x)
      # minka: adjust loc for labels
      if(!is.null(labels)) {
        # technically, strheight doesn't work until we've made axes
        spc <- max(strheight(labels,cex=cex))
        loc.y <- rev((2.1+spc)*((1:md)-1))
      }
      loc <- expand.grid(loc.x, loc.y)[1:n.loc, ]
    }
    else {
      if (is.numeric(locations) && length(locations) == 2) {
        ## all stars around the same origin
        loc <- cbind(rep(locations[1],n.loc),
                     rep(locations[2],n.loc))
        if(!missing(labels) && n.loc > 1)
          warning("labels don't make sense for a single location")
        else labels <- NULL
      }
      else {
        if (is.data.frame(locations))
          locations <- data.matrix(locations)
        if (!is.matrix(locations) || ncol(locations) != 2)
          stop("locations must be a 2-column matrix.")
        loc <- .Alias(locations)
        if (n.loc != nrow(loc))
          stop("number of rows of locations and x must be equal.")
      }
    }

    ## Angles start at zero and pace around the circle counter
    ## clock-wise in equal increments.
    angles <-
      if(full)
        seq(0, 2*pi, length=n.seg+1)[-(n.seg+1)]
    else if (draw.segments)
      seq(0, pi, length=n.seg+1)[-(n.seg+1)]
    else
      seq(0, pi, length=n.seg)

    if (length(angles) != n.seg)
      stop("length(angles) must be the same as ncol(x)")

    ## Missing values are treated as 0
    x[is.na(x)] <- 0

    if (scale) {
      # minka: use full range of [0.1,1]
      x <- sweep(x,2,apply(x,2,min),FUN="-")
      x <- sweep(x,2,apply(x,2,max)/0.9,FUN="/")
      x <- sweep(x,2,0.1,FUN="+")
      ## Columns of 0s will put NAs in x, next line gets rid of them
      x[is.na(x)] <- 0
      # minka: sqrt because we are depicting area
      x = sqrt(x)
    }

    x <- x * len

    # minka: place key automatically at top
    if(!is.null(key.loc) && is.na(key.loc)) {
      key.loc <- c(mean(loc[[1]]), max(loc[[2]]) + 2.5)
    }
    if(is.null(key.labels)) key.labels <- dimnames(x)[[2]]

    if(is.null(xlim)) {
      xlim <- range(loc[,1] + len, loc[,1] - len)
      if(!is.null(key.loc)) {
        w <- max(strwidth(key.labels,cex=cex))
        xlim[1] <- min(xlim[1],key.loc[1]-1-w)
        xlim[2] <- max(xlim[2],key.loc[1]+1+w)
      }
    }
    if(is.null(ylim)) {
      ylim <- range(loc[,2] + len, loc[,2] - len)
      if(!is.null(labels)) {
        spc <- max(strheight(labels,cex=cex))
        ylim[1] <- min(loc[,2] - 1 - spc)
      }
      if(!is.null(key.loc)) {
        h <- max(strheight(key.labels,cex=cex))
        ylim[1] <- min(ylim[1],key.loc[2]-1-h)
        ylim[2] <- max(ylim[2],key.loc[2]+1+h)
      }
    }

    ## The asp argument keeps everything square
    plot(0, type="n", ..., xlim=xlim, ylim=ylim,
         xlab="", ylab="", asp = 1, axes = axes)

    # minka: abbreviate labels if necessary
    if(!is.null(labels) && abbrev) {
      labels <- abbreviate.width(labels,2,cex=cex)
    }

    if ( draw.segments ) {
      ## for each location, draw a segment diagram
      for ( i in 1:n.loc ) {
        poly.x <- NA
        poly.y <- NA
        start.x <- x[i,] * cos( angles ) + loc[i,1]
        start.y <- x[i,] * sin( angles ) + loc[i,2]

        ### FIXME : we can do without the following inner loop !

        for (j in 1:n.seg) {
          poly.x <- c(poly.x,loc[i,1],start.x[j])
          poly.y <- c(poly.y,loc[i,2],start.y[j])

          next.angle <-
            if(j < n.seg)
              angles[j+1]
          else (if(full) 360 else 180) * deg

          k <- seq(from = angles[j], to = next.angle, by = deg)
          poly.x <- c(poly.x, x[i,j] * cos( k ) + loc[i,1], NA)
          poly.y <- c(poly.y, x[i,j] * sin( k ) + loc[i,2], NA)
        }
        polygon(poly.x, poly.y, col = seg.colors, lwd=lwd)
        if (!is.null(labels))
          text(loc[i,1], loc[i,2] - if(full)max(x) else 0.1 * max(x),
               labels[i], cex=cex, adj=c(0.5,1))
      }
    } # Segment diagrams

    else { # Draw stars instead
      for ( i in 1:n.loc ) {
        temp.x <- x[i,] * cos( angles ) + loc[i,1]
        temp.y <- x[i,] * sin( angles ) + loc[i,2]
        if (radius)
          segments(rep(loc[i,1],n.seg),
                   rep(loc[i,2],n.seg),
                   temp.x, temp.y, lwd=lwd)
        lines(c(temp.x, temp.x[1]),
              c(temp.y, temp.y[1]), lwd=lwd)
        if (!is.null(labels))
          text(loc[i,1], loc[i,2] - if(full)max(x) else 0.1 * max(x),
               labels[i], cex=cex, adj=c(0.5,1))
      }
    }

    if ( !is.null(key.loc) ) { ## Draw unit key

      ## allow drawing outside plot region (inside figure region):
      op <- par(xpd = key.xpd); on.exit(par(op))
      key.r <- len * 0.8
      key.x.coord <- cos( angles ) * key.r + key.loc[1]
      key.y.coord <- sin( angles ) * key.r + key.loc[2]
      if ( draw.segments ) {
        key.x <- NA
        key.y <- NA
        for (j in 1:n.seg){
          key.x <- c(key.x,key.loc[1],key.x.coord[j])
          key.y <- c(key.y,key.loc[2],key.y.coord[j])
          k <- angles[j] + deg
          next.angle <-
            if (j < n.seg) angles[j+1]
          else (if(full) 360 else 180) * deg

          while (k < next.angle) {
            key.x <- c(key.x, key.r * cos( k ) + key.loc[1])
            key.y <- c(key.y, key.r * sin( k ) + key.loc[2])
            k <- k + deg
          }
          key.x <- c(key.x, key.r * cos( next.angle ) + key.loc[1], NA)
          key.y <- c(key.y, key.r * sin( next.angle ) + key.loc[2], NA)
        }
        polygon(key.x, key.y, col = seg.colors, lwd=lwd)
      }
      else { # draw unit star
        if ( radius )
          segments(rep(key.loc[1],n.seg), rep(key.loc[2],n.seg),
                   key.x.coord, key.y.coord, lwd=lwd)
        lines(c(key.x.coord, key.x.coord[1]),
              c(key.y.coord, key.y.coord[1]), lwd=lwd)
      }

      lab.angl <- angles +
        if(draw.segments) (angles[2] - angles[1]) / 2 else 0

      label.r <- len*0.9
      label.x <- cos( lab.angl ) * label.r + key.loc[1]
      label.y <- sin( lab.angl ) * label.r + key.loc[2]

      ##-- FIXME : Do the following witout loop !
      for (k in 1:n.seg) {
        text.adj <-
          c(## horizontal
            if      (lab.angl[k] < 90*deg || lab.angl[k] > 270*deg) 0
            else if (lab.angl[k] > 90*deg && lab.angl[k] < 270*deg) 1
            else 0.5,
            ## vertical
            if (lab.angl[k] <= 90*deg) (1 - lab.angl[k] / (90*deg)) /2
            else if (lab.angl[k] <= 270*deg)
              (lab.angl[k] - 90*deg) / (180*deg)
            else ## lab.angl[k] > 270*deg
              1 - (lab.angl[k] - 270*deg) / (180*deg)
          )

        text(label.x[k], label.y[k],
             labels= key.labels[k], cex = cex, adj = text.adj)
      }
    } # Unit key is drawn and labelled
    invisible()
  }

qp <- function(A,b,tol=1e-10) {
  # returns the vector x which minimizes x'x subject to Ax=b.
  #library(MASS)
  # standardize the problem
  s = sqrt(rowSums(A * A))
  A = A/rep.col(s,ncol(A))
  b = b/s

  d = ncol(A)
  x = array(0,c(d,1))
  ws = c()
  for(iter in 1:100) {
    # find most violated constraint
    s = (A %*% x) - b
    # constraints in the working set must be satisfied
    # this is to avoid problems with numerical inaccuracy
    s[ws] = 0
    i = which.min(s)
    if(s[i] >= -tol) break
    if(i %in% ws) {
      print(s[i])
      stop("double violated constraint")
    }
    ws = c(ws,i)
    while(T) {
      # solve equality-constrained problem
      Aw = A[ws,,drop=F]
      bw = b[ws,drop=F]
      if(T) {
        y = solve(Aw %*% t(Aw), bw)
      } else {
        y = ginv(Aw %*% t(Aw)) %*% bw
      }
      x = t(Aw) %*% y
      if(any(y<=eps)) {
        # remove y<0 from working set
        ws = ws[y>eps]
      } else break
    }
    #cat("working set:",ws,"\n",fill=60)
  }
  x
}

# not used
collision = function(x,y,w,h) {
  # returns T iff the first box collides with any of other boxes
  n = length(x)
  for(i in 2:n) {
    if((abs(x[1]-x[i]) < (w[1]+w[i])/2) &&
       (abs(y[1]-y[i]) < (h[1]+h[i])/2)) return(T)
  }
  F
}

find.gap <- function(base=0,h0,y,h,symmetric=T) {
  # returns the smallest position of a box of size h0 which will not collide
  # with boxes of size h centered on y.
  # sort y
  ord = order(y)
  y = y[ord]
  h = h[ord]
  n = length(y)
  # compute all gaps
  dy = diff(y)-(h[1:(n-1)]+h[2:n])/2
  # dy must be positive since there are no collisions among y's
  i = which(dy > h0)
  # upper end
  y.new = max(y[n]+(h[n]+h0)/2,base)
  if(length(i) > 0) {
    # this isn't quite optimal
    # lower bound
    y.new = c(y.new, y[i]+(h[i]+h0)/2)
    i.next = pmin(i+1,n)
    # upper bound
    y.new = c(y.new, y[i.next]-(h[i.next]+h0)/2)
  }
  if(symmetric) {
    # lower end
    y.new = c(y.new, y[1]-(h[1]+h0)/2)
  }
  if(!symmetric) {
    y.new = y.new[y.new >= base]
    if(length(y.new) == 0) stop("empty")
  }
  y.new[which.min(abs(y.new - base))]
}

move.collisions.vertical <- function(x,w,h,base=0,symmetric=T) {
  # find the minimal vertical shifts so that boxes of width w and height h
  # centered on x do not overlap.
  n = length(x)
  if(length(w) == 1) w = rep(w,n)
  if(length(h) == 1) h = rep(h,n)
  if(length(base) == 1) base = rep(base,n)
  y = numeric(n)
  step = median(h)/4
  # loop in reverse order
  y[n] = base[n]
  for(i in (n-1):1) {
    # find nearby boxes (those for which overlap is possible)
    is = (i+1):n
    j = is[abs(x[i]-x[is]) < (w[i]+w[is])/2]
    # now only need to avoid vertical overlap
    if(length(j) == 0) y[i] = base[i]
    else y[i] = find.gap(base[i],h[i],y[j],h[j],symmetric)
    if(y[i] < base[i] && !symmetric) {
      print(y[i])
      print(base[i])
      print(h[i])
      print(y[j])
      print(h[j])
      stop("")
    }
  }
  y
}

move.collisions <- function(x,w) {
  # move x's the minimal amount so that boxes of width w centered on x
  # do not overlap.
  # units of x and w must match (preferably inches).
  if(length(x) == 1) return(x)
  if(F) {
    # jitter x
    x = x + 2*eps*runif(length(x))
    x = x + runif(length(x))
  }
  ord = order(x)
  x = x[ord]
  w = w[ord]

  # construct constraint matrix
  # constraints are (x[j]+dx[j])-(x[i]+dx[i]) > (w[j]+w[i])/2
  # or dx[j]-dx[i] > (w[j]+w[i])/2 -(x[j]-x[i])
  # for all j>i
  n = length(x)
  nc = n*(n-1)/2
  A = array(0,c(nc,n))
  b = array(0,c(nc,1))
  index = 1
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      A[index,j] = 1
      A[index,i] = -1
      b[index] = (w[j]+w[i])/2 - (x[j]-x[i])
      index = index + 1
    }
  }
  colnames(A) = names(x)
  x = x + drop(qp(A,b))
  # restore original order
  x[ord] = x
  x
}

absdist <- function(x,y=x) {
  abs(rep.row(y,length(x)) - rep.col(x,length(y)))
}
clip <- function(x) {
  (x > 0)*x
}
line.of.sight <- function(bbox,bbox.seg,dx,dy,offset) {
  n = length(dx)
  e = rep(0,n)
  for(i in 1:n) {
    bbox.i = bbox[setdiff(1:n,i),]
    s = bbox.i - rep.row(bbox.seg[i,],nrow(bbox.i))
    s[,c(1,3)] = -s[,c(1,3)]
    m1 = apply(s,1,min)
    bbox.i = bbox.i[(m1 > 0),,drop=F]
    if(nrow(bbox.i) > 0) {
      #print(nrow(bbox.i))
      dm = array(c(dy[i],0,-dx[i],0, 0,dy[i],-dx[i],0,
                   dy[i],0,0,-dx[i], 0,dy[i],0,-dx[i]),c(4,4))
      s = (bbox.i %*% dm) - offset[i]
      m1 = apply(s,1,min)
      m2 = apply(s,1,max)
      e[i] = sum((m1*m2 < 0)*pmin(-m1,m2)^2)
    }
  }
  sum(e)
}

move.collisions2 <- function(x,y,w,h,cost=rep(1,length(x)),los=F) {
  # move (x,y) the minimal amount so that boxes of width w and height h
  # centered on (x,y) do not overlap.
  # units of (x,y) and (w,h) must match (preferably inches).
  # cost is a vector specifying how immovable each box is.
  # cost[i] = 0 means box i will not be moved.
  # returns the new (x,y) values, as a matrix.
  ok = !is.na(x) & !is.na(y) & !is.na(w) & !is.na(h)
  x0=x;y0=y;w0=w;h0=h;cost0=cost
  x=x[ok];y=y[ok];w=w[ok];h=h[ok];cost=cost[ok]
  n = length(x)
  cost2 = c(cost,cost)
  n2 = 2*n
  i = c((1:n2)[cost2 > 0],(1:n2)[cost2 == 0])
  ir = 1:n2
  ir[i] = 1:n2
  nz = sum(cost2 == 0)
  wij = (rep.row(w,n)+rep.col(w,n))/2
  diag(wij) = 0
  hij = (rep.row(h,n)+rep.col(h,n))/2
  diag(hij) = 0

  cost.fcn <- function(dxy) {
    # returns the "badness" of a given c(dx,dy) vector
    dxy = c(dxy,rep(0,nz))[ir]
    dx = dxy[1:n]
    dy = dxy[n+(1:n)]
    x = x+dx
    y = y+dy
    # want total movement to be small
    e = dx*dx+dy*dy
    # moving beyond the original box is especially bad
    e = e/100 + clip(abs(dx) - w/2)^2 + clip(abs(dy) - h/2)^2
    e = sum(cost*e)
    # total overlap
    eo = clip(wij-absdist(x))*clip(hij-absdist(y))
    eo = sum(eo^2)
    el = 0
    if(los) {
      # want to have an unblocked "line of sight" to the original position
      bbox = cbind(x-w/2,x+w/2,y-h/2,y+h/2)
      offset = x*dy - y*dx
      bbox.seg = cbind(x-dx,x,y-dy,y)
      bbox.seg = cbind(pmax(bbox.seg[,1],bbox.seg[,2]),
                       pmin(bbox.seg[,1],bbox.seg[,2]),
                       pmax(bbox.seg[,3],bbox.seg[,4]),
                       pmin(bbox.seg[,3],bbox.seg[,4]))
      el = line.of.sight(bbox,bbox.seg,dx,dy,offset)
    }
    #cat(e," ",eo," ",el,"\n")
    e + 100*eo + 100*el
  }
  dxy = rep(0,n2-nz)
  if(F) {
    res = optim(dxy,cost.fcn,control=list(maxit=10000))
    dxy = res$par
    cat("optim code: ",res$convergence,"\n")
    cat("iters: ",res$counts,"\n")
    cat("fcn value: ",res$value,"\n")
  }
  if(T) {
    res = nlm(cost.fcn,dxy)
    dxy = res$estimate
    cat("nlm code: ",res$code,"\n")
    cat("iters: ",res$iterations,"\n")
    cat("fcn value: ",res$minimum,"\n")
  }
  # see if it is really optimal
  #dxy[1] = 0
  #print(cost.fcn(dxy))
  dxy = c(dxy,rep(0,nz))[ir]
  dx = dxy[1:n]
  dy = dxy[n+(1:n)]
  x = x+dx
  y = y+dy
  x0[ok]=x;y0[ok]=y;x=x0;y=y0
  cbind(x,y)
}


panel.hist <- function(x, ...)
{
  # helper function for pairs(), from examples
  # usage: pairs(x,diag.panel=panel.hist)
  # put histograms on the diagonal
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

# default color for partial residuals
formals(termplot)$col.res <- 1
# default lowess span
#formals(panel.smooth)$span <- 4/5
#formals(lowess)$f <- 3/4

# returns a model frame where the responses have been replaced by their
# residuals under the fit
residual.frame <- function(object,data) {
  if(missing(data)) {
    x <- my.model.frame(object)
    resp <- response.var(x)
    if(inherits(object,"glm")) {
      p <- object$fitted.value
      x[[resp]] <- (object$y - p)/p/(1-p)
    }
    else x[[resp]] <- residuals(object)
  }
  else {
    x <- data
    resp <- response.var(object)
    if(inherits(object,"glm") || inherits(object,"tree")) {
      x[[resp]] = misclassResiduals(object,data)
    } else {
      # expand p to include NAs
      p <- predict(object,x)[rownames(x)]
      names(p) <- NULL
      x[[resp]] <- x[[resp]] - p
      #attr(x,"terms") <- terms(object)
    }
    x = my.model.frame(formula(paste(resp,"~.")),x)
  }
  x
}

partial.residuals <- function(object,omit,condensed=T) {
  resp <- response.var(object)
  pred <- predictor.vars(object)
  if(condensed) {
    # remove the predictor and refit the model
    x <- list()
    for(i in 1:length(omit)) {
      fit <- my.update(object,paste(".~. -",omit[i]))
      #cat(pred[i],coef(fit),"\n")
      x[[i]] <- residuals(fit)
    }
    names(x) <- omit
  } else {
    # do not refit the model
    x <- data.frame(residuals(object,type="partial"),check.names=F)[omit]
    r = residuals(object)
    for(i in 1:length(x)) {
      names(x[[i]]) = names(r)
    }
  }
  x
}

partial.residual.frame <- function(object,omit,...) {
  x = model.frame(object)
  resp <- response.var(x)
  x[[resp]] = partial.residuals(object,omit,...)[[1]]
  pred = predictor.vars(formula(paste(resp,"~",omit)))
  x[c(pred,resp)]
}


expand.frame <- function(data,v) {
  # also see expand.model.frame
  if(length(v) == 0) return(data)
  vs = strsplit(v,":")
  for(i in 1:length(v)) {
    x = data[vs[[i]]]
    if(any(sapply(x,is.factor))) {
      as.contrasts = function(f) if(is.factor(f)) as.numeric(f==levels(f)[2]) else f
      x = apply.df(x,as.contrasts)
    }
    data[[v[[i]]]] <- apply(x,1,prod)
  }
  data
}


#' Plot predictors versus response
#'
#' Makes a matrix of pairwise scatterplots with lowess-type trend lines.
#'
#' @author{Tom Minka}
#'
#' @aliases predict_plot predict_plot.data.frame predict_plot.formula
#' @param formula a formula specifying the response and predictor variables
#' @param data,x a data frame with at least two columns
#' @param partial a model from which to compute partial residuals (used
#'                                        by \code{\link{predict_plot.lm}}).
#' @param mcol,mlwd If plotting partial residuals of an \code{lm},
#'     the color and width of the model predictions.
#' @param layout a vector \code{c(rows,cols)} specifying the desired
#'     layout of panels.  Otherwise chosen automatically based on the size
#'     of the plotting window.
#' @param highlight a logical vector specifying which predictors to highlight.
#' @param se If \code{TRUE}, show standard errors in linecharts.
#' @param scol,slwd color and width of trend lines.
#' @param span,degree,family parameters for the trend line (see \code{loess}).
#' @param rtype how a factor response should be handled when drawing a
#'     trend line.
#' @param identify.pred A character vector of predictor names for which
#'     to interactively \code{\link{identify}} points.  If \code{TRUE},
#'     done for all predictors.
#' @param mar margins within each panel
#' @param xaxt,yaxt arguments to \code{par}
#' @param col plotting color for symbols
#' @param asp Aspect ratio for each panel.  If \code{"auto"}, the aspect
#'     ratio is chosen automatically based on the trend line and
#'     \code{\link{auto.aspect}}.
#' @param given,given.lab,nlevels,pretty,key,bg,color.palette,pch.palette used for conditioning plots.
#' @param main,xlab,ylab axis labels.
#' @param ... extra arguments passed to \code{predict_plot.data.frame} or
#'     \code{\link{plot}}.
#' @return
#'   If the predictor is numeric, makes a scatterplot with loess line on top.
#'   If the predictor is a factor, makes a \code{\link{linechart}}.
#'
#' @seealso \code{\link{loess}}, \code{\link{model.plot}}
#' @examples
#' data(Cars)
#' predict_plot(Price ~ ., CarsT)
#' fit <- lm(Price ~ ., CarsT)
#' predict_plot(Price~ ., CarsT, partial=fit)
#' # same thing using predict_plot.lm
#' predict_plot(fit, partial = TRUE)
#' @export
predict_plot <- function(...) UseMethod("predict_plot")

#' @export
predict_plot.default <- function(x,y,labels,xlab,ylab,...) {
  frame = data.frame(x,y)
  if(!missing(labels)) rownames(frame) = labels
  if(missing(xlab)) xlab = deparse(substitute(x))
  if(missing(ylab)) ylab = deparse(substitute(y))
  predict_plot.data.frame(frame,xlab=xlab,ylab=ylab,...)
}

std.error = function(x) sd(x)/sqrt(length(x))

#' @export
predict_plot.data.frame <- function(x,layout=NULL,partial=NULL,condensed=T,
                                    rtype=c("prob","logit","probit"),
                                    highlight,se=T,
                                    identify.pred=F,
                                    mcol=0,mlwd=2,
                                    scol="red",slwd=2,
                                    span=2/3,degree=2,family="symmetric",
                                    mar=NA,bg=par("bg"),
                                    xaxt=par("xaxt"),yaxt=par("yaxt"),
                                    col=1,asp=NULL,
                                    given=NULL,given.lab=NULL,nlevels=2,
                                    pretty=T,key=!is.null(given),
                                    color.palette=default.colors,
                                    pch.palette=default.pch,
                                    main=NULL,xlab,ylab=NULL,...) {
  rtype <- match.arg(rtype)
  resp <- response.var(x)
  if(missing(degree) && is.factor(x[[resp]])) degree = 0
  pred <- predictor.terms(x)
  x <- expand.frame(x,setdiff(pred,predictor.vars(x)))
  if(!is.null(partial)) {
    if(missing(highlight)) highlight <- predictor.terms(partial)
    if(!is.null(given.lab)) {
      if(given.lab %in% predictor.vars(partial)) {
        partial = my.update(partial,paste(".~.-",given.lab))
      }
    }
    pres <- partial.residuals(partial,pred,condensed=condensed)
    if(is.null(ylab)) ylab = paste("Partial",resp)
  }
  else {
    if(missing(highlight)) highlight <- NULL
    if(is.null(ylab)) ylab = resp
  }
  if(!is.null(given)) {
    # key==1 means separate panel, key==2 means on each panel
    if(key==2 && is.null(main)) main = given.lab
    b = NULL
    if(!is.factor(given)) {
      # make given a factor
      b <- break.quantile(given,nlevels,pretty=pretty)
      given <- rcut(given,b)
    } else nlevels = length(levels(given))
    if(is.function(color.palette)) color.palette = color.palette(nlevels)
    lev = 1:nlevels
    level.color = color.palette[((lev-1) %% length(color.palette))+1]
    level.pch <- pch.palette[(((lev-1) %/% length(color.palette)) %% length(pch.palette))+1]
  }
  if(is.null(layout)) layout <- auto.layout(length(pred) + (key==1))
  opar <- par(bg=bg)
  if(any(layout != c(1,1))) {
    opar = c(opar,par(mfrow=layout))
  }
  if(!is.null(mar)) {
    if(is.na(mar)) {
      mar = auto.mar(main=main,xlab=if(missing(xlab)) "y" else xlab,
                     ylab=if(missing(ylab)) "y" else ylab,xaxt=xaxt,yaxt=yaxt,key=(key==2))
    }
    opar <- c(opar,par(mar=mar))
  }
  on.exit(par(opar))
  for(i in pred) {
    if(!is.null(partial)) {
      # sometimes pres is longer than x, because of NAs
      x[[resp]] <- pres[[i]][rownames(x)]
      if(is.null(x[[resp]])) stop(paste("partial of",i,"not found"))
    }
    if(is.factor(x[[i]]) && !is.ordered(x[[i]])) {
      x[[i]] <- sort.levels(x[[i]],as.numeric(x[[resp]]))
    }
    # plot points
    if(is.null(given)) {
      # no given
      if(!is.na(scol)) {
        # smooth curve
        #k <- is.finite(x[,i]) & is.finite(x[[resp]])
        # keep warnings about infinite points
        k <- !is.na(x[,i]) & !is.na(x[[resp]])
        fit = NULL
        if(scol != 0) {
          if(is.factor(x[[resp]])) {
            # factor response
            if(length(levels(x[[resp]]))==2 && !is.factor(x[[i]])) {
              if(rtype=="prob") {
                fit = loprob(x[k,i], x[k,resp], degree=degree, span=span)
                if(F) {
                  # piecewise constant approx
                  library(tree)
                  frame = data.frame(x=x[k,i],y=x[k,resp])
                  tr = tree(y~x,frame)
                  partition.tree(tr,add=T,force=F)
                }
              } else {
                fit <- loprob(x[k,i], x[k,resp])
                p <- fit$y-1
                p <- switch(rtype,logit=log(p/(1-p)),probit=qnorm(p))
                fit$y <- p+1.5
                plot(x[[i]],2*as.numeric(x[[resp]])-3,xlab=i,ylab=rtype)
              }
            }
          } else {
            # numeric response
            if(!is.factor(x[[i]])) {
              if(!all(is.finite(x[k,i])))
                stop(paste(i,"has infinite values.\nCheck your transformation."))
              if(!all(is.finite(x[k,resp])))
                stop(paste(resp,"has infinite values.\nCheck your transformation."))
              if(length(k) > 1500) {
                # save time on dense plots
                fit = lowess(x[k,i], x[k,resp], f=span)
              } else {
                fit = lowess2(x[k,i], x[k,resp], f=span, degree=degree, family=family)
              }
            }
          }
        }
      }
      real.asp = if(identical(asp,"auto")) auto.aspect(fit) else asp
      if(is.factor(x[[resp]])) {
        # factor response
        if(rtype=="prob") {
          if(is.factor(x[[i]])) {
            # should sort if the factor is not ordered
            #mosaicplot(table(x[[i]], x[[resp]]), xlab=i, ylab=ylab, ...)
            tab = table(x[[resp]],x[[i]])
            names(dimnames(tab)) = c(resp,i)
            linechart(row.probs(tab,se=T), xlab=" ",...)
          } else {
            plot(x[[i]], x[[resp]], xlab="", ylab=ylab, ...)
          }
        }
      } else if(is.factor(x[[i]])) {
        # factor predictor, numeric response
        rt = rtable(formula(paste(resp,"~",i)),x)
        se = rtable(formula(paste(resp,"~",i)),x,fun=std.error) * 1.64
        linechart(rt,se,labels=NULL,ylab=ylab,xlab=" ",...)
      } else {
        # numeric response and predictor
        plot(x[[i]], x[[resp]], xlab="", ylab=ylab, main=main,
             asp=real.asp, col=col, ...)
        if(!is.null(partial) && inherits(partial,"lm")) {
          a <- coef(partial)
          #a0 = a["(Intercept)"]
          #if(is.null(a0)) a0 = 0
          a0 = -mean(x[[i]])*a[i]
          if(mcol != 0 && !is.na(a[i])) abline(a0,a[i],col=mcol,lwd=mlwd)
        }
      }
      if(missing(xlab)) {
        if(i %in% highlight) title(xlab=i,font.lab=2)
        else title(xlab=i)
      } else title(xlab=xlab)
      if(!is.null(fit)) lines(fit,col=scol,lwd=slwd)
      # identify
      if((identify.pred == T) || (i %in% identify.pred)) {
        identify(x[k,i],x[k,resp],labels=rownames(x)[k])
      }
    } else if(is.factor(x[[i]])) {
      # given
      fmla = formula(paste(resp,"~",i,"+",given.lab))
      x[[given.lab]] = given
      rt = rtable(fmla,x)
      if(se) {
        rt.se = rtable(fmla,x, fun=std.error) * 1.64
        rt.se[is.na(rt.se)] = sd(x[[resp]])
      } else rt.se = NULL
      linechart(t(rt),t(rt.se),xlab=" ",
                xaxt=xaxt,yaxt=yaxt,row.labels=NULL,...)
      font.lab = if(i %in% highlight) 2 else 1
      title(xlab=i,font.lab=font.lab)
    } else {
      # condition on the given variable
      plot.new()
      box()
      title(main=main,ylab=ylab,line=if(key==1) NA else 1)
      if(i %in% highlight) title(xlab=i,font.lab=2)
      else title(xlab=i)
      if(is.factor(x[[i]])) {
        # setup window for factor
        xlim <- length(levels(x[[i]]))
        plot.window(xlim=c(0.5,xlim+0.5),ylim=range(x[[resp]]))
        if(xaxt != "n") axis(1,1:xlim,labels=levels(x[[i]]))
        if(yaxt != "n") axis(2)
        cat(paste("jittering",i,"\n"))
      } else {
        # setup window for numeric
        plot.window(xlim=range(x[[i]],na.rm=T),ylim=range(x[[resp]],na.rm=T))
        if(xaxt != "n") axis(1)
        if(yaxt != "n") axis(2)
      }
      lev <- levels(given)
      for(g in 1:length(lev)) {
        color <- level.color[g]
        pch <- level.pch[g]
        val <- lev[g]
        k <- !is.na(x[,i]) & !is.na(x[,resp])
        k = k & (given == val)
        if(is.factor(x[[i]])) {
          jitter <- (runif(length(x[k,i]))-0.5)/5
          points(as.numeric(x[k,i])+jitter, x[k,resp], col=color,
                 pch=pch, ...)
        } else {
          points(x[k,i], x[k,resp], col=color, pch=pch, ...)
        }
        if(is.factor(x[[resp]])) {
          fit = loprob(x[k,i], x[k,resp])
        } else {
          fit = lowess2(x[k,i], x[k,resp], f=span, degree=degree)
        }
        lines(fit,col=color,lwd=slwd)
      }
    }
  }
  # show legend for the given variable, as a separate panel
  # this is better than putting a color key on each plot
  if(!is.null(given)) {
    if(key==1) {
      plot.new()
      if(!is.null(given.lab)) {
        font.lab = if(given.lab %in% highlight) 2 else 1
        title(xlab=given.lab,font.lab=font.lab)
      }
      lev = levels(given)
      y <- cumsum(strheight(lev)+0.02)
      for(i in 1:length(lev)) {
        color <- color.palette[i]
        val <- lev[i]
        text(0.5,0.5+y[i]-y[length(y)]/2,val,col=color,adj=0.5)
        if(length(unique(level.pch))>1)
          points(0.6+strwidth(val)/2,0.5+y[i],pch=level.pch[i])
      }
    } else if(key==2) {
      i = 1:length(levels(given))
      if(is.null(b))
        color.key(color.palette[i],level.pch[i],levels(given))
      else
        color.key(color.palette[i],level.pch[i],breaks=b)
    }
  }
}

#' Plot predictors versus residuals.
#'
#' Makes a matrix of pairwise scatterplots with lowess-type trend lines.
#'
#' Partial residuals are computed by fitting a new model with the
#' predictor removed, which is different from \code{residuals(type="partial")}.
#' @author Tom Minka
#' @param object the output of \code{lm}.
#' @param data a data frame to use instead of \code{model.frame(object)}.
#' @param partial If \code{TRUE}, plot partial residuals instead of residuals.
#' @param ylab axis label.
#' @param ... extra arguments for \code{\link{predict_plot.data.frame}}.
#' @return
#'   A plot similar to \code{\link{predict_plot}}, but where the vertical
#'   axis is residuals.  These plots can be used to judge which predictors
#'   should be added to the model.
#' @seealso \code{\link{predict_plot}}
#' @examples
#' # see the examples for predict_plot
#' @export
predict_plot.lm <- function(object,data,partial=F,ylab=NULL,...) {
  if(!partial) {
    if(is.null(ylab)) ylab = paste("Residual",response.var(object))
    if(missing(data)) {
      res <- residual.frame(object)
    } else {
      res <- residual.frame(object,data)
    }
    cat("plotting residuals\n")
    predict_plot.data.frame(res,highlight=predictor.terms(object),ylab=ylab,...)
  } else {
    if(missing(data)) data <- model.frame(object)
    else data = my.model.frame(formula(paste(response.var(object),"~.")),data)
    cat("plotting partial residuals\n")
    predict_plot.data.frame(data,partial=object,ylab=ylab,...)
  }
}

given.vars <- function(formula) {
  rhs <- formula[[3]]
  if(is.call(rhs) && (deparse(rhs[[1]]) == "|")) {
    # remove givens from formula
    deparse(rhs[[3]])
  } else NULL
}

remove.given <- function(formula) {
  rhs <- formula[[3]]
  formula[[3]] <- rhs[[2]]
  formula
}

#' @exportS3Method
predict_plot.formula <- function(formula,data=parent.frame(),...) {
  # formula has givens?
  given = given.vars(formula)
  if(!is.null(given)) {
    formula = remove.given(formula)
    if(is.environment(data)) g <- get(given,env=data)
    else g <- data[[given]]
    if(is.null(g))
      stop(paste("variable \"",given,"\" not found",sep=""))
    return(predict_plot.formula(formula,data,
                                given=g,given.lab=given,...))
  }
  if(F) {
    expr <- match.call(expand = F)
    expr$... <- NULL
    expr$na.action <- na.pass
    expr[[1]] <- as.name("model.frame.default")
    x <- eval(expr, parent.frame())
  } else {
    # formula has its own environment
    x <- model.frame.default(formula,data,na.action=na.pass)
  }
  predict_plot.data.frame(x,...)
}


# only used once in lab 10
# expand.cross <- function(object) {
#   resp <- response.var(object)
#   pred <- predictor.vars(object)
#   # quadratic terms only
#   pred2 <- c()
#   for(i in pred) {
#     for(j in pred) {
#       if(match(i,pred) < match(j,pred)) pred2 = c(pred2,paste(i,j,sep=":"))
#     }
#   }
#   pred = c(pred,pred2)
#   formula(paste(resp,"~",paste(pred,collapse="+")))
# }
# step.up <- function(object,scope=expand.cross(object)) {
#   step(object,list(upper=scope,lower=formula(object)))
# }


#' Contour plot matrix
#'
#' @examples
#' fit <- lm(Fertility ~ ., data=swiss)
#' interact.plot(fit)
#' @export
interact.plot <- function(object, ...) UseMethod("interact.plot")

#' @export
interact.plot.lm <- function(object,data,partial=F,main=NULL,...) {
  if(!partial) {
    if(missing(data)) {
      res <- residual.frame(object)
    } else {
      res <- residual.frame(object,data)
    }
    cat("plotting residuals\n")
    if(is.null(main)) main = "Residuals"
    interact.plot.data.frame(res,highlight=predictor.terms(object),main=main,...)
  } else {
    if(missing(data)) data <- model.frame(object)
    else data = my.model.frame(formula(paste(response.var(object),"~.")),data)
    cat("plotting partial residuals\n")
    interact.plot.data.frame(data,partial=object,main=main,...)
  }
}

#' @export
interact.plot.data.frame <- function(x,ypred,partial=NULL,highlight,span=0.75,
                                     scol="red",slwd=2,type=c("*",":"),
                                     xaxt="n",yaxt="n",se=T,
                                     bg=par("bg"),nlev=8,main=NULL,...) {
  # type=":" means partials remove cross term only
  # type="*" means partials remove everything
  #library(modreg)
  type = match.arg(type)
  resp <- response.var(x)
  pred = predictor.vars(x)
  if(missing(ypred)) ypred = pred
  if(missing(bg) & any(sapply(x[pred],is.numeric))) bg = grey(0.65)
  n = length(pred)
  mar=c(0.05,0.05,0,0.05)
  if(all(ypred == pred)) layout = c(n,n)
  else layout = auto.layout(n*length(ypred))
  opar = par(mfrow=layout,mar=mar,bg=bg)
  on.exit(par(opar))
  if(!is.null(partial)) {
    if(missing(highlight)) highlight <- predictor.terms(partial)
    if(is.null(main)) {
      main = paste("Partial residuals\n",
                   if(type == ":") "cross-term" else "both predictors","absent",sep="\n")
    }
  }
  else if(missing(highlight)) highlight <- NULL
  for(py in ypred) {
    for(px in pred) {
      if(px == py) {
        plot.new()
        box()
        y = 0.5
        if(is.factor(x[[py]])) {
          auto.legend(levels(x[[py]]),default.colors(nlevels(x[[py]])))
          #y = 0.5
        }
        text(0.5,y,py,font=if(py %in% highlight) 2 else 1)
      } else if(!is.factor(x[[py]]) && (match(py,pred) < match(px,pred))) {
        plot.new()
        if(match(py,pred) == 1 && match(px,pred) == 2)
          text(0.5,0.5,main)
      } else {
        # do this to convert the names so formula(y) will work
        y = data.frame(x[c(px,py,resp)])
        if(!is.null(partial)) {
          predxy = paste(px,type,py,sep="")
          fit = my.update(partial,paste(".~.-",predxy))
          y[[resp]] = residuals(fit)
        }
        mar = c(0.05,0.05,0,0.05)
        if(is.factor(y[[px]]) && is.factor(y[[py]])) {
          rt = rtable(formula(y),y)
          xlab = if(px %in% ypred) "" else px
          if(se) {
            rt.se = rtable(formula(y),y,fun=std.error) * 1.64
            rt.se[is.na(rt.se)] = sd(x[[resp]])
          } else rt.se = NULL
          linechart(t(rt),t(rt.se),xlab=xlab,ylab="",yaxt=yaxt,labels=NULL,...)
        } else if(is.factor(y[[px]])) {
          plot.new()
        } else if(is.factor(y[[py]])) {
          fmla = formula(paste(resp,"~",px,"|",py))
          predict_plot(fmla,y,xlab="",ylab="",xaxt=xaxt,yaxt=yaxt,key=F)
        } else if(is.factor(y[[resp]])) {
          color.plot(y,
                     xlab="",ylab="",mar=mar,key=F,axes=F,n=nlev,...)
        } else {
          color.plot(loess(formula(y),y,span=span),
                     xlab="",ylab="",zlab="",mar=mar,key=F,axes=F,n=nlev,...)
        }
        predxy = paste(px,":",py,sep="")
        predyx = paste(py,":",px,sep="")
        if(any(c(predxy,predyx) %in% highlight)) {
          box(col="red")
        } else {
          box()
        }
      }
    }
  }
}


if(!exists("my.factor.scope")) old.factor.scope = factor.scope
my.factor.scope <- function(factor, scope)
{
  drop <- scope$drop
  add <- scope$add

  if(length(factor) && !is.null(drop)) {# have base model
    nmdrop <- colnames(drop)
    facs <- factor
    if(length(drop)) {
      nmfac <- colnames(factor)
      where <- match(nmdrop, nmfac, 0)
      if(any(!where)) stop("lower scope is not included in model")
      facs <- factor[, -where, drop = FALSE]
      nmdrop <- nmfac[-where]
    } else nmdrop <- colnames(factor)
    # minka: allow non-hierarchy
  } else nmdrop <- character(0)

  if(!length(add)) nmadd <- character(0)
  else {
    nmfac <- colnames(factor)
    nmadd <- colnames(add)
    if(!is.null(nmfac)) {
      where <- match(nmfac, nmadd, 0)
      if(any(!where)) stop("upper scope does not include model")
      nmadd <- nmadd[-where]
      add <- add[, -where, drop = FALSE]
    }
    # minka: allow non-hierarchy
  }
  list(drop = nmdrop, add = nmadd)
}
#replaceInNamespace("factor.scope",my.factor.scope)
#replaceInNamespace("factor.scope",old.factor.scope)

