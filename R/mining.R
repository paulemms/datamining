# Routines for controlling aspect ratio

# taken from library(lattice)
banking <- function (dx, dy) {
  if (is.na(dx)) NA
  else {
    if (is.list(dx)) {
      dy <- dx[[2]]
      dx <- dx[[1]]
    }
    if (length(dx) != length(dy))
      stop("Non matching lengths")
    id <- dx != 0 & dy != 0
    if (any(id)) {
      r <- abs(dx[id]/dy[id])
      median(r)
    }
    else 1
  }
}

# from Splus
banking <- function(dx, dy, iter = 100, tolerance = 0.5)
{
  dx <- abs(dx)
  dy <- abs(dy)
  if(sum(is.na(dx)) != 0) {
    warning("infinite in x, no banking")
    dx[is.na(dx)] <- 0
  }
  if(sum(is.na(dy)) != 0) {
    warning("infinity in y, no banking")
    dy[is.na(dy)] <- 0
  }
  if(sum(dx) == 0 || sum(dy) == 0) return(NA)
  good <- dx > 0
  dx <- dx[good]
  dy <- dy[good]
  # minka: fixed bug in original (dy/dx)
  mad <- median(dx/dy)
  # initial guess
  ar <- if(mad > 0) mad else 1
  radians.to.angle <- 180/pi
  for(i in 1:iter) {
    distances <- sqrt(dx^2 + (ar * dy)^2)
    orientations <- atan2(ar * dy, dx)
    avg <- (radians.to.angle * sum(orientations * distances))/
      sum(distances)
    if(abs(45 - avg) < tolerance)
      break
    ar <- ar * (45/avg)
  }
  ar
}

banking.fcn <- function(dx,dy,asp) {
  # want result to be 45
  i = (dx != 0) & (dy != 0)
  dx = dx[i]
  dy = dy[i]
  len = sqrt(dx^2 + (asp*dy)^2)
  radians.to.angle <- 180/pi
  orient = abs(atan(asp*dy,dx)*radians.to.angle)
  #hist(orient)
  sum(orient*len)/sum(len)
}

# this function is broken


#' Choose aspect ratio
#'
#' Uses Cleveland's "banking to 45" rule to determine an optimal aspect ratio
#' for a plot.
#'
#' @param x,y numeric vectors of coordinates defining a continuous curve, or
#' multiple curves delimited by \code{NA}s.  Alternatively, \code{x} can be a
#' list of two vectors.
#' @return An aspect ratio, suitable for the \code{asp} parameter to
#' \code{\link{plot}} or \code{\link{plot.window}}.
#' @author Tom Minka
#' @export
#' @examples
#'
#' data(lynx)
#' plot(lynx)
#' plot(lynx,asp=auto.aspect(time(lynx),lynx))
#'
auto.aspect <- function(x,y) {
  # returns an aspect ratio in user coordinates, suitable for asp=
  # x and y are vectors
  # NA specifies a break in the curve, just as in lines()
  if(is.list(x)) {
    y = x[[2]]
    x = x[[1]]
  }
  dx = diff(x); dy = diff(y)
  i <- !is.na(dx) & !is.na(dy)
  dx <- dx[i]; dy <- dy[i]
  asp = banking(dx,dy)
  #print(banking.fcn(dx,dy,asp))
  asp
}

# this version doesn't preserve margins properly
set.aspect <- function(aspect,margins=T) {
  aspect <- c(1, aspect)
  din <- par("din")
  pin <- aspect*min(din/aspect)
  if(margins) {
    # mai is (bottom,left,top,right)
    m <- din-pin
    mai <- c(sum(par("mai")[c(2,4)]), sum(par("mai")[c(1,3)]))
    pin <- din - pmax(m, mai)
  }
  par(pin=pin)
}
reset.aspect <- function() {
  mai <- c(sum(par("mai")[c(2,4)]), sum(par("mai")[c(1,3)]))
  par(pin=par("din")-mai)
}

# Sets graphics parameters so the next plot has aspect ratio asp.
# It is okay to run set.aspect repeatedly without resetting the graphics
# parameters.
#
# To reset:
#   reset.aspect()
# or:
#   opar <- set.aspect(a); par(opar)
#
# In R, it is also possible to force a given aspect ratio using
#   layout(1,width=1,height=aspect,respect=T)
# and reset with
#   layout(1)
# However, the aspect ratio it gives is not quite correct.
#
# Test:
#   set.aspect(2); plot(runif(100)); reset.aspect()
#   layout(1,width=1,height=2,respect=T); plot(runif(100)); layout(1)
#
set.aspect <- function(aspect,margins=T) {
  fin <- par("fin")
  a <- aspect*fin[1]/fin[2]
  # plt is (left,right,bottom,top)
  if(a > 1) plt <- c(0.5-1/a/2,0.5+1/a/2,0,1)
  else plt <- c(0,1,0.5-a/2,0.5+a/2)
  if(margins) {
    # add room for margins
    # mai is (bottom,left,top,right)
    mplt <- par("mai")[c(2,4,1,3)]/c(fin[1],fin[1],fin[2],fin[2])
    mplt[c(2,4)] <- 1-mplt[c(2,4)]
    # compose plt and mplt
    plt[1:2] <- mplt[1] + (mplt[2]-mplt[1])*plt[1:2]
    plt[3:4] <- mplt[3] + (mplt[4]-mplt[3])*plt[3:4]
  }
  par(plt=plt)
}
reset.aspect <- function() {
  if(F) {
    fin <- par("fin")
    mplt <- par("mai")[c(2,4,1,3)]/c(fin[1],fin[1],fin[2],fin[2])
    mplt[c(2,4)] <- 1-mplt[c(2,4)]
    par(plt=mplt)
  } else {
    # simpler way
    par(mar=par("mar"))
  }
}


# histograms with confidence intervals and bin merging

bhist.stats <- function(x,b) {
  if(length(b)==1) {
    r <- range(x,na.rm=T)
    b <- seq(0,1,by=1/b)*diff(r) + r[1]
  }
  if(is.null(b)) {
    h <- hist(x, probability=TRUE, plot=FALSE)
  } else {
    h <- hist(x, probability=TRUE, plot=FALSE, breaks=b)
    #print("Using given breaks")
  }
  b <- h$breaks
  len <- length(b)
  wx <- b[2:len] - b[1:len-1]
  xm <- h$mids
  nx <- h$counts + wx/sum(wx)*8
  n <- sum(nx)
  p <- nx/n
  list(p=p,n=n,nx=nx,xm=xm,wx=wx,h=h,b=b)
}

# b is starting number of bins
# n is final number of bins
# note: large b can cause problems, because of the greediness of merging


#' Merge histogram bins
#'
#' Quantize a variable by merging similar histogram bins.
#'
#' The desired number of bins is achieved by successively merging the two most
#' similar histogram bins.  The distance between bins of height (f1,f2) and
#' width (w1,w2) is measured according to the chi-square statistic
#' \deqn{w1*(f1-f)^2/f + w2*(f2-f)^2/f} where f is the height of the merged
#' bin: \deqn{f = (f1*w1 + f2*w2)/(w1 + w2)}
#'
#' @param x a numerical vector
#' @param b the starting number of bins, or a vector of starting break
#' locations. If NULL, chosen automatically by \code{\link{hist}}.
#' @param n the desired number of bins.
#' @return A vector of bin breaks, suitable for use in \code{\link{hist}},
#' \code{\link{bhist}}, or \code{\link{cut}}. Two plots are shown: a
#' \code{\link{bhist}} using the returned bin breaks, and a merging trace.  The
#' trace shows, for each merge, the chi-square distance of the bins which were
#' merged.  This is useful for determining the appropriate number of bins.  An
#' interesting number of bins is one that directly precedes a sudden jump in
#' the chi-square distance.
#' @author Tom Minka
#' @export
#' @examples
#'
#' x <- c(rnorm(100,-2,0.5),rnorm(100,2,0.5))
#' b <- seq(-4,4,by=0.25)
#' merge_hist(x,b,10)
#' # according to the merging trace, n=5 and n=11 are most interesting.
#'
#' x <- runif(1000)
#' b <- seq(0,1,by=0.05)
#' merge_hist(x,b,10)
#' # according to the merging trace, n=6 and n=9 are most interesting.
#' # because the data is uniform, there should only be one bin,
#' # but chance deviations in density prevent this.
#' # a multiple comparisons correction in merge_hist may fix this.
#'
merge_hist <- function(x,b=NULL,n=b,trace=T) {
  ss <- c()
  bn <- NULL
  repeat {
    bh <- bhist.stats(x,b)
    p <- bh$p; wx <- bh$wx
    b <- bh$b
    if(is.null(bn)) bn <- b
    # compute scores
    len <- length(p)
    p1 <- p[1:len-1]
    w1 <- wx[1:len-1]
    f1 <- p1/w1
    p2 <- p[2:len]
    w2 <- wx[2:len]
    f2 <- p2/w2
    fs <- (p1+p2)/(w1+w2)
    # s is difference in expected log-prob of future data
    #s <- p1*log(f1/fs) + p2*log(f2/fs)
    # needs a multiple-comparisons correction, as in CHAID
    # chi-square statistic
    s <- w1*(f1-fs)^2/fs + w2*(f2-fs)^2/fs
    #s <- p1*p2
    #s <- (p1 - p2)^2
    #s <- p1+p2
    ss[len-1] <- min(s)
    j <- which.min(s)

    b <- not(b, j+1);
    # bn is the one they want
    if(length(b) == n+1) bn <- b
    if(length(b) <= 3) break
  }
  b <- bn
  if(trace) {
    split.screen(c(2,1))
    plot(1:length(ss),ss^(0.25),xlab="bins",ylab="chi-sq")
    #plot(2:length(ss),-diff(ss^(0.25)),xlab="bins",ylab="chi-sq")
    screen(2)
  }
  bhist(x,b)
  title(paste(length(b)-1,"bins"))
  if(trace) close.screen(all = TRUE)
  b
}



#' Histogram with confidence intervals
#'
#' Same as \code{\link{hist}} except a confidence interval is drawn around each
#' bin height.
#'
#' The width of the interval for height p is
#' \code{sqrt(p*(1-p)/n)*exp(-1/6/p/n)}.
#'
#' @param x a numerical vector
#' @param b the number of bins, or a vector of break locations.  If NULL,
#' chosen automatically by \code{\link{hist}}.
#' @return None.
#' @author Tom Minka
#' @export
#' @examples
#'
#' x <- c(rnorm(100,-2,0.5), rnorm(100,2,0.5))
#' b <- seq(-4,4,by=0.25)
#' bhist(x,b)
#'
bhist <- function(x,b=NULL,z=1.64,xlab="",main="") {
  bh <- bhist.stats(x,b)
  p <- bh$p; nx <- bh$nx; n <- bh$n; xm <- bh$xm; wx <- bh$wx; h <- bh$h

  if(FALSE) {
    # exact interval
    # z=1 means q=0.1586
    f1 <- qbeta(q, nx, n-nx)/wx
    f2 <- qbeta(1-q, nx, n-nx)/wx
  } else {
    # approx interval
    s <- sqrt(p*(1-p)/n)*exp(-1/6/nx)
    f1 <- (p - z*s)/wx
    f2 <- (p + z*s)/wx
  }
  f <- p/wx
  h$counts <- nx
  h$intensities <- f
  h$density <- f
  plot(h,col="bisque",freq=FALSE,xlab=xlab,main=main)
  points(xm, f, col="blue", pch=18)
  arrows(xm, f1, xm, f2,
	 code = 3, col = "green", angle = 75, length = .1)
  if(FALSE) {
    p <- exp(digamma(nx)-digamma(n))
    s <- sqrt(1/nx + 1/n)
    s <- sqrt(trigamma(nx)+trigamma(n))
    f1 <- exp(log(p) - z*s)/wx
    f2 <- exp(log(p) + z*s)/wx
  }
  if(FALSE) {
    # compare intervals
    s <- sqrt(p*(1-p)/n)*exp(-1/6/nx)
    f1 <- (p - s)/wx
    f2 <- (p + s)/wx
    arrows(xm+0.1, f1, xm+0.1, f2,
	   code = 3, col = "red", angle = 75, length = .1)
  }
}
# routines for quantization and one-dimensional clustering
# Tom Minka 1/3/02

# reduce data to n points, by local averaging
# use quantiles so size of each bin is equal
squash <- function(x,n=1000) {
  b <- floor(length(x)/n)
  rb <- rep(b,n)
  # spread the overflow broadly over the range
  overflow <- length(x) - n*b
  rb <- rb + c(0,diff(trunc(overflow/n*(1:n))))
  # f is 1 1 1 2 2 2 3 3 3 ...
  f <- rep(1:n, rb)
  f <- f[1:length(x)]
  x <- sort(x)
  tapply(x,f,mean)
}


#' Draw break locations
#'
#' Draws lines to indicate break locations.
#' @param b a numeric vector
#' @return Uses \code{\link{segments}} to draw thick blue vertical lines on top
#'   of the current plot, at the given locations.
#' @author Tom Minka
#' @seealso
#'   \code{\link{plot.hclust.breaks}},
#'   \code{\link{plot.segments.ts}},
#'   \code{\link{break_kmeans}}
#' @export
plot_breaks <- function(b) {
  r <- par("usr")[3:4]
  segments(b,rep(r[1],length(b)), b,rep(r[2],length(b)), col="blue", lwd=2)
}



#' Histogram with hierarchical cluster breaks
#'
#' A representation of a one-dimensional hierarchical clustering
#'
#'
#' @param h an \code{hclust} object
#' @param x the numerical vector that was clustered to produce \code{h}
#' @param k a vector of the cluster cardinalities to plot
#' @return A histogram of \code{x} is shown with blue lines cutting the x-axis.
#' The tallest lines correspond to divisions made at the top of the hierarchy.
#' By reading top to bottom, you can see how each cluster is subdivided. This
#' can be far more illuminating than a plot of the hierarchy as a tree.
#' @author Tom Minka
#' @exportS3Method
#' @seealso \code{\link{boxplot.hclust}}, \code{\link{ward}},
#' \code{\link{break_ward}}
hist.hclust <- function(hc,x,k=2:5) {
  hb <- break.equal(x,40)
  h <- hist(x,hb,freq=FALSE,col="bisque",main="",xlab="")
  top <- par("usr")[4]
  # for cutree
  for(i in 1:length(k)) {
    q <- cutree(hc,k[i])
    b <- break.cut(x,q)
    scale <- (length(k)-i+1)/length(k)
    segments(b,rep(0,length(b)), b,rep(top*scale,length(b)), col="blue",lwd=3)
  }
}

# evenly-spaced breaks
break.equal <- function(x,n=2) {
  r <- range(x,na.rm=T)
  seq(0,1,by=1/n)*diff(r) + r[1]
}

# this is "equal-frequency" binning


#' Equal-count breaks of a dataset
#'
#' Computes a set of breaks using the equal-count algorithm.
#'
#'
#' @author Tom Minka
break.quantile <- function(x,n=2,plot=F,pretty=F) {
  if(is.list(x)) {
    b = lapply(x, function(y) break.quantile(y,n=n,pretty=pretty))
    return(b)
  }
  probs = seq(0,1,len=n+1)
  x = unique(na.omit(x))
  b <- as.numeric(quantile(x, probs))
  if(pretty) {
    for(i in 1:length(b)) {
      for(d in 1:4) {
        b.new = signif(b[i],d)
        p.new = mean(x <= b.new)
        if(probs[i] == 0 || probs[i] == 1) {
          ok = (p.new == probs[i])
        } else {
          ok = (abs(p.new - probs[i]) < 0.1/n)
        }
        #cat("i=",i,"d=",d,"ok=",ok,"p.new=",p.new,"probs[i]=",probs[i],"\n");
        if(ok) {
          b[i] = b.new
          break
        }
      }
    }
  }
  if(plot) {
    hb <- break.equal(x,40)
    h <- hist(x,hb,freq=FALSE,col="bisque",main="",xlab="")
    plot_breaks(b)
  }
  b
}

# convert a vector of cluster centers into a break vector
break.centers <- function(x,m) {
  m <- sort(m)
  n <- length(m)
  b <- (m[1:n-1]+m[2:n])/2
  r <- range(x)
  c(r[1]-eps,b,r[2])
}

#' @export
break_kmeans <- function(x,n=2,plot=T) {
  km <- kmeans(x,n)
  m <- as.numeric(km$centers)
  ss <- sum(km$withinss)
  if(plot) cat("sum of squares =", format(ss), "\n")
  b <- break.centers(x,m)
  if(plot) {
    hb <- break.equal(x,40)
    h <- hist(x,hb,freq=FALSE,col="bisque",main="",xlab="")
    plot_breaks(b)
  }
  b
}

# divide the range of x into n pieces, based on the largest spaces
# returns a vector of breaks, suitable for input to "hist" or "cut"
# tends to overdivide low-density regions
break.diff <- function(x,n=2) {
  n <- n-1
  len <- length(x)
  x <- sort(x)
  diff <- x[2:len] - x[1:len-1]
  i <- order(-diff)
  mid <- (x[2:len] + x[1:len-1])/2
  b <- sort(mid[i[1:n]])

  split.screen(c(2,1))
  split.screen(c(1,2),screen=1)
  screen(4)
  par(mar=c(2,4,1,1))
  last <- n+10
  plot(diff[i[1:last]],ylab="diff")
  m <- (diff[i[n]] + diff[i[n+1]])/2
  #arrows(1,m, last,m, col="black",code=0,lty=2)
  segments(n+0.5, diff[i[last]], n+0.5, diff[i[1]], col="black",lty=2)

  screen(2)
  par(mar=c(2,4,1,1))
  plot(x)
  segments(0,b, len+1,b, col="blue",lty=2)

  screen(3)
  par(mar=c(2,4,1,1))
  hist(x,30,col="bisque",main="")
  segments(b,rep(0,n), b,rep(100,n), col="blue")
  close.screen(all = TRUE)

  c(x[1]-eps,b,x[len])
}

# convert a cut of x into breaks
# f is a factor which divides x into disjoint groups
# this is only possible for certain kinds of cuts
break.cut <- function(x,f) {
  low <- sort(as.numeric(tapply(x,f,min)))
  high <- sort(as.numeric(tapply(x,f,max)))
  n <- length(high)
  b <- (low[2:n]+high[1:n-1])/2
  c(low[1]-eps,b,high[n])
}

scatter <- function(x) {
  if(is.null(dim(x))) sum((x - mean(x))^2)
  else sum((x - rep.row(colMeans(x),nrow(x)))^2)
}
sum.of.squares <- function(x,q) {
  # x is matrix, q is factor
  q = as.factor(q)
  n = table(q)
  f = indicators.factor(q)
  avg = scale.rows(f,1/n)
  x = as.matrix(x)
  m = avg %*% x
  # m is levels by dims
  x = x - (t(f) %*% m)
  sum(x*x)
}

#' Create a hierarchy by Ward's method
#'
#' Produces a hierarchical clustering of one-dimensional data via Ward's method.
#'
#' @param x a numerical vector, or a list of vectors.
#' @param n if x is a vector of cluster means, n is the size of each cluster.
#' @param s if x is a vector of cluster means, s is the sum of squares in each
#'     cluster.  only needed if \code{same.var=F}.
#' @param sortx if \code{sortx=F}, only clusters which are adjacent in \code{x}
#'     can be merged.  Used by \code{\link{break.ts}}.
#' @param same.var if \code{same.var=T}, clusters are assumed to have the same
#'   true variance, otherwise not.  This affects the cost function for merging.
#' @details Repeatedly merges clusters in order to minimize the clustering cost.
#'   By default, it is the same as \code{hclust(method="ward")}.
#'   If \code{same.var=T}, the cost is the sum of squares:
#'     \deqn{sum_c sum_{i in c} (x_i - m_c)^2}
#'   The incremental cost of merging clusters ca and cb is
#'   \deqn{(n_a*n_b)/(n_a+n_b)*(m_a - m_b)^2}
#'   It prefers to merge clusters which are small and have similar means.
#'
#'   If \code{same.var=F}, the cost is the sum of log-variances:
#'     \deqn{sum_c n_c*log(1/n_c*sum_{i in c} (x_i - m_c)^2)}
#'   It prefers to merge clusters which are small, have similar means,
#'   and have similar variances.
#'
#'   If \code{x} is a list of vectors, each vector is assumed to be a
#'   cluster.  \code{n} and \code{s} are computed for each cluster and
#'   \code{x} is replaced by the cluster means.
#'   Thus you can say \code{ward(split(x,f))} to cluster the data for different
#'   factors.
#'
#' @return The same type of object returned by \code{\link{hclust}}.
#' @author Tom Minka
#' @section Bugs: Because of the adjacency constraint used in implementation,
#'   the clustering that results
#'   from \code{sortx=T} and \code{same.var=F} may occasionally be suboptimal.
#' @seealso
#'   \code{\link{hclust}},
#'   \code{\link{plot_hclust_trace}},
#'   \code{\link{hist.hclust}},
#'   \code{\link{boxplot.hclust}},
#'   \code{\link{break_ward}},
#'   \code{\link{break.ts}},
#'   \code{\link{merge_factor}}
#' @examples
#' x <- c(rnorm(700,-2,1.5),rnorm(300,3,0.5))
#' hc <- ward(x)
#' opar <- par(mfrow=c(2,1))
#' # use dev.new() in RStudio
#' plot_hclust_trace(hc)
#' hist(hc,x)
#' par(opar)
#'
#' x <- c(rnorm(700,-2,0.5),rnorm(1000,2.5,1.5),rnorm(500,7,0.1))
#' hc <- ward(x)
#' opar <- par(mfrow=c(2,1))
#' plot_hclust_trace(hc)
#' hist(hc,x)
#' par(opar)
#'
#' data(OrchardSprays)
#' x <- OrchardSprays$decrease
#' f <- factor(OrchardSprays$treatment)
#' # shuffle levels
#' #lev <- levels(OrchardSprays$treatment)
#' #f <- factor(OrchardSprays$treatment,levels=sample(lev))
#' hc <- ward(split(x,f))
#' # is equivalent to:
#' #n <- tapply(x,f,length)
#' #m <- tapply(x,f,mean)
#' #s <- tapply(x,f,var)*n
#' #hc <- ward(m,n,s)
#' boxplot(hc,split(x,f))
#' @export
ward <- function(x,n=rep(1,length(x)),s=rep(1,length(x)),
                 sortx=TRUE,same.var=T) {
  if(is.list(x)) {
    n <- sapply(x,length)
    s <- sapply(x,scatter)
    #if(!same.var) s <- s + 0.1
    x <- sapply(x,mean)
  }
  # don't sort and it becomes time-series clustering
  if(sortx) {
    ord <- order(x)
    x <- x[ord]
    n <- n[ord]
    s <- s[ord]
  } else {
    ord <- 1:length(x)
  }
  height <- c()
  label <- -(1:length(x))
  label <- label[ord]
  merge <- matrix(0,length(x)-1,2)
  if(length(x) == 1) ks <- c()
  else ks <- 1:(length(x)-1)
  for(k in ks) {
    len <- length(x)
    # here we exploit the fact that the data is one-dimensional.
    # only clusters which are adjacent in sorted order need be considered
    # for merging.
    x1 <- x[1:len-1]
    x2 <- x[2:len]
    n1 <- n[1:len-1]
    n2 <- n[2:len]
    # merging cost
    cost <- n1*n2/(n1+n2)*(x1-x2)^2
    if(!same.var) {
      s1 <- s[1:len-1]
      s2 <- s[2:len]
      cost <- (n1+n2)*log((s1+s2+cost)/(n1+n2)) - n1*log(s1/n1) - n2*log(s2/n2)
    }
    j <- which.min(cost)

    if(!same.var) {
      s[j] <- s[j] + s[j+1] +
              n[j]*n[j+1]/(n[j]+n[j+1])*(x[j]-x[j+1])^2
      s <- not(s,j+1)
    }

    x[j] <- (x[j]*n[j] + x[j+1]*n[j+1])/(n[j] + n[j+1])
    n[j] <- n[j] + n[j+1]
    x <- not(x,j+1)
    n <- not(n,j+1)

    height[k] <- min(cost)
    merge[k,1] <- label[j]
    merge[k,2] <- label[j+1]
    label[j] <- k
    label <- not(label,j+1)
  }
  hc <- list(merge=merge,height=cumsum(height),order=ord,labels=NULL,
             method="ward")
  class(hc) <- "hclust"
  hc
}


#' Plot a merging trace
#'
#' Plots the cost of successive merges in a hierarchical clustering
#' @param h an \code{hclust} object
#' @param ka vector of the cluster cardinalities to plot
#' @return
#'   The trace shows, for each merge, the
#'   distance of the clusters which were merged.  This is useful for determining
#'   the appropriate number of clusters.  An interesting number of clusters is one
#'   that directly precedes a sudden jump in distance.
#' @author Tom Minka
#' @seealso
#'   \code{\link{ward}}, \code{\link{break_ward}}
#' @export
plot_hclust_trace <- function(h,k=1:10) {
  g <- c(rev(h$height),0)
  k <- k[k < length(g)]
  #g <- sqrt(g)
  if(h$method != "single") {
    g <- -diff(g)
  }
  plot(k,g[k],ylab="merging cost",xlab="clusters")
}



#' Quantize by clustering
#'
#' Quantize a one-dimensional variable by calling a clustering routine.
#'
#' These are convenience routines which simply call the appropriate clustering
#' routine (\code{\link{ward}}, \code{\link{hclust}}, or \code{\link{kmeans}}),
#' convert the output to a break vector, and make plots.
#'
#' @aliases break_ward break_kmeans break_hclust
#' @param x a numerical vector
#' @param n the desired number of bins
#' @param method argument given to \code{\link{hclust}}
#' @param plot If TRUE, a histogram with break lines is plotted
#' (\code{\link{hist.hclust}} or \code{\link{plot_breaks}}). For
#' \code{break_ward} and \code{break_hclust}, also shows a merging trace
#' (\code{\link{plot_hclust_trace}}).
#' @return A break vector, suitable for use in \code{\link{cut}} or
#' \code{\link{hist}}.
#' @author Tom Minka
#' @export
#' @examples
#'
#' x <- c(rnorm(700,-2,1.5),rnorm(300,3,0.5))
#' break_ward(x,2)
#' break_hclust(x,2,method="complete")
#' break_kmeans(x,2)
#'
#' x <- c(rnorm(700,-2,0.5),rnorm(1000,2.5,1.5),rnorm(500,7,0.1))
#' break_ward(x,3)
#' break_hclust(x,3,method="complete")
#' break_kmeans(x,3)
#'
break_ward <- function(x,n=2,plot=T) {
  h <- ward(x)
  q <- cutree(h,n)
  ss <- sum(tapply(x,q,scatter))
  if(plot) cat("sum of squares =", format(ss), "\n")
  b <- break.cut(x,q)

  if(plot) {
    split.screen(c(2,1))
    plot_hclust_trace(h)
    screen(2)
    hist(h,x,k=2:n)
    close.screen(all=TRUE)
  }
  b
}

# "ward" is best
# "ave","centroid" are good
# "single", "complete","median","mcq" are bad
#' @export
break_hclust <- function(x,n=2,method="ward",plot=T) {
  h <- hclust(dist(x)^2,method)
  q <- cutree(h,n)
  ss <- sum(tapply(x,q,scatter))
  cat("sum of squares =", format(ss), "\n")
  b <- break.cut(x,q)

  if(plot) {
    split.screen(c(2,1))
    plot_hclust_trace(h)
    screen(2)
    hist(h,x,k=2:n)
    close.screen(all=TRUE)
  }
  b
}

name.clusters <- function(f,y) {
  # given cluster memberships f and a factor y, returns the majority level
  # in each cluster
  tab = table(f,y)
  #tab = scale.cols(tab,1/colSums(tab))
  #tab = scale.rows(tab,1/rowSums(tab))
  make.unique(levels(y)[apply(tab,1,which.max)])
}

ward.kmeans <- function(x,k=2) {
  hc <- hclust(dist(x)^2,method="ward")
  f <- factor(cutree(hc,k=k))
  m = prototypes(x,f)
  cl = kmeans(x,m,iter.max=100)
  f = factor(cl$cluster)
  names(f) = rownames(x)
  ss = sum.of.squares(x,f)
  cat("sum of squares =", format(ss), "\n")
  f
}

retry.kmeans <- function(x,k=2,n=100) {
  best.ss = Inf
  for(iter in 1:n) {
    cl = kmeans(x,k,iter.max=100)
    f = factor(cl$cluster)
    ss = sum.of.squares(x,f)
    if(ss < best.ss) {
      best.ss = ss
      best.f = f
    }
  }
  cat("sum of squares =", format(best.ss), "\n")
  f = best.f
  names(f) = rownames(x)
  f
}

kernel.kmeans <- function(x,k=2,s=1) {
  # what if used 1/rank?
  # is there an explanation like similarity templates?
  d <- distances(x)^2
  s = s*sqrt(mean(diag(cov(x))))
  d = d/s
  x2 <- exp(-d)
  x2 = div_by_sum(x2)
  x2 = as.data.frame(x2)
  hc <- hclust(dist(x2)^2,method="ward")
  f <- factor(cutree(hc,k=k))
  m = prototypes(x2,f)
  cl = kmeans(x2,m,iter.max=100)
  f = factor(cl$cluster)
  ss = sum.of.squares(x2,f)
  cat("sum of squares =", format(ss), "\n")
  f
}

cluster.errors <- function(y,f) {
  # returns the number of pairs which are grouped differently in y than f
  n = length(y)
  tab = table(y,f)
  err = sum(rowSums(tab)^2) + sum(colSums(tab)^2) - 2*sum(tab^2)
  #err/2/n
  err/2
}
cluster.errors <- function(y,f) {
  tab = table(y,f)
  err = sum(apply(tab,2,function(s) sum(s)-max(s)))
  err
}

