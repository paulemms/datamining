# time series changepoints

# convert a cut of x into a changepoint plot


#' Plot time-series segments
#'
#' Shows how the mean level of a time-series varies between
#'   segments.
#' @param x a numeric vector or \code{ts} object
#' @param b a numeric vector of break times
#' @return Plots the time series as a line, then draws a blue horizontal line
#'   through each segment, at the mean value of the segment.
#' @author Tom Minka
#' @seealso
#'   \code{\link{break.ts}},
#'   \code{\link{plot_breaks}}
#' @export
plot.segments.ts <- function(x,b,lwd=2,...) {
  plot(x,type="l",xlab="time",ylab="level")
  f <- cut(time(x),b,include.lowest=T)
  m <- tapply(as.numeric(x),f,mean)
  low <- tapply(time(x),f,min)
  high <- tapply(time(x),f,max)
  segments(low, m, high, m, col="blue",lwd=lwd,...)
}




#' Change-point analysis by clustering
#'
#' Divide a time-series into homogenous segments.
#'
#' Calls \code{\link{ward}} with \code{sortx=F} to cluster the series into
#' segments.  Only the marginal distribution of data is used; temporal
#' smoothness, for example, is ignored.
#'
#' @param x a numerical vector or \code{ts} object
#' @param n the desired number of segments
#' @param trace If TRUE, shows a merging trace via
#' \code{\link{plot_hclust_trace}}
#' @param same.var argument passed to \code{\link{ward}}
#' @return A vector of time breaks.  The breaks are also plotted visually via
#' \code{\link{plot.segments.ts}}.
#' @author Tom Minka
#' @seealso \code{\link{plot.segments.ts}}, \code{\link{plot_breaks}}
#' @examples
#'
#' data(LakeHuron)
#' # single major change
#' break.ts(LakeHuron,2)
#' # merging trace suggests n=6 is also interesting:
#' break.ts(LakeHuron,6)
#' # interesting oscillation
#'
#' data(treering)
#' break.ts(treering[1:500],9,same=T)
#' break.ts(treering[1:100],7,same=T)
#' # interesting multiscale structure
#'
#' x <- c(rnorm(100),rnorm(300)*3,rnorm(200)*2)
#' b <- break.ts(x,3,same=F)
#' plot(x,type="l")
#' plot_breaks(b)
#'
#' @export
break.ts <- function(x,n=2,trace=T,same.var=T,...) {
  h <- ward(as.numeric(x),sortx=F,same.var=same.var)
  q <- cutree(h,n)
  b <- break.cut(time(x),q)
  #b <- c(0,which(diff(q)>0),length(x)+1)
  ss <- sum(tapply(x,q,scatter))
  cat("sum of squares =", format(ss), "\n")

  if(length(h$height) < 2) trace <- F
  if(trace) {
    split.screen(c(2,1))
    plot_hclust_trace(h)
    screen(2)
  }
  plot.segments.ts(x,b,...)
  #points(x)
  if(trace) close.screen(all=TRUE)
  b
}


#' Boxplot with hierarchical cluster breaks
#'
#' A representation of a hierarchical clustering of predefined groups
#'
#'
#' @param h an \code{hclust} object
#' @param x the list of vectors that was clustered to produce \code{h}
#' (typically via \code{\link{ward}})
#' @param k a vector of the cluster cardinalities to plot
#' @param ... arguments passed to \code{boxplot}
#' @return A boxplot of \code{x} is shown with blue lines cutting the x-axis.
#' The tallest lines correspond to divisions made at the top of the hierarchy.
#' By reading top to bottom, you can see how each cluster is subdivided. This
#' can be far more illuminating than a plot of the hierarchy as a tree.
#' @author Tom Minka
#' @exportS3Method
#' @seealso \code{\link{hist.hclust}}, \code{\link{ward}},
#' \code{\link{break_ward}}
boxplot.hclust <- function(hc,x,k=2:5,col="bisque",...) {
  x <- x[hc$order]
  boxplot(x,col=col,...)
  r <- range(x)
  opar <- par(lwd=3)
  on.exit(par(opar))
  for(i in 1:length(k)) {
    if(k[i] > length(hc$height)) q <- 1:(length(hc$height)+1)
    else q <- cutree(hc,k[i])[hc$order]
    b <- which(diff(q)!=0)+0.5
    b <- c(0,b,length(x)+1)
    scale <- (length(k)-i+1)/length(k)
    segments(b,rep(r[1],length(b)), b,rep(r[2]*scale,length(b)), col="blue")
  }
}

# breaks a factor into nbins bins in order to preserve the prediction of x


#' Merge factor levels
#'
#' Merges factor levels with similar response distributions (assumed normal).
#'
#' Calls \code{ward(split(x,f))} to get a tree, cuts the tree, and constructs a
#' new factor.  The tree is shown via \code{\link{boxplot.hclust}}.
#'
#' @param f a factor
#' @param x a numerical vector, same length as \code{f}
#' @param n the desired number of factor levels
#' @param same.var argument passed to \code{\link{ward}}
#' @param trace If TRUE, a merging trace is plotted
#' (\code{\link{plot_hclust_trace}})
#' @param xlab,ylab axis labels.  If NA, taken from f and x arguments.
#' @return A new factor, same length as \code{f}, but with \code{n} levels.
#' @author Tom Minka
#' @export
#' @examples
#'
#' n <- 20
#' x <- c(rnorm(n)+1, rnorm(n)+2, rnorm(n)*4+2)
#' f <- gl(3,n)
#' levels(f) <- c("a","b","c")
#' merge_factor(f,x,2,same.var=T)
#' merge_factor(f,x,2,same.var=F)
#'
#' # an ordered factor
#' data(va.deaths)
#' merge_factor(va.deaths$Age,va.deaths$Rate,2)
#'
merge_factor <- function(f,x,n,same.var=T,trace=T,xlab=NA,ylab=NA) {
  if(is.na(xlab)) xlab <- deparse(substitute(f))
  if(is.na(ylab)) ylab <- deparse(substitute(x))
  ord <- is.ordered(f)
  if(!ord) f <- sort_levels(f,x,fun=mean)
  # drop unused categories
  f <- factor(f)
  h <- ward(split(x,f),sortx=!ord,same.var=same.var)
  if(n > length(h$height)) q <- 1:(length(h$height)+1)
  else q <- cutree(h,n)
  nam <- tapply(levels(f),q,function(x) merge.names(x,i=ord))
  of <- f
  levels(f) <- nam[q]
  if(same.var) {
    ss <- sum(tapply(x,f,scatter))
    cat("sum of squares =", format(ss), "\n")
  } else {
    n <- tapply(x,f,length)
    ss <- sum(n*log(tapply(x,f,scatter)/n))
    cat("sum of log-variance =", format(ss), "\n")
  }
  if(trace) {
    split.screen(c(2,1))
    on.exit(close.screen(all=TRUE))
    plot_hclust_trace(h)
    screen(2)
  }
  boxplot(h,split(x,of),xlab=xlab,ylab=ylab)
  f
}



