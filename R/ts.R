# time series changepoints

# convert a cut of x into a changepoint plot
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
#' Apply Ward's method to find changepoints in a time-series.
#' Divide a time-series into homogenous segments.
#' @details
#' Calls \code{\link{ward}} with \code{sortx=F} to cluster the series into
#' segments.  Only the marginal distribution of data is used;
#' temporal smoothness, for example, is ignored.
#' @param x a numerical vector or \code{ts} object
#' @param n the desired number of segments
#' @param trace If TRUE, shows a merging trace via
#'  \code{\link{plot.hclust.trace}}
#' @param same.var argument passed to \code{\link{ward}}
#' @param ...
#'
#' @return A vector of time breaks. The breaks are also plotted visually via
#' \code{\link{plot.segments.ts}}
#' @seealso
#'  \code{\link{plot.segments.ts}}, \code{\link{plot.breaks}}
#' @export
#'
#' @examples
#' library(ts)
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
#' plot.breaks(b)
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
    plot.hclust.trace(h)
    screen(2)
  }
  plot.segments.ts(x,b,...)
  #points(x)
  if(trace) close.screen(all=TRUE)
  b
}

#############################################################################

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
merge.factor <- function(f,x,n,same.var=T,trace=T,xlab=NA,ylab=NA) {
  if(is.na(xlab)) xlab <- deparse(substitute(f))
  if(is.na(ylab)) ylab <- deparse(substitute(x))
  ord <- is.ordered(f)
  if(!ord) f <- sort.levels(f,x,fun=mean)
  # drop unused categories
  f <- factor(f)
  h <- ward(split(x,f),sortx=!ord,same.var=same.var)
  if(n > length(h$height)) q <- 1:(length(h$height)+1)
  else q <- cutree(h,n)
  nam <- tapply(levels(f),q,function(x) merge.names(x,ordered=ord))
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
    plot.hclust.trace(h)
    screen(2)
  }
  boxplot.hclust(h,split(x,of),xlab=xlab,ylab=ylab)
  f
}



