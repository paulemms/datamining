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
    orientations <- atan(ar * dy, dx)
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
# Tom Minka 9/19/01

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
merge.hist <- function(x,b=NULL,n=b,trace=T) {
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
  plot.histogram(h,col="bisque",freq=FALSE,xlab=xlab,main=main)
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

plot.breaks <- function(b) {
  r <- par("usr")[3:4]
  segments(b,rep(r[1],length(b)), b,rep(r[2],length(b)), col="blue", lwd=2)
}

hist.hclust <- function(hc,x,k=2:5) {
  hb <- break.equal(x,40)
  h <- hist(x,hb,freq=FALSE,col="bisque",main="",xlab="")
  top <- par("usr")[4]
  # for cutree
  library(mva)
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
    plot.breaks(b)
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

break.kmeans <- function(x,n=2,plot=T) {
  library(mva)
  km <- kmeans(x,n)
  m <- as.numeric(km$centers)
  ss <- sum(km$withinss)
  if(plot) cat("sum of squares =", format(ss), "\n")
  b <- break.centers(x,m)
  if(plot) {
    hb <- break.equal(x,40)
    h <- hist(x,hb,freq=FALSE,col="bisque",main="",xlab="")
    plot.breaks(b)
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

# Create a hierarchy by Ward's method
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

plot.hclust.trace <- function(h,k=1:10) {
  g <- c(rev(h$height),0)
  k <- k[k < length(g)]
  #g <- sqrt(g)
  if(h$method != "single") {
    g <- -diff(g)
  }
  plot(k,g[k],ylab="merging cost",xlab="clusters")
}

break.ward <- function(x,n=2,plot=T) {
  h <- ward(x)
  library(mva)
  q <- cutree(h,n)
  ss <- sum(tapply(x,q,scatter))
  if(plot) cat("sum of squares =", format(ss), "\n")
  b <- break.cut(x,q)

  if(plot) {
    split.screen(c(2,1))
    plot.hclust.trace(h)
    screen(2)
    hist.hclust(h,x,k=2:n)
    close.screen(all=TRUE)
  }
  b
}

# "ward" is best
# "ave","centroid" are good
# "single", "complete","median","mcq" are bad
break.hclust <- function(x,n=2,method="ward",plot=T) {
  library(mva)
  h <- hclust(dist(x)^2,method)
  q <- cutree(h,n)
  ss <- sum(tapply(x,q,scatter))
  cat("sum of squares =", format(ss), "\n")
  b <- break.cut(x,q)

  if(plot) {
    split.screen(c(2,1))
    plot.hclust.trace(h)
    screen(2)
    hist.hclust(h,x,k=2:n)
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
  x2 = div.by.sum(x2)
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

#############################################################################
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

# Apply Ward's method to find changepoints in a time-series
break.ts <- function(x,n=2,trace=T,same.var=T,...) {
  h <- ward(as.numeric(x),sortx=F,same.var=same.var)
  library(mva)
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
  library(mva)
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
  library(mva)
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
  i = (is.na(x) | is.na(y) | is.na(labels))
  x = x[!i]; y = y[!i]; labels = labels[!i]
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
  hull[,ix] <- (hull[,ix] - adj[1])*rep.col(lab.w,4)
  hull[,iy] <- (hull[,iy] - adj[2])*rep.col(lab.h,4)
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

color.plot <- function(object, ...) UseMethod("color.plot")
color.plot.formula <- function(formula,data=parent.frame(),...) {
  x <- model.frame.default(formula,data,na.action=na.pass)
  # put response at end (my preferred ordering)
  x = cbind(x[-1],x[1])
  color.plot.data.frame(x,...)
}
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
    w = strwidth(labels,units="inch",cex=cex)
    h = strheight(labels,units="inch",cex=cex)
    n = length(x)
    s = inch.symbol(pch=pch,cex=cex) + 0.1
    w2 = c(w,rep(s,n))
    h2 = c(h,rep(s,n))
    x2 = c(x,x)
    y2 = c(y,y)
    cost = c(rep(1,n),rep(0,n))
    r = move.collisions2(inchx(x2),inchy(y2),w2,h2,cost)
    points(x,y,cex=cex,pch=pch,...)
    r = r[1:n,]
    x = xinch(r[,1])
    y = yinch(r[,2])
  }
  text(x,y,labels=labels,cex=cex,adj=adj,srt=srt,...)
}
# test: text.plot(runif(10),runif(10),rep("a",10))

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

##############################################################################

rgb2name <- function(r,g,b) {
  if(missing(g)) {
    r <- as.matrix(r)
    if(ncol(r) != 3) {
      if(nrow(r) == 3) r <- t(r)
      else stop("matrix must have 3 rows or 3 columns")
    }
  } else {
    r <- array(c(r,g,b),c(length(r),3))
  }
  if(max(r) > 1) r <- r/255
  cc <- colors()
  # remove duplicates
  #cc <- setdiff(cc, c("gray100","grey100"))
  crgb <- col2rgb(cc)/255
  # distance matrix
  d <- sqdist(r,t(crgb))
  cc[apply(d, 1, which.min)]
}

col2hsv <- function(col) {
  rgb2hsv(col2rgb(col))
}

rgb2hsv <- function(r,g,b) {
  matrix.out = F
  if(missing(g)) {
    if(is.list(r)) {
      b = r[[3]]
      g = r[[2]]
      r = r[[1]]
    } else {
      matrix.out = T
      if(nrow(r) != 3) {
        if(ncol(r) == 3) r <- t(r)
        else stop("matrix must have 3 rows or 3 columns")
      }
      b <- r[3,]
      g <- r[2,]
      r <- r[1,]
    }
  }
  n <- length(r)
  lim <- array(0,c(2,n))
  lim[1,] <- pmin(r,g,b)
  lim[2,] <- pmax(r,g,b)
  rang <- lim[2,]-lim[1,]
  v <- lim[2,]
  s <- v
  i <- (lim[2,] > 0)
  s[!i] <- 0
  s[i] <- rang[i]/lim[2,i]
  h <- v
  h[s == 0] <- 0
  i <- (s > 0) & (r == lim[2,])
  h[i] <- (g[i] - b[i])/rang[i]/6
  i <- (s > 0) & (g == lim[2,])
  h[i] <- (2 + (b[i]-r[i])/rang[i])/6
  i <- (s > 0) & (b == lim[2,])
  h[i] <- (4 + (r[i]-g[i])/rang[i])/6
  i <- (h < 0)
  h[i] <- h[i] + 1
  if(matrix.out)
    t(array(c(h,s,v),c(n,3),list(NULL,c("hue","saturation","value"))))
  else
    data.frame(hue=h,saturation=s,value=v)
}

test.col2hsv <- function() {
  h <- col2hsv(palette())
  h[3,] = h[3,]/255
  print(palette())
  print(sapply(1:ncol(h), function(i) hsv(h[1,i],h[2,i],h[3,i])))
}

hsv2cone <- function(h,s,v) {
  if(missing(s)) {
    if(nrow(h) != 3) {
      if(ncol(h) == 3) h <- t(h)
      else stop("matrix must have 3 rows or 3 columns")
    }
    v <- h[3,]
    s <- h[2,]
    h <- h[1,]
  }
  x <- s*v*cos(h*2*pi)
  y <- s*v*sin(h*2*pi)
  z <- v
  n <- length(x)
  t(array(c(x,y,z),c(n,3),list(NULL,c("x","y","z"))))
}
cone2hsv <- function(x,y,z) {
  if(missing(y)) return(cone2hsv(x[1,],x[2,],x[3,]))
  v <- z
  h <- (atan2(y,x)/2/pi + 1) %% 1
  s <- sqrt(x^2 + y^2)/v
  i <- which(s > 1 | is.na(s))
  if(length(i) > 0) s[i] <- 0
  r <- rbind(h,s,v)
  rownames(r) <- c("hue","saturation","value")
  r
}

# this is only for routines used by minka library

round.partition <- function(x) {
  # rounds the elements of x, while preserving the sum
  diff(c(0,round(cumsum(x))))
}
rand.partition <- function(n) {
  # sample from hypergeometric
  # n is vector of desired table counts
  sample(rep(1:length(n),n))
}

clean <- function(x,threshold=0.5) {
  mostly.na = function(x) (mean(is.na(x))>threshold)
  bad = sapply(x,mostly.na)
  # remove variables with many NAs
  x = x[!bad]
  # remove cases with any remaining NAs
  x = na.omit(x)
  # recompute factor levels
  for(v in colnames(x)) {
    if(is.factor(x[,v])) x[,v] = factor(x[,v])
  }
  x
}

na.dummy <- function(x,value=666) {
  # replace NAs with a dummy value
  if(is.data.frame(x)) return(apply.df(x,na.dummy))
  x[which(is.na(x))] = value
  x
}

rank.stable <- function(x) {
  # handles ties better than "rank"
  n <- length(x)
  r <- rep(0,n)
  r[order(x)] <- 1:n
  r
}

shuffle <- function(x) {
  if(length(dim(x)) == 2) {
    x[order(runif(nrow(x))),order(runif(ncol(x)))]
  } else {
    x[order(runif(length(x)))]
  }
}

named.list <- function(...,names.) {
  lst <- list(...)
  names(lst) = names.
  lst
}
named.numeric <- function(names.) {
  x <- numeric(length(names.))
  names(x) = names.
  x
}
as.named.vector <- function(a) {
  if(is.matrix(a) || is.data.frame(a)) {
    if(ncol(a) == 1) {
      x = a[,1]
      names(x) = rownames(a)
    } else {
      x = as.vector(as.matrix(a))
      #names(x) = outer(rownames(a),colnames(a),FUN=function(a,b) paste(a,b,sep="."))
      names(x) = outer(rownames(a),colnames(a),FUN=paste)
    }
    x
  } else if(is.list(a)) {
    unlist(a)
  } else stop("don't know what to do")
}

# companions to xinch,yinch
inchx = function(x,warn.log=T) {
  if (warn.log && par("xlog"))
    warning("x log scale:  inchx() is non-sense")
  x / diff(par("usr")[1:2]) * par("pin")[1]
}
inchy = function (y, warn.log = TRUE)
{
  if (warn.log && par("ylog"))
    warning("y log scale:  inchy() is non-sense")
  y / diff(par("usr")[3:4]) * par("pin")[2]
}

inch.symbol <- function(pch=par("pch"),cex=par("cex")) {
  # returns the height in inches of the given plotting symbol
  # with current graphics settings.
  # this is for pch="H"
  inch = par("csi")*cex/1.5
  if(pch <= 20) inch = inch/2
  inch
}

format.bytes = function(x) {
  # returns byte sizes in pretty form
  # example: format.bytes(c(2,1030,2e6))
  # based on code by Fang Chen
  round.to.half = function(x) round(x*2)/2
  s = as.character(x)
  names(s) = names(x)
  kilo = 2^10; mega = 2^20
  i = (x >= mega)
  s[i] = paste(round.to.half(x[i]/mega),"M",sep="")
  i = (x >= kilo) & (x < mega)
  s[i] = paste(round(x[i]/kilo),"k",sep="")
  s
}
ls.sizes <- function(name=parent.frame(), pretty = T, top = 10, ...) {
  # lists the `top' biggest objects in the given environment
  # suggested by Fang Chen
  x <- sapply(ls(name,...),function(x) object.size(get(x,name)))
  if(top > length(x)) top = length(x);
  x = rev(rev(sort(x))[1:top])
  if(pretty) x = format.bytes(x)
  x <- as.data.frame(x)
  colnames(x) <- "bytes"
  x
}

# smallest number such that 1+eps > 1
eps <- 1e-15

# returns the elements of x not indexed by i
not <- function(x,i) {
  if(is.numeric(i)) return(x[-i])
  if(is.character(i)) return(x[setdiff(names(x),i)])
  if(is.logical(i)) return(x[!i])
}

sort.levels <- function(f,x,fun=median) {
  if(length(x) == length(f)) x <- tapply(x,f,fun)
  if(length(x) != length(levels(f))) stop("wrong number of values to sort by")
  factor(f,levels=levels(f)[order(x)])
}

# n is a vector of names
# i is a vector of indices
# merges names n[i] into a composite name
# returns a shorter vector of names
merge.names <- function(n,i=1:length(n)) {
  if(is.dimOrdered(n)) {
    first <- n[i[1]]
    last <- n[i[length(i)]]
    # what kind of range notation should we use?
    if(F) {
      split <- "\\.\\."; sep <- ".."
    } else {
      split <- "-"; sep <- "-"
    }
    s1 <- unlist(strsplit(first,split))
    s2 <- unlist(strsplit(last,split))
    # work around a bug in strsplit
    if(length(s2) == 1) s2 <- c("",s2)
    n[i] <- paste(s1[[1]], s2[[2]], sep=sep)
  } else {
    n[i] <- paste(n[i],collapse=",")
  }
  a <- attributes(n)
  n <- n[-i[-1]]
  attributes(n) <- a
  n
}

which2 <- function(x) {
  i <- which(x)-1
  j <- floor(i/nrow(x))
  i <- i - j*nrow(x)
  cbind(i+1,j+1)
}

which.max2 <- function(x) {
  i <- which.max(x)-1
  j <- floor(i/nrow(x))
  i <- i - j*nrow(x)
  c(i+1,j+1)
}

which.min2 <- function(x) {
  i <- which.min(x)-1
  j <- floor(i/nrow(x))
  i <- i - j*nrow(x)
  c(i+1,j+1)
}


zip <- function(a,b,fun="*") {
  fun <- match.fun(fun)
  r <- lapply(1:length(a), function(i) fun(a[[i]],b[[i]]))
  names(r) <- names(a)
  r
}

accumulate <- function(lst,fun="+") {
  fun <- match.fun(fun)
  r <- lst[[1]]
  for(i in 2:length(lst)) {
    #r <- r + lst[[i]]
    r <- fun(r,lst[[i]])
  }
  r
}

# returns split(x,f) where f is a random factor
# if p is a single number then f has two factors in the proportion (p,1-p)
rsplit <- function(x,p) {
  if(is.data.frame(x)) n <- nrow(x)
  else n <- length(x)
  n1 <- ceiling(p*n)
  i <- c(rep(1,n1),rep(2,n-n1))
  split(x,sample(i))
}

is.dimOrdered <- function(d) {
  # dimensions are ordered by default
  a <- attr(d,"ordered")
  is.null(a) || a
}
dimOrdered <- function(x) {
  dn = dimnames(x)
  # never return NULL
  if(is.null(dn)) dn = dim(x)
  sapply(dn,is.dimOrdered)
}
"dimOrdered<-" <- function(x,value) {
  d = 1:length(dim(x))
  if(length(value) == 1) value = rep(value,length(d))
  for(i in d) {
    attr(dimnames(x)[[i]],"ordered") <- value[i]
  }
  x
}

# convert from interval to range notation
interval2range <- function(x) {
  x <- unlist(strsplit(x,","))
  x[1] <- substr(x[1],2,nchar(x[1]))
  x[1] <- format(as.numeric(x[1]))
  if(x[1] == "-Inf") x[1] <- "."
  x[2] <- substr(x[2],1,nchar(x[2])-1)
  x[2] <- format(as.numeric(x[2]))
  if(x[2] == "Inf") x[2] <- "."
  paste(x[1],x[2],sep="-")
}

auto.signif <- function(x,tol=0.1) {
  # rounds x to minimum number of digits to achieve specified tolerance.
  # i.e. if rounded value is y then abs(y-x)/abs(x) <= tol.
  # also see break.quantile(pretty=T)
  for(i in 1:length(x)) {
    if(x[i] != 0) {
      for(d in 1:10) {
        y = signif(x[i],d)
        if(abs(y - x[i])/abs(x[i]) <= tol) break
      }
      x[i] = y
    }
  }
  x
}

auto.format <- function(x) {
  # use the minimal number of digits to make x's unique
  # similar to abbrev
  for(digits in 0:6) {
    s = formatC(x,digits=digits,format="f")
    if(all(duplicated(s) == duplicated(x))) break
  }
  s
}

# used by color.plot and predict.plot given
rcut <- function(x,b) {
  f <- cut(x,b,include.lowest=T)
  levels(f) <- sapply(levels(f),interval2range)
  f
}

gradient <- function(x) {
  n <- length(x)
  g <- (diff(x[1:(n-1)])+diff(x[2:n]))/2
  g = c(x[2]-x[1],g,x[n]-x[n-1])
  names(g) = names(x)
  g
}

scale.range <- function(x) {
  (x-min(x,na.rm=T))/(diff(range(x,na.rm=T))+eps)
}

cut.quantile <- function(x,n=3,suffix="") {
  #if(missing(suffix)) suffix = deparse(substitute(x))
  #x[is.na(x)] <- 0
  b <- quantile(unique(x),probs=seq(0,1,len=n+1),na.rm=T)
  if(length(suffix) == 1) {
    if(n == 5) labels <- c("low","med.low","med","med.high","high")
    else if(n == 4) labels <- c("low","med.low","med.high","high")
    else if(n == 3) labels <- c("low","med","high")
    else if(n == 2) labels <- c("low","high")
    else labels <- paste("cut",1:n,sep="")
    labels = paste(labels,suffix)
  } else labels = suffix
  cut(x,b,include=T,labels=labels)
}

indicators.factor <- function(y) {
  # convert a factor into a matrix of indicators
  # result is level by case
  # works if y contains NAs
  r <- array(FALSE,c(length(levels(y)),length(y)),list(levels(y),names(y)))
  for(lev in levels(y)) r[lev,y==lev] <- TRUE
  r
}

reorder.rows <- function(y,i) {
  # indexes a matrix while preserving attributes
  if(is.na(nrow(y))) {
    y[i]
  } else {
    a <- attributes(y); y <- y[i,]
    if(!is.null(a$dimnames[[1]])) a$dimnames[[1]] <- a$dimnames[[1]][i]
    attributes(y) <- a
    y
  }
}
reorder.cols <- function(y,i) {
  # indexes a matrix while preserving attributes
  if(is.na(ncol(y))) {
    y[i]
  } else {
    a <- attributes(y); y <- y[,i]
    if(!is.null(a$dimnames[[2]])) a$dimnames[[2]] <- a$dimnames[[2]][i]
    attributes(y) <- a
    y
  }
}

##############################################################################
# Graphics tools

auto.layout <- function(len,asp=1) {
  # asp is the aspect ratio of one cell
  layout <- c(1,len)
  din <- par("din")
  asp <- din[2]/din[1]/asp
  if(asp > 1) {
    layout[2] <- round(sqrt(len/asp))
    layout[1] <- ceiling(len/layout[2])
  } else {
    layout[1] = round(sqrt(len*asp))
    layout[2] = ceiling(len/layout[1])
  }
  layout
}

auto.legend <- function(labels, col=1, lty, pch, text.col, ...) {
  # barplot could use an auto.legend
  usr <- par("usr")
  top <- usr[4]
  left <- usr[2]-max(strwidth(labels))-xinch(0.04)
  y <- cumsum(strheight(labels)+yinch(0.04))
  if(missing(lty) && missing(pch)) {
    text.only <- T
    if(missing(text.col)) text.col <- col
  }
  for(i in 1:length(labels)) {
    text(left,top-y[i],labels[i],col=text.col[i],adj=0,...)
  }
}

auto.mar <- function(main="",sub="",xlab="x",ylab="y",axes=T,
                     xaxt=if(axes) par("xaxt") else "n",
                     yaxt=if(axes) par("yaxt") else "n",
                     las=par("las"),mar.scale=NULL,key=F,...) {
  # mar is c(bottom,left,top,right)
  mar = c(4.5,4,0,0.1)
  if(!is.null(main) && main != "") mar[3] = mar[3] + 1
  if(key) mar[3] = mar[3] + 1
  if(is.null(xlab) || xlab == "") {
    if(xaxt == "n") mar[1] = 0.05
    else mar[1] = 2.5
  }
  if(is.null(ylab) || ylab == "") {
    if(yaxt == "n") mar[2] = 0
    else mar[2] = 2
  }
  if(is.null(mar.scale)) mar.scale = 2
  if(las == 3 || las == 2) mar[1] = mar[1]*mar.scale
  if(las == 1 || las == 2) mar[2] = mar[2]*mar.scale
  if(!is.null(sub) && sub != "") {
    mar[1] = mar[1] + 2
  }
  mar
}

##############################################################################
# matlab compatibility

rep.mat <- function(x,n,m) {
  array(rep(x,n*m),c(length(x)*n,m))
}
rep.col <- function(x,n) {
  rep.mat(x,1,n)
}
rep.row <- function(x,n) {
  t(rep.col(x,n))
}

eye <- function(n) {
  diag(n)
}
norm <- function(x) {
  sqrt(sum(x^2))
}
# behaves like matlab "sum"
dim.sum <- function(x,dim) {
  apply(x,setdiff(1:length(dim(x)),dim),sum)
}
scale.cols <- function(x,s) {
  x * rep.row(s,nrow(x))
}
scale.rows <- function(x,s) {
  x * rep.col(s,ncol(x))
}

##############################################################################
# data.frames

formula.with.data <- function(fmla, data)
{
  # put data into the environment of formula
  # from Brian Ripley
  if(identical(as.character(fmla[[3]]),".")) {
    # dot must be expanded
    resp = response.var(fmla)
    pred = setdiff(names(data),resp)
    fmla = formula(paste(resp,"~",paste(pred,collapse="+")))
  }
  env <- new.env()
  for(i in names(data)) assign(i, data[[i]], envir=env)
  environment(fmla) <- env
  fmla
}

# sort the rows of a data frame (df) according to column f
# f can also be a vector with length(f)==nrow(df)
sort.data.frame <- function(df,f=ncol(df),...) {
  if(length(f) == 1) f <- df[[f]]
  else if(length(f) != nrow(df))
    stop("length of sort vector doesn't match nrows")
  df[order(f,...),,drop=F]
}

sort.cells <- function(x) {
  sort.data.frame(as.data.frame(x))
}

empty.data.frame <- function(col.names) {
  if(missing(col.names)) y <- NULL
  else {
    y <- sapply(col.names,function(x) NULL)
  }
  structure(y,class="data.frame")
}

rbind.extend <- function(df,df2) {
  # like rbind, but allows different numbers of columns (NAs inserted)
  x <- list()
  col.names <- unique(c(names(df),names(df2)))
  for(j in 1:length(col.names)) {
    k <- col.names[j]
    # must check in names first to avoid abbreviation match
    if(k %in% names(df)) v1 <- df[[k]]
    else v1 <- factor(rep(NA,nrow(df)))
    if(k %in% names(df2)) v2 <- df2[[k]]
    else v2 <- factor(rep(NA,nrow(df2)))
    # requires cfactor
    # force dispatch on second argument too
    x[[j]] <- cfactor(v1,v2)
  }
  names(x) <- col.names
  row.names <- make.unique(c(rownames(df),rownames(df2)))
  attr(x,"row.names") <- row.names
  class(x) <- "data.frame"
  x
}
rbind.extend <- function(df,df2) {
  if(nrow(df) == 0) return(df2)
  if(nrow(df2) == 0) return(df)
  not.1 <- setdiff(names(df2),names(df))
  not.2 <- setdiff(names(df),names(df2))
  if(length(not.1) > 0) df[not.1] <- rep(NA,nrow(df))
  if(length(not.2) > 0) df2[not.2] <- rep(NA,nrow(df2))
  col.names <- unique(c(names(df),names(df2)))
  x <- rbind(df[col.names],df2[col.names])
  row.names <- make.unique(c(rownames(df),rownames(df2)))
  attr(x,"row.names") <- row.names
  x
}

cbind.extend <- function(df,df2) {
  # like cbind.data.frame, but pads with NA until rownames match
  if(is.null(df) || nrow(df) == 0) return(df2)
  if(is.null(df2) || nrow(df2) == 0) return(df)
  not.1 = setdiff(rownames(df2),rownames(df))
  not.2 = setdiff(rownames(df),rownames(df2))
  if(length(not.1) > 0) df[not.1,] = NA
  if(length(not.2) > 0) df2[not.2,] = NA
  r = cbind(df,df2[rownames(df),])
  names(r) = c(names(df),names(df2))
  r
}

# returns the name of the response variable
# works for formulas, model frames, and fitted models
response.var <- function(object) {
  if(is.null(object)) return(NULL)
  if(inherits(object, "terms")) {
    a <- attributes(object)
    if(!a$response) return(character(0))
    return(as.character(a$variables[2]))
  }
  if(is.null(attr(object,"terms"))) {
    if(is.data.frame(object)) {
      # shortcut to avoid make.names
      return(names(object)[length(object)])
    }
    if(is.table(object)) return(response.var(terms(object)))
    #if(is.list(object) || is.vector(object) || is.array(object)) return(NULL)
  }
  response.var(terms(object))
}

# returns the variables out of which the predictors are constructed
predictor.vars <- function(object) {
  if(inherits(object, "terms")) {
    a <- attributes(object)
    pred = rownames(a$factors)[apply(a$factors,1,any)]
    # skip cross-products
    #pred = a$term.labels[a$order == 1]
    # drop I() terms
    #pred = pred[substring(pred,1,2) != "I("]
    return(pred)
  }
  if(inherits(object,"data.frame") && is.null(attr(object,"terms"))) {
    # shortcut to avoid make.names
    return(names(object)[-length(object)])
  }
  predictor.vars(terms(object))
}
# returns all terms on the rhs, including higher-order terms
predictor.terms <- function(object) {
  attr(terms(object),"term.labels")
}

as.data.frame.col <- function(x,n="V1") {
  frame = as.data.frame(as.matrix(x))
  names(frame) = n
  frame
}
as.data.frame.row <- function(x,row.name="") {
  extra.args <- list(check.names=F,row.names=row.name)
  do.call("data.frame",append(as.list(x),extra.args))
}

nocheck.data.frame <- function(...,row.names=NULL) {
  # must be done this way to prevent check.names (bug in data.frame)
  # example: data.frame(list("a b"=3),check.names=F)
  do.call("data.frame",append(as.list(...),list(row.names=row.names,check.names=F)))
}
apply.df <- function(x,fun) {
  y = lapply(x,fun)
  if(length(y[[1]]) == length(x[[1]]))
    nocheck.data.frame(y,row.names=rownames(x),check.names=F)
  else
    nocheck.data.frame(y,check.names=F)
}

my.model.frame <- function(...) {
  # bugfix for model.frame - puts response at end
  x = model.frame(...)
  a = attr(x,"terms")
  x = nocheck.data.frame(c(x[-1],x[1]),check.names=F,row.names=rownames(x))
  attr(x,"terms") = a
  x
}

set.response.var <- function(fmla,v) {
  # v is character string
  vars = c(response.var(fmla),predictor.terms(fmla))
  i = grep(v,vars)
  if(length(i)>0) vars = vars[-i[1]]
  formula(paste(v,paste(vars,collapse="+"),sep="~"))
}

# concatenate factors (should be built in)
# unfortunately "sort" requires the original c(), which drops labels
cfactor <- function(f1,f2=NULL) {
  if(length(f1) == 0) return(f2)
  else if(length(f2) == 0) return(f1)
  if(!is.factor(f1) || !is.factor(f2)) return(c(f1,f2))
  all.levs <- unique(c(levels(f1),levels(f2)))
  factor(c(as.character(f1),as.character(f2)),levels=all.levs)
}

# should be built into R
sum.data.frame <- function(x,...) sapply(x,sum,...)

##############################################################################
# bug fixes

which.environment <- function(x,...) {
  env = environment(...)
  while(!is.null(env) && !exists(x,envir=env,inherits=F)) {
    env = parent.env(env)
  }
  env
}

setInNamespace <- function(subx,x,ns=environment(x)) {
  if(!exists("asNamespace")) return(assign(subx,x,env=ns))
  # from fixInNamespace
  if (bindingIsLocked(subx, ns)) {
    unlockBinding(subx, ns)
    assign(subx, x, env = ns)
    w <- options("warn")
    on.exit(options(w))
    options(warn = -1)
    lockBinding(subx, ns)
  }
  else assign(subx, x, env = ns)
  if (isNamespace(ns) && !isBaseNamespace(ns)) {
    S3 <- getNamespaceInfo(ns, "S3methods")
    if (!length(S3))
      return(invisible(NULL))
    S3names <- sapply(S3, function(x) x[[3]])
    if (subx %in% S3names) {
      i <- match(subx, S3names)
      genfun <- get(S3[[i]][[1]])
      defenv <- if (typeof(genfun) == "closure")
        environment(genfun)
      else .BaseNamespaceEnv
      S3Table <- get(".__S3MethodsTable__.", envir = defenv)
      if (exists(subx, envir = S3Table, inherits = FALSE))
        assign(subx, x, S3Table)
    }
  }
}
replaceInNamespace <- function(name,value) {
  if(exists("asNamespace")) {
    environment(value) = environment(get(name))
    setInNamespace(name,value)
    setInNamespace(name,value,which.environment(name))
  } else {
    assign(name,value,env=.GlobalEnv)
  }
}

aggregate.data.frame <- function(x, by, FUN, ...) {
    if(!is.data.frame(x))
        x <- as.data.frame(x)
    if(!is.list(by))
        stop("`by' must be a list")
    if(is.null(names(by)))
        names(by) <- paste("Group", seq(along = by), sep = ".")
    else {
        nam <- names(by)
        ind <- which(nchar(nam) == 0)
        names(by)[ind] <- paste("Group", ind, sep = ".")
    }
    y <- lapply(x, tapply, by, FUN, ..., simplify = FALSE)
    if(any(sapply(unlist(y, recursive = FALSE), length) > 1))
        stop("`FUN' must always return a scalar")
    z <- y[[1]]
    d <- dim(z)
    w <- list()
    for (i in seq(along = d)) {
        j <- rep(rep(seq(1 : d[i]),
                     prod(d[seq(length = i - 1)]) * rep(1, d[i])),
                 prod(d[seq(from = i + 1, length = length(d) - i)]))
        #w <- cbind(w, dimnames(z)[[i]][j])
        # minka: preserve original levels
        f <- factor(dimnames(z)[[i]][j],levels=levels(by[[i]]))
        w[[names(by)[i]]] = f
    }
    # minka: not sure what this is for
    #w <- w[which(!unlist(lapply(z, is.null))), ]
    y <- data.frame(w, lapply(y, unlist, use.names = FALSE))
    #names(y) <- c(names(by), names(x))
    y
}
replaceInNamespace("aggregate.data.frame",aggregate.data.frame)

# desiderata for make.unique:
# 1. if A is unique, make.unique(c(A,B)) preserves A
# 2. make.unique(c(A,B)) = make.unique(c(make.unique(A),B))

# internal version
# does not satisfy desideratum #2
# make.unique(c("a","a","a")) != make.unique(c(make.unique(c("a","a")),"a"))
make.unique <- function(names) {
  # names is character vector
  if(!is.character(names)) stop("names must be a character vector")
  while(any(dups <- duplicated(names))) {
    names[dups] <- paste(names[dups], seq(length = sum(dups)), sep = "")
  }
  names
}
# satifies both desiderata
# make.unique(c("a","a","a.2","a")) ==
# make.unique(c(make.unique(c("a","a")),"a.2","a"))
make.unique <- function(names,sep=".",start=2) {
  # names is character vector
  if(!is.character(names)) stop("names must be a character vector")
  repeat {
    dups <- which(duplicated(names))
    if(length(dups) == 0) break
    # loop duplicates
    for(j in dups) {
      # for each duplicate, find the lowest value of cnt which makes it
      # different from previous names.
      i <- start
      repeat {
        newnam <- paste(names[j],i,sep=sep)
        # compare to previous elements only
        if(!any(is.element(newnam,names[1:(j-1)]))) break
        i <- i + 1
      }
      names[j] <- newnam
    }
    # repeat in case new duplicates have been created (see examples)
  }
  names
}

make.names <- function(names, unique=F) {
  names <- .Internal(make.names(as.character(names)))
  # minka: change keywords
  i <- is.element(names, c("for","while","repeat","if","else","function"))
  if(any(i)) names[i] <- paste(names[i],".",sep="")
  if(unique) names <- make.unique(names)
  names
}
replaceInNamespace("make.names",make.names)

# does NOT throw an error if no terms
terms.default <- function(x,...) x$terms
terms.data.frame <- function(x,env=parent.frame(),...) {
  fmla <- attr(x,"terms")
  if(is.null(fmla)) {
    # minka: assume the last column is the response
    #nm <- make.names(names(x))
    nm = names(x)
    if(length(nm) > 1) {
      lhs <- nm[length(nm)]
      rhs <- nm[-length(nm)]
    }
    else {
      lhs <- NULL
      rhs <- nm
    }
    fmla <- terms(formula(paste(lhs,"~",paste(rhs,collapse="+")),env=env,...))
  }
  fmla
}
formula.data.frame <- function(x,env=parent.frame(),...) {
  formula(terms(x,env=env),...)
}
terms.table <- function(x,...) {
  terms.data.frame(as.data.frame(x),...)
}

formula.default <- function (x,env=parent.frame(), ...)
{
    if (!is.null(x$formula))		eval(x$formula)
    else if (!is.null(x$call$formula))	eval(x$call$formula)
    # minka: always return formula, not terms
    else if (!is.null(x$terms))		formula.terms(x$terms)
    else if (!is.null(attr(x, "formula"))) attr(x, "formula")
    else {form<-switch(mode(x),
		NULL = structure(NULL, class = "formula"),
		character = formula(eval(parse(text = x)[[1]])),
		call = eval(x), stop("invalid formula"))
        environment(form)<-env
        form
    }
}

# Ripley says the original update is not broken (just non-modular)
# related problem:
# http://maths.newcastle.edu.au/~rking/R/help/02b/1318.html
my.update <-
    function (object, formula., ..., evaluate = TRUE)
{
    call <- object$call
    if (is.null(call))
	stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
	call$formula <- update.formula(formula(object), formula.)
    if(length(extras) > 0) {
	existing <- !is.na(match(names(extras), names(call)))
	## do these individually to allow NULL to remove entries.
	for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
	if(any(!existing)) {
	    call <- c(as.list(call), extras[!existing])
	    call <- as.call(call)
	}
    }
    if(evaluate) {
      # minka: use environment of formula instead of parent.frame
      # see the man page for formula
      env<-environment(call$formula)
      if (is.null(env)) env<-parent.frame()
      eval(call,env)
    }
    else call
}

plot.default <- function(x, y=NULL, type="p", xlim=NULL, ylim=NULL,
			 log="", main=NULL, sub=NULL, xlab=NULL, ylab=NULL,
			 ann=par("ann"), axes=TRUE, frame.plot=axes,
			 panel.first=NULL, panel.last=NULL,
			 col=par("col"), bg=NA, pch=par("pch"),
			 cex = 1, lty=par("lty"), lab=par("lab"),
			 lwd=par("lwd"), asp=NA, ...)
{
    xlabel <- if (!missing(x)) deparse(substitute(x))
    ylabel <- if (!missing(y)) deparse(substitute(y))
    # minka
    # treat character as factor
    if(mode(y) == "character") y <- factor(y)
    if(mode(x) == "character") x <- factor(x)
    # this works even if x and y are factors
    xy <- xy.coords(x, y, xlabel, ylabel, log)
    if(all(is.na(xy$x))) xy$x <- 1:length(x)
    xlab <- if (is.null(xlab)) xy$xlab else xlab
    ylab <- if (is.null(ylab)) xy$ylab else ylab
    xlim <- if (is.null(xlim)) range(xy$x[is.finite(xy$x)]) else xlim
    ylim <- if (is.null(ylim)) range(xy$y[is.finite(xy$y)]) else ylim
    plot.new()
    plot.window(xlim, ylim, log, asp, ...)
    panel.first
    plot.xy(xy, type, col=col, pch=pch, cex=cex, bg=bg, lty=lty, lwd=lwd, ...)
    panel.last
    if (axes) {
        # minka
        if(is.factor(x)) {
          lev <- levels(x)
          axis(1, at=1:length(lev), labels=lev, ...)
        } else axis(1, ...)
        # minka
        if(is.factor(y)) {
          lev <- levels(y)
          axis(2, at=1:length(lev), labels=lev, ...)
        } else axis(2, ...)
    }
    if (frame.plot)
	box(...)
    if (ann)
	title(main=main, sub=sub, xlab=xlab, ylab=ylab, ...)
    invisible()
}

model.frame.multinom <- function(object) {
  oc <- object$call
  oc[[1]] <- NULL
  do.call("model.frame",as.list(oc))
}

model.frame.glm <-
function (object, data, na.action, ...)
{
    if (is.null(object$model)) {
        fcall <- object$call
        fcall$method <- "model.frame"
        #fcall[[1]] <- as.name("glm")
        # minka
        if (!is.null(object$terms))
            env <- environment(object$terms)
        else
            env <- environment(fcall$formula)
        if (is.null(env))
            env <- parent.frame()
        eval(fcall, env)
    }
    else object$model
}
model.frame.nnet <- model.frame.glm
# missing from modreg
model.frame.loess <- model.frame.glm

#############################################################################
# bugfix

my.loess <-
function(formula, data=NULL, weights, subset, na.action, model = FALSE,
	 span = 0.75, enp.target, degree = 2, parametric = FALSE,
	 drop.square = FALSE, normalize = TRUE,
	 family = c("gaussian", "symmetric"),
	 method = c("loess", "model.frame"),
	 control = loess.control(...), ...)
{
    family <- match.arg(family)
    method <- match.arg(method)
    mt <- terms(formula, data = data)
    mf <- match.call(expand.dots=FALSE)
    mf$model <- mf$span <- mf$enp.target <- mf$degree <-
      mf$parametric <- mf$drop.square <- mf$normalize <- mf$family <-
        mf$control <- mf$... <- NULL
    # minka: bugfix
    mf$method <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (match.arg(method) == "model.frame") return(mf)
    na.act <- attr(mf, "na.action")
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    if(is.null(w)) w <- rep(1, length(y))
    nmx <- as.character(attr(mt, "variables"))[-(1:2)]
    x <- mf[, nmx, drop=FALSE]
    if(any(sapply(x, is.factor))) stop("predictors must all be numeric")
    x <- as.matrix(x)
    D <- ncol(x)
    nmx <- colnames(x)
    names(nmx) <- nmx
    drop.square <- match(nmx, nmx[drop.square], 0) > 0
    parametric <- match(nmx, nmx[parametric], 0) > 0
    if(!match(degree, 0:2, 0)) stop("degree must be 0, 1 or 2")
    iterations <- if(family=="gaussian") 1 else control$iterations
    if(!missing(enp.target))
	if(!missing(span))
	    warning("both span and enp.target specified: span will be used")
	else {				# White book p.321
	    tau <- switch(degree+1, 1, D+1, (D+1)*(D+2)/2) - sum(drop.square)
	    span <- 1.2 * tau/enp.target
	}
    fit <- simpleLoess(y, x, w, span, degree, parametric, drop.square,
		       normalize, control$statistics, control$surface,
		       control$cell, iterations, control$trace.hat)
    fit$call <- match.call()
    fit$terms <- mt
    fit$xnames <- nmx
    fit$x <- x
    fit$y <- y
    fit$weights <- w
    if(model) fit$model <- mf
    if(!is.null(na.act)) fit$na.action <- na.act
    fit
}
# minka: robust by default
formals(my.loess)$family = "symmetric"

my.predict.loess <- function(object, newdata = NULL, se = FALSE, ...)
{
    if(!inherits(object, "loess"))
	stop("First argument must be a loess object")
    if(is.null(newdata) & (se == FALSE)) return(fitted(object))

    if(is.null(newdata)) newx <- object$x
    else {
      vars <- as.character(attr(delete.response(terms(object)),
                                "variables"))[-1]
      if(length(vars) > 1 || NCOL(newdata) > 1) {
        if(any(!match(vars, colnames(newdata), FALSE)))
          # minka
          newdata = model.frame.default(delete.response(terms(object)),newdata)
        else
          newdata = newdata[,vars,drop=F]
      }
      newx = as.matrix(newdata)
    }
    res <- predLoess(object$y, object$x, newx, object$s, object$weights,
		     object$pars$robust, object$pars$span, object$pars$degree,
		     object$pars$normalize, object$pars$parametric,
		     object$pars$drop.square, object$pars$surface,
		     object$pars$cell, object$pars$family,
		     object$kd, object$divisor, se=se)
    if(se)
	res$df <- object$one.delta^2/object$two.delta
    res
}

library(modreg)
replaceInNamespace("loess",my.loess)
environment(my.predict.loess) = environment(loess)
setInNamespace("predict.loess",my.predict.loess,environment(loess))

##############################################################################
# bug fix

# minka: added `key' and `main' arguments
my.filled.contour <-
function (x = seq(0, 1, len = nrow(z)),
          y = seq(0, 1, len = ncol(z)),
          z,
          xlim = range(x, finite=TRUE),
          ylim = range(y, finite=TRUE),
          zlim = range(z, finite=TRUE),
          levels = pretty(zlim, nlevels), nlevels = 20,
          color.palette = cm.colors,
          col = color.palette(length(levels) - 1),
          plot.title, plot.axes, key.title, key.axes, key=T,
          asp = NA, xaxs="i", yaxs="i", las = 1, axes = TRUE, main="", ...)
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq(0, 1, len = nrow(z))
            }
        }
        else stop("no `z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
        stop("increasing x and y values expected")

    if(key) {
      mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
      on.exit(par(par.orig))
      par(las = las)

      w <- (3 + mar.orig[2]) * par('csi') * 2.54
      layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
      ## Plot the `plot key' (scale):
      mar <- mar.orig
      mar[4] <- mar[2]
      mar[2] <- 1
      par(mar = mar)
      plot.new()
      plot.window(xlim=c(0,1), ylim=range(levels), xaxs="i", yaxs="i")
      rect(0, levels[-length(levels)], 1, levels[-1], col = col)
      if (missing(key.axes)) {
        if (axes)
          axis(4)
      }
      else key.axes
      box()
      if (!missing(key.title))
	key.title
    }
    ## Plot contour-image::
    if(key) {
      mar <- mar.orig
      mar[4] <- 1
      par(mar=mar)
    }
    plot.new()
    plot.window(xlim, ylim, "", xaxs=xaxs, yaxs=yaxs, asp=asp)

    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
        stop("no proper `z' matrix specified")
    if (!is.double(z))
        storage.mode(z) <- "double"
    .Internal(filledcontour(as.double(x),
                            as.double(y),
                            z,
                            as.double(levels),
                            col = col))
    if (axes) {
      if (missing(plot.axes)) {
        # minka: no title here
        axis(1)
        axis(2)
      } else plot.axes
      if (missing(plot.title))
        title(main=main,...)
      else
	plot.title
    }
    box()
    invisible()
}
replaceInNamespace("filled.contour",my.filled.contour)

##############################################################################

# Routines for contingency tables
# Tom Minka

barplot.table <- function(m,col=default.colors.w(nrow(m)),
                          xlab,main,
                          ylab=response.var(m),legend.text=T,...) {
  nam <- names(dimnames(m))
  opar <- par(mar=auto.mar(main=nam[1]))
  on.exit(par(opar))
  if(missing(xlab)) xlab = nam[2]
  if(missing(main)) main = nam[1]
  barplot.default(m,legend.text=legend.text,col=col,...,
                  xlab=xlab,ylab=ylab,main=main)
}

flipud <- function(m) {
  m[nrow(m):1,]
}
fliplr <- function(m) {
  m[,ncol(m):1]
}

image.table <- function(m,col=YR.colors(64),bg=gray(0.5),mar,mar.scale=NULL,
                        ...) {
  m <- flipud(m)
  b = break.quantile(as.vector(m),length(col))
  #b = break.equal(as.vector(m),length(col))
  if(is.null(colnames(m))) colnames(m) = 1:ncol(m)
  if(is.null(rownames(m))) rownames(m) = 1:nrow(m)
  nam <- names(dimnames(m))
  if(missing(mar)) mar = auto.mar(xlab=nam[2],ylab=nam[1],mar.scale=mar.scale,...)
  opar = par(bg=bg,mar=mar)
  on.exit(par(opar))
  image.default(t(m),axes=F,col=col,breaks=b,...)
  box()
  # axis 1 is horizontal
  axis(1,at=seq(0,1,len=ncol(m)),labels=colnames(m),...)
  # axis 2 is vertical
  axis(2,at=seq(0,1,len=nrow(m)),labels=rownames(m),...)
  title(xlab=nam[2],ylab=nam[1])
}

if(!exists('pie.default')) {
  pie.default = pie
  pie = function(x,...) UseMethod("pie")
}

# draws a pie for each column
pie.table <- function(m,layout,col=default.colors.w(nrow(m)),...) {
  if(length(dim(m)) == 1) return(pie.default(m,col=col,...))
  if(missing(layout)) layout <- auto.layout(ncol(m))
  opar <- par(mfrow=layout,mar=c(1,1,1,0)*2)
  on.exit(par(opar))
  nam <- names(dimnames(m))
  for(i in 1:ncol(m)) {
    pie(m[,i],col=col,...)
    title(paste(nam[1],"\n",nam[2],"=",colnames(m)[i]))
  }
}

##############################################################################

# only for 2D tables right now
sort.table <- function(x,k=NULL) {
  if(!is.null(k)) {
    if(k == 1) {
      # sort rows only
      n <- margin.table(x,k)+1e-15
      v2 <- 1:ncol(x)
      v1 <- drop(x %*% v2)/n
      return(x[order(v1),])
    }
    if(k == 2) {
      # sort cols only
      n <- margin.table(x,k)+1e-15
      v1 <- 1:nrow(x)
      v2 <- drop(v1 %*% x)/n
      return(x[,order(v2)])
    }
  }

  n1 <- margin.table(x,1)
  n2 <- margin.table(x,2)
  n <- sum(n1)
  if(n == 0) return(x)

  # iterate
  v2 <- rnorm(ncol(x))
  iter <- 1
  repeat {
    old.v2 <- v2
    # update rows
    # add eps in case there are zeros
    v1 <- drop(x %*% v2)/(n1+eps)
    # shift to zero mean
    v1 <- v1 - sum(v1*n1)/n
    # scale to unit sd
    v1 <- v1/sqrt(sum(v1^2*n1)/n)
    # update cols
    v2 <- drop(v1 %*% x)/(n2+eps)
    if(max(abs(v2 - old.v2)) < 1e-8) break
    iter <- iter + 1
  }
  # remove sign ambiguity
  nz <- which(n1 > 0)[1]
  s <- sign(v1[nz])
  v1 <- v1*s
  v2 <- v2*s

  # empty rows/cols get placed at end
  z1 <- which(n1 == 0)
  z2 <- which(n2 == 0)
  v1[z1] <- 1e10
  v2[z2] <- 1e10
  i1 <- order(v1)
  i2 <- order(v2)
  y <- x[i1,i2]
  class(y) <- class(x)
  y
  #list(v1=v1,v2=v2)
}

indep.fit <- function(x) {
  ## Compute the expected cell counts under the independence model.
  loglin(x, as.list(1:length(dim(x))), fit = TRUE, print = FALSE)$fit
}

# returns top k associations in table x
# uses sort.cells, sort.data.frame, rbind.extend
mine.associations <- function(x,top=10,targets=NULL,z=1) {
  if(!inherits(x,"class")) x <- as.table(x)
  dm <- dim(x)
  nd <- length(dm)
  dn <- attr(dimnames(x),"names")
  if(is.null(targets)) targets <- 1:nd
  else if(is.character(targets)) {
    targets <- pmatch(targets,dn)
    if(any(is.na(targets))) stop("unrecognized target")
  }
  not.targets <- setdiff(1:nd, targets)

  # create a data frame with 0 rows
  res <- empty.data.frame(c("Lift",dn))
  # loop rows
  for(i in targets) {
    predictors <- c(not.targets, targets[targets > i])
    for(j in predictors) {
      y <- margin.table(x, c(i,j))
      df <- sort.cells((y - z*sqrt(y))/(indep.fit(y)+eps))
      names(df)[3] <- "Lift"
      # take top k lifts
      if(nrow(df) > top) df <- df[nrow(df)+1-(top:1),]
      # extend cols
      for(v in 1:nd) {
        if(v != i && v != j) df[[dn[v]]] <- factor(rep(NA,nrow(df)))
      }
      res <- rbind.extend(res,df)
      # ensure unique rownames
      #rownames(res) <- nrow(res):1
    }
  }
  res <- sort.data.frame(res,1)
  if(nrow(res) > top) res <- res[nrow(res)+1-(top:1),]
  rownames(res) <- nrow(res):1
  res
}

##############################################################################
# merging

# Returns Pearson's chi-square statistic for testing whether two samples come
# from the same population.
# n1 and n2 are count vectors (same length)
homogeneity <- function(n1,n2) {
  # most time is spent in this function so it must be efficient
  if(sum(n1)==0 || sum(n2)==0) return(0)
  n12 <- n1+n2
  e1 <- n12*sum(n1)/sum(n12)
  e2 <- n12*sum(n2)/sum(n12)
  # chisq approx
  sum((n1 - e1)^2/(e1+eps)) + sum((n2 - e2)^2/(e2+eps))
}

js.divergence <- function(n1,n2) {
  if(sum(n1)==0 || sum(n2)==0) return(0)
  n12 <- n1+n2

  p1 <- n1/(sum(n1)+eps)
  p2 <- n2/(sum(n2)+eps)
  p12 <- n12/(sum(n12)+eps)

  i <- (p1 > 0)
  g <- sum(p1[i]*log(p1[i]/p12[i]))*sum(n1)
  i <- (p2 > 0)
  g <- g + sum(p2[i]*log(p2[i]/p12[i]))*sum(n2)
  g
}

# returns a matrix of merging costs for the values along dimension v
merge.table.cost <- function(x,v,cost=homogeneity) {
  ordered <- dimOrdered(x)[[v]]
  # make v the first dimension
  p <- 1:length(dim(x))
  p[1] <- v
  p[v] <- 1
  x <- aperm(x,p)
  d <- dim(x)
  # collapse all other dimensions
  dim(x) <- c(d[1], prod(d[2:length(d)]))

  if(ordered) {
    # compute cost between all adjacent values
    s <- 0
    for(i in 1:(d[1]-1)) {
      s[i] <- cost(x[i,],x[i+1,])
    }
  } else {
    # compute cost between all value pairs
    s <- array(0,c(d[1],d[1]))
    for(i in 1:nrow(s)) {
      for(j in 1:ncol(s)) {
        if(i == j) s[i,j] <- Inf
        else s[i,j] <- cost(x[i,],x[j,])
      }
    }
  }
  s
}

test.merge.table.cost <- function(cost=homogeneity) {
  x <- array(c(0,0,1,0,1,1),c(3,2))
  i <- rep(1:3,rep(2,3))
  x <- x[i,]
  merge.table.cost(x,1,cost)
}

# x is a table
# returns a new table with values i and j along dimension v merged
merge.table.cells <- function(x,v,i,j) {
  p <- 1:length(dim(x))
  p[1] <- v
  p[v] <- 1
  x <- aperm(x,p)
  dn <- dimnames(x)
  d <- dim(x)
  dim(x) <- c(d[1], prod(d[2:length(d)]))
  x[i,] <- x[i,]+x[j,]
  x <- x[-j,]
  d[1] <- d[1] - 1
  dim(x) <- d
  dn[[1]] <- merge.names(dn[[1]],c(i,j))
  dimnames(x) <- dn
  # aperm removes the class attribute
  x <- aperm(x,p)
  class(x) <- "table"
  x
}

# cost is the merging cost function to use (default homogeneity)
# trace is the max number of costs to report in the merging trace
#
# if trees are desired, returns list(table=x,trees=t)
# x is the merged table
# t is a list of length(ds) trees, or one tree if length(ds)==1
merge.table <- function(x, bins=rep(2,length(ds)),
                        ds=1:length(dim(x)),
                        cost=homogeneity, trace=12) {
  if(length(bins) != length(ds)) {
    stop("must specify the desired number of bins for each dimension being merged")
  }
  dn <- attr(dimnames(x),"names")
  if(is.character(ds)) {
    # convert ds to integers
    if(is.null(dn)) stop("table doesn't have category names")
    ds <- pmatch(ds, dn)
    if(any(is.na(ds))) stop("unrecognized variable")
  }
  d <- dim(x)
  if(F) {
    # initialize trees
    hs <- list()
    for(k in 1:length(ds)) {
      dk <- d[ds[k]]
      hs[[k]] <- list(merge=matrix(0,dk-1,2),height=numeric(0),
                      order=1:dk,labels=-(1:dk),var=dn[k])
    }
  }

  # which.dim[i] is the dimension which was merged at step i
  which.dim <- numeric(0)
  costs <- numeric(0)
  nbins <- character(0)
  # main loop
  dim.cost <- array(0,c(length(ds),1))
  dim.merge <- list()
  iter <- 1
  # has desired table been found?
  desired <- F
  repeat {
    d <- dim(x)
    # score all dimensions
    for(k in 1:length(ds)) {
      dk <- d[ds[k]]
      # can we make this dimension any smaller?
      if(dk > bins[k]) {
        s <- merge.table.cost(x,ds[k],cost)
        if(is.null(dim(s))) {
          # ordered case
          # s is a vector
          dim.cost[k] <- min(s)
          i <- which.min(s)
          dim.merge[[k]] <- c(i,i+1)
        } else {
          # unordered case
          # s is a matrix
          dim.cost[k] <- min(s)
          dim.merge[[k]] <- which.min2(s)
        }
        #cat("dim",ds[k],":",dim.cost[k],"\n")
      } else {
        dim.cost[k] <- Inf
      }
    }
    k <- which.min(dim.cost)
    if(length(k) == 0) {
      if(!desired) {
        # desired table size has been reached.
        desired <- T
        cat("total cost =",sum(costs),"\n")
        if(trace) {
          # if a full trace is wanted, keep going until 1 bin.
          desired.x <- x
          bins <- rep(1,length(ds))
        } else break
      } else {
        x <- desired.x
        break
      }
    } else {
      # merge
      i <- dim.merge[[k]][1]
      j <- dim.merge[[k]][2]
      #print(paste("merging value", i, "of dimension",ds[k]))
      if(!desired) {
        cat("merging",dn[ds[k]],"=",dimnames(x)[[ds[k]]][i],
            "and",dimnames(x)[[ds[k]]][j],"\n")
      }
      x <- merge.table.cells(x,ds[k],i,j)

      # record
      which.dim[iter] <- k
      costs[iter] <- dim.cost[k]
      nbins[iter] <- paste(dim(x)[ds],collapse="x")
      if(F) {
        len <- length(hs[[k]]$height)
        hs[[k]]$height[len+1] <- dim.cost[k]
        hs[[k]]$merge[len+1,1] <- hs[[k]]$labels[i]
        hs[[k]]$merge[len+1,2] <- hs[[k]]$labels[j]
        hs[[k]]$labels[i] <- len+1
        hs[[k]]$labels <- hs[[k]]$labels[-j]
      }
      iter <- iter + 1
    }
  }
  #cat("merge order:", dn[ds[which.dim]], "\n")
  if(F) {
    for(k in 1:length(hc)) {
      hs[[k]]$height <- cumsum(hs[[k]]$height)
      hs[[k]]$labels <- NULL
    }
  }
  if(trace) {
    # show trace of merging cost
    if(length(costs) > trace) {
      costs <- rev(costs)
      costs <- costs[1:trace]
      costs <- rev(costs)
      nbins <- rev(nbins)
      nbins <- nbins[1:trace]
      nbins <- rev(nbins)
    }
    par(mai=c(1,1,0.5,0.25))
    dotchart(costs,labels=nbins,xlab="merging cost",ylab="nbins")
  }
  x
}
# Geometry functions
# rest are in maps

# result is d[i,j] = squared Euclidean distance betw x[i,] and y[j,]
sqdist <- function(x,y=x,A) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  if(missing(A)) {
    xmag <- rowSums(x * x)
    ymag <- rowSums(y * y)
    d = rep.row(ymag, nrow(x)) + rep.col(xmag, nrow(y)) - 2*(x %*% t(y))
  } else {
    xA = x %*% A
    yA = y %*% A
    xmag <- rowSums(x * xA)
    ymag <- rowSums(y * yA)
    d = rep.row(ymag, nrow(x)) + rep.col(xmag, nrow(y)) - 2*(x %*% t(yA))
  }
  if(identical(x,y)) diag(d) = 0
  # fix numeric errors
  d[which(d < 0)] = 0
  d
}

nearest.points <- function(a,b) {
  # for each point (row) in matrix A, returns the index of the closest point
  # in matrix B, as well as the distance to that closest point.
  # if B is omitted, it is assumed to be the same as A, but where a
  # point cannot match itself.
  b.which = numeric(nrow(a))
  b.dist = numeric(nrow(a))

  if(missing(b)) block.a = ceiling(100000/nrow(a))
  else block.a = ceiling(100000/nrow(b))
  i = 1
  while(i <= nrow(a)) {
    ind = i:min(i+block.a,nrow(a))
    if(missing(b)) {
      d = sqdist(a[ind,],a)
      d[(ind-1)*length(ind) + ind] = Inf
    } else {
      d = sqdist(a[ind,],b)
    }
    b.which[ind] = apply(d,1,which.min)
    b.dist[ind] = d[ind + (b.which[ind]-1)*nrow(d)]
    i = i + block.a
  }
  list(which=b.which,dist=b.dist)
}

distances <- function(x,fun) {
  if(missing(fun)) sqrt(sqdist(x))
  else {
    n <- nrow(x)
    d <- array(NA,c(n,n),list(rownames(x),rownames(x)))
    for(i in 1:n) {
      for(j in 1:n) {
        d[i,j] <- fun(x[i,],x[j,])
      }
    }
    d
  }
}

as.matrix.polygon <- function(p) {
  if(is.list(p) && !is.data.frame(p)) p <- cbind(p$x,p$y)
  p
}
as.list.polygon <- function(p) {
  if(is.matrix(p)) p <- data.frame(p,c("x","y"))
  p
}

# samples n points uniformly in p using rejection from bounding box
# result is a matrix
polygon.sample <- function(p,n) {
  r <- polygon.range(p)
  x <- NULL
  total.ok <- 0.9
  total.tried <- 1
  while(n > 0) {
    # m is the number of points to sample in order to have n in the polygon
    pct <- total.ok/total.tried
    m <- n + qnbinom(0.9, n, pct)
    xt <- cbind(runif(m,r$x[1],r$x[2]), runif(m,r$y[1],r$y[2]))

    # add the points which fell in the polygon to our list
    ok <- which(in.polygon(p,xt))
    n.ok <- length(ok)
    if(length(ok) > n) ok <- ok[1:n]
    x <- rbind(x,xt[ok,])
    #cat("got",n.ok,"points using pct=",pct,"\n")

    # update percentage and try again
    n <- n - n.ok
    total.ok <- total.ok + n.ok
    total.tried <- total.tried + nrow(xt)
  }
  x
}

# returns T if x is inside the polygon with given vertices
# poly is a matrix where each row is a vertex
# polygon is closed automatically
# or a list of x and y (converted to a matrix)
in.polygon <- function(poly,x) {
  poly <- as.matrix.polygon(poly)
  if(is.list(x) && !is.data.frame(x)) x <- cbind(x$x,x$y)
  if(is.vector(x)) in.polygon1(poly,x)
  num.above <- rep(0,nrow(x))
  start <- 1
  end <- nrow(poly)
  for(i in 1:nrow(poly)) {
    x1 <- poly[i,1]
    if(is.na(x1)) start <- i+1
    else {
      if(i == nrow(poly) || is.na(poly[i+1,1])) { i2 <- start }
      else { i2 <- i+1 }
      x2 <- poly[i2,1]
      y1 <- poly[i,2]
      y2 <- poly[i2,2]
      xmax <- max(x1,x2)
      xmin <- min(x1,x2)
      # xi indexes the points within the line extent
      xi <- which((x[,1]>=xmin) & (x[,1]<=xmax))
      # a is the height of the line at x[xi,1]
      a <- (x[xi,1]-x1)*(y2-y1)/(x2-x1) + y1
      vert <- (x1==x2)
      a[vert] <- y1[vert]
      # j indexes the points below the line
      j <- (a>x[xi,2])
      num.above[xi[j]] <- num.above[xi[j]] + 1
    }
  }
  (num.above %% 2) == 1
}
in.polygon1 <- function(poly,x) {
  x1 <- poly[,1]
  x2 <- c(poly[2:nrow(poly),1],poly[1,1])
  y1 <- poly[,2]
  y2 <- c(poly[2:nrow(poly),2],poly[1,2])
  breaks <- which(is.na(x1))
  if(length(breaks) > 0) {
    starts <- c(1,breaks+1)
    ends <- c(breaks-1,nrow(poly))
    x2[ends] <- x1[starts]; y2[ends] <- y1[starts]
    x1 <- x1[-breaks]; x2 <- x2[-breaks]; y1 <- y1[-breaks]; y2 <- y2[-breaks]
  }
  vert <- which(x1==x2)
  # a is the height of the lines at x[1]
  a <- (x[1]-x1)*(y2-y1)/(x2-x1) + y1
  a[vert] <- y1[vert]
  xmin <- pmin(x1,x2)
  xmax <- pmax(x1,x2)
  # i indexes the lines above x
  i <- (a>x[2]) & (x[1] >= xmin) & (x[1] <= xmax)
  (sum(i) %% 2) == 1
}

test.in.polygon <- function() {
  n <- 1000
  x <- cbind(runif(n),runif(n))
  plot(x)
  p <- cbind(runif(7),runif(7))
  p <- p[chull(p),]
  # close the polygon
  p <- rbind(p,p[1,])
  if(T) {
    p2 <- cbind(runif(7),runif(7))
    p2 <- p2[chull(p2),]
    # close the polygon
    p2 <- rbind(p2,p2[1,])
    p = rbind(p,cbind(NA,NA),p2)
  }
  points(p,col=2)
  polygon(p,border=2)
  if(F) {
    for(i in 1:nrow(x)) {
      if(in.polygon1(p,x[i,])) points(x[i,,drop=F],col=3)
    }
  } else {
    points(x[in.polygon(p,x),],col=3)
  }
}

# Routines for manipulating directed graphs
# Tom Minka 10/16/01

# returns a graph with n nodes and no edges
# n is an integer or a vector of node names
graph <- function(n) {
  if(length(n) == 1) n <- 1:n
  # a graph is a list of vectors
  g <- sapply(n, function(x) numeric(0))
  names(g) <- as.character(n)
  g
}

# adds an edge from i to j
# returns new graph
add.edge <- function(g,i,j) {
  # i and j can be names or numbers
  if(is.character(j)) j <- pmatch(j,names(g))
  # node vector contains node indices, not names
  # this allows node names to be changed by the user via names(g)
  g[[i]] <- union(g[[i]],j)
  g
}

has.edge <- function(g,i,j) {
  if(is.character(j)) j <- pmatch(j,names(g))
  j %in% g[[i]]
}

from <- function(g,i) {
  # i is a name or index
  # result is indices
  g[[i]]
}
edges.to <- function(g,j) {
  if(is.character(j)) j <- pmatch(j,names(g))
  which(sapply(g,function(x) (j %in% x)))
}

is.leaf <- function(g,i) {
  length(g[[i]]) == 0
}

leaves <- function(g,i=1) {
  # returns a vector of leaf indices
  child <- from(g,i)
  if(length(child) == 0) {
    i
  } else {
    r <- c()
    for(j in child) {
      r <- c(r,leaves(g,j))
    }
    r
  }
}

# returns a directed graph with edges from parents to children
as.graph.tree <- function(tr) {
  frame <- tr$frame
  g <- graph(rownames(frame))
  p <- parents.tree(tr)
  for(i in 1:nrow(frame)) {
    if(!is.na(p[i])) g <- add.edge(g,p[i],i)
  }
  g
}
plot.graph.tree <- function(tr,digits=2,tol=0.1,abbrev=F,sizes=T,extra.labels=NULL,...) {
  g <- as.graph.tree(tr)
  frame <- tr$frame
  labels <- as.character(frame$var)
  names(labels) = rownames(frame)
  i <- (labels == "<leaf>")
  #labels[i] <- paste(response.var(tr),"=",format(auto.signif(frame$yval[i],tol=tol)))
  y = frame$yval
  if(is.numeric(y)) y = signif(y,digits)
  labels[i] <- paste(response.var(tr),"=",as.character(y[i]))
  if(!is.null(frame$yprob)) {
    f = indicators.factor(frame$yval[i])
    p = rowSums(frame$yprob[i,] * t(f))
    labels[i] = paste(labels[i]," (",as.character(signif(100*p,2)),"%)",sep="")
  }
  if(sizes) {
    # size of each leaf
    labels[i] <- paste(paste("(",as.character(frame$n[i]),")",sep=""),
                       labels[i])
  }
  if(!is.null(extra.labels)) {
    i = names(extra.labels)
    labels[i] = paste(labels[i],extra.labels)
  }
  #labels[i] <- "?"
  #labels[1:length(labels)] <- "?"
  p <- parents.tree(tr)
  left <- is.left.child(tr)
  xlevels <- attr(tr,"xlevels")
  if(abbrev) xlevels <- lapply(xlevels, abbreviate, minlength=abbrev)
  # clean up the split labels
  for(i in 1:nrow(frame)) {
    splits <- frame$splits[i,]
    splits1 <- substring(splits,1,1)
    splits2 <- substring(splits,2)
    if(splits1[1] == "<") {
      # numeric split
      x <- format(auto.signif(as.numeric(splits2),tol=tol))
      frame$splits[i,] <- paste(splits1,x,sep="")
    } else if(splits[1] != "") {
      # categorical split
      lev <- xlevels[[labels[[i]]]]
      for(j in 1:2) {
        nc <- nchar(splits2[j])
        sh <- substring(splits2[j],1:nc,1:nc)
        frame$splits[i,j] <- paste(lev[match(sh,letters)],collapse=",")
      }
    }
  }
  branch.labels <- c()
  for(i in 1:length(g)) {
    if(!is.na(p[i])) {
      branch.labels[[i]] <-
        if(left[i]) as.character(frame$splits[p[i],1])
        else as.character(frame$splits[p[i],2])
    }
  }
  #branch.labels[1:length(branch.labels)] <- "?"
  plot.graph(g,labels=labels,branch.labels=branch.labels,...)
}

root.of.graph <- function(g) {
  # the root is the first node with no parents
  # result is an index, not a name
  setdiff(1:length(g),accumulate(g,union))[1]
}
sort.graph <- function(g) {
  depths <- node.depths(g,root.of.graph(g))
  ord <- order(depths)
  n <- length(g)
  ranks <- rank.stable(depths)
  for(i in 1:length(g)) {
    g[[i]] <- ranks[g[[i]]]
  }
  g[ord]
}

as.graph.hclust <- function(hc) {
  n.leaves <- length(hc$labels)
  n.internal <- nrow(hc$merge)
  g <- graph(n.internal + n.leaves)
  i.internal <- 1:n.internal
  i.leaves <- n.internal + (1:n.leaves)
  names(g)[i.leaves] <- hc$labels
  for(i in 1:nrow(hc$merge)) {
    for(k in 1:2) {
      j <- hc$merge[i,k]
      if(j < 0) j <- hc$labels[-j]
      else j <- as.character(j)
      g <- add.edge(g,i,j)
    }
  }
  sort.graph(g)
}
plot.graph.hclust <- function(hc,col=par("fg"),...) {
  g = as.graph.hclust(hc)
  labels = rep("",length(g))
  i.leaves = leaves(g)
  labels[i.leaves] = names(g)[i.leaves]
  if(length(col) > 1) {
    if(!is.null(names(col))) col = col[labels]
    else stop("col vector must have names")
  }
  plot.graph(g,labels,col=col,...)
}

node.depths <- function(g,i=1,depth=0,depths=c()) {
  child <- from(g,i)
  depths[[i]] <- depth
  for(j in child) {
    depths <- node.depths(g,j,depth+1,depths)
  }
  depths
}

bbox.graph <- function(g,i=1,heights=NULL,labels=NULL,branch.labels=NULL,
                       space=c(1,1),inch.usr,
                       cex=par("cex"),srt=90,orient=90,bbox=list(),...) {
  # returns c(w.usr,w.inch,h.usr,h.inch)
  # "usr" is flexible space, while "inch" is not
  # the meaning of "height" depends on the orientation of the plot
  # space is in user coordinates
  orient.h <- ((orient == 0) || (orient == 180))
  srt.h <- (srt == 0)
  if(missing(inch.usr)) {
    pin <- par("pin")
    usr <- par("usr")
    inch.usr <- c(pin[1]/diff(usr[1:2]),pin[2]/diff(usr[3:4]))
    if(orient.h) inch.usr <- inch.usr[2:1]
  }
  lab <- labels[[i]]
  if(!is.null(lab)) {
    if(orient.h == srt.h) {
      if(orient < 135) lab <- paste(lab,"")
      else lab <- paste("",lab)
    }
    if(srt.h == orient.h) {
      lab.width <- strheight(lab, units="inches", cex=cex)
      lab.height <- strwidth(lab, units="inches", cex=cex)
    } else {
      lab.width <- strwidth(lab, units="inches", cex=cex)
      lab.height <- strheight(lab, units="inches", cex=cex)
    }
  }
  child <- from(g,i)
  if(length(child) == 0) {
    # leaf
    w.usr <- 0
    w.inch <- 0
    h.usr <- 0
    h.inch <- 0
    if(!is.null(lab)) {
      w.inch <- max(w.inch,lab.width)
      h.inch <- max(h.inch,lab.height)
    }
  } else {
    # not leaf
    w.usr <- 0
    w.inch <- 0
    h.usr <- 0
    h.inch <- 0
    height <- 0
    if(!is.null(lab)) {
      h.inch <- lab.height
      height <- h.usr*inch.usr[2]+h.inch
    }
    for(j in child) {
      bbox <- bbox.graph(g,j,heights,labels,branch.labels,space,inch.usr,
                         cex,srt,orient,bbox)
      w.usr <- w.usr + bbox[[j]][1]
      w.inch <- w.inch + bbox[[j]][2]
      h <- c(bbox[[j]][3]+space[2],bbox[[j]][4])
      if(h[1]*inch.usr[2]+h[2] > height) {
        max.height <- h
        height <- h[1]*inch.usr[2]+h[2]
      }
    }
    w.usr <- w.usr + space[1]*(length(child)-1)
    h.usr <- max.height[1]
    h.inch <- max.height[2]
  }
  bbox[[i]] <- c(w.usr,w.inch,h.usr,h.inch)
  bbox
}

label.adj <- function(orient,srt) {
  if(orient == 270) {
    if(srt == 90) adj <- c(0,0.5)
    else adj <- c(0.5,-0.2)
  } else if(orient == 90) {
    if(srt == 90) adj <- c(1,0.5)
    else adj <- c(0.5,1)
  } else if(orient == 180) {
    if(srt == 90) adj <- c(0.5,1.2)
    else adj <- c(0,0.5)
  } else {
    if(srt == 90) adj <- c(0.5,-0.2)
    else adj <- c(1,0.5)
  }
  adj
}
draw.graph <- function(g,i=1,x=0,y=0,
                       heights=NULL,labels=NULL,branch.labels=NULL,
                       space=c(1,1),
                       usr=par("usr"),inch.usr,
                       cex=par("cex"),srt=90,orient=90,bbox,
                       col.labels=NULL,...) {
  orient.h <- ((orient == 0) || (orient == 180))
  srt.h <- (srt == 0)
  if(missing(inch.usr)) {
    pin <- par("pin")
    inch.usr <- c(pin[1]/diff(usr[1:2]),pin[2]/diff(usr[3:4]))
    if(orient.h) inch.usr <- inch.usr[2:1]
  }
  if(missing(bbox)) {
    print(1)
    bbox <- bbox.graph(g,i,heights,labels,branch.labels,
                       space,inch.usr,cex,srt,orient)
  }
  child <- from(g,i)
  lab <- labels[[i]]
  if(!is.null(lab)) {
    if(orient.h == srt.h) {
      if(orient < 135) lab <- paste(lab,"")
      else lab <- paste("",lab)
    }
    if(srt.h == orient.h) {
      lab.width <- strheight(lab, units="inches", cex=cex)
      lab.height <- strwidth(lab, units="inches", cex=cex)
    } else {
      lab.width <- strwidth(lab, units="inches", cex=cex)
      lab.height <- strheight(lab, units="inches", cex=cex)
    }
  }
  # label adjustment
  adj <- label.adj(orient,srt)
  # branch label adjustment
  branch.adj <- label.adj(if(orient.h) 270 else 0, srt)
  if(length(child) == 0) {
    # leaf
    if(!is.null(lab)) {
      col = par("col")
      if(!is.null(col.labels)) col = col.labels[i]
      if(orient.h) {
        text(y,x,lab,cex=cex,srt=srt,adj=adj,col=col,...)
      } else {
        text(x,y,lab,cex=cex,srt=srt,adj=adj,col=col,...)
      }
    }
  } else {
    # not leaf
    # xc is position of child
    w <- (bbox[[i]][1]*inch.usr[1] + bbox[[i]][2])/abs(inch.usr[1])
    xc <- x - w/2
    yc <- y - space[2]
    if(!is.null(lab)) {
      # draw node label
      if(length(child) == 2) {
        j <- child[1]
        w1 <- (bbox[[j]][1]*inch.usr[1] + bbox[[j]][2])/abs(inch.usr[1])
        x.lab <- xc + w1/2 + w/4 + space[1]/4
        # check for collision between node label and child
        if(lab.height > space[2]*abs(inch.usr[2])) {
          # place between the children
          x.lab <- xc + w1 + space[1]/2
        }
        col = par("col")
        if(!is.null(col.labels)) col = col.labels[i]
        if(orient.h) {
          text(y,x.lab,lab,cex=cex,srt=srt,adj=adj,col=col,...)
        } else {
          text(x.lab,y,lab,cex=cex,srt=srt,adj=adj,col=col,...)
        }
      } else stop("don't know how to label with >2 children")
    }
    for(j in child) {
      w <- (bbox[[j]][1]*inch.usr[1] + bbox[[j]][2])/abs(inch.usr[1])
      xc <- xc + w/2
      # draw line to child
      if(orient.h) {
        segments(y,x,y,xc,...)
        segments(y,xc,yc,xc,...)
      } else {
        segments(x,y,xc,y,...)
        segments(xc,y,xc,yc,...)
      }
      # draw branch label
      lab <- branch.labels[[j]]
      if(!is.null(lab)) {
        y.lab <- (y+yc)/2
        if(orient.h) {
          text(y.lab,xc,lab,cex=cex,srt=srt,adj=branch.adj,...)
        } else {
          text(xc,y.lab,lab,cex=cex,srt=srt,adj=branch.adj,...)
        }
      }
      #cat("node",i,"child",j,"at yc=",yc,"\n")
      draw.graph(g,j,xc,yc,heights,labels,branch.labels,space=space,
                 usr,inch.usr,cex,srt,orient,bbox,col.labels,...)
      xc <- xc + w/2 + space[1]
    }
  }
}

plot.graph <- function(g,labels=NULL,branch.labels=NULL,orient=180,srt=0,
                       main="",col.labels=NULL,mar=NULL,...) {
  if(!match(srt,c(0,90))) stop("srt must be 0 or 90")
  if(!match(orient,c(0,90,180,270))) stop("orient must be 0, 90, 180, or 270")
  if(is.null(mar)) mar = c(0.05,0,if(main=="") 0 else 1,0.05)
  opar <- par(mar=mar)
  on.exit(par(opar))
  plot.new()
  pin <- par("pin")
  if(orient == 0 || orient == 180) pin <- pin[2:1]
  bbox <- bbox.graph(g,labels=labels,branch.labels=branch.labels,
                     orient=orient,srt=srt,...)
  #print(bbox)
  w.usr <- bbox[[1]][1]
  w.inch <- bbox[[1]][2]
  h.usr <- bbox[[1]][3]
  h.inch <- bbox[[1]][4]
  w.frac <- w.inch/pin[1]
  w <- w.usr/(1 - w.frac)
  h.frac <- h.inch/pin[2]
  h <- h.usr/(1 - h.frac)
  if(orient == 270) {
    usr <- c(-w/2,w/2,0,-h)
  } else if(orient == 90) {
    usr <- c(-w/2,w/2,-h,0)
  } else if(orient == 180) {
    usr <- c(0,-h,-w/2,w/2)
  } else {
    usr <- c(-h,0,-w/2,w/2)
  }
  plot.window(usr[1:2],usr[3:4])
  title(main)
  space <- c(sign(w),sign(h))
  draw.graph(g,labels=labels,branch.labels=branch.labels,
             bbox=bbox,orient=orient,srt=srt,space=space,
             usr=usr,col.labels=col.labels,...)
  #rect(usr[1],usr[3],usr[2],usr[4])
}

##############################################################################
# knn

# use these to preprocess categorical predictors for knn
code.factor <- function(x,y) {
  if(!is.factor(x)) return(x)
  if(is.factor(y)) {
    tab = table(x,y)
    p = tab[,2]/(tab[,1]+tab[,2])
  } else {
    p = as.vector(as.matrix(prototypes(y,x)))
    #p = p - mean(p)
    names(p) = levels(x)
  }
  print(round(sort(p),2))
  p[x]
}
code.factors <- function(x,type=c("indicator","effect")) {
  # returns a new data frame where all categorical predictors are numerically
  # coded.
  resp = response.var(x)
  type = match.arg(type)
  if(type == "indicator") {
    # ordered factors are made integer, others made indicator
    x = apply.df(x,function(s) if(is.ordered(s)) as.numeric(s) else s)
    m = model.matrix(terms(x),x)[,-1,drop=F]
    cbind(as.data.frame(m),x[resp])
  } else {
    y = x[[resp]]
    pred = predictor.vars(x)
    for(k in pred) {
      if(is.factor(x[[k]])) cat(k,":\n")
      x[[k]] = code.factor(x[[k]],y)
    }
    x
  }
}

set.sublist <- function(lst,i,v) {
  if(i == 1) c(v,lst[-1])
  else if(i == length(lst)) c(lst[-i],v)
  else c(lst[1:(i-1)],v,lst[(i+1):length(lst)])
}
make.ordered.factors <- function(x,exclude=NULL) {
  # convert unordered factors in x with >2 levels into multiple
  # binary factors.
  is.unordered = function(f){is.factor(f) && !is.ordered(f) && nlevels(f)>2}
  vs = names(which(sapply(x,is.unordered)))
  vs = setdiff(vs,exclude)
  x = as.list(x)
  for(v in vs) {
    y = indicators.factor(x[[v]])[-1,,drop=F]
    rownames(y) = paste(v,rownames(y),sep="")
    y = as.list(as.data.frame(t(y)))
    y = lapply(y,factor.logical)
    i = which(names(x) == v)
    x = set.sublist(x,i,y)
  }
  as.data.frame(x)
}

make.numeric.data.frame <- function(data,exclude=NULL,warn=F) {
  i = !sapply(data,is.numeric)
  i[exclude] = FALSE
  if(any(i)) {
    if(warn) warning("converting factors to numeric")
    data[i] = apply.df(data[i],as.numeric)
  }
  data
}

knn.model <- function(formula,data=NULL,k=1) {
  # works for Splus too
  library(class)
  if(is.null(data)) data <- formula
  else data <- model.frame(formula,data)
  pred <- predictor.vars(data)
  data[pred] = make.numeric.data.frame(data[pred],warn=T)
  object <- list(terms=terms(data),data=data,k=k)
  class(object) <- "knn"
  object$scale <- sd(data[pred])
  object$scale[object$scale == 0] = 1
  resp <- response.var(data)
  p <- table(data[[resp]])
  object$tie.breaker <- levels(data[[resp]])[which.max(p)]
  object
}

print.knn <- function(object) {
  cat("nearest neighbor classifier, k =",format(object$k),"\n")
  pred <- predictor.vars(object)
  cat(nrow(object$data), "examples in", length(pred),"dimensions\n")
}

model.frame.knn <- function(object) {
  object$data
}



predict.knn <- function(object,test,k,type=c("class","vector","response"),...) {
  train <- object$data
  resp <- response.var(object)
  pred <- predictor.vars(object)
  test[pred] = make.numeric.data.frame(test[pred])
  s <- object$scale
  # subtracting mean is irrelevant
  x <- scale(train[pred],center=F,scale=s)
  xt <- scale(test[pred],center=F,scale=s)
  if(missing(k)) k <- object$k
  type <- match.arg(type)
  r <- knn(x,xt,train[[resp]],k=k,l=floor(k/2+1),prob=(type!="class"))
  r[is.na(r)] <- object$tie.breaker
  if(type == "vector") {
    p <- attr(r,"prob")
    lev <- levels(train[[resp]])
    v <- array(0,c(length(r),length(lev)),list(rownames(test),lev))
    for(j in 1:length(lev)) {
      i <- (r == lev[j])
      v[i,j] <- p[i]
      v[i,-j] <- 1-p[i]
    }
    r <- v
  } else if(type == "response") {
    p = attr(r,"prob")
    lev <- levels(train[[resp]])
    v = p
    i = which(r == lev[1])
    v[i] = 1-p[i]
    names(v) = rownames(test)
    r = v
  }
  r
}

misclass.knn <- function(object,data,rate=F) {
  if(missing(data)) data <- model.frame(object)
  resp <- response.var(object)
  r = sum(predict(object,data) != data[[resp]])
  if(rate) r/nrow(data) else r
}
deviance.knn <- function(object,data,rate=F) {
  if(missing(data)) data <- model.frame(object)
  resp <- response.var(object)
  y = data[[resp]]
  p = predict(object,data,type="vector")
  dev = 0
  for(lev in levels(y)) {
    dev = dev + sum(log(p[y==lev,lev]))
  }
  dev = -2*dev
  if(rate) dev/nrow(data) else dev
}

# should use deviance, not misclass
best.k.knn <- function(object,ks=1:50) {
  # find k with best cross-validation performance
  # uses leave-one-out only
  train <- object$data
  resp <- response.var(object)
  pred <- predictor.vars(object)
  s <- object$scale
  x <- scale(train[pred],center=F,scale=s)
  max.k = max(table(train[[resp]]))
  ks = ks[ks <= max.k]
  r <- c()
  for(i in 1:length(ks)) {
    y <- knn.cv(x,train[[resp]],k=ks[i],l=floor(ks[i]/2+1))
    y[is.na(y)] <- object$tie.breaker
    r[i] <- sum(y != train[[resp]])
  }
  plot(ks,r,xlab="number of neighbors (k)",ylab="misclass",type="o")
  title("Leave-one-out cross-validation")
  i <- min(which(r==min(r)))
  cat("best k is",ks[i],"\n")
  object$k <- ks[i]
  object
}
test.k.knn <- function(object,test,ks=1:50) {
  # find k with fewest test errors
  train <- object$data
  resp <- response.var(object)
  pred <- predictor.vars(object)
  test[pred] = make.numeric.data.frame(test[pred])
  s <- object$scale
  x <- scale(train[pred],center=F,scale=s)
  xt = scale(test[pred],center=F,scale=s)
  max.k = max(table(train[[resp]]))
  ks = ks[ks <= max.k]
  r = c()
  for(i in 1:length(ks)) {
    y = knn(x,xt,train[[resp]],k=ks[i],l=floor(ks[i]/2+1))
    y[is.na(y)] <- object$tie.breaker
    r[i] = sum(y != test[[resp]])
  }
  plot(ks,r,xlab="k",ylab="misclass",type="o")
  i <- which.min(r)
  cat("best k on test is",ks[i],"\n")
  object$k <- ks[i]
  object
  names(r) = ks
  invisible(r)
}

# score predictors via leave-one-out
# uses sort.data.frame
drop1.knn <- function(object) {
  train <- object$data
  resp <- response.var(object)
  pred <- predictor.vars(object)
  s <- object$scale
  x <- scale(train[pred],center=F,scale=s)
  x <- data.frame(x)
  yt <- train[[resp]]

  # score original model
  y <- knn.cv(x,yt,k=object$k,l=floor(object$k/2+1))
  y[is.na(y)] <- object$tie.breaker
  e <- sum(y != yt)
  res <- data.frame(k=object$k,misclass=e,row.names="<none>")
  # score each reduced model
  ks = 1:20
  for(i in 1:length(pred)) {
    e = c()
    for(j in 1:length(ks)) {
      k = ks[j]
      y <- knn.cv(not(x,pred[i]),yt,k=k,l=floor(k/2+1))
      y[is.na(y)] <- object$tie.breaker
      e[j] <- sum(y != yt)
    }
    k = ks[which.min(e)]
    e = min(e)
    res <- rbind(res, data.frame(k=k,misclass=e,row.names=paste("-",pred[i])))
  }
  sort.data.frame(res)
}

# remove predictors, using drop1.knn
step.knn <- function(object) {
  resp <- response.var(object)
  repeat {
    a <- drop1.knn(object)
    i <- which(a$misclass == min(a$misclass))
    i <- i[rownames(a)[i] != "<none>"]
    if(length(i) == 0) {
      print(a)
      break
    }
    i = i[1]
    pred <- predictor.vars(object)
    omit <- substring(rownames(a)[i], 3)
    cat("dropping",omit,"\n")
    pred <- pred[!(pred %in% omit)]
    object$data <- not(object$data, omit)
    object$scale <- not(object$scale, omit)
    object$k = a$k[i]
    object$terms <- terms(formula(paste(resp,"~",paste(pred,collapse="+"))))
  }
  object
}

reduce.knn <- function(object) {
  train <- object$data
  resp <- response.var(object)
  pred <- predictor.vars(object)
  s <- object$scale
  x <- scale(train[pred],center=F,scale=s)
  y <- train[[resp]]
  if(object$k == 1) {
    if(F) keep <- condense(x,y,trace=F)
    else keep <- 1:nrow(x)
    keep <- reduce.nn(x,keep,y)
  } else {
    keep <- multiedit(x,y,object$k)
  }
  object$data <- object$data[keep,]
  if(object$k > nrow(object$data)) {
    object$k <- nrow(object$data)-1
    cat("changing k to",object$k,"\n")
  }
  object
}

proto.knn <- function(object,np=1) {
  # reduce the examples to np cluster centers in each class
  library(mva)
  train <- object$data
  resp <- response.var(object)
  pred <- predictor.vars(object)
  s <- object$scale
  x <- scale(train[pred],center=F,scale=s)
  y <- train[[resp]]
  nc = nlevels(y)
  if(length(np) == 1) {
    np = rep(np,nc)
    if(F) {
      p = table(y)
      p = p/min(p)
      np = round(np*p)
      print(np)
    }
  }
  if(is.null(names(np))) names(np) = levels(y)
  r = list()
  for(lev in levels(y)) {
    x.lev = x[y==lev,]
    if(np[[lev]] == 1) {
      r[[lev]] = as.data.frame.row(mean.data.frame(x.lev),"1")
    } else if(np[[lev]] == nrow(x.lev)) {
      r[[lev]] = x.lev
    } else {
      cl = kmeans(x.lev,np[[lev]])
      r[[lev]] = as.data.frame(cl$centers)
    }
    rownames(r[[lev]]) = paste(lev,rownames(r[[lev]]),sep=".")
  }
  x = scale(rbind(r[[1]],r[[2]]),center=F,scale=1/s)
  y = rep(levels(y),times=sapply(r,nrow))
  object$data = cbind(y,x)
  names(object$data)[1] = resp
  object
}

gcl.model <- function(formula,data,equal.mean=F,
                      type=c("full","diag","sphere")) {
  # Gaussian classifier
  type = match.arg(type)
  if(is.null(data)) data <- formula
  else data <- model.frame(formula,data)
  pred = predictor.vars(data)
  data[pred] = make.numeric.data.frame(data[pred],warn=T)
  y = data[[response.var(data)]]
  m = v = p = bias = list()
  for(lev in levels(y)) {
    x = data[y==lev,pred,drop=F]
    m[[lev]] = if(equal.mean) mean(data[pred]) else mean(x)
    if(type == "full") {
      v[[lev]] = t(chol(ml.cov(x)))
    } else if(type == "sphere") {
      v[[lev]] = diag(ncol(x))*sqrt(mean(diag(cov(x))))
    } else {
      v[[lev]] = diag(sqrt(diag(cov(x))))
    }
    p[[lev]] = mean(y==lev)
    bias[[lev]] = -2*log(p[[lev]]) - sum(log(diag(v[[lev]])))
  }
  object = list(m=m,v=v,p=p,bias=bias,terms=terms(formula),call=match.call())
  class(object) = "gcl"
  object
}
predict.gcl <- function(object,test,type=c("class","response"),se=F) {
  pred <- predictor.vars(object)
  test[pred] = make.numeric.data.frame(test[pred])
  classes = names(object$m)
  e = array(0,c(nrow(test),length(classes)),list(rownames(test),classes))
  for(lev in classes) {
    xt = test[pred] - rep.row(object$m[[lev]],nrow(test[pred]))
    xt = forwardsolve(object$v[[lev]],t(xt))
    e[,lev] = t(colSums(xt*xt)) + object$bias[[lev]]
  }
  type = match.arg(type)
  if(type == "class") classes[apply(e,1,which.min)]
  else if(type == "response") 1/(1+exp(-0.5*(e[,1]-e[,2])))
}
misclass.gcl = misclass.knn
# Useful functions for making maps

locate.bbox <- function() {
  # use the mouse to pick a bounding box
  coord = locator(2,type="p",pch="+",col="red",cex=3)
  bbox = c(sort(coord$x),sort(coord$y))
  rect(bbox[1],bbox[3],bbox[2],bbox[4],col=NA,border="red")
  bbox
}

# convert census COB file into map object
# if centers=T, returns a data frame with columns x and y
# note: names may have duplicates (multi-polygon tracts)
read.cob <- function(pfile,pafile=NULL,centers=F) {
  s <- scan(pfile,"",quiet=T)
  # remove last two ENDs
  s <- s[-length(s)]
  s <- s[-length(s)]
  breaks <- which(s == "END")
  index <- s[c(1,breaks+1)]
  if(is.null(pafile)) {
    nam <- index
  } else {
    spa <- scan(pafile,"",quiet=T)
    nam <- spa[seq(2,length(spa),by=2)]
    names(nam) <- spa[seq(1,length(spa),by=2)]
    # order names to match the p file
    # tracts with bogus index will have name "NA"
    nam <- nam[index]
    if(centers) {
      # make names unique
      dup <- which(duplicated(nam))
      for(i in dup) {
        if(!is.na(nam[i])) {
          j <- which(nam == nam[i])
          nam[j] <- paste(nam[j],as.character(1:length(j)),sep=":")
        }
      }
    }
  }
  # remove internal names
  names(nam) = NULL
  if(centers) {
    # keep only first two coordinates of each block
    x <- as.numeric(s[c(2,breaks+2)])
    y <- as.numeric(s[c(3,breaks+3)])
    i <- !is.na(nam)
    data.frame(x=x[i],y=y[i],row.names=nam[i])
  } else {
    # change END to NA,NA
    s[c(breaks,breaks+1)] <- "NA"
    # remove first two coordinates of each block
    i <- rep(T,length(s))
    i[1:3] <- F
    i[c(breaks+2, breaks+3)] <- F
    s <- as.numeric(s[i])
    list(x=s[seq(1,length(s),by=2)],y=s[seq(2,length(s),by=2)],names=nam)
  }
}
# convert a FIPS code into a comma-separated region name
fips.split = function(x) {
  sapply(x, function(s) paste(substring(s,1,5),substring(s,6),sep=","))
}

# choropleth map
map.vector <- function(m,x,main="",nlevels=NULL,key=T,warn=T,
                   col=color.palette,color.palette=YlGnBu.colors,
                   bg=gray(0.5),breaks,mar,...) {
  if(is.data.frame(x)) {
    df <- x; x <- df[[1]]; names(x) <- rownames(df)
    if(missing(main)) main <- colnames(df)[1]
  }
  if(!warn && length(setdiff(names(x),m$names)) > 0) {
    warning("some regions have data but do not appear on the map")
  }
  #print(setdiff(m$names,names(x)))
  j <- match.map(m,names(x),warn=warn)
  x <- as.numeric(x)
  # determine breaks,nlevels
  if(missing(breaks)) {
    if(is.null(nlevels)) {
      if(is.function(col)) nlevels = 5
      else nlevels = length(col)
    }
    breaks <- break.quantile(x,nlevels,pretty=T)
  } else nlevels = length(breaks)-1
  # breaks,nlevels are now known
  if(is.function(col)) col = col(nlevels)
  else if(nlevels != length(col)) warning("some colors not used")
  a <- cut(x,breaks,include=T)
  if(missing(mar)) {
    mar = c(0,0,0,0.1)
    mar[3] = if(key) 2 else 1
  }
  opar = par(mar=mar)
  on.exit(par(opar))
  map(m,fill=T,color=col[as.numeric(a)[j]],res=0,bg=bg,mar=mar,...)
  if(key) {
    if(nlevels > 50) color.key(col)
    else color.key(col,breaks=breaks)
    title(main,line=1)
  }
  else title(main)
  invisible(breaks)
}

mosaicplot <- function(x, ...) UseMethod("mosaicplot")

### Changes by MM:
## - NULL instead of NA for default arguments, etc  [R / S convention]
## - plotting at end; cosmetic
## - mosaic.cell():
### Changes by KH:
##   Shading of boxes to visualize deviations from independence by
##   displaying sign and magnitude of the standardized residuals.
### Changes by minka:
##   Correct label placement and margins.

# Changes for Splus:
# mosaic.cell must be external, needs shade,rev,rot arguments.
# main cannot be NULL
# adj cannot be a vector
# draw border of polygon separately, to get line style.

mosaicplot.default <-
function(X, main = "", xlab = NULL, ylab = NULL, sort = NULL, space = NULL,
         dir = NULL, color = FALSE, shade = FALSE, margin = NULL,
         type = c("pearson", "deviance", "FT"), rev.y = F, cex=1,
         rotate = F, ...)
{
    # allocate space for the first column of X, then recurse.
    # relies on (1,1000) user coordinates.
    # X is a matrix containing indices and data.
    # space[1] is the percentage spacing.
    # maxdim[1] is the cardinality.
    # label[[1]] are the values.
    # (x1,y1)(x2,y2) are the corners of the plotting region.
    # label.x is the position of row labels
    # label.y is the position of column labels
    mosaic.cell <- function(X, x1, y1, x2, y2,
                            space, dir, color, label.x, label.y,
                            maxdim, currlev, label, shade)
    {
        # p is the column containing counts
        p <- ncol(X) - 2
        if (dir[1] == "v") {            # split here on the X-axis.
            xdim <- maxdim[1]
            # XP are the column widths
            XP <- rep(0, xdim)
            for (i in 1:xdim) {
                XP[i] <- sum(X[X[,1]==i,p]) / sum(X[,p])
            }
            white <- space[1] * (x2 - x1) / max(1, xdim-1)
            # (x.l,x.r)[i] is the drawing region for column i
            x.l <- x1
            x.r <- x1 + (1 - space[1]) * XP[1] * (x2 - x1)
            if (xdim > 1) {
                for (i in 2:xdim) {
                    x.l <- c(x.l, x.r[i-1] + white)
                    x.r <- c(x.r, x.r[i-1] + white +
                             (1 - space[1]) * XP[i] * (x2 - x1))
                }
            }
            if (label.y > 0) {
                this.lab <-
                    if (is.null(label[[1]][1])) {
                        paste(rep(as.character(currlev),
                                  length(currlev)),
                              as.character(1:xdim), sep=".")
                    } else label[[1]]
                # minka
                text(x= x.l + (x.r - x.l) / 2,
                     y= label.y,
                     srt=0, adj=c(.5,1), cex=cex, this.lab)
                h <- max(strheight(this.lab,cex=cex))
		label.y <- label.y - h
            }
            if (p > 2) {          # recursive call.
                for (i in 1:xdim) {
                    if (XP[i] > 0) {
                        mosaic.cell(as.matrix(X[X[,1]==i, -1]),
                                    x.l[i], y1, x.r[i], y2,
                                    space[-1],dir[-1],
                                    color, (i==1)*label.x, label.y,
                                    maxdim[-1],
                                    currlev+1, label[-1], shade)
                    } else {
                        segments(rep(x.l[i],3), y1+(y2-y1)*c(0,2,4)/5,
                                 rep(x.l[i],3), y1+(y2-y1)*c(1,3,5)/5)
                    }
                }
            } else { # ncol(X) <= 1 : final split polygon and segments.
                for (i in 1:xdim) {
                    if (XP[i] > 0) {
                        polygon(c(x.l[i], x.r[i], x.r[i], x.l[i]),
                                c(y1, y1, y2, y2),
                                lty = if(shade) X[i, p+1] else 1,
                                col = if(shade) {
                                    color[X[i, p+2]]
                                } else color[i])
			if(F) {
                          segments(c(rep(x.l[i],3),x.r[i]),
                                   c(y1,y1,y2,y2),
                                   c(x.r[i],x.l[i],x.r[i],x.r[i]),
                                   c(y1,y2,y2,y1))
			}
                    } else {
                        segments(rep(x.l[i],3), y1+(y2-y1)*c(0,2,4)/5,
                                 rep(x.l[i],3), y1+(y2-y1)*c(1,3,5)/5)
                    }
                }
            }
        } else { ## dir[1] - "horizontal" : split here on the Y-axis.
            ydim <- maxdim[1]
	    # YP are the row heights
            YP <- rep(0, ydim)
            for (j in 1:ydim) {
                YP[j] <- sum(X[X[,1]==j,p]) / sum(X[,p])
            }
            if(rev.y) YP <- rev(YP)
            white <- space[1] * (y2 - y1) / max(1, ydim - 1)
            # (y.b,y.t)[i] is the drawing region for row i
            y.b <- y2 - (1 - space[1]) * YP[1] * (y2 - y1)
            y.t <- y2
            if (ydim > 1) {
                for (j in 2:ydim) {
                    y.b <- c(y.b, y.b[j-1] - white -
                             (1 - space[1]) * YP[j] * (y2 - y1))
                    y.t <- c(y.t, y.b[j-1] - white)
                }
            }
            if (label.x > 0) {
                this.lab <-
                    if (is.null(label[[1]][1])) {
                        paste(rep(as.character(currlev),
                                  length(currlev)),
                              as.character(1:ydim), sep=".")
                    } else label[[1]]
                # minka
                if(rev.y) this.lab <- rev(this.lab)
                if(rotate) {
                  text(x= label.x,
                       y= y.b + (y.t - y.b) / 2,
                       srt=90, adj=c(.5,1), cex=cex, this.lab)
                  h <- max(strheight(this.lab,cex=cex,units="inch"))
                  h = h/par("pin")[1]*diff(par("usr")[1:2])
                  label.x <- label.x + h
                } else {
                  w <- max(strwidth(this.lab,cex=cex))
                  label.x <- label.x + w
                  text(x= label.x,
                       y= y.b + (y.t - y.b) / 2,
                       srt=0, adj=1, cex=cex, this.lab)
                }
                #if(label.x > left.limit) label.x <- 0
            }
            if (p > 2) {          # recursive call.
                for (j in 1:ydim) {
                    if (YP[j] > 0) {
                        mosaic.cell(as.matrix(X[X[,1]==j,2:(p+2)]),
                                    x1, y.b[j], x2, y.t[j],
                                    space[-1], dir[-1], color,
                                    label.x, (j==1)*label.y,
                                    maxdim[-1],
                                    currlev+1, label[-1], shade)
                    } else {
                        segments(x1+(x2-x1)*c(0,2,4)/5, rep(y.b[j],3),
                                 x1+(x2-x1)*c(1,3,5)/5, rep(y.b[j],3))
                    }
                }
            } else {  # ncol(X) <= 1: final split polygon and segments.
                for (j in 1:ydim) {
                    if (YP[j] > 0) {
                        polygon(c(x1,x2,x2,x1),
                                c(y.b[j],y.b[j],y.t[j],y.t[j]),
                                lty = if(shade) X[j, p+1] else 1,
                                col = if(shade) {
                                    color[X[j, p+2]]
                                } else color[j])
                        if(F) {
                          segments(c(x1,x1,x1,x2),
                                   c(y.b[j],y.b[j],y.t[j],y.t[j]),
                                   c(x2,x1,x2,x2),
                                   c(y.b[j],y.t[j],y.t[j],y.b[j]))
			}
                    } else {
                        segments(x1+(x2-x1)*c(0,2,4)/5, rep(y.b[j],3),
                                 x1+(x2-x1)*c(1,3,5)/5, rep(y.b[j],3))
                    }
                }
            }
        }
    }

    ##-- Begin main function
    if(is.null(dim(X)))
        X <- as.array(X)
    else if(is.data.frame(X))
        X <- data.matrix(X)
    dimd <- length(dX <- dim(X))
    if(dimd == 0 || any(dX == 0))
        stop("`X' must not have 0 dimensionality")
    ##-- Set up `Ind' matrix : to contain indices and data
    Ind <- 1:dX[1]
    if(dimd > 1) {
        Ind <- rep(Ind, prod(dX[2:dimd]))
        for (i in 2:dimd) {
            Ind <- cbind(Ind,
                         c(matrix(1:dX[i], byrow=TRUE,
                                  nr = prod(dX[1:(i-1)]),
                                  nc = prod(dX[i:dimd]))))
        }
    }
    Ind <- cbind(Ind, c(X))
    ## Ok, now the columns of `Ind' are the cell indices (which could
    ## also have been created by `expand.grid()' and the corresponding
    ## cell counts.  We add two more columns for dealing with *EXTENDED*
    ## mosaic plots which are produced unless `shade' is FALSE, which
    ## currently is the default.  These columns have NAs for the simple
    ## case.  Otherwise, they specify the line type (1 for positive and
    ## 2 for negative residuals) and color (by giving the index in the
    ## color vector which ranges from the ``most negative'' to the
    ## ``most positive'' residuals.
    if(is.logical(shade) && !shade) {
        Ind <- cbind(Ind, NA, NA)
    }
    else {
        if(is.logical(shade))
          shades <- c(2, 4)
        else {
          if(any(shade <= 0) || length(shade) > 5)
            stop("invalid shade specification")
          shades <- sort(shade)
          shade <- T
        }
        breaks <- c(-Inf, - rev(shades), 0, shades, Inf)
        if(T) {
          color <- c(hsv(0,               # red
                         s = seq(1, to = 0, length = length(shades) + 1)),
                     hsv(4/6,             # blue
                         s = seq(0, to = 1, length = length(shades) + 1)))
        } else {
          # for S
          color <- c(6, 6, 0,    # red
                     0, 5, 5)    # blue
        }
        if(is.null(margin))
            margin <- as.list(1:dimd)
        ## Fit the loglinear model.
        E <- loglin(X, margin, fit = TRUE, print = FALSE)$fit
        ## Compute the residuals.
        type <- match.arg(type)
        residuals <-
            switch(type,
                   pearson = (X - E) / sqrt(E),
                   deviance = {
                       tmp <- 2 * (X * log(ifelse(X==0, 1, X/E)) - (X-E))
                       tmp <- sqrt(pmax(tmp, 0))
                       ifelse(X > E, tmp, -tmp)
                   },
                   FT = sqrt(X) + sqrt(X + 1) - sqrt(4 * E + 1))
        ## And add the information to the data matrix.
        Ind <- cbind(Ind,
                     c(1 + (residuals < 0)),
                     as.numeric(cut(residuals, breaks)))
    }

    ## The next four may all be NULL:
    label <- dimnames(X)
    nam.dn <- names(label)
    if(is.null(space)) space <- rep(10, length=dimd)
    if(is.null(dir)) dir <- rep(c("v","h"), length=dimd)
    if(is.null(xlab)) xlab <- paste(nam.dn[dir=="v"],collapse=" / ")
    if(is.null(ylab)) ylab <- paste(nam.dn[dir=="h"],collapse=" / ")

    if (!is.null(sort)) {
        if(length(sort) != dimd)
            stop("length(sort) doesn't conform to dim(X)")
        ## Sort columns.
        Ind[,1:dimd] <- Ind[,sort]
        space <- space[sort]
        dir <- dir[sort]
        label <- label[sort]
    }

    ncolors <- length(tabulate(Ind[,dimd]))
    if(!shade && ((is.null(color) || length(color) != ncolors))) {
        color <- if (is.null(color) || !color[1])
            rep(0, ncolors)
        else
            2:(ncolors+1)
    }

    ##-- Plotting
    # leave space for labels, but not tick marks
    opar <- par(mar = c(1.1,1,0,0.05), mgp = c(0, 0, 0))
    #opar <- par(mar = c(0.1,0,0,0.05), mgp = c(1,1,0))
    on.exit(par(opar))
    plot.new()
    if(!shade) {
        # make the limits (1,1000) exactly
        par(usr=c(1,1000,1,1000))
    }
    else {
        ## This code is extremely ugly, and certainly can be improved.
        ## In the case of extended displays, we also need to provide a
        ## legend for the shading and outline patterns.  The code works
        ## o.k. with integer breaks in `shade'; rounding to two 2 digits
        ## will not be good enough if `shade' has length 5.
        pin <- par("pin")
        rtxt <- "Standardized\nResiduals:"
        ## Compute cex so that the rotated legend text does not take up
        ## more than 1/12 of the of the plot region horizontally and not
        ## more than 1/4 vertically.
        rtxtCex <- min(1,
                       pin[1] / (strheight(rtxt, units = "i") * 12),
                       pin[2] / (strwidth (rtxt, units = "i") / 4))
        rtxtWidth <- 0.1                # unconditionally ...
        ## We put the legend to the right of the third axis.
        # this is same as par(usr=)
        par(usr=c(1, 1000 * (1.1 + rtxtWidth), 1, 1000))
        rtxtHeight <-
            strwidth(rtxt, units = "i", cex = rtxtCex) / pin[2]
        text(1000 * (1.05 + 0.5 * rtxtWidth), 0, labels = rtxt,
             adj = c(0, 0.25), srt = 90, cex = rtxtCex)
        ## `len' is the number of positive or negative intervals of
        ## residuals (so overall, there are `2 * len')
        len <- length(shade) + 1
        ## `bh' is the height of each box in the legend (including the
        ## separating whitespace
        bh <- 0.95 * (0.95 - rtxtHeight) / (2 * len)
        x.l <- 1000 * 1.05
        x.r <- 1000 * (1.05 + 0.7 * rtxtWidth)
        y.t <- 1000 * rev(seq(from = 0.95, by = - bh, length = 2 * len))
        y.b <- y.t - 1000 * 0.8 * bh
        ltype <- c(rep(2, len), rep(1, len))
        for(i in 1 : (2 * len)) {
            polygon(c(x.l, x.r, x.r, x.l),
                    c(y.b[i], y.b[i], y.t[i], y.t[i]),
                    col = color[i],
                    lty = ltype[i])
        }
        brks <- round(breaks, 2)
        y.m <- y.b + 1000 * 0.4 * bh
        text(1000 * (1.05 + rtxtWidth), y.m,
             c(paste("<", brks[2], sep = ""),
               paste(brks[2 : (2 * len - 1)],
                     brks[3 : (2 * len)],
                     sep = ":"),
               paste(">", brks[2 * len], sep = "")),
             srt = 90, cex = cex)
    }

    if (!is.null(main) || !is.null(xlab) || !is.null(ylab))
        title(main, xlab=xlab, ylab=ylab, ...)

    # compute margins
    top.margin = right.margin = 0
    for(k in which(dir == "v")) {
      top.margin = top.margin + max(strheight(label[[k]],cex=cex))
      # right.margin must accommodate half of the last label
      # (in the worst case)
      w = strwidth(label[[k]],cex=cex)
      right.margin = max(right.margin,w[length(w)]/2)
    }
    # 0.05 in between text and plot
    top.margin = top.margin + 0.05/par("pin")[2]*diff(par("usr")[3:4])
    left.margin = 0
    for(k in which(dir == "h")) {
      if(rotate) {
        h = max(strheight(label[[k]],cex=cex,units="inch"))
        h = h/par("pin")[1]*diff(par("usr")[1:2])
        left.margin = left.margin + max(h)
      } else {
        w = strwidth(label[[k]],cex=cex)
        left.margin = left.margin + max(w)
      }
    }
    left.margin = left.margin + 0.05/par("pin")[1]*diff(par("usr")[1:2])

    # x2 is not exactly 1000 to accommodate labels which run off the end
    mosaic.cell(Ind,
                x1=left.margin, y1=1,
                x2=1000-right.margin, y2=1000-top.margin,
                space/100, dir,
                color, par("usr")[1], par("usr")[4],
                maxdim= apply(as.matrix(Ind[,1:dimd]), 2, max),
                currlev= 1, label, shade)

}
# extra routines for multivariate analysis
# Tom Minka 12/3/01

eigen2.top <- function(A,B,n=1,x,tol=1e-10) {
  d <- nrow(A)
  if(missing(x)) x <- array(rnorm(d),c(d,1))
  U <- chol(B)
  for(iter in 1:100) {
    x0 <- x
    Ax0 <- A %*% x0
    x <- backsolve(U,Ax0,transpose=T)
    x <- backsolve(U,x)
    x <- x/sqrt(sum(x*x))
    if(max(abs(x - x0)/max(abs(x))) < tol) break
  }
  e <- c()
  for(i in 1:ncol(x)) {
    xAx <- t(x[,i]) %*% A %*% x[,i]
    xBx <- t(x[,i]) %*% B %*% x[,i]
    e[i] <- xAx/xBx
  }
  list(values=e,vectors=x)
}
solve.chol <- function(a,b) {
  ch <- chol(a)
  backsolve(ch,backsolve(ch,b,transpose=T))
}

# returns x which solves a*x = b*x*eigenvalues
# k is desired number of eigenvectors (optional)
eigen2 <- function(a,b,orth=NULL,k) {
  # solve(b,a) is same as inv(b)*a
  #iba <- solve.chol(b,a)
  iba <- solve(b,a)
  if(is.null(orth)) e <- La.eigen(iba,symmetric=F)
  else {
    if(is.vector(orth)) orth <- array(orth,c(length(orth),1))
    n <- nrow(orth)
    orth.k <- ncol(orth)
    orth.m <- eye(n) - (orth %*% t(orth))
    x <- orth.m %*% iba %*% orth.m
    e <- La.eigen(x,symmetric=F)
    # remove zero eigenvalues from result
    e$vectors <- e$vectors[,1:(ncol(x)-orth.k),drop=F]
    e$values <- e$values[1:(ncol(x)-orth.k)]
  }
  if(!missing(k)) {
    e$vectors <- e$vectors[,1:k,drop=F]
    e$values <- e$values[1:k,drop=F]
  }
  # remove trivial imaginary parts
  e$vectors <- Re(e$vectors)
  e$values <- Re(e$values)
  e
}

ml.cov <- function(x) {
  n <- nrow(x)
  v <- cov(x)
  if(n == 1) v <- array(0,dim(v))
  # avoid singularity
  condition.matrix(v,n)
}
condition.matrix <- function(x,n=100) {
  #a <- max(eigen(v,only=T)$values)
  a = mean(diag(x))
  (x*(n-1) + eye(ncol(x))*a*1e-5)/n
}

#############################################################################

project <- function(x,w) {
  pred <- rownames(w)
  xw <- data.frame(data.matrix(x[pred]) %*% w)
  other <- setdiff(colnames(x),pred)
  if(length(other) > 0) cbind(xw, x[other])
  else xw
}
slice <- function(x,w,b,k=1/4) {
  # returns a subset of x containing k percent of the points
  # where w'x is near b
  px = project(x,w)[,1]
  k = ceiling(k*length(px))
  i = order(abs(px - b))[1:k]
  x[i,]
}

standardize.projection <- function(w) {
  # break symmetries by making early dimensions positive
  s = (nrow(w):1) %*% w
  for(j in 1:ncol(w)) {
    if(s[j] < 0) w[,j] <- -w[,j]
  }
  #print(t(w)%*%w)
  w
}



pca <- function(x,k=1,...) {
  x <- as.data.frame(x)
  x <- data.matrix(x[sapply(x,is.numeric)])
  s <- svd(x)
  cat("R-squared =",format(sum(s$d[1:k]^2)/sum(s$d^2)),"\n")
  w <- s$v[,1:k]
  if(k == 1) dim(w) <- c(length(w),1)
  rownames(w) <- colnames(x)
  colnames(w) <- paste("h",1:ncol(w),sep="")
  standardize.projection(w)
}


pca2 <- function(x,k=1) {
  # PCA vis Roweis's EM algorithm
  # takes O(dnk) time
  x <- as.data.frame(x)
  x <- data.matrix(x[sapply(x,is.numeric)])
  x <- t(x)
  d <- nrow(x)
  n <- ncol(x)
  w <- array(rnorm(d*k),c(d,k))
  # EM
  for(iter in 1:100) {
    old.w <- w
    if(d >= n) {
      # h is k by n
      # flops is kdk + kdn + kkn + dnk + knk + kkn = O(dnk)
      h <- solve.chol(t(w) %*% w, t(w) %*% x)
      w <- x %*% t(solve.chol(h %*% t(h), h))
    } else {
      # flops is kdk + kdn + kkd + knd + knk + kkd = O(dnk)
      h <- solve.chol(t(w) %*% w, t(w)) %*% x
      w <- t(solve.chol(h %*% t(h), h %*% t(x)))
    }
    if(max(abs(w - old.w)) < 1e-5) break
  }
  if(iter == 100) warning("not enough iterations")
  # postprocessing
  w <- qr.Q(qr(w))
  wx <- t(w) %*% x
  s <- svd(wx)
  cat("R-squared =",format(sum(s$d^2)/sum(x*x)),"\n")
  w <- w %*% s$u
  rownames(w) <- rownames(x)
  colnames(w) <- paste("h",1:ncol(w),sep="")
  standardize.projection(w)
}

projection <- function(x,y=NULL,k=1,type=c("mv","m","v","nn"),...) {
  if(is.null(y)) {
    resp <- response.var(x)
    y <- x[[resp]]
    pred <- predictor.terms(x)
    x <- x[pred]
  }
  type <- match.arg(type)
  if(type == "nn") return(projection.nn(x,y,k,...))
  if(ncol(x) > 500) stop("data has too many dimensions")
  if(all(is.na(y))) {
    # unsupervised case
    #if(type == "m") return(pca(x,k=k))
    s = projection.stats.unsup(x,type,...)
    va <- s[[1]]
    vs <- s[[2]]
    ns <- s[[3]]
  } else if(!is.factor(y)) {
    # regression projection
    s <- projection.stats.reg(x,y,...)
    va <- s[[1]]
    vs <- s[[2]]
    ns <- s[[3]]
  } else {
    xs <- split(x,y)
    vs <- lapply(xs,ml.cov)
    ns <- lapply(xs,nrow)
    va <- ml.cov(x)
  }
  w = projection.from.stats(va,vs,ns,k,type,...)
  rownames(w) <- names(x)
  colnames(w) <- paste("h",1:ncol(w),sep="")
  standardize.projection(w)
}
projection.from.stats <- function(va,vs,ns,k,type=c("mv","m","v"),...) {
  fun <- switch(type, m=projection.m.stats,
                      v=projection.cov.stats,
                      mv=projection.mv.stats)
  fun(va,vs,ns,k,...)
}
slicing <- function(x,y=NULL,...) {
  # best slice is worst projection
  if(is.null(y)) {
    resp <- response.var(x)
    y <- x[[resp]]
    pred <- predictor.terms(x)
    x <- x[pred]
  }
  k = ncol(x)
  w = projection(x,y,k=k,...)
  w[,k,drop=F]
}
optimize.slice <- function(x,y,...) {
  # returns the setting of b which gives best separation
  # x must already be projected
  if(missing(y)) {
    resp <- response.var(x)
    y <- x[[resp]]
    pred <- predictor.terms(x)
    x <- x[pred]
  }
  optimize.slice.reg(x,y,...)
}

# uses linear regression to get projection
projection.lm <- function(x,y,k=1,given=NULL) {
  vxy <- cov(x,y)
  vxx <- cov(x)
  e <- eigen2(vxy %*% t(vxy),vxx,given)
  w <- e$vectors[,1:k]
  if(k == 1) dim(w) <- c(length(w),1)
  if(is.null(given)) w
  else cbind(given,w)
}

projection.stats.unsup <- function(x,type,span,nmax=40,...) {
  d = ncol(x)
  if(missing(span)) span = 2*d
  subset = 1:nrow(x)
  if(nrow(x) > nmax) subset = sample(1:nrow(x),nmax)
  va = ml.cov(x[subset,])
  vs = list()
  ns = list()
  if(type == "m") {
    for(i in 1:length(subset)) {
      vs[[i]] = diag(d)
      ns[[i]] = 1
    }
  } else {
    for(i in 1:length(subset)) {
      # learn a covariance around x[i,]
      v = diag(d)
      for(iter in 1:100) {
        # generally takes 10 iters
        dis = sqdist(x[subset[i],],x,solve(v))
        j = order(dis)[1:(span+1)]
        if(iter > 1 && all(j == old.j)) break
        old.j = j
        v = ml.cov(x[j,])
      }
      vs[[i]] = v
      ns[[i]] = length(j)
    }
  }
  list(va,vs,ns)
}

midpoints <- function(x) (x[-length(x)]+x[-1])/2

# returns stats for a projection which separates x values leading to
# different y's
# treats it as a classification problem with n classes
# it is invariant to transformations of y
projection.stats.reg <- function(x,y,span=NULL,means=F,res=2,...) {
  ny = length(unique(y))
  # estimate covariance matrix about each point
  if(is.null(span)) {
    span <- min(c(5*ncol(x),ny/2))
    #cat("projection.reg: using nearest ",span," responses\n")
  } else {
    span <- span*ny
  }
  vs <- list()
  ms = list()
  # ys is the "classes" of y to be separated
  # res determines how many classes will be used
  ys = midpoints(break.quantile(y,ny/span*res))
  #print(length(ys))
  #ys = y
  ns = list()
  for(i in seq(ys)) {
    d <- abs(y - ys[i])
    ord <- order(d)
    my.span = span
    while(my.span < length(d)) {
      # include all points whose distance is within that of last neighbor
      j <- which(d <= d[ord[my.span]] + eps)
      #cat(i,"has",length(j),"neighbors\n")
      vs[[i]] <- ml.cov(x[j,,drop=F])
      if(sum(diag(vs[[i]])) > 0) break
      my.span = my.span + 1
    }
    ns[[i]] = my.span
    if(means) ms[[i]] = colMeans(x[j,,drop=F])
  }
  va <- ml.cov(x)
  if(means) list(va,vs,ns,ms)
  else list(va,vs,ns)
}
optimize.slice.reg <- function(x,y,...) {
  s = projection.stats.reg(x,y,means=T,...)
  va <- s[[1]]
  vs <- s[[2]]
  ns <- s[[3]]
  ms = s[[4]]
  optimize.slice.stats(va,vs,ns,ms)
}
optimize.slice.stats <- function(va,vs,ns,ms) {
  ms = as.numeric(ms)
  vs = as.numeric(vs)
  a = as.numeric(ns)*log(vs/va)/vs
  b = mean(ms)
  # this iteration isn't quite right, but seems to work
  repeat {
    old.b = b
    g = a*dnorm(b,ms,sqrt(vs))
    b = sum(g*ms)/sum(g)
    if(abs(b - old.b) < 1e-6) break
  }
  b
}

# returns a projection which separates the class means
projection.m.stats <- function(va,vs,ns,k=1,given=NULL,...) {
  n <- accumulate(ns)
  vw <- accumulate(zip(vs,ns))/n
  e <- eigen2(va,vw,given)
  i <- rev(order(e$values))[1:k]
  w <- e$vectors[,i]
  if(k == 1) dim(w) <- c(length(w),1)
  if(is.null(given)) w
  else cbind(given,w)
}

# returns a projection which separates the class covariances
projection.cov.stats <- function(va,vs,ns,k=1,given=NULL,...) {
  if(length(vs) > 2) {
    n <- accumulate(ns)
    va <- accumulate(zip(vs,ns))/n
    return(projection.mv.stats(va,vs,ns,k,given))
  }
  v1 <- vs[[1]]
  v2 <- vs[[2]]
  n <- accumulate(ns)
  a <- ns[[1]]/n

  if(T) {
    tim <- system.time(e <- eigen2(v1,v2,given))[3]
    #cat("time to compute",nrow(v1),"eigenvectors:",tim,"\n")
  } else {
    tim1 <- system.time(e1 <- eigen2.top(v1,v2,k))[3]
    tim2 <- system.time(e2 <- eigen2.top(v2,v1,k))[3]
    #cat("time to compute",2*k,"eigenvectors:",tim1+tim2,"\n")
    e <- list(values=c(e1$values,1/e2$values),vectors=cbind(e1$vectors,e2$vectors))
  }
  lambda <- abs(e$values)
  J <- a*log(lambda) - log(a*lambda + 1-a)
  # keep eigenvectors with smallest J
  i <- order(J)[1:k]
  w <- e$vectors[,i]
  if(k == 1) dim(w) <- c(length(w),1)
  if(is.null(given)) w
  else cbind(given,w)
}

score.mv <- function(w,va,vs,ns) {
  n <- accumulate(ns)
  J <- logdet(t(w) %*% va %*% w)
  for(i in 1:length(vs)) {
    J <- J - ns[[i]]/n*logdet(t(w) %*% vs[[i]] %*% w)
  }
  J
}
score.m <- function(w,va,vs,ns) {
  n <- accumulate(ns)
  J <- logdet(t(w) %*% va %*% w)
  avg.v <- array(0,dim(va))
  for(i in 1:length(vs)) {
    avg.v <- avg.v + (ns[[i]]/n)*vs[[i]]
  }
  J - logdet(t(w) %*% avg.v %*% w)
}
vec <- function(x) {
  if(is.list(x)) array(unlist(x),c(prod(dim(x[[1]])),length(x)))
  else array(as.vector(x),c(length(x),1))
}

# returns a projection which separates the class means and covariances
projection.mv.stats <- function(va,vs,ns,k=1,given=NULL,joint=F,...) {
  n <- accumulate(ns)
  if(!joint && k > 1) {
    # initialize with greedy solution
    # find projections one at a time
    w <- projection.mv.stats(va,vs,ns,k=1,given,...)
    for(i in 1:(k-1)) {
      w <- projection.mv.stats(va,vs,ns,k=1,given=w,...)
    }
    return(w)
  } else {
    # initialize with best of m,v
    w1 <- projection.m.stats(va,vs,ns,k=k,given)
    w1 <- w1[,ncol(w1)+1-(k:1),drop=F]
    J1 <- score.mv(w1,va,vs,ns)
    if(length(vs) == 2) {
      w2 <- projection.cov.stats(va,vs,ns,k=k,given)
      w2 <- w2[,ncol(w2)+1-(k:1),drop=F]
      J2 <- score.mv(w2,va,vs,ns)
    } else J2 <- -Inf
    if(J2 > J1) w <- w2
    else w <- w1
  }
  ns <- unlist(ns)
  d <- nrow(va)
  if(k == 1) vec.vs <- vec(vs)
  # EM optimization
  for(iter in 1:100) {
    oldw <- w
    vb <- array(0,dim(va))
    if(k > 1) {
      # slow
      for(i in 1:length(vs)) {
        wv <- det(t(w) %*% vs[[i]] %*% w)^(1/k)
        vb <- vb + (ns[[i]]/n/wv)*vs[[i]]
      }
    } else {
      # fast for k == 1
      wv <- t(w %x% w) %*% vec.vs
      vb <- vec.vs %*% vec(ns/n/wv)
      dim(vb) <- c(d,d)
    }
    if(F) {
      # check that objective always increases
      sw <- if(is.null(given)) w else cbind(given,w)
      print(score.mv(sw,va,vs,ns))
    }
    w <- eigen2(va,vb,given,k=k)$vector
    # force vectors to be orthogonal
    if(k > 1) w <- qr.Q(qr(w))
    if(max(abs(w - oldw)) < 1e-5) break
  }
  #cat("projection.mv: ",iter," EM iterations\n")
  if(is.null(given)) w
  else cbind(given,w)
}

projection.glm <- function(fit,k=2,...) {
  x <- model.frame(fit)
  fw <- coef(fit)[-1]
  fw <- fw/norm(fw)
  projection(x,given=fw,k=k-1,...)
}
color.plot.project.glm <- function(fit,data,type="mv",col=3,coef=F,lwd=2,...) {
  x <- model.frame(fit)
  fw <- coef(fit)[-1]
  a <- coef(fit)[1]/norm(fw)
  fw <- fw/norm(fw)
  w <- projection(x,k=1,given=fw,type=type)
  a = a*w[1,1]/fw[1]
  px = project(if(missing(data)) x else data, w)
  color.plot(px,...)
  abline(v=-a,col=col,lwd=lwd)
  invisible(w)
}

plot.axes <- function(w,col=2,origin=NULL,keep=NULL,top=NULL,
                      cex=par("cex"),labels,...) {
  # labels is to prevent it going to text()
  if(is.null(keep)) keep <- 0.2
  usr <- par("usr")
  # is the origin in the plot?
  if(is.null(origin)) {
    origin <- (usr[1] < 0) && (usr[2] > 0) && (usr[3] < 0) && (usr[4] > 0)
  }
  if(origin) m <- c(0,0)
  else m <- c((usr[2]+usr[1])/2, (usr[4]+usr[3])/2)

  width <- strwidth(rownames(w),cex=cex)/2
  height <- strheight(rownames(w),cex=cex)/2
  # "a" is scale factor to place text at arrow tips
  a <- pmin(width/abs(w[,1]), height/abs(w[,2]))
  # txt.w is the offset of the text from the arrow tip
  txt.w <- w * rep.col(a,2)

  xlim <- usr[1:2]-m[1]
  # make room on left
  xlim[1] <- xlim[1] + diff(xlim)/par("pin")[1]*0.02
  xscale <- c()
  # find biggest xscale so that
  #   xscale*w[,1] + txt.w[,1] - width > xlim[1]  (w < 0)
  #   xscale*w[,1] + txt.w[,1] + width < xlim[2]  (w > 0)
  # ignoring constraints which cannot be satisfied
  i <- which(w[,1] < 0)
  if(length(i) > 0) {
    xscale[1] <- min((xlim[1] + width[i] - txt.w[i,1])/w[i,1])
  } else xscale[1] = Inf
  i <- which(w[,1] > 0)
  if(length(i) > 0) {
    xscale[2] <- min((xlim[2] - width[i] - txt.w[i,1])/w[i,1])
  } else xscale[2] = Inf
  # if one of the sets is empty, the scale will be NA
  xscale <- min(xscale,na.rm=T)

  ylim <- usr[3:4]-m[2]
  # make room for color key
  ylim[2] <- ylim[2] - diff(ylim)/par("pin")[2]*0.06
  yscale <- c()
  # find yscale so that
  #   yscale*w[,2] + txt.w[,2] - height > ylim[1]  (w < 0)
  i <- (w[,2] < 0)
  yscale[1] <- min((ylim[1] + height[i] - txt.w[i,2])/w[i,2])
  i <- (w[,2] > 0)
  yscale[2] <- min((ylim[2] - height[i] - txt.w[i,2])/w[i,2])
  yscale <- min(yscale,na.rm=T)

  # scale equally
  xscale <- min(c(xscale,yscale))
  if(xscale < 0) {
    e = match.call()
    e$origin = F
    return(eval(e))
    stop("cannot fit axis labels on the plot.  try origin=F")
  }
  yscale <- xscale
  w <- scale(w,center=F,scale=1/c(xscale,yscale))
  if(!is.null(top)) {
    # only plot longest arrows
    wn = apply(w,1,norm)
    i = rep(F,nrow(w))
    i[rev(order(wn))[seq(top)]] = T
  } else {
    # only plot arrows of significant length
    wi = w
    wi[,1] = inchx(w[,1])
    wi[,2] = inchy(w[,2])
    wn = apply(wi,1,norm)
    i <- (wn > keep)
  }
  if(any(!i)) cat("omitting arrow for",rownames(w)[!i],"\n")
  if(!any(i)) { warning("all arrows omitted"); return() }
  len <- min(par("pin"))/20
  arrows(m[1], m[2], m[1]+w[i,1], m[2]+w[i,2], col=col, len=len)
  w <- scale(w,center=-m,scale=F)
  w <- w + txt.w
  # note: rownames(w[i,]) does not work if i picks only one
  text(w[i,1],w[i,2],labels=rownames(w)[i],col=col,cex=cex,...)
}

#############################################################################

projection.nn <- function(x,y=NULL,k=1,...) {
  if(is.null(y)) {
    resp <- response.var(x)
    y <- x[[resp]]
    pred <- predictor.terms(x)
    x <- x[pred]
  }
  # this is for y=NA
  j = nearest.points(x)$which
  # center only
  x = scale(x,scale=F)
  x = as.matrix(x)
  xj = x[j,]
  A = t(x) %*% xj
  A = A + t(A)
  xx = t(x) %*% x
  A = A + xx - (t(xj) %*% xj)
  # avoid singularity
  A = A + eye(nrow(A))*mean(diag(A))*1e-5
  s = La.eigen(A,symm=T)
  w = s$vectors[,1:k,drop=F]
  n = nrow(x)
  for(iter in 1:10) {
    old.w = w
    A = array(0,c(ncol(x),ncol(x)))
    d = (xj %*% w) - (x %*% w)
    e = rowMeans(d*d)
    s.sum = numeric(n)
    # e[i] is the distance between projected i and j
    for(ik in 1:nrow(x)) {
      d = rep.row(x[ik,] %*% w, nrow(x)) - (x %*% w)
      # ek[i] is the distance between projected i and k
      ek = rowMeans(d*d)
      # want difference to be large
      ek = ek - e
      # s[i] is weight for term (i,k)
      s = 1/(1 + exp(ek))
      xki = rep.row(x[ik,],n) - x
      xki = scale.rows(xki,s)
      A = A + (t(xki) %*% xki)
      s.sum = s.sum + s
    }
    #print(s.sum)
    cat("pts with most weight: ",order(s.sum),"\n")
    xji = xj - x
    xji = scale.rows(xji,s.sum)
    A = A - (t(xji) %*% xji)
    g = A %*% w

    if(max(abs(w - old.w)) < 1e-4) break
  }
  print(iter)
  rownames(w) <- colnames(x)
  colnames(w) <- paste("h",1:ncol(w),sep="")
  standardize.projection(w)
}

#############################################################################

identify.data.frame <- function(x,n=1,index=F,...) {
  resp = response.var(x)
  pred = predictor.vars(x)[1]
  i = identify(x[,pred],x[,resp],labels=rownames(x),n=n,...)
  if(index) i else rownames(x)[i]
}
identify.formula <- function(fmla,data=parent.frame(),...) {
  x = model.frame(fmla,data)
  identify.data.frame(x,...)
}

scale.data.frame <- function(x,center=T,scale=T,exclude.response=F,...) {
  # scales the numeric columns, leaving rest alone
  a = attr(x,"terms")
  i <- sapply(x,is.numeric)
  if(exclude.response) {
    resp = response.var(x)
    j = which(colnames(x) == resp)
    i[j] = F
  }
  if(length(center) == 1) center = rep(center,ncol(x))
  if(length(scale) == 1) scale = rep(scale,ncol(x))
  for(j in which(i)) {
    x[,j] = scale.default(x[,j],center=center[j],scale=scale[j],...)
  }
  #x[i] = apply.df(x[i],function(y) scale(y,...))
  #attr(x,"terms") = a
  x
}
hist.data.frame <- function(x,breaks,layout,col="tomato",...) {
  if(missing(layout)) layout <- auto.layout(length(x))
  opar <- par(mfrow=layout)
  on.exit(par(opar))
  for(i in names(x)) {
    if(is.numeric(x[[i]])) {
      if(missing(breaks)) breaks <- max(2,nclass.Sturges(x[[i]]))
      hist(x[[i]],xlab=i,main="",breaks=breaks,col=col,...)
    } else {
      barplot(table(x[[i]]),xlab=i,legend=F,...)
    }
  }
}

#############################################################################

factor.dichotomy <- function(f,lev,label=NULL) {
  # lev can be a vector of levels
  if(is.null(label)) label = paste(lev,collapse=",")
  labels <- c(label,paste("NOT",label))
  fd <- rep(labels[2],length(f))
  if(T) {
    fd[f %in% lev] <- labels[1]
  } else {
    for(one.lev in lev) {
      fd[f == one.lev] <- labels[1]
    }
  }
  fd <- factor(fd,levels=labels)
  names(fd) <- names(f)
  fd
}
factor.isolate <- function(f,i) {
  # add a new level containing only i
  lev = levels(f)
  nam <- names(f)
  f <- as.character(f)
  names(f) <- nam
  f[i] <- i
  lev = intersect(lev,unique(f))
  factor(f,levels=c(lev,i))
}


separate.level <- function(x,lev,type="mv",layout,keep=NULL,identify.flag=F,axes=F,...) {
  # one vs rest
  # only makes sense if length(lev) > 2
  resp <- response.var(x)
  y <- x[[resp]]
  if(missing(lev)) lev <- levels(y)
  if(missing(layout)) layout <- auto.layout(length(lev))
  # must do this or mfrow will clobber cex (bug in par?)
  wpar <- par(cex=par("cex"))
  opar <- par(mfrow=layout)
  par(wpar)
  on.exit(par(opar,wpar))
  for(k in 1:length(lev)) {
    x[[resp]] <- factor.dichotomy(y,lev[k])
    w <- projection(x,k=2,type=type)
    px <- project(x,w)
    color.plot(px,axes=axes,...)
    plot.axes(w,origin=F,keep=keep,col=3)
    if(identify.flag) {
      h <- identify(px,plot=F)
      if(length(h) > 0) text(px[h,1],px[h,2],rownames(px)[h],cex=0.9)
    }
  }
}

separate.pairwise <- function(x,level,type="mv",layout,identify.flag=F,axes=F,...) {
  # level vs each
  resp <- response.var(x)
  y <- x[[resp]]
  lev <- setdiff(levels(y),level)
  if(missing(layout)) layout <- auto.layout(length(lev))
  opar <- par(mfrow=layout)
  on.exit(par(opar))
  for(k in 1:length(lev)) {
    i <- (y == lev[k]) | (y == level)
    xk <- x[i,]
    xk[[resp]] <- factor(xk[[resp]],levels=c(level,lev[k]))
    w <- projection(xk,k=2,type=type)
    px <- project(xk,w)
    color.plot(px,axes=axes,...)
    plot.axes(w,origin=F,col=3)
    if(identify.flag) {
      h <- identify(px,plot=F)
      if(length(h) > 0) text(px[h,1],px[h,2],rownames(px)[h],cex=0.9)
    }
  }
}

logdet <- function(x) {
  2*sum(log(diag(chol(x))))
}

separation <- function(x,type="mv") {
  x <- as.data.frame(x)
  y <- x[[response.var(x)]]
  xp <- x[predictor.vars(x)]
  xs <- split(xp,y)
  vs <- lapply(xs,ml.cov)
  ns <- lapply(xs,nrow)
  n <- nrow(x)
  if(T) {
    va <- ml.cov(xp)
    J <- logdet(va)
    for(i in 1:length(ns)) {
      J <- J - ns[[i]]/n*logdet(vs[[i]])
    }
  } else {
    # loop pairs
    r <- nchoosek(length(xs),2)
    J <- 0
    for(i in 1:length(r)) {
      ri <- r[[i]]
      r1 <- ri[1]
      r2 <- ri[2]
      n1 <- ns[[r1]]
      n2 <- ns[[r2]]
      v1 <- vs[[r1]]
      v2 <- vs[[r2]]
      if(type == "m") {
        v1 <- (n1*v1 + n2*v2)/(n1+n2)
        v2 <- v1
      }
      v12 <- ml.cov(rbind(xs[[r1]],xs[[r2]]))
      Ji <- n1*logdet(v1) + n2*logdet(v2) - (n1+n2)*logdet(v12)
      J <- J - exp(Ji/2/(n1+n2))
    }
  }
  J
}

nchoosek <- function(n,k) {
  if(k == 1) return(as.list(1:n))
  r <- list()
  for(i in k:n) {
    r <- c(r,lapply(nchoosek(i-1,k-1), function(x) c(i,x)))
  }
  r
}
simple.projection <- function(d,i) {
  if(length(d) == 1) {
    d.names <- as.character(1:d)
  } else {
    d.names <- d
    d <- length(d)
  }
  w <- array(0,c(length(d.names),length(i)),list(d.names,i))
  for(j in 1:ncol(w)) w[i[j],j] <- 1
  w
}
score.simple.stats <- function(va,vs,ns,k=1,type="mv",...) {
  if(type == "v") {
    n <- accumulate(ns)
    va <- accumulate(zip(vs,ns))/n
  }
  pred <- colnames(va)
  f <- nchoosek(length(pred),k)
  r <- empty.data.frame(c(paste("Var",1:k,sep=""),"Score"))
  for(i in 1:length(f)) {
    p <- pred[f[[i]]]
    w <- simple.projection(pred,p)
    ri <- c(as.list(p), list(score.mv(w,va,vs,ns)))
    names(ri) <- names(r)
    r <- rbind.extend(r,data.frame(ri))
  }
  sort.data.frame(r)
}
# returns a frame of scores for each feature combination
score.features <- function(x,y,type=c("mv","m","v","r","lm"),...) {
  x <- as.data.frame(x)
  if(missing(y)) {
    resp <- response.var(x)
    y <- x[[resp]]
    pred <- predictor.vars(x)
    x <- x[pred]
  }
  if(missing(type)) {
    if(is.factor(y)) type <- "mv"
    else type <- "r"
  } else type <- match.arg(type)
  if(type == "r") {
    s <- projection.stats.reg(x,y,...)
    va <- s[[1]]
    vs <- s[[2]]
    ns <- s[[3]]
  } else {
    xs <- split(x,y)
    vs <- lapply(xs,ml.cov)
    ns <- lapply(xs,nrow)
    va <- ml.cov(x)
  }
  score.simple.stats(va,vs,ns,type=type,...)
}
top.features <- function(x,r,m=4,layout,type="mv",...) {
  #r <- score.features(x,k)
  k <- ncol(r)-1
  pred <- predictor.vars(x)
  if(missing(layout)) layout <- auto.layout(m)
  opar <- par(mfrow=layout)
  on.exit(par(opar))
  r <- r[nrow(r)+1-(1:m),]
  for(i in 1:m) {
    names <- as.vector(as.matrix(r[i,1:k]))
    w <- simple.projection(pred,names)
    xw <- project(x,w)
    if(k > 2) {
      w <- projection(xw,k=2,type=type)
      xw <- project(xw,w)
    }
    color.plot(xw,...)
    if(k > 2) plot.axes(w)
  }
}


#############################################################################

outliers <- function(x,p=0.05,names=T) {
  x <- x[sapply(x,is.numeric)]
  if(ncol(x) >= nrow(x)) {
    warning("more variables than cases")
    v <- diag(sd(x)^2)
  } else v <- cov(x)
  D2 <- mahalanobis(x,mean(x),v)
  pv <- 1-pchisq(D2,df=ncol(x))
  r <- data.frame(D2=D2,"p.value"=pv,row.names=rownames(x))
  r <- sort.data.frame(r)
  i <- (r$p.value < p)
  if(T) {
    #qqplot(qchisq(ppoints(length(D2)), df=ncol(x)), D2)
    if(any(i)) print(r[i,])
  }
  nam <- rownames(r)[i]
  if(names) nam else match(nam,rownames(x))

}

separate <- function(x,i,label=NULL,axes=F,...) {
  # i can be a vector of integers or row names
  if(is.integer(i)) {
    f <- factor.dichotomy(seq(nrow(x)),i,label=label)
  } else {
    f <- factor.dichotomy(rownames(x),i,label=label)
  }
  sx <- scale(as.data.frame(x))
  w <- projection(sx,f,k=2,type="m")
  mar = par("mar")
  on.exit(par(mar=mar))
  color.plot(project(sx,w),f,axes=axes,...)
  plot.axes(w,col=3,origin=F,...)
  # clean up after creating sx
  rm(sx)
  gc()
  sort(w[,1])
}

residuals.pca <- function(x,w) {
  pred <- rownames(w)
  y <- x[pred]
  m <- apply(y,2,mean)
  y <- sweep(y,2,m)
  y <- data.frame((data.matrix(y) %*% w) %*% t(w))
  y <- sweep(y,2,-m)
  x[pred] <- x[pred] - y
  x
}

#############################################################################

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
  library(mva)
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
apply.order <- function(x,ord) {
  x = x[ord[[1]],]
  x = x[,ord[[2]]]
  x
}

data.image <- function(x,las=1,reorder=T,scale=T,...) {
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
    x = apply.order(x,ord)
  }
  image.table(as.matrix(x),las=las,...)
  invisible(ord)
}

star.plot <- function(x,draw.segments=T,proto=T,reorder=T,scale=T,...) {
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
    ord = reorder.svd.data.frame(x,which(reorder))
    x = apply.order(x,ord)
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

  labeled.curves(x,y,xlab="",ylab="",yaxt=yaxt,...)
  i
}

parallel.plot <- function(x,yscale=c("linear","range","none"),flipy=F,
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

  labeled.curves(x,t(y),group=group,xlab="",ylab="",yaxt=yaxt,...)
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

#############################################################################

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
qp = function(A,b,tol=1e-10) {
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
  cbind(x,y)
}
# modified from rw1030/library/base/R/base/read.ftable
read.array <- function(file, sep = "", quote = "\"", skip = 0)
{
  z <- count.fields(file, sep, quote, skip)
  if(z[2] != max(z)) stop("unknown array representation")
  if(z[1] == 1) {
    ## Variable names and levels.  File looks like
    ##                                cvar.nam
    ## rvar.1.nam   ... rvar.k.nam    cvar.lev.1 ... cvar.lev.l
    ## rvar.1.lev.1 ... rvar.k.lev.1  ...        ... ...
    ##
    # read the header
    n.col.vars <- 1
    col.vars <- vector("list", length = n.col.vars)
    s <- scan(file, what = "", sep = sep, quote = quote,
              nlines = 2, skip = skip, quiet = TRUE)
    names(col.vars) <- s[1]
    s <- s[-1]
    n.row.vars <- z[max(which(z == max(z)))] - z[length(z)] + 1
    row.vars <- vector("list", length = n.row.vars)
    i <- 1:n.row.vars
    names(row.vars) <- s[i]
    col.vars[[1]] <- s[-i]
    z <- z[3:length(z)]
    skip <- skip + 2
  } else if(z[1] == min(z[2:length(z)])-1) {
    ## Variable levels, no names.  File looks like
    ##                                cvar.lev.1 ... cvar.lev.l
    ## rvar.1.lev.1 ... rvar.k.lev.1  ...        ... ...
    n.col.vars <- 1
    col.vars <- vector("list", length = n.col.vars)
    s <- scan(file, what = "", sep = sep, quote = quote,
              nlines = 1, skip = skip, quiet = TRUE)
    col.vars[[1]] <- s
    n.row.vars <- z[max(which(z == max(z)))] - z[length(z)] + 1
    row.vars <- vector("list", length = n.row.vars)
    z <- z[2:length(z)]
    skip <- skip + 1
  }
  else if(all(z == max(z))) {
    ## Just a block of numbers.
    values <- scan(file, sep = sep, quote = quote, quiet = TRUE, skip = skip)
    dim(values) <- c(z[1],length(z))
    return(t(values))
  } else {
    stop("unknown array representation")
  }
  # read the rest of the file
  s <- scan(file, what = "", sep = sep, quote = quote, quiet = TRUE,
            skip = skip)
  # select the data portion
  is.row.lab <- rep(rep(c(TRUE, FALSE), length(z)),
                    c(rbind(z - min(z) + 1, min(z) - 1)))
  values <- as.numeric(s[!is.row.lab])
  # label portion
  tmp <- s[is.row.lab]
  len <- length(tmp)
  p <- 1
  n <- integer(n.row.vars)
  for(k in seq(from = 1, to = n.row.vars)) {
    # n[k] is the cardinality of variable k
    n[k] <- sum(z >= max(z) - k + 1) / p
    # minka: added this line
    p <- n[k]
    # i indexes the lines
    i <- seq(from = 1, to = len, by = len / n[k])
    row.vars[[k]] <- unique(tmp[i])
    tmp <- tmp[seq(from = 2, to = len / n[k])]
    len <- length(tmp)
  }
  # minka: transposed from R version
  dim(values) <- c(sapply(col.vars, length),sapply(rev(row.vars), length))
  dimnames(values) <- c(col.vars,rev(row.vars))
  aperm(values)
}

# from rw1030/library/base/R/base/write.ftable
# modified to accept arrays
write.array <- function(x, file = "", quote = TRUE,
                        digits = getOption("digits"))
{
  ox <- x
  charQuote <- function(s) {
    if(length(s) > 1) sapply(s,charQuote)
    else {
      if(!is.null(s) && length(grep(" ",s))>0)
        paste("\"", s, "\"", sep = "")
      else s
    }
  }
  makeLabels <- function(lst) {
      lens <- sapply(lst, length)
      cplensU <- c(1, cumprod(lens))
      cplensD <- rev(c(1, cumprod(rev(lens))))
      y <- NULL
      for (i in rev(seq(along = lst))) {
          ind <- 1 + seq(from = 0, to = lens[i] - 1) * cplensD[i + 1]
          tmp <- character(length = cplensD[i])
          tmp[ind] <- charQuote(lst[[i]])
          y <- cbind(rep(tmp, times = cplensU[i]), y)
      }
      y
  }
  dn <- dimnames(x)
  if(is.null(dn)) {
    x <- format(unclass(x),digits=digits)
    cat(t(x), file = file, sep = c(rep(" ", ncol(x) - 1), "\n"))
    return(invisible(ox))
  }
  xrv <- dn[1:(length(dn)-1)]
  xcv <- dn[length(dn)]
  LABS <- rbind(matrix("", nr = length(xcv), nc = length(xrv)),
                charQuote(names(xrv)),
                makeLabels(xrv))
  x <- aperm(x)
  dim(x) <- c(prod(sapply(xcv,length)),prod(sapply(xrv,length)))
  x <- t(x)
  DATA <- rbind(t(makeLabels(xcv)),
                format(unclass(x), digits = digits))
  if(!is.null(names(xcv))) {
    DATA <- rbind(c(charQuote(names(xcv)),rep("",ncol(x)-1)),DATA)
  }
  x <- cbind(apply(LABS, 2, format, justify = "left"),
             apply(DATA, 2, format, justify = "right"))
  cat(t(x), file = file, sep = c(rep(" ", ncol(x) - 1), "\n"))
  invisible(ox)
}
# extra routines for linear regression
# Tom Minka

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

##############################################################################

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
predict.plot.default <- function(x,y,labels,xlab,ylab,...) {
  frame = data.frame(x,y)
  if(!missing(labels)) rownames(frame) = labels
  if(missing(xlab)) xlab = deparse(substitute(x))
  if(missing(ylab)) ylab = deparse(substitute(y))
  predict.plot.data.frame(frame,xlab=xlab,ylab=ylab,...)
}
std.error = function(x) sd(x)/sqrt(length(x))
predict.plot.data.frame <- function(x,layout=NULL,partial=NULL,condensed=T,
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
predict.plot.lm <- function(object,data,partial=F,ylab=NULL,...) {
  if(!partial) {
    if(is.null(ylab)) ylab = paste("Residual",response.var(object))
    if(missing(data)) {
      res <- residual.frame(object)
    } else {
      res <- residual.frame(object,data)
    }
    cat("plotting residuals\n")
    predict.plot.data.frame(res,highlight=predictor.terms(object),ylab=ylab,...)
  } else {
    if(missing(data)) data <- model.frame(object)
    else data = my.model.frame(formula(paste(response.var(object),"~.")),data)
    cat("plotting partial residuals\n")
    predict.plot.data.frame(data,partial=object,ylab=ylab,...)
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
predict.plot.formula <- function(formula,data=parent.frame(),...) {
  # formula has givens?
  given = given.vars(formula)
  if(!is.null(given)) {
    formula = remove.given(formula)
    if(is.environment(data)) g <- get(given,env=data)
    else g <- data[[given]]
    if(is.null(g))
      stop(paste("variable \"",given,"\" not found",sep=""))
    return(predict.plot.formula(formula,data,
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
  predict.plot.data.frame(x,...)
}
predict.plot <- function(...) UseMethod("predict.plot")

##############################################################################

expand.cross <- function(object) {
  resp <- response.var(object)
  pred <- predictor.vars(object)
  # quadratic terms only
  pred2 <- c()
  for(i in pred) {
    for(j in pred) {
      if(match(i,pred) < match(j,pred)) pred2 = c(pred2,paste(i,j,sep=":"))
    }
  }
  pred = c(pred,pred2)
  formula(paste(resp,"~",paste(pred,collapse="+")))
}
step.up <- function(object,scope=expand.cross(object)) {
  step(object,list(upper=scope,lower=formula(object)))
}

interact.plot.data.frame <- function(x,ypred,partial=NULL,highlight,span=0.75,
                                     scol="red",slwd=2,type=c("*",":"),
                                     xaxt="n",yaxt="n",se=T,
                                     bg=par("bg"),nlev=8,main=NULL,...) {
  # type=":" means partials remove cross term only
  # type="*" means partials remove everything
  library(modreg)
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
          predict.plot(fmla,y,xlab="",ylab="",xaxt=xaxt,yaxt=yaxt,key=F)
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
interact.plot <- function(object, ...) UseMethod("interact.plot")

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

##############################################################################
# glm stuff

predict.frame <- function(object) {
  # returns a frame relating predicted values to actual values
  # for glms, the predicted value is the value before transformation
  resp = response.var(object)
  frame = data.frame(z=predict(object),response=model.frame(object)[[resp]])
  names(frame)[2] = resp
  frame
}
model.plot <- function(object,data,pred=predictor.terms(data),
                       highlight,layout=NULL,
                       partial=F,smooth=F,se=F,span=2/3,
                       col="green",scol=2,lwd=2,add=F,lty=1,
                       bg=par("bg"),key=T,type="p",...) {
  # pred=predictor.vars(data)?
  # this only makes predictions at the data, while plot.loess does more
  if(missing(data)) data <- model.frame(object)
  resp = response.var(object)
  data <- expand.frame(data,setdiff(pred,predictor.vars(data)))
  if(inherits(object,"tree")) {
    p = predict(object,data,type="vector")[,2]
  } else {
    p <- predict(object,data,type="response",se=se)
  }
  if(is.factor(data[[resp]])) {
    if(!se) p = p+1
    else p$fit = p$fit+1
  }
  if(add) {
    # must have one predictor term
    if(length(pred) > 1) stop("model has more than one predictor")
    x <- data[[pred]]
    i <- order(x)
    x <- sort(x)
    if(!se) {
      lines(x,p[i],col=col,lwd=lwd,lty=lty,...)
    } else {
      lines(x,p$fit[i],col=col,lwd=lwd,lty=lty,...)
      lines(x,p$fit[i]+p$se[i],col=col,lty=2,lwd=lwd,...)
      lines(x,p$fit[i]-p$se[i],col=col,lty=2,lwd=lwd,...)
    }
    return(invisible(p))
  }
  if(is.null(layout)) layout <- auto.layout(length(pred))
  opar <- par(bg=bg,mar=auto.mar(main="f"))
  if(any(layout != c(1,1))) {
    opar = c(opar,par(mfrow=layout))
  }
  on.exit(par(opar))
  if(missing(highlight)) highlight = predictor.terms(object)
  for(v in pred) {
    #print(data[[v]])
    plot.default(data[[v]],data[[resp]],xlab="",ylab=resp,type=type,...)
    font.lab = if(v %in% highlight) 2 else 1
    title(xlab=v,font.lab=font.lab)
    x = data[[v]]
    if(length(predictor.terms(object)) == 1 &&
       all(pred==predictor.terms(object))) {
      i <- order(x)
      lines(x[i],p[i],col=col,lwd=lwd,lty=lty,...)
    } else {
      if(partial) {
        fit <- my.update(object,paste(".~. -",v))
        if(v %in% predictor.terms(fit)) warning("update didn't work")
        p = predict(fit,data,type="response",se=se)
        if(!se) p = p+1
        else p$fit = p$fit+1
      }
      # doesn't work
      #i <- order(x)
      #lines(x[i],p[i],col=col,lwd=lwd,lty=lty,...)

      #fit = smooth(p~x)
      #fit = loprob(x,p,span=span)
      # unfortunately lowess uses degree=1
      fit = lowess(x,p,f=span,iter=0)
      lines(fit,col=col,lwd=lwd,lty=lty,...)
    }
    if(smooth) {
      if(is.factor(data[[resp]])) fit = loprob(x,data[[resp]],span=span)
      else fit = lowess(x,data[[resp]],f=span)
      lines(fit,col=scol,lwd=lwd,lty=lty,...)
      if(key) color.key(c(col,scol),labels=c("predicted","actual"))
    }
  }
}

predictClass <- function(object,data,p=0.5) {
  resp <- response.var(object)
  if(inherits(object,"glm")) {
    r = predict(object,data,type="response")
    factor.logical(r>p,labels=levels(data[[resp]]))
  } else if(p == 0.5) {
    predict(object,data,type="class")
  } else {
    r = predict(object,data,type="vector")[,2]
    factor.logical(r>p,labels=levels(data[[resp]]))
  }
}
misclassResiduals <- function(object,data,p=0.5) {
  r = predictClass(object,data,p=p)
  y = data[[response.var(object)]]
  factor.logical(r != y, labels=c("Correct","Error"))
}

misclass.default <- function(fit,data=NULL,rate=F) {
  if(is.null(data)) data = model.frame(fit)
  resp <- response.var(fit)
  r = predictClass(fit,data)
  r = sum(r != data[[resp]])
  if(rate) r/nrow(data) else r
}
misclass <- function(object, ...) UseMethod("misclass")

formula.nnet = formula.lm

deviance.glm <- function(object,data,rate=F) {
  if(missing(data)) {
    dev = object$deviance
    if(rate) dev = dev/nrow(model.frame(fit))
    return(dev)
  }
  resp <- response.var(object)
  if(object$family$family == "binomial") {
    truth <- (as.numeric(data[[resp]]) == 2)
    p <- predict(object,data,type="response")
    p[p == 0] <- 1e-3
    p[!truth] <- 1-p[!truth]
    dev = -2*sum(log(p))
  } else {
    stop("family not handled")
  }
  if(rate) dev/nrow(data) else dev
}

confusion.default <- function(object,data=NULL,p=0.5) {
  if(is.null(data)) data = model.frame(object)
  resp <- response.var(object)
  r = predictClass(object,data,p=p)
  table(truth=data[[resp]],predicted=r)
}
confusion <- function(object, ...) UseMethod("confusion")

factor.logical <- function(x,labels=c("No","Yes")) {
  f <- factor(x,levels=c(F,T))
  levels(f) <- labels
  f
}


# expands a model formula to include all squares and cross-products of the
# predictors
expand.quadratic <- function(fmla,cross=T) {
  resp <- response.var(fmla)
  pred <- predictor.vars(fmla)
  pred2 <- sapply(pred,function(s) paste("I(",s,"^2)",sep=""))
  s <- paste(resp,"~",paste(pred,collapse="+"),
             "+",paste(pred2,collapse="+"))
  len <- length(pred)
  if(cross && len > 1) {
    pred3 <- c()
    for(i in 1:(len-1)) {
      for(j in (i+1):len) {
        pred3 <- c(pred3, paste(pred[i],":",pred[j],sep=""))
      }
    }
    s <- paste(s,"+",paste(pred3,collapse="+"))
  }
  formula(s)
}

logistic <- function(formula,data,...) {
  #environment(formula) = sys.frame(sys.nframe())
  control = glm.control(maxit=30,trace=F)
  # can't use control else update will fail
  glm(formula.with.data(formula,data),family=binomial,...)
  #glm(formula,data,family=binomial,control=control,...)
}

my.logistic <- function(formula,data,lambda=NULL,...) {
  data = model.frame(formula,data)
  pred = predictor.terms(formula)
  resp = response.var(formula)
  edata = expand.frame(data,setdiff(pred,colnames(data)))
  x = as.matrix(edata[pred])
  x = cbind("(Intercept)"=numeric(nrow(x))+1,x)
  y = 2*as.numeric(data[[resp]])-3
  w = logistic.fit(scale.rows(x,y),lambda=lambda)
  object = list(coefficients=w,model=data,levels=levels(data[[resp]]),
    terms=terms(data))
  class(object) = "logistic"
  object
}
print.logistic <- function(object) {
  cat("Logistic regression on ")
  dput(object$levels)
  print(object$coefficients)
}
coef.logistic = coef.glm
model.frame.logistic = model.frame.glm
predict.logistic <- function(object,data=NULL,
                             type=c("class","response","vector"),se=F) {
  pred = predictor.terms(object)
  if(is.null(data)) data = model.frame(object)
  else {
    no.resp = formula(paste("~",paste(pred,collapse="+")))
    data = model.frame(no.resp,data)
  }
  edata = expand.frame(data,setdiff(pred,colnames(data)))
  x = as.matrix(edata[pred])
  w = t(t(object$coefficients[-1]))
  z = (x %*% w) + object$coefficients[1]
  type = match.arg(type)
  s = 1/(1+exp(-z))
  if(type == "response") s
  else if(type == "class") object$levels[(s>0.5)+1]
  else if(type == "vector") {
    r = array(0,c(nrow(x),2),list(rownames(data),object$levels))
    names(dimnames(r)) = c("case",response.var(object))
    r[,1] = 1-s
    r[,2] = s
    r
  }
}

deviance.logistic <- function(object,data=NULL,rate=F) {
  if(is.null(data)) data = model.frame(object)
  resp <- response.var(object)
  truth = as.numeric(data[[resp]])
  p = predict(object,data,type="vector")
  p[p == 0] = 1e-3
  i = (1:nrow(p)) + (truth-1)*nrow(p)
  dev = -2*sum(log(p[i]))
  if(rate) dev/length(truth) else dev
}

logistic.fit <- function(x,lambda=NULL) {
  # each row of x is a data point
  d = ncol(x)
  if(is.null(lambda)) lambda = 1e-4
  w = array(0,c(1,d))
  colnames(w) = colnames(x)
  for(iter in 1:100) {
    old.w = w
    # s1 is column
    s1 = 1/(1+exp(x %*% t(w)))
    a = s1*(1-s1)
    # g is row
    g = (t(s1) %*% x) - lambda*w
    h = (t(x) %*% scale.rows(x,a)) + lambda*diag(d)
    w = w + t(solve(h,t(g)))
    #cat("loglik =",sum(log(1-s1)),"\n")
    #print(drop(w))
    if(max(abs(w - old.w)) < 1e-6) break
  }
  drop(w)
}

logistic.fit.bohning <- function(x,lambda=1e-4) {
  # each row of x is a data point
  d = ncol(x)
  h = (t(x) %*% x)/4 + lambda*diag(d)
  # r is upper triangular
  r = chol(h)
  w = array(0,c(1,d))
  for(iter in 1:10) {
    old.w = w
    # s1 is column
    s1 = 1/(1+exp(x %*% t(w)))
    # g is row
    g = (t(s1) %*% x) - lambda*w
    # u is column
    u = backsolve(r,forwardsolve(t(r),t(g)))
    # line search along u
    ug = drop(g %*% u)
    ux = x %*% u
    a = s1*(1-s1)
    uhu = sum((ux^2)*a) + lambda*sum(u*u)
    w = w + (ug/uhu)*t(u)
    print(drop(w))
    if(max(abs(w - old.w)) < 1e-6) break
  }
  colnames(w) = colnames(x)
  w
}
# Functions for document and image retrieval
# By Tom Minka

# Functions for lab 1

strip.text <- function(txt) {
  # remove apostrophes
  txt <- gsub("'","",txt)
  # convert to lowercase
  txt <- tolower(txt)
  # change other non-alphanumeric characters to spaces
  # "man 7 regex" for info
  txt <- gsub("[^a-z0-9]"," ",txt)
  # change digits to #
  txt <- gsub("[0-9]+","#",txt)
  # split and make one array
  txt <- unlist(strsplit(txt," "))
  # remove empty words
  txt <- txt[txt != ""]
  txt
}
read.doc <- function(fname,remove.header=T) {
  txt <- readLines(fname)
  if(remove.header) {
    # remove header
    i <- which(txt == "")
    if(length(i) > 0) {
      txt <- txt[-(1:i[1])]
    }
  }
  strip.text(txt)
}

remove.singletons <- function(x) {
  # remove infrequent words
  count <- colSums(x>0)
  x[,(count > 1)]
}

idf.weight <- function(x) {
  # IDF weighting
  doc.freq <- colSums(x>0)
  doc.freq[doc.freq == 0] <- 1
  w <- log(nrow(x)/doc.freq)
  scale.cols(x,w)
}
idf.weight2 <- function(x) {
  # Joachims' IDF
  doc.freq <- colMeans(div.by.sum(x))
  w <- sqrt(1/doc.freq)
  scale.cols(x,w)
}

div.by.sum <- function(x) {
  scale.rows(x,1/(rowSums(x)+1e-16))
}
div.by.euc.length <- function(x) {
  scale.rows(x,1/sqrt(rowSums(x^2)+1e-16))
}


if(T) {
remove.singletons.ragged <- function(x) {
  # x is a list of vectors (a ragged array)
  col.names <- c()
  for(i in 1:length(x)) {
    col.names <- c(col.names, names(x[[i]]))
  }
  count <- table(col.names)
  for(i in 1:length(x)) {
    not.single <- (count[names(x[[i]])] > 1)
    x[[i]] <- x[[i]][not.single]
  }
  x
}
standardize.ragged <- function(x) {
  # x is a list of vectors (a ragged array)
  # standardize each vector to have same length and ordering, by adding NAs
  col.names <- c()
  for(i in 1:length(x)) {
    # this is faster than unique(c(col.names, names(x[[i]])))
    col.names <- c(col.names, setdiff(names(x[[i]]),col.names))
  }
  for(i in 1:length(x)) {
    #not.in <- setdiff(col.names,names(x[[i]]))
    #x[[i]][not.in] <- NA
    #is.na(x[[i]][not.in]) <- T
    # this automatically returns NA for missing names
    x[[i]] <- x[[i]][col.names]
    names(x[[i]]) <- col.names
  }
  x
}
}

mds <- function(d,y=NULL,labels=T,k=2,main="Multidimensional scaling",...) {
  library(MASS)
  for(iter in 1:2) {
    px = try(as.data.frame(sammon(d,k=k,trace=F)$points), silent=TRUE)
    if(!inherits(px,"try-error")) break
    d = d + 1e-4*array(runif(prod(dim(d))),dim(d))
    diag(d) = 0
  }
  if(inherits(px,"try-error")) stop(px)
  if(k == 2) {
    if(is.null(y)) {
      if(is.null(labels))
        plot(px,asp=1,main=main,...)
      else
        text.plot(px,labels,asp=1,main=main,...)
    } else {
      color.plot(px,y,labels=labels,asp=1,zlab=main,...)
    }
  } else if(k == 3) {
    if(is.null(y)) {
      vrml.plot3d(formula(px),px,labels=labels,...)
    } else {
      vrml.color.plot(px,y,labels=labels,...)
    }
  }
  invisible(px)
}

# Functions for lab 2

as.oldpix <- function(pix) {
  d = attr(pix,"size")
  x = array(0,c(d,3))
  x[,,1] = attr(pix,"red")
  x[,,2] = attr(pix,"green")
  x[,,3] = attr(pix,"blue")
  x
}

flatten.pixmap <- function(x) {
  if(inherits(x,"pixmapRGB")) x = as.oldpix(x)
  len <- prod(dim(x)[1:2])
  # convert pixmap to array with cols r,g,b
  r <- array(0,c(len,3),list(NULL,c("red","green","blue")))
  for(i in 1:3) {
    r[,i] <- as.vector(x[,,i])
  }
  as.data.frame(r)
}

plot.pixels <- function(x,k=2,as.hsv=T) {
  col <- rgb(x[[1]],x[[2]],x[[3]])
  if(as.hsv) x = rgb2hsv(x)
  if(k==3)
    vrml.plot3d(formula(x),x,col=col,bg=grey(0.5))
  else
    plot3d(formula(x),x,color=col,bg=grey(0.5))
}

quantize.cube <- function(x,k) {
  # quantize to integer
  # x has three columns
  # first col varies fastest
  r <- 0
  d <- 1
  for(i in 1:ncol(x)) {
    r <- r + as.integer(x[,i]*k/1.00001)*d
    d <- d*k
  }
  r
}

closest <- function(d) {
  # for each column, return the name of the row with smallest distance
  if(dim(d)[1] == dim(d)[2] && all(diag(d)) == 0) {
    diag(d) <- Inf
  }
  r <- rownames(d)[apply(d,2,which.min)]
  names(r) <- colnames(d)
  r
}

if(T) {
doc.page <- function(index.file,base="",cols=5) {
  # make a web page with document names
  # output is index.file with extension changed to html
  base = path.expand(base)
  tab <- read.table(index.file,header=T)
  col <- 1
  s <- strsplit(index.file,"\\.")
  fname <- paste(s[[1]][1],"html",sep=".")
  cat("wrote",fname,"\n")
  con <- file(fname,"w")
  cat("<html><body><table>\n",file=con)
  for(i in 1:nrow(tab)) {
    name <- as.character(tab[i,1])
    fname <- as.character(tab[i,3])
    fname = file.path(base,fname)
    if(col == 1) cat("<tr align=center>\n",file=con)
    cat("<td>",file=con)
    cat("<a href=\"",fname,"\">",name,"</a>\n",sep="",file=con)
    cat("</td>",file=con)
    if(col == cols) { cat("</tr>\n",file=con); col <- 0 }
    col <- col + 1
  }
  cat("</table></body></html>\n",file=con)
  close(con)
}

image.page <- function(index.file,base="",cols=5) {
  # make a web page with images
  base = path.expand(base)
  tab <- read.table(index.file,header=T)
  col <- 1
  s <- strsplit(index.file,"\\.")
  fname <- paste(s[[1]][1],"html",sep=".")
  cat("wrote",fname,"\n")
  con <- file(fname,"w")
  cat("<html><body><table>\n",file=con)
  for(i in 1:nrow(tab)) {
    name <- as.character(tab[i,1])
    fname <- as.character(tab[i,2])
    s <- strsplit(fname,"\\.")
    fname <- paste(s[[1]][1],"jpg",sep=".")
    if(nchar(base) > 0) fname = file.path(base,fname)
    if(col == 1) cat("<tr align=center>\n",file=con)
    cat("<td>",file=con)
    cat("<img src=\"",fname,"\"><br>",name,"\n",sep="",file=con)
    cat("</td>",file=con)
    if(col == cols) { cat("</tr>\n",file=con); col <- 0 }
    col <- col + 1
  }
  cat("</table></body></html>\n",file=con)
  close(con)
}

read.images <- function(index.file,k=4) {
  # read images and convert to matrix of color counts
  library(pixmap)
  tab <- read.table(index.file,header=T)
  if(T) {
    midpoints <- function(x) (x[-1]+x[-length(x)])/2
    col <- midpoints(seq(0,1,len=k+1))
    col <- as.matrix(expand.grid(col,col,col))
    labels <- rgb2name(col2rgb(hsv(col[,1],col[,2],col[,3]))/255)
    labels <- make.unique(labels)
    #if(any(duplicated(labels))) stop("duplicate color names")
  } else {
    labels <- 0:(k^3-1)
  }
  imgs <- array(0,c(nrow(tab),k^3))
  rownames(imgs) <- as.character(tab[,1])
  for(i in 1:nrow(tab)) {
    fname <- as.character(tab[i,2])
    print(fname)
    img <- flatten.pixmap(read.pnm(fname))
    img <- rgb2hsv(img)
    #img <- hsv2cone(img)
    q <- quantize.cube(img,k)
    x <- table(factor(q,levels=0:(k^3-1),labels=labels))
    imgs[i,] <- x
    colnames(imgs) <- names(x)
  }
  imgs
}
}

# Functions for lab 3

subtable <- function(x,i,j) {
  if(!inherits(x,"table")) x <- as.table(as.matrix(x))
  if(is.character(i)) i <- pmatch(i,rownames(x))
  if(is.character(j)) j <- pmatch(j,colnames(x))
  rs <- sum(x[i,])
  cs <- sum(x[,j])
  n <- sum(x)
  r <- array(c(x[i,j],cs-x[i,j],rs-x[i,j],n-cs-rs+x[i,j]),c(2,2))
  rownames(r) <- c(rownames(x)[i], paste("not",rownames(x)[i]))
  colnames(r) <- c(colnames(x)[j], paste("not",colnames(x)[j]))
  names(dimnames(r)) <- names(dimnames(x))
  r
}

entropy <- function(x) {
  x = as.vector(x)
  x = x/sum(x)
  -sum(x*log2(x+eps))
}
as.factor.data.frame <- function(x) {
  for(i in 1:length(x)) if(!is.factor(x[[i]])) x[[i]] = factor(x[[i]])
  x
}
tabulate.matrix <- function(x) {
  if(is.logical(x)) {
    s = colSums(x)
    rbind("FALSE"=nrow(x)-s,"TRUE"=s)
  } else {
    nb = max(x)
    apply(x,2,function(f) tabulate(f,nb=nb))
  }
}
col.entropy <- function(x) {
  # x is matrix of counts
  # entropy of each column of x
  x = x+eps
  s = colSums(x)
  -colSums(x*log2(x))/s +log2(s)
}
information <- function(x,y=NULL,actual=NULL,smooth=0) {
  if(is.null(y)) {
    # x is a matrix of counts
    x = x + smooth
    if(is.null(actual)) {
      # average information
      if(F) {
        p = colSums(x)
        p = p/sum(p)
        entropy(rowSums(x)) - sum(p*apply(x,2,entropy))
      } else {
        entropy(rowSums(x))+entropy(colSums(x))-entropy(x)
      }
    } else {
      # conditional information
      entropy(rowSums(x)) - entropy(x[,actual])
    }
  } else {
    # x is a matrix of factors, y is a factor
    # return the information in each table(x[,i],y)
    p = table(y) + 2*smooth
    p = p/sum(p)
    if(is.null(actual)) {
      h = col.entropy(tabulate.matrix(x)+smooth)
      for(lev in levels(y)) {
        h = h - col.entropy(tabulate.matrix(x[y==lev,])+smooth)*p[lev]
      }
      h
    } else {
      tab = array(0,c(nlevels(y),ncol(x)),list(levels(y),colnames(x)))
      for(lev in levels(y)) {
        tab[lev,] = tabulate.matrix(x[y==lev,])[actual,]
      }
      tab = tab + smooth
      entropy(p) - col.entropy(tab)
    }
  }
}
all.subsets <- function(n) {
  if(n == 1) subsets = list(1)
  else if(n == 2) subsets = list(1,2,1:2)
  else if(n == 3) subsets = list(1,2,3,1:2,2:3,c(1,3),1:3)
  else {
    subsets = list()
    for(k in 1:n) {
      subsets = c(subsets,nchoosek(n,k))
    }
  }
  subsets
}
interaction.information <- function(x,y=NULL,actual=NULL) {
  if(is.null(y)) {
    if(!is.null(actual)) {
      return(interaction.information(x[,,actual]) -
             interaction.information(margin.table(x,1:2)))
    }
    nd = if(is.null(dim(x))) 1 else length(dim(x))
    subsets = all.subsets(nd)
    total = 0
    for(s in subsets) {
      sig = (-1)^(length(s)-nd+1)
      total = total + entropy(margin.table(x,s))*sig
    }
    total
  } else {
    # x is a matrix of factors, y is a factor
    # return the interaction in each table(x[,i],x[,j],y)
  }
}
interaction.variance <- function(x) {
  n = dim(x)[length(dim(x))]
  a = c()
  for(i in 1:n) {
    a[i] = interaction.information(x,actual=i)
  }
  w = margin.table(x,length(dim(x)))
  w = w/sum(w)
  m = sum(w*a)
  print(m)
  sum(w*(a-m)^2)
}
information.graph <- function(tab,main="",color=grey(0.7)) {
  tab = as.table(tab)
  nd = if(is.null(dim(tab))) 1 else length(dim(tab))
  subsets = all.subsets(nd)
  theta = seq(0,2*pi,len=nd+1)[-1]
  x = cbind(cos(theta),sin(theta))*4
  if(nd > 3) {
    x = x + array(rnorm(prod(dim(x))),dim(x))
  }
  mar=auto.mar(main=main,axes=F,xlab="",ylab="")
  opar = par(mar=mar)
  on.exit(par(opar))
  plot.new()
  plot.window(c(-5,5),c(-5,5),asp=1)
  for(s in rev(subsets)) {
    mid = colMeans(x[s,,drop=F])
    for(i in s) {
      segments(mid[1],mid[2],x[i,1],x[i,2])
    }
    h = interaction.information(margin.table(tab,s))
    col = if(h < 0) color else "white"
    draw.circle(mid[1],mid[2],sqrt(abs(h)),col=col,fill=T)
  }
  labels = names(dimnames(tab))
  text(x[,1],x[,2],labels)
}

score.binary <- function(x,y,type=c("info","a.info"),a=1) {
  if(!inherits(x,"table")) x <- as.table(as.matrix(x))
  y = table(as.numeric(y))
  x = x + a
  y = y + 2*a
  n = sum(y)
  type <- match.arg(type)
  if(type == "info") {
    not.x = rep.col(y,ncol(x)) - x
    r = (x/n)*log2(x/n) + (not.x/n)*log2(not.x/n)
    r = colSums(r)
    r = r - sum((y/n)*log2(y/n))
    x = colSums(x)
    not.x = n - x
    r - (x/n)*log2(x/n) - (not.x/n)*log2(not.x/n)
  } else if(type == "a.info") {
    col.total = rep.row(colSums(x),nrow(x))
    r = (x/col.total)*log2(x/col.total)
    r = colSums(r)
    r = r - sum((y/n)*log2(y/n))
    r
  }
}

score.words <- function(x,type=c("info","bayes","chisq",
                            "a.info","odds","lift","overlap"),
                        a=1,z=1) {
  # info,bayes,chisq are very similar
  # a.info,odds are similar
  # input is contingency table (matrix of counts)
  if(!inherits(x,"table")) x <- as.table(as.matrix(x))
  x = x + a
  fit <- indep.fit(x)
  type <- match.arg(type)
  if(type == "info") {
    row.total <- rep.col(rowSums(x),ncol(x))
    not.x = row.total - x
    n = sum(as.vector(x))
    r = (x/n)*log2(x/n) + (not.x/n)*log2(not.x/n)
    r = r - (row.total/n)*log2(row.total/n)
    r = colSums(r)
    x = colSums(x)
    not.x = n - x
    r - (x/n)*log2(x/n) - (not.x/n)*log2(not.x/n)
  } else if(type == "a.info") {
    row.total <- rep.col(rowSums(x),ncol(x))
    col.total = rep.row(colSums(x),nrow(x))
    n = sum(as.vector(x))
    r = (x/col.total)*log2(x/col.total)
    r = r - (row.total/n)*log2(row.total/n)
    r = colSums(r)
    r
  } else if(type == "bayes") {
    # extra smoothing
    x = x + 1
    row.total <- rep.col(rowSums(x),ncol(x))
    not.x = row.total - x
    r = lgamma(x) + lgamma(not.x) - lgamma(row.total)
    r = colSums(r)
    x = colSums(x)
    row.total = colSums(row.total)
    not.x = row.total - x
    r - (lgamma(x)+lgamma(not.x)-lgamma(row.total))
  } else if(type == "chisq") {
    r <- (x - fit)^2/fit
    total.words <- rep.col(rowSums(x),ncol(x))
    r <- r + (x - fit)^2/(total.words - fit)
    colSums(r)
  } else if(type == "lift") {
    r <- (x - sqrt(x))/fit
    #colSums(-log(abs(r)+1e-4))
    #colSums(r^2)
    apply(r,2,max)
  } else if(type == "odds") {
    n <- sum(x)
    cs <- rep.row(colSums(x),nrow(x))
    rs <- rep.col(rowSums(x),ncol(x))
    r <- log(x/(rs-x)/(cs-x)*(n-rs-cs+x))
    v <- 1/x + 1/(rs-x) + 1/(cs-x) + 1/(n-rs-cs+x)
    #colSums((r - z*sqrt(v))^2)
    # subtract the expected sum of squares under independence
    colSums(r^2 - z*v)
    #colSums(abs(r)) - z*sqrt(colSums(v))
  } else if(type == "overlap") {
    cs <- colSums(x)
    r <- apply(x,2,max)/cs
    v <- r*(1-r)/cs
    r-z*sqrt(v)
  }
}

odds.ratios <- function(x) {
  # returns the log odds-ratio for all cells in the table
  n <- sum(as.vector(x))
  cs <- rep.row(colSums(x),nrow(x))
  rs <- rep.col(rowSums(x),ncol(x))
  log(x/(rs-x)/(cs-x)*(n-rs-cs+x))
  #log(x/(cs-x))
}
odds.ratios.se <- function(x) {
  # returns the standard error of each log odds-ratio
  n <- sum(x)
  cs <- rep.row(colSums(x),nrow(x))
  rs <- rep.col(rowSums(x),ncol(x))
  sqrt(1/x + 1/(rs-x) + 1/(cs-x) + 1/(n-rs-cs+x))
  #sqrt(1/x + 1/(cs-x))
}

discriminate.node <- function(g,i=1,x) {
  k <- names(g)[leaves(g,i)]
  if(length(k) == 0) return("")
  not.k <- setdiff(rownames(x),k)
  if(length(not.k) == 0) return("")
  if(!is.character(i)) i <- names(g)[i]
  not.i <- paste("NOT",i)
  xp <- array(0,c(2,ncol(x)),list(c(i,not.i),colnames(x)))
  xp[1,] <- colSums(x[k,,drop=F])
  xp[2,] <- colSums(x[not.k,,drop=F])
  s <- score.features(xp,type="odds")
  f <- colnames(x)[which.max(s)]
  #print(subtable(xp,1,f))
  f
}

discriminate.from.sibling <- function(g,i=1,x) {
  leaf1 <- names(g)[leaves(g,i)]
  pa = edges.to(g,i)
  if(length(pa) == 0) return("")
  sibling = setdiff(from(g,pa),i)
  if(length(sibling) == 0) return("")
  leaf2 <- names(g)[leaves(g,sibling)]
  xp <- array(0,c(2,ncol(x)),list(1:2,colnames(x)))
  xp[1,] = colSums(x[leaf1,,drop=F])
  xp[2,] = colSums(x[leaf2,,drop=F])
  xp <- xp + 1e-10*rep.col(rowSums(xp),ncol(xp))
  s = odds.ratios(xp)[1,]
  f <- colnames(x)[which.max(s)]
  if(T) {
    cat(f," ",xp[,f],"\n")
    ord = rev(order(s))[1:5]
    for(i in 1:length(ord)) {
      j = colnames(x)[ord[i]]
      cat(" ",j," ",xp[,j]," ",s[j],"\n")
    }
  }
  f
}

join.scores <- function(xp,verbose=F,...) {
  # xp is matrix with three rows: child1,child2,other
  #xp <- xp + 1e-3*rep.col(rowSums(xp),ncol(xp))
  #xp <- xp + 1e-2*rep.col(rowSums(xp),ncol(xp))
  xp <- xp + 1e-2
  if(T) {
    # unlike score.words, this measure is asymmetric
    # the word must favor the node
    #s12 <- abs(odds.ratios(xp[1:2,])[1,])
    #s12.se <- odds.ratios.se(xp[1:2,])[1,]
    s13 <- odds.ratios(xp[c(1,3),])[1,]
    s23 <- odds.ratios(xp[2:3,])[1,]
    s13.se <- odds.ratios.se(xp[c(1,3),])[1,]
    s23.se <- odds.ratios.se(xp[2:3,])[1,]
    if(T) {
      # hedging
      a = 1
      s13 = s13 - s13.se*a
      s23 = s23 - s23.se*a
    }
    if(F) {
      s13 = s13/s13.se
      s23 = s23/s23.se
    }
    #s <- pmin(s12/s13,s12/s23)
    #s = s12 - pmin(s13,s23)
    #s = s12+s12.se - pmin(s13+s13.se,s23+s23.se)
  } else {
    type = "a.info"
    a = 0
    #print(subtable(xp[c(1,3),],1,"the"))
    s13 = score.words(xp[c(1,3),],type=type,a=a)
    s23 = score.words(xp[c(2,3),],type=type,a=a)
  }
  # hard min or softmin?
  if(F) s = pmin(s13,s23)
  else if(T) s = -log(exp(-s13)+exp(-s23))
  else {
    # rank combination
    s = rank(s13)+rank(s23)
  }
  if(verbose) {
    # debugging
    f <- colnames(xp)[which.max(s)]
    cat(f," ",xp[,f],"\n")
    ord = rev(order(s))[1:5]
    for(i in 1:length(ord)) {
      j = colnames(xp)[ord[i]]
      cat(" ",j," ",xp[,j]," ",s[j],s13[j],s23[j],"\n")
    }
  }
  s
}

discriminate.node <- function(g,i=1,x,...) {
  child <- from(g,i)
  if(length(child) < 2) return("")
  leaf1 <- names(g)[leaves(g,child[1])]
  leaf2 <- names(g)[leaves(g,child[2])]
  other <- setdiff(rownames(x),c(leaf1,leaf2))
  if(length(other) == 0) return("")
  xp <- array(0,c(3,ncol(x)),list(1:3,colnames(x)))
  xp[1,] <- colSums(x[leaf1,,drop=F])
  xp[2,] <- colSums(x[leaf2,,drop=F])
  xp[3,] <- colSums(x[other,,drop=F])
  if(F) {
    pa = edges.to(g,i)
    if(length(pa) == 0) return("")
    sibling = setdiff(from(g,pa),i)
    if(length(sibling) == 0) return("")
    sib.leaves = names(g)[leaves(g,sibling)]
    xp.sib = colSums(x[sib.leaves,,drop=F])
    xp[3,] = xp[3,]/1 + xp.sib
  }
  s = join.scores(xp,...)
  colnames(x)[which.max(s)]
}

plot.join <- function(hc,doc,cex=0.8,...) {
  g <- as.graph.hclust(hc)
  labels <- names(g)
  for(i in 1:length(g)) {
    if(!is.leaf(g,i)) {
      labels[i] <- discriminate.node(g,i,doc,...)
    }
  }
  if(T) {
    n = length(labels)
    col = rep(1,n)
    for(i in 1:n) {
      lab = strsplit(labels[i],"\\.")[[1]][1]
      res = try(col2rgb(lab),silent=T)
      if(!inherits(res,"try-error")) col[i] = lab
    }
    plot.graph(g,labels=labels,orient=180,srt=0,cex=cex,col=col,...)
  }
  else plot.graph(g,labels=labels,orient=180,srt=0,cex=cex,...)
}
# functions for response tables
# Tom Minka 11/8/01

# given an aov fit, plot the effects for each predictor
# as well as a boxplot of the residuals
# you can do an F-test visually by comparing the spread of the effects
# with the spread of the residuals
effects.plot <- function(object,se=F) {
  mt <- model.tables(object,se=se)
  vars <- names(mt$tables)
  nvar <- length(vars)

  opar <- par(mar=c(2.5,4,0,0.1))
  on.exit(par(opar))
  plot.new()
  ylim <- c(min(sapply(mt$tables,min)), max(sapply(mt$tables,max)))
  plot.window(xlim=c(0.5,nvar+1+0.5),ylim=ylim)
  axis(1,1:(nvar+1), labels=c(vars,"residuals"))
  axis(2)
  title(ylab=response.var(object))

  if(se) {
    if(!is.numeric(se)) se <- 1.96
    p <- (1-2*pnorm(-se))*100
    cat("Displaying ",format(p,digits=2),"% confidence intervals\n",sep="")
  }
  for(k in 1:nvar) {
    eff <- mt$tables[[k]]
    text(k, eff, names(eff))
    if(se) {
      jitter <- runif(length(eff))*0.1
      arrows(k+jitter, eff-se*mt$se[[k]], k+jitter, eff+se*mt$se[[k]],
             code = 3, col = "green", angle = 75, length = .1)
    }
  }

  res <- residuals(object)
  x <- model.frame(object)
  res <- tapply(res, x[vars], mean)
  boxplot(res, at=nvar+1, add=T)
}

#############################################################################

as.rtable <- function(a,resp="Response",check.names=T) {
  if(check.names) resp = make.names(resp)
  nam = names(dimnames(a))
  if(is.null(nam)) stop("Dimensions must have names")
  fmla <- paste(resp,"~",
                paste(make.names(nam),collapse="+"))
  attr(a,"terms") <- terms(formula(fmla))
  #attr(a,"ordered") <- dimOrdered(a)
  class(a) <- c("rtable","table")
  a
}

print.rtable <- function(rt,...) {
  cat(response.var(rt),"\n")
  # strip attributes and print as a matrix
  attributes(rt) <- attributes(rt)[c("dim","dimnames")]
  print(rt)
}

terms.rtable <- function(rt) {
  attr(rt,"terms")
}

as.data.frame.rtable <- function(rt,...) {
  df <- as.data.frame.table(rt,...)
  # preserve ordering of factors
  i = names(which(dimOrdered(rt)))
  #if(length(i) > 0) df[i] = apply.df(df[i],as.ordered)
  len <- length(names(df))
  names(df)[len] <- response.var(rt)
  df
}
as.matrix.rtable <- function(rt,...) {
  x = as.matrix(unclass(rt))
  if(length(dim(rt)) == 1) {
    # why does as.matrix drop the name??
    # as.matrix doesn't understand this data type
    names(dimnames(x))[[1]] = names(dimnames(rt))
    t(x)
  } else x
}

row.probs <- function(x,smooth=0.5,se=F) {
  # converts ctable to rtable
  x = x + smooth
  col.sum = rep.row(colSums(x),nrow(x))
  nam = names(dimnames(x))
  if(is.null(nam) || nam[1] == "") nam = "row"
  p = as.rtable(x/col.sum,paste("Probability of",nam[1]))
  if(se != F) {
    if(se == T) se = 1.64
    p.se = sqrt(p*(1-p)/col.sum) * se
    list(p=p,se=p.se)
  } else p
}
log.counts <- function(x,smooth=0.5,se=F) {
  x = x+smooth
  r = as.rtable(log(x),"log(count)")
  if(se != F) {
    if(se == T) se = 1.64
    se = 1/sqrt(x) * se
    list(r=r,se=se)
  } else r
}
logit.rtable = function(x,smooth=0.5,se=F) {
  resp = response.var(x)
  pred = predictor.vars(x)
  lev = levels(x[[resp]])
  i = (x[[resp]]==lev[1])
  y1 = table(x[i,pred])+smooth
  y2 = table(x[!i,pred])+smooth
  r = log(y2/y1)
  r = as.rtable(r,paste("Logit of",resp))
  if(se != F) {
    if(se == T) se = 1.64
    se = (1/sqrt(y1)+1/sqrt(y2)) * se
    list(r=r,se=se)
  } else r
}
pclass.rtable = function(x,pclass=2,smooth=0.5,se=F) {
  resp = response.var(x)
  pred = predictor.vars(x)
  lev = levels(x[[resp]])
  i = (x[[resp]]==lev[pclass])
  y1 = table(x[i,pred])+smooth
  y2 = table(x[!i,pred])+smooth
  n = y1+y2
  p = y1/(y1+y2)
  p = as.rtable(p,paste("Probability of",resp,"=",lev[pclass]))
  if(se != F) {
    if(se == T) se = 1.64
    se = sqrt(p*(1-p)/n) * se
    list(p=p,se=se)
  } else p
}



aov.rtable <- function(rt,med=F) {
  frame <- as.data.frame.rtable(rt)
  if(med) {
    p <- medpolish(rt,trace.iter=F)
    return(structure(terms=terms(rt), model=frame, class="aov"))
  }
  aov(terms(rt), frame)
}
aov.residual <- function(rt,...) {
  rt - rtable(aov.rtable(rt,...))
}

# returns a table of responses, stratified according to the terms in fmla
# and aggregated according to fun.
rtable.terms <- function(fmla, x, fun = mean) {
  resp <- response.var(fmla)
  pred <- predictor.vars(fmla)
  rt <- tapply(x[[resp]], x[pred], fun)
  if(F && (length(dim(rt)) == 1)) {
    # force it to be a 1-row table
    dim(rt) = c(1,length(rt))
    dimnames(rt) = named.list(NULL,levels(x[[pred]]),names.=c("",pred))
  }
  class(rt) <- "rtable"
  attr(rt,"terms") <- fmla
  dimOrdered(rt) <- sapply(x[pred],is.ordered)
  rt
}
# model.frame is a data.frame with a "terms" attribute describing a model
rtable.data.frame <- function(object, ...) {
  rtable(terms(object), object, ...)
}
rtable.aov <- function(object, ...) {
  x <- model.frame(object)
  resp <- response.var(x)
  pred <- predictor.vars(x)
  y <- expand.grid(lapply(x[pred],levels))
  y[[resp]] <- predict(object,y)
  y <- model.frame(terms(x),y)
  rtable(y, ...)
}
# legal input:
# rtable(y~f1+f2,x)
# rtable(y~.,x)
# rtable(x$y ~ x$f1 + x$f2)
rtable.formula <- function(formula, data, ...) {
  # convert "|" into "+"
  rhs = formula[[3]]
  if(is.call(rhs) && (deparse(rhs[[1]]) == "|")) {
    rhs[[1]] = as.name("+")
    formula[[3]] = rhs
    model = model.frame.default(formula,data)
  } else {
    expr <- match.call(expand = F)
    expr$... <- NULL
    expr[[1]] <- as.name("model.frame.default")
    model <- eval(expr, parent.frame())
  }
  return(rtable(model, ...))
}
rtable <- function(object, ...) UseMethod("rtable")

# plots rows of x as traces
# similar to parallel.plot
# replacement for dotchart
# y is a named vector or matrix
linechart <- function(y,se=NULL,xlab=NULL,ylab,effects=F,med=F,
                      xscale=c("equal","linear","none"),...) {
  if(is.list(y) && missing(se)) {
    # y is list(y,se)
    se = y[[2]]
    y = y[[1]]
  }
  if(is.vector(y)) {
    y = t(as.matrix(y))
    if(!is.null(se)) se = t(as.matrix(se))
  }
  if(!is.matrix(y)) {
    y = as.matrix(y)
    if(!is.null(se)) se = as.matrix(se)
  }
  #if(missing(ylab)) ylab <- deparse(substitute(y))
  row.var <- names(dimnames(y))[1]
  col.var <- names(dimnames(y))[2]
  if(is.null(xlab)) if(!is.null(col.var) && !is.na(col.var)) xlab = col.var
  if(F) {
    # remove all-NA rows and columns
    y = y[!apply(y,1,function(z) all(is.na(z))),]
    y = y[,!apply(y,2,function(z) all(is.na(z)))]
  }

  xscale <- match.arg(xscale)
  if(effects) {
    if(missing(ylab)) ylab <- paste(row.var, "effect on", response.var(y))
    # choose x to make the profiles linear
    # fails if there are missing values
    if(med) {
      library(eda)
      fit <- medpolish(y,trace.iter=F)
      col.centers <- fit$col + fit$overall
      resid <- fit$residuals
    } else {
      col.centers <- colMeans(y,na.rm=T)
      row.effects <- rowMeans(y,na.rm=T) - mean(as.vector(y))
      resid <- y - outer(row.effects, col.centers, "+")
    }
    # subtract column centers
    y <- y - rep.row(col.centers,nrow(y))
  } else {
    if(missing(ylab)) ylab <- response.var(y)
    # sort columns by center
    #v = col.centers
    if(length(dim(y)) == 2)
      resid <- y - rep.col(rowMeans(y),ncol(y))
    else
      resid = y
  }
  if(any(is.na(resid))) {
    if(xscale != "none" && !dimOrdered(y)[2])
      cat("Warning: NAs prevent automatic ordering\n")
    xscale = "none"
  }
  if(xscale == "none") {
    x = 1:ncol(y)
  } else {
    s <- svd(resid)
    x <- s$v[,1]
  }
  # x is the desired x positions
  # remove flip ambiguity
  # choose the ordering which has maximum inner product with 1:ncol
  if(sum((1:length(x))*x) < sum((length(x):1)*x)) x = -x
  i = rank.stable(x)
  in.order = (!dimOrdered(y)[2] || all(i == 1:length(i)))
  if(in.order) {
    x <- switch(xscale,linear=x,equal=i,none=i)
  } else {
    x <- 1:ncol(y)
  }
  i <- order(x)
  y <- reorder.cols(y,i)
  x <- x[i]
  if(!is.null(se)) {
    se = reorder.cols(se,i)
    # kludge for call below
    se = t(se)
  }
  names(x) = colnames(y)
  labeled.curves(x,t(y),se,xlab=xlab,ylab=ylab,...)
}

labeled.curves <- function(x,y,se=NULL,labels,xlab=NULL,ylab,type="o",
                           group,
                           color.palette=default.colors(6),col,
                           lty.palette=1:6,lty,lwd=2,jitter=0,
                           legend.=NULL,move=T,
                           cex=par("cex"),horizontal=F,
                           mar,bg=par("bg"),ylim,main="",
                           xaxt=par("xaxt"),yaxt=par("yaxt"),...) {
  # must include (xaxt,yaxt) because they are swapped for horizontal
  # x is vector
  # y is a list or matrix of columns
  if(horizontal && missing(legend.)) legend. = 1
  if(horizontal) x = -x
  nam = names(dimnames(y))
  if(is.null(names(x))) names(x) = rownames(y)
  if(is.matrix(y)) y = as.data.frame.matrix(y)
  n = length(y)
  #if(is.data.frame(y)) y = as.list(y)
  if(is.null(xlab)) {
    if(!is.null(nam)) xlab = nam[1]
    else xlab = ""
  }
  if(missing(ylab)) {
    if(!is.null(nam)) ylab = nam[2]
    else ylab = ""
  }
  if(missing(mar)) {
    mar = if(horizontal)
      auto.mar(main=main,xlab=ylab,ylab=xlab,xaxt=yaxt,yaxt=xaxt,...)
    else
      auto.mar(main=main,xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
  }
  opar <- par(mar=mar,bg=bg)
  on.exit(par(opar))
  plot.new()
  xlim <- range(x)
  if(missing(labels)) labels = names(y)
  if(!is.null(labels) && is.null(legend.)) {
    # make space for labels
    xspc = 0.1
    w <- (max(strwidth(labels,units="inches",cex=cex))+xspc)/par("pin")[1]
    xlim[2] <- xlim[2] + diff(xlim)*w/(1-w)
  }
  if(!is.null(se)) {
    if(is.matrix(se)) se = as.data.frame.matrix(se)
    if(any(dim(se) != dim(y))) stop("se not same size as y")
    if(missing(ylim)) {
      vec.y = unlist(y)
      vec.se = unlist(se)
      # need y in case se is NA
      ylim = range(c(vec.y,vec.y-vec.se,vec.y+vec.se),finite=T)
    }
    # make room for arrow heads
    w = 0.2/par("pin")[if(horizontal) 2 else 1]
    dw = diff(xlim)*w/(1-w)
    xlim[1] = xlim[1] - dw/2
    xlim[2] = xlim[2] + dw/2
  } else {
    if(missing(ylim)) ylim = range(unlist(y),finite=T)
  }
  if(horizontal) {
    plot.window(xlim=ylim,ylim=xlim)
    if(is.null(names(x))) {
      axis(2,yaxt=xaxt,...)
      grid(nx=0)
    } else {
      #cat("columns are",names(x),"\n")
      axis(2,x, labels=names(x),yaxt=xaxt,...)
      # grid
      abline(h=seq(1,length(x)),lty="dotted",col="lightgray")
    }
    axis(1,xaxt=yaxt,...)
  }
  else {
    plot.window(xlim=xlim,ylim=ylim)
    if(is.null(names(x))) {
      axis(1,xaxt=xaxt,...)
      grid(ny=0)
    } else {
      #cat("columns are",names(x),"\n")
      axis(1,x, labels=names(x),xaxt=xaxt,...)
      # grid
      abline(v=x,lty="dotted",col="lightgray")
    }
    axis(2,yaxt=yaxt,...)
  }
  box()
  title(main=main)
  if(horizontal) title.las(ylab=xlab,xlab=ylab,...)
  else title.las(xlab=xlab,ylab=ylab,...)

  if(missing(group)) group = 1:n
  if(missing(col)) {
    if(mode(color.palette) == "function")
      color.palette <- color.palette(n)
    col <- color.palette[((group-1) %% length(color.palette))+1]
  }
  if(missing(lty)) {
    lty <- lty.palette[(((group-1) %/% length(color.palette)) %% length(lty.palette))+1]
  }

  jitter = xinch(jitter)
  jitter = jitter*((1:n)-0.5-n/2)
  #if(length(x) > 1) jitter = jitter*min(diff(x))
  if(!is.null(labels) && is.null(legend.)) {
    xspc = xinch(xspc)
    txt.x = txt.y = c()
  }
  for(i in 1:n) {
    if(horizontal) {
      lines(y[[i]],x,col=col[i],lty=lty[i],type=type,lwd=lwd, ...)
      if(!is.null(se)) {
        # arrows do not have lty
        arrows(y[[i]]-se[[i]],x+jitter[i],y[[i]]+se[[i]],x+jitter[i],
               len=0.1,angle=75,code=3,col=col[i],lty=lty[i],lwd=lwd,xpd=T)
      }
    } else {
      lines(x, y[[i]],col=col[i],lty=lty[i],type=type,lwd=lwd, ...)
      if(lty[i] == 7 && type != "p") {
        lines(x, y[[i]],col=8,lty=lty[i],type="l",lwd=1, ...)
      }
      if(!is.null(se)) {
        # arrows do not have lty
        arrows(x+jitter[i],y[[i]]-se[[i]],x+jitter[i],y[[i]]+se[[i]],
               len=0.1,angle=75,code=3,col=col[i],lwd=lwd,xpd=T)
      }
    }
    if(!is.null(labels) && is.null(legend.)) {
      # find the rightmost non-NA column in this row
      j <- rev(which(!is.na(y[[i]])))[1]
      txt.x[i] = x[j]+xspc
      txt.y[i] = y[[i]][j]
    }
  }
  if(!is.null(labels) && is.null(legend.)) {
    # move text vertically to avoid collisions
    h = strheight(labels,units="inch",cex=cex)
    names(txt.y) = labels
    if(move) txt.y = yinch(move.collisions(inchy(txt.y),h))
    text(txt.x, txt.y, labels, col=col, adj=0, cex=cex, ...)
  }
  if(!is.null(labels) && !is.null(legend.)) {
    #auto.legend(labels,col=col,lty=lty,text.col=col)
    my.legend(legend.,labels,lty=lty,col=col,lwd=lwd)
  }
}

my.legend <- function(pos=c(1,1),...) {
  if(length(pos) == 1) pos = c(pos,pos)
  xlim = par("usr")[1:2]
  ylim = par("usr")[3:4]
  x = xlim[1] + diff(xlim)*pos[1]
  y = ylim[1] + diff(ylim)*pos[2]
  legend(x,y,...,xjust=pos[1],yjust=pos[2])
}

title.las <- function(xlab=NULL,ylab=NULL,las=par("las"),...) {
  if(las == 0) {
    title(xlab=xlab,ylab=ylab,...)
  } else if(las == 1) {
    title(xlab=xlab,...)
    title(ylab=ylab,line=5,...)
  } else if(las == 2) {
    title(xlab=xlab,line=5.5,...)
    title(ylab=ylab,line=5,...)
  } else if(las == 3) {
    title(ylab=ylab,...)
    title(xlab=xlab,line=5.5,...)
  }
}
# Extensions and uses of loess

smooth <- function(x,y) {
  library(modreg)
  if(length(unique(x)) >= 4) {
    fit = smooth.spline(x,y)
    predict(fit)
  } else {
    lowess(x,y)
  }
}
plot.loess <- function(object,xlim,col=2,lwd=2,res=200,asp=NULL,add=F,...) {
  x <- model.frame(object)
  pred <- predictor.vars(x)
  resp <- response.var(x)
  if(missing(xlim)) xlim <- range(x[[pred]])
  if(add) {
    xlim[1] = max(xlim[1], par("usr")[1])
    xlim[2] = min(xlim[2], par("usr")[2])
  }
  xt <- seq(xlim[1],xlim[2],len=res)
  y <- predict(object,xt)
  real.asp = if(identical(asp,"auto")) auto.aspect(xt,y) else asp
  if(!add)
    plot(x[[pred]],x[[resp]],xlab=pred,ylab=resp,asp=real.asp,...)
  lines(xt,y,col=col,lwd=lwd,...)
}

if(T) {
smooth <- function(...) {
  library(modreg)
  expr <- match.call()
  expr[[1]] <- as.name("loess")
  object <- eval(expr, parent.frame())
  #object = loess(...)
  class(object) = c("smooth","loess")
  object
}
print.smooth = get("print.loess",environment(loess))
plot.smooth = plot.loess
predict.smooth = get("predict.loess",environment(loess))
} else {
smooth <- function(fmla,data=parent.frame(),span=2/3,...) {
  given = given.vars(fmla)
  if(!is.null(given)) {
    fmla2 = remove.given(fmla)
    g = data[[given]]
    if(!is.factor(g)) stop("given must be a factor")
    children = list()
    for(lev in levels(g)) {
      children[[lev]] = smooth(fmla2,data[g == lev,],span,...)
    }
    object = list(model=data,terms=terms(fmla),given=given,children=children)
    class(object) = "multi.smooth"
    return(object)
  }
  x = model.frame(fmla,data)
  pred = predictor.vars(x)[1]
  resp = response.var(x)
  xy = lowess(x[[pred]],x[[resp]],f=span,...)
  names(xy$x) = rownames(x)
  names(xy$y) = rownames(x)
  res = x[[resp]] - xy$y
  orig.x = as.named.vector(data[pred])
  object = c(xy,list(model=x,terms=terms(x),residuals=res,orig.x=orig.x))
  class(object) = "smooth"
  object
}
print.multi.smooth <- function(object) {
  cat("Multi-smooth object: ")
  print(formula(object))
  cat(object$given,"values:", names(object$children), "\n")
  #print(object$children)
}
print.smooth <- function(object) {
  cat("Smooth object: ")
  print(formula(object))
}
predict.multi.smooth <- function(object,frame) {
  if(missing(frame)) frame = object$model
  g = frame[,object$given]
  r = named.numeric(rownames(frame))
  for(lev in levels(g)) {
    i = (g == lev)
    r[i] = predict.smooth(object$children[[lev]], frame[i,])
  }
  r
}
predict.smooth <- function(object,x=object$orig.x) {
  if(is.data.frame(x)) x = x[[predictor.vars(object)[1]]]
  y = approx(object$x,object$y,x)$y
  if(is.null(names(y))) names(y) = names(x)
  y
}
}
model.frame.smooth = model.frame.loess

lowess2 <- function(x,y,f=2/3,res=200,...) {
  if(length(unique(x)) < length(x)/4) return(lowess(x,y,f=f))
  object = smooth(y~x,data.frame(x,y),span=f,...)
  xlim <- range(x)
  xo <- seq(xlim[1],xlim[2],len=res)
  yo <- predict(object,xo)
  list(x=xo,y=yo)
}

loprob <- function(x,y,degree=3,...) {
  # simulates lowess for binary response
  # returns list(x,y)
  lev <- levels(y)
  if(length(lev) != 2) stop("y must have two levels")
  if(!is.factor(x)) {
    from <- min(x)
    to <- max(x)
  }
  if(T) {
    my.lowess = function(x,y,span=2/3,...) lowess(x,y,f=span,iter=0,...)
    my.lowess(x,y,...)
  } else if(T) {
    # smoothing spline
    library(modreg)
    y <- as.numeric(y)-1
    r <- seq(from,to,length=50)
    if(degree == 3 && length(unique(x)) >= 4) {
      fit = smooth.spline(x,y)
      p <- predict(fit,r)$y
    } else {
      fit <- loess(y~x,degree=degree,...)
      p <- predict(fit,r)
    }
    list(x=r,y=p+1)
  } else if(F) {
    # density estimate - too bumpy
    densf <- function(x) {
      if(is.factor(x)) {
        list(x=levels(x),y=tabulate(x)/length(x))
      } else {
        density(x,from=from,to=to,adjust=1.5)
      }
    }
    dens <- lapply(split(x,y),densf)
    p <- sum(y==levels(y)[2])/length(y)
    p <- p*dens[[2]]$y/((1-p)*dens[[1]]$y + p*dens[[2]]$y)
    list(x=dens[[1]]$x,y=p+1)
  } else if(F) {
    # hill density estimate - too noisy
    xr <- seq(from,to,length=50)
    densf <- function(x) hill.density(x,xr)
    dens <- lapply(split(x,y),densf)
    p <- dens[[2]]/(dens[[1]] + dens[[2]])
    list(x=xr,y=p+1)
  } else {
    # gam smooth - too slow
    library(mgcv)
    y <- as.numeric(y)-1
    fit <- gam(y~s(x),family=binomial)
    r <- seq(from,to,length=50)
    p <- predict(fit,data.frame(x=r),type="response")
    list(x=r,y=p+1)
  }
}

##############################################################################
# 3D plots

# x,y,r are in user coordinates
draw.circle <- function(x,y,r,res=64,fill=F,...) {
  theta <- seq(0,2*pi,len=res)
  if(length(r) == 1) r <- rep(r,length(x))
  xt <- cos(theta)
  yt <- sin(theta)
  for(i in 1:length(x)) {
    if(fill)
      polygon(x[i]+xt*r[i],y[i]+yt*r[i],...)
    else
      lines(x[i]+xt*r[i],y[i]+yt*r[i],...)
  }
}
draw.vane <- function(x,y,a,r) {
  xlim <- par("usr")[1:2]
  ylim <- par("usr")[3:4]
  xd <- r*cos(a)*diff(xlim)/par("pin")[1]
  yd <- r*sin(a)*diff(ylim)/par("pin")[2]
  segments(x-xd/2,y-yd/2,x+xd/2,y+yd/2)
}
draw.brick <- function(x,y,w,h) {
  xlim <- par("usr")[1:2]
  ylim <- par("usr")[3:4]
  w <- w*diff(xlim)/par("pin")[1]
  h <- h*diff(ylim)/par("pin")[2]
  # these drawing directions are important to get completed corners
  segments(x-w/2,y+h/2,x+w/2,y+h/2)
  segments(x+w/2,y-h/2,x-w/2,y-h/2)
  segments(x-w/2,y-h/2,x-w/2,y+h/2)
  segments(x+w/2,y+h/2,x+w/2,y-h/2)
}
bubbles <- function(object, ...) UseMethod("bubbles")
bubbles.formula <- function(formula,data=parent.frame(),...,xlab,ylab,zlab) {
  x <- model.frame.default(formula,data,na.action=na.pass)
  resp = response.var(x)
  pred = predictor.vars(x)
  if(missing(xlab)) xlab = pred[1]
  if(missing(ylab)) ylab = pred[2]
  if(missing(zlab)) zlab = resp
  bubbles.default(x[,pred[1]],x[,pred[2]],x[,resp],
                  xlab=xlab,ylab=ylab,zlab=zlab,...)
}
bubbles.default <- function(x,y,z,data=parent.frame(),
                    zlim,zex=1,xlab,ylab,zlab,...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(missing(zlab)) zlab <- deparse(substitute(z))
  if(is.null(data)) data <- parent.frame()
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  z <- eval(substitute(z),data)
  opar <- par(mar=auto.mar(main=zlab,xlab=xlab,ylab=ylab))
  on.exit(par(opar))

  plot(x,y,type="n",xlab=xlab,ylab=ylab,main=zlab,...)
  if(!missing(zlim)) z <- z[(z >= zlim[1]) & (z <= zlim[2])]
  z = rank(z)
  z <- scale.range(z)
  # avoid zero-size bubbles
  z <- (z + 0.025)/1.025
  if(T){
    # encode z as area
    z <- sqrt(z)
  }
  z <- z*zex
  for(i in 1:length(x)) {
    points(x[i],y[i],pch=1,cex=z[i]*3,...)
  }
}
vanes <- function(x,y,z,data=parent.frame(),w=NULL,
                  zlim,wlim,wex=0.1,xlab,ylab,zlab,...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(missing(zlab)) zlab <- deparse(substitute(z))
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  z <- eval(substitute(z),data)
  opar <- par(mar=auto.mar(main=zlab,xlab=xlab,ylab=ylab))
  on.exit(par(opar))

  plot(x,y,type="n",xlab=xlab,ylab=ylab,main=zlab,...)
  if(!missing(zlim)) z <- z[(z >= zlim[1]) & (z <= zlim[2])]
  if(!missing(wlim)) w <- w[(w >= wlim[1]) & (w <= wlim[2])]
  # map z to (0,pi/2)
  z = rank(z)
  z <- scale.range(z)*pi/2
  if(is.null(w)) w <- 1
  if(length(w) == 1) w <- rep(w,length(z))*wex
  else w <- (w - min(w,na.rm=T))/diff(range(w,na.rm=T))*wex
  draw.vane(x,y,z,w)
}
bricks <- function(x,y,w,data=parent.frame(),h=NULL,
                   wex=0.1,hex=0.1,xlab,ylab,wlab,...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(missing(wlab)) wlab <- deparse(substitute(w))
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  w <- eval(substitute(w),data)
  opar <- par(mar=auto.mar(main=zlab,xlab=xlab,ylab=ylab))
  on.exit(par(opar))

  plot(x,y,type="n",xlab=xlab,ylab=ylab,main=wlab,...)
  w <- scale.range(w)*wex
  if(is.null(h)) h <- 0.5
  if(length(h) == 1) h <- rep(h,length(x))*hex
  else h <- scale.range(h)*hex
  draw.brick(x,y,w,h)
}
plot3d <- function(object, ...) UseMethod("plot3d")
plot3d.formula <- function(formula,data=parent.frame(),...,xlab,ylab,zlab) {
  x <- model.frame.default(formula,data,na.action=na.pass)
  resp = response.var(x)
  pred = predictor.vars(x)
  if(missing(xlab)) xlab = pred[1]
  if(missing(ylab)) ylab = pred[2]
  if(missing(zlab)) zlab = resp
  plot3d.default(x[,pred[1]],x[,pred[2]],x[,resp],
                  xlab=xlab,ylab=ylab,zlab=zlab,...)
}
plot3d.default <- function(x,y,z,data=parent.frame(),xlab,ylab,zlab,angle=60,
                           bg=par("bg"),...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(missing(zlab)) zlab <- deparse(substitute(z))
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  z <- eval(substitute(z),data)
  library(scatterplot3d)
  opar = par(bg=bg)
  on.exit(par(opar))
  scatterplot3d(x,y,z,xlab=xlab,ylab=ylab,zlab=zlab,angle=angle,...)
}

formals(persp.default)$theta = 30
formals(persp.default)$shade = 0.65
surface = function(object,...) UseMethod("surface")
surface.loess <- function(object,res=20,xlim,ylim,clip=T,...) {
  if(length(res) == 1) res <- rep(res,2)
  x <- model.frame(object)
  resp <- response.var(object)
  pred <- predictor.vars(object)
  if(missing(xlim)) xlim <- range(x[[pred[1]]])
  x1 <- seq(xlim[1],xlim[2],length=res[1])
  if(missing(ylim)) ylim <- range(x[[pred[2]]])
  x2 <- seq(ylim[1],ylim[2],length=res[2])
  xt <- expand.grid(x1,x2)
  names(xt) <- pred[1:2]
  z <- predict(object,xt)
  if(identical(clip,T)) {
    i <- chull(x[pred])
    clip <- x[pred][i,]
  }
  if(!identical(clip,F)) {
    is.na(z[!in.polygon(clip,x[pred])]) <- T
  }
  xt[,resp] = z
  surface.rtable(rtable(xt),xlim=xlim,ylim=ylim,...)
}
surface.rtable <- function(zm,xlim,ylim,zlim,asp,xlab,ylab,zlab,...) {
  # clip specifies a polygon in (x,y) over which the surface is to be defined.
  # possible values are F (no clipping), T (clip to the convex hull of the
  # data), a list(x,y), or a matrix with two columns.
  z = as.vector(zm)
  x = as.numeric(rownames(zm))
  y = as.numeric(colnames(zm))
  if(!missing(zlim)) {
    z[z < zlim[1]] <- zlim[1]
    z[z > zlim[2]] <- zlim[2]
  }
  if(missing(xlim)) xlim = range(x)
  if(missing(ylim)) ylim = range(y)
  if(!missing(asp)) {
    # should be included in surface()
    # fake asp for routines that don't support it
    s = diff(xlim)/diff(ylim)/asp
    if(s >= 1)
      ylim = (ylim-mean(ylim))*s + mean(ylim)
    else
      xlim = (xlim-mean(xlim))/s + mean(xlim)
  }
  pred = predictor.vars(zm)
  if(missing(xlab)) xlab <- pred[1]
  if(missing(ylab)) ylab <- pred[2]
  if(missing(zlab)) zlab <- response.var(zm)
  opar <- par(mar=c(1,0,0,0))
  on.exit(par(opar))
  persp(x,y,zm,
        xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,zlab=zlab,...)
}

image.loess <- function(object,res=40,xlim,ylim,clip=T,...) {
  if(length(res) == 1) res <- rep(res,2)
  x <- model.frame(object)
  resp <- response.var(object)
  pred <- predictor.vars(object)
  if(missing(xlim)) xlim <- range(x[[pred[1]]])
  x1 <- seq(xlim[1],xlim[2],length=res[1])
  if(missing(ylim)) ylim <- range(x[[pred[2]]])
  x2 <- seq(ylim[1],ylim[2],length=res[2])
  xt <- expand.grid(x1,x2)
  names(xt) <- pred[1:2]
  z <- predict(object,xt)
  if(identical(clip,T)) {
    i <- chull(x[pred])
    clip <- x[pred][i,]
  }
  xt[,resp] = z
  image.data.frame(xt,clip=clip,...)
}
image.formula <- function(formula,data=parent.frame(),...) {
  x <- model.frame.default(formula,data,na.action=na.pass)
  image.data.frame(x,...)
}
image.data.frame <- function(data,res=200,
                         clip=F,xlim,ylim,xaxs="r",yaxs="r",mar=par("mar"),
               xlab,ylab,zlab,main=zlab,nlevels=128,breaks,col,bg=grey(0.5),
               color.palette=YR.colors,asp=NULL,axes=T,...) {
  # asp and axes are there to be overridden by formals
  library(akima)
  pred = predictor.vars(data)
  resp = response.var(data)
  x = data[,pred[1]]
  y = data[,pred[2]]
  z = data[,resp]
  if(missing(xlab)) xlab <- pred[1]
  if(missing(ylab)) ylab <- pred[2]
  if(missing(zlab)) zlab <- resp

  if(length(res) == 1) res = rep(res,2)
  if(missing(xlim)) xlim <- range(x)
  xo <- seq(xlim[1],xlim[2],len=res[1])
  if(missing(ylim)) ylim <- range(y)
  yo <- seq(ylim[1],ylim[2],len=res[2])
  p <- interp(x,y,z,xo,yo)
  #p <- interp.new(x,y,z,xo,yo)
  if(!identical(clip,F)) {
    i <- !in.polygon(clip,expand.grid(xo,yo))
    dim(i) <- c(length(xo),length(yo))
    p$z[i] <- NA
  }
  if(missing(breaks)) {
    if(missing(col)) {
      col <- color.palette(nlevels)
    } else {
      nlevels = length(col)
    }
    breaks <- break.quantile(z,nlevels,pretty=T)
  } else {
    if(missing(col)) {
      nlevels = length(breaks)-1
      col = color.palette(nlevels)
    }
    # clip at the breaks
    zlim = range(breaks)
    p$z[p$z < zlim[1]] <- zlim[1]
    p$z[p$z > zlim[2]] <- zlim[2]
  }
  opar <- par(bg=bg,mar=mar)
  on.exit(par(opar))
  image.default(p,xlab=xlab,ylab=ylab,main=main,breaks=breaks,col=col,xaxs=xaxs,yaxs=yaxs,xlim=xlim,ylim=ylim,asp=asp,axes=axes,...)
}

slices <- function(object,...) UseMethod("slices")
slices.formula <- function(fmla,data=parent.frame(),xlab,ylab,zlab,...) {
  given = given.vars(fmla)
  if(is.null(given)) {
    #stop("formula must have conditioning (|)")
    pred = predictor.vars(fmla)
    fmla = formula(paste(pred[2],"~",pred[1],"|",response.var(fmla)))
    given = given.vars(fmla)
  }
  fmla = remove.given(fmla)
  pred = predictor.vars(fmla)[1]
  resp = response.var(fmla)
  frame = model.frame.default(formula(paste(resp,"~",pred,"+",given)),data)
  given.vec = strsplit(given," \\+ ")[[1]]
  x <- frame[,pred]
  y <- frame[,resp]
  z <- frame[,given.vec]
  if(missing(xlab)) xlab <- pred
  if(missing(ylab)) ylab <- resp
  if(missing(zlab)) zlab <- given
  slices.default(x,y,z,xlab=xlab,ylab=ylab,zlab=zlab,...)
}
slices.default <- function(x,y,z,data=parent.frame(),nlevels=4,
                   layout,xlab,ylab,zlab,xlim,ylim,...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(missing(zlab)) zlab <- deparse(substitute(z))
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  z <- eval(substitute(z),data)

  # because of unique, this is different from equal.count
  b <- break.quantile(z,nlevels)
  # should be done by my.cut
  if(is.list(z)) {
    zf = list()
    for(i in 1:length(z)) {
      zf[[i]] = cut(z[[i]],b[[i]],include=T)
    }
    zf = do.call("interaction",zf)
  } else {
    zf <- cut(z,b,include=T)
  }

  if(missing(layout)) layout = auto.layout(length(levels(zf)))
  opar <- par(mfrow=layout)
  on.exit(par(opar))

  if(missing(xlim)) xlim = range(x)
  if(missing(ylim)) ylim = range(y)

  for(lev in levels(zf)) {
    i <- which(zf == lev)
    main = paste(zlab,"=",lev)
    # limits must be consistent
    predict.plot.default(x[i],y[i],xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,...)
  }
}

project.slices <- function(x,y=NULL,nlevels=4,aligned=T,type="mv",...) {
  sw = slicing(x,y,type=type,...)
  print(auto.signif(sw))
  z = project(x,sw)[,1]
  f = cut.quantile(z,nlevels)
  opar = par(mfrow=auto.layout(nlevels))
  on.exit(par(opar))
  if(aligned) {
    w = projection(x,y,k=2,given=sw,type=type,...)
    w = w[,-1]
    px = project(x,w)
    xlim = range(px[,1])
    ylim = range(px[,2])
  } else {
    xlim = ylim = NULL
  }
  for(lev in levels(f)) {
    i = which(f == lev)
    y.lev = y[i]
    if(!aligned) {
      x.lev = x[i,]
      w = projection(x.lev,y.lev,k=2,type=type,...)
      #w = pca(x.lev,k=2,...)
      #w = pca(x.lev,y.lev,k=2,...)
      px.lev = project(x.lev,w)
    } else {
      px.lev = px[i,]
    }
    main = lev
    if(!is.null(y.lev) && all(is.na(y.lev)))
      plot(px.lev,main=main,xlim=xlim,ylim=ylim,...)
    else {
      color.plot(px.lev,zlab=main,xlim=xlim,ylim=ylim,...)
      #color.plot(px,y.lev,...)
    }
    plot.axes(w,...)
  }
}

slices2 <- function(x,y,z,data=parent.frame(),nlevels=4,
                   digits=2,xlab,ylab,zlab,...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(missing(zlab)) zlab <- deparse(substitute(z))
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  z <- eval(substitute(z),data)

  opar <- par(mar=c(0,0,0,0))
  # lattice changes the mar
  library(lattice)
  par(opar)
  trellis.par.set("add.line",list(col=2,lwd=2))
  #trellis.settings$background$col <- 0
  #trellis.settings$plot.symbol$col <- 1

  df <- data.frame(x,y,z)
  df$z <- equal.count(z,nlevels,overlap=0)
  nam <- sapply(c(xlab,ylab,zlab),make.names)
  names(df) <- nam
  fmla <- formula(paste(nam[2],"~",nam[1],"|",nam[3]))
  xyplot(fmla, data=df,
         panel = function(x, y) {
           panel.grid(v=2); panel.xyplot(x, y);
           panel.loess(x,y)
         })
}

# extra routines for using regression trees
# requires: tree package, cluster.r
# Tom Minka 10/26/01

# plot tree partitions along with data colored by response
color.plot.tree <- function(tr,data,col=2,lwd=2,add=F,...) {
  if(missing(data)) data <- model.frame(tr)
  resp <- response.var(tr)
  # reduce to the variables used in tr
  pred <- levels(factor(tr$frame$var[tr$frame$var != "<leaf>"]))
  fmla <- formula(paste(resp,"~",paste(pred,collapse="+")))
  data <- model.frame(fmla,data)
  if(!add) color.plot.data.frame(data,...)
  # use the par options set by color.plot.data.frame
  partition.tree(tr,add=T,ordvars=pred,col=col,lwd=lwd,...)
}

# returns the parent index of every node
parents.tree <- function(tr) {
  i <- as.numeric(rownames(tr$frame))
  match(i %/% 2, i)
}
is.left.child <- function(tr,i) {
  if(missing(i)) i <- 1:nrow(tr$frame)
  if(is.numeric(i)) i <- rownames(tr$frame)[i]
  (as.numeric(i) %% 2 == 0)
}

# converts tr$where into a factor of node names
# (actually letters are more readable)
factor.tree <- function(tr,nodes=F) {
  f <- factor(rownames(tr$frame)[tr$where])
  if(!nodes) levels(f) <- toupper(letters[1:nlevels(f)])
  names(f) <- names(tr$where)
  f
}

# returns a vector of node indices under node
# the input is a node name, but output is node number
descendants.tree <- function(tr,node) {
  p <- parents.tree(tr)
  node <- match(node,rownames(tr$frame))
  i <- 1:length(p)
  r <- c()
  repeat {
    r <- c(r,which(i==node))
    i <- p[i]
    if(all(is.na(i) | (i<node))) break
  }
  r
}

# returns a vector of case numbers under node
cases.tree <- function(tr,node) {
  which(tr$where %in% descendants.tree(tr,node))
}

path.tree <- function(tr,node) {
  lab <- labels.tree(tr)
  p <- parents.tree(tr)
  node <- match(node,rownames(tr$frame))
  s <- ""
  repeat {
    if(nchar(s) == 0) s <- lab[node]
    else s <- paste(lab[node],s,sep="; ")
    node <- p[node]
    if(node == 1) break
  }
  s
}

best.size.tree <- function(tr, fold=10, misclass=F) {
  # wrapper for cv.tree
  # returns the tree with best size (actually best cost-complexity param)
  m <- model.frame(tr)
  n <- length(m[[1]])
  if(T) {
    # do the cv in a more conventional way
    rand = rand.partition(round.partition(rep(n,fold)/fold))
    if(misclass) cv <- cv.tree(tr,rand,FUN=prune.misclass,K=fold)
    else cv <- cv.tree(tr,rand,K=fold)
  } else {
    if(misclass) cv <- cv.tree(tr,FUN=prune.misclass,K=fold)
    else cv <- cv.tree(tr,K=fold)
  }
  plot(cv,type="o")
  title(paste(fold,"fold cross-validation"))
  cv$dev <- rev(cv$dev)
  cv$size <- rev(cv$size)
  i <- which.min(cv$dev)
  cat("best size is",cv$size[i],"leaves\n")
  tr = if(misclass) prune.misclass(tr,best=cv$size[i])
       else prune.tree(tr,best=cv$size[i])
  # necessary for model.frame.tree to work
  tr$model = m
  tr
}

#############################################################################
# bug fixes

my.tree <- function(formula,data,...) {
  tr = tree(formula,data,...)
  tr$model = model.frame(formula,data)
  tr
}

# from library/tree/R/tree
my.print.tree <-
    function(x, pretty = 0, spaces = 2, digits = getOption("digits")-3, ...)
{
    if(!inherits(x, "tree")) stop("Not legitimate tree")
    is.prob <- !is.null(ylevels <- attr(x, "ylevels"))
    # minka: omit deviance
    if(is.prob) cat("node), split, n, yval, (yprob)\n")
    else cat("node), split, n, yval\n")
    cat("      * denotes terminal node\n\n")
    frame <- x$frame
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32), collapse = "")
                                     #32 is the maximal depth
    if(length(node) > 1) {
        indent <- substring(indent, 1, spaces * seq(depth))
        indent <- paste(c("", indent[depth]), format(node), ")", sep = "")
    } else
    indent <- paste(format(node), ")", sep = "")
    if(is.prob) {
        yval <- paste(as.character(frame$yval), " (", sep = "")
        yprob <- format(frame$yprob, digits = digits)
        for(i in seq(ylevels)) yval <- paste(yval, yprob[, i])
        yval <- paste(yval, ")")
    } else
    yval <- format(signif(frame$yval, digits = digits))
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels.tree(x, pretty = pretty)
    # minka: omit deviance
    z <- paste(indent, z, format(round(frame$n,2)),
               yval, term)
    cat(z, sep = "\n")
    invisible(x)
}

# partition.tree must require cont predictors because split sets on a
# factor can vary
# minka: added digits as argument
my.partition.tree <- function(tree, label = "yval", add = FALSE, ordvars,
                           digits=2, force = T, main=NULL, ...)
{
    ptXlines <- function(x, v, xrange, xcoord = NULL, ycoord = NULL, tvar, i = 1)
    {
        if(v[i] == "<leaf>") {
            y1 <- (xrange[1] + xrange[3])/2
            y2 <- (xrange[2] + xrange[4])/2
            return(xcoord, ycoord = c(ycoord, y1, y2), i = i)
        }
        if(v[i] == tvar[1]) {
            xcoord <- c(xcoord, x[i], xrange[2], x[i], xrange[4])
            # minka: indent the partitions
            #xrange[c(2,4)] = xrange[c(2,4)] + c(1,-1)*inchy(0.05)
            xr <- xrange
            xr[3] <- x[i]
            ll2 <- Recall(x, v, xr, xcoord, ycoord, tvar, i + 1)
            xr <- xrange
            xr[1] <- x[i]
            return(Recall(x, v, xr, ll2$xcoord, ll2$ycoord, tvar, ll2$i + 1))
        } else if(v[i] == tvar[2]) {
            xcoord <- c(xcoord, xrange[1], x[i], xrange[3], x[i])
            # minka
            #xrange[c(1,3)] = xrange[c(1,3)] + c(1,-1)*inchx(0.05)
            xr <- xrange
            xr[4] <- x[i]
            ll2 <- Recall(x, v, xr, xcoord, ycoord, tvar, i + 1)
            xr <- xrange
            xr[2] <- x[i]
            return(Recall(x, v, xr, ll2$xcoord, ll2$ycoord, tvar, ll2$i + 1))
        }
        else stop("Wrong variable numbers in tree.")
    }
    if(inherits(tree, "singlenode")) {
      warning("Singlenode tree has no partitions")
      return(invisible(NULL))
    }
    if(!inherits(tree, "tree")) stop("Not legitimate tree")
    frame <- tree$frame
    leaves <- frame$var == "<leaf>"
    var <- unique(as.character(frame$var[!leaves]))
    if(length(var) > 2 || length(var) < 1)
        stop("Tree can only have one or two predictors")
    nlevels <- sapply(xlevels <- attr(tree, "xlevels"), length)
    if(any(nlevels[var] > 0))
        stop("Tree can only have continuous predictors")
    x <- rep(NA, length(leaves))
    x[!leaves] <- as.double(substring(frame$splits[!leaves, "cutleft"], 2, 100))
    m <- model.frame(tree)
    # minka: treat 1 var as if var 2 was yvar
    if(force && (length(var) == 1)) {
      var[2] <- response.var(tree)
    }
    if(length(var) == 1) {
        ## one x variable
        x <- sort(c(range(m[[var]]), x[!leaves]))
        if(is.null(attr(tree, "ylevels"))) y <- frame$yval[leaves]
        # minka: ??
        else y <- frame$yprob[leaves, 2] + 1
        # minka: must extend y to match x
        y <- c(y,y[length(y)])
        if(add) lines(x, y, type = "s", ...)
        else {
          # minka: use response.var
          resp <- response.var(m)
          xo <- m[[resp]]
          ylim = NULL
          #if(is.factor(xo)) ylim <- c(0,1) else ylim <- range(xo)
          plot(x, y, ylab = resp, xlab = var, type = "s", ylim = ylim,
               xaxs = "i", ...)
        }
        invisible(list(x = x, y = y))
    } else {
        ## two x variables
        if(!missing(ordvars)) {
            ind <- match(var, ordvars)
            if(any(is.na(ind))) stop("unmatched names in vars")
            var <- ordvars[sort(ind)]
        }
        lab <- frame$yval[leaves]
        resp = response.var(tree)
        if(is.null(frame$yprob)) lab <- format(signif(lab, digits))
        else if(match(label, attr(tree, "ylevels"), nomatch = 0)) {
          lab <- format(signif(frame$yprob[leaves, label], digits))
          # minka: useful title
          if(is.null(main)) main = paste("Pr(",resp," = ",label,")",sep="")
        }
        # minka: useful title
        if(is.null(main)) main = resp
        rx <- range(m[[var[1]]])
        rx <- rx + c(-0.025, 0.025) * diff(rx)
        rz <- range(as.numeric(m[[var[2]]]))
        rz <- rz + c(-0.025, 0.025) * diff(rz)
        xrange <- c(rx, rz)[c(1, 3, 2, 4)]
        xcoord <- NULL                  # x1lo, x2lo, x1hi, x2hi
        ycoord <- NULL                  # y1, y2
        xy <- ptXlines(x, frame$var, xrange, xcoord, ycoord, var)
        xx <- matrix(xy$xcoord, nrow = 4)
        yy <- matrix(xy$ycoord, nrow = 2)
        if(!add)
            plot(rx, rz, xlab = var[1], ylab = var[2], type = "n",
                 xaxs = "i", yaxs = "i", main = main, ...)
        # minka: use ...
        segments(xx[1,  ], xx[2,  ], xx[3,  ], xx[4,  ], ...)
        text(yy[1,  ], yy[2,  ], as.character(lab), ...)
    }
}

# minka: added newdata,rate options
my.misclass.tree <- function(tree, newdata, rate = F, detail = FALSE)
{
    if(!inherits(tree, "tree"))
        stop("Not legitimate tree")
    if(is.null(attr(tree, "ylevels")))
        stop("Misclassification error rate is appropriate for factor responses only")
    if(!missing(newdata)) {
      resp <- response.var(tree)
      r = sum(predict.tree(tree, newdata, type="class") != newdata[[resp]])
      return(if(rate) r/nrow(newdata) else r)
    }
    if(is.null(y <- tree$y))
        y <- model.extract(model.frame(tree), "response")
    if(is.null(wts <- tree$weights))
        wts <- model.weights(model.frame(tree))
    if(is.null(wts)) wts <- rep(1, length(y))
    frame <- tree$frame
    if(detail) {
        which <- descendants(as.numeric(row.names(frame)))
        tmp <- as.vector((which[, tree$where] *
                          outer(frame$yval, y, "!=")) %*% wts)
        names(tmp) <- row.names(tree$frame)
        r = tmp
    }
    else r = sum(wts*(frame$yval[tree$where] != y))
    if(rate) r/length(y) else r
}

# minka: added newdata,rate options
my.deviance.tree <- function(object, newdata, rate = F, detail = FALSE, ...)
{
    if(!inherits(object, "tree")) stop("Not legitimate tree")
    if(!missing(newdata)) {
      yvar <- response.var(object)
      truth <- as.numeric(newdata[[yvar]])
      p <- predict.tree(object, newdata)
      # correction used in prune.tree
      p[p == 0] <- 1e-3
      i <- (1:nrow(p)) + (truth-1)*nrow(p)
      dev = -2*sum(log(p[i]))
      return(if(rate) dev/length(truth) else dev)
    }
    frame <- object$frame
    if(detail) dev = frame$dev
    else dev = sum(frame$dev[frame$var == "<leaf>"])
    if(rate) dev/sum(frame$var == "<leaf>") else dev
}

my.cv.tree <- function(object, rand, FUN = prune.tree, K = 10, ...)
{
    if(!inherits(object, "tree")) stop("Not legitimate tree")
    m <- model.frame(object)
    extras <- match.call(expand.dots = FALSE)$...
    FUN <- deparse(substitute(FUN))
    init <- do.call(FUN, c(list(object), extras))
    if(missing(rand)) rand <- sample(K, length(m[[1]]), replace = TRUE)
    cvdev <- 0
    for(i in unique(rand)) {
        # minka: must build a tree as big as original
      # can use object$call to see what controls were used
        tlearn <- tree(model = m[rand != i,  , drop = FALSE], mindev=-1)
        plearn <- do.call(FUN, c(list(tlearn, newdata =
                                      m[rand ==i, , drop = FALSE],
                                      k = init$k), extras))
        cvdev <- cvdev + plearn$dev
    }
    init$dev <- cvdev
    init
}

my.model.frame.tree <- function(formula, ...)
{
    m <- formula$model
    if(!is.null(m)) return(m)
    oc <- formula$call
    if(substring(deparse(oc[[1]]), 1, 7) == "predict") {
        m <- eval(oc$newdata)
        if(is.null(attr(m, "terms"))) {
            object <- eval(oc$object)
            m <- model.frame(object$terms, m, na.pass)
        }
        return(m)
    }
    while(deparse(oc[[1]]) != "tree")  oc <- eval(oc[[2]])$call
    oc$subset <- names(formula$where)
    oc$method <- "model.frame"
    # minka: use the environment of terms, not the current env
    env = environment(formula$terms)
    if (is.null(env)) env<-parent.frame()
    eval(oc,env)
}

# must do this before minka for bugfixes to apply
#library(tree)
if("tree" %in% .packages()) {
replaceInNamespace("print.tree",my.print.tree)
replaceInNamespace("partition.tree",my.partition.tree)
replaceInNamespace("misclass.tree",my.misclass.tree)
replaceInNamespace("deviance.tree",my.deviance.tree)
replaceInNamespace("cv.tree",my.cv.tree)
replaceInNamespace("model.frame.tree",my.model.frame.tree)
}
# bugfix for ts library
library(ts,warn=F)
as.data.frame.ts <- function(x,optional) {
  nam <- deparse(substitute(x))
  df <- as.data.frame.vector(as.vector(x))
  names(df) <- nam
  cbind(data.frame(time=as.vector(time(x))),df)
}
my.ccf <- function(x, y, lag.max = NULL,
                type = c("correlation", "covariance"),
                plot = TRUE, na.action = na.fail,
                   ylab=paste("Cross",type,sep="-"), ...)
{
    type <- match.arg(type)
    if(is.matrix(x) || is.matrix(y))
        stop("univariate time series only")
    X <- na.action(ts.union(as.ts(x), as.ts(y)))
    colnames(X) <- c(deparse(substitute(x)), deparse(substitute(y)))
    acf.out <- acf(X, lag.max = lag.max, plot = FALSE, type = type)
    lag <- c(rev(acf.out$lag[-1,2,1]), acf.out$lag[,1,2])
    y   <- c(rev(acf.out$acf[-1,2,1]), acf.out$acf[,1,2])
    acf.out$acf <- array(y, dim=c(length(y),1,1))
    acf.out$lag <- array(lag, dim=c(length(y),1,1))
    acf.out$snames <- paste(acf.out$snames, collapse = " & ")
    if (plot) {
      # minka: added ylab
      plot.acf(acf.out, ylab=ylab, ...)
      return(invisible(acf.out))
    } else return(acf.out)
}
replaceInNamespace("ccf",my.ccf)
# workaround bug in formals()
env = environment(acf)
formals(acf)$na.action = na.pass
environment(acf) = env
formals(ccf)$na.action = na.pass
environment(ccf) = env

##############################################################################
# new functions
# also see cluster.r

# lag.default should really be lag.ts
library(ts)
lag.ts = get("lag.default",asNamespace("ts"))
lag.default = function(x,k=1) {
  c(rep(NA,k),x[1:(length(x)-k)])
}
make.ar <- function(n,a,e,start) {
  k = length(a)
  if(missing(start)) {
    x = e*rnorm(k)/rnorm(k)
  } else {
    x = start
  }
  for(i in (k+1):n) {
    x[i] = sum(a*x[(i-1):(i-k)]) + e*rnorm(1)/rnorm(1)
  }
  x
}

growth.chart <- function(y,ref=1,ylab,...) {
  pct.change = function(x,r) (x-r)/r*100
  abs.change = function(x,r) (x-r)
  y = t(y)
  if(is.character(ref)) ref = pmatch(ref,rownames(y))
  nam = names(dimnames(y))
  if(missing(ylab)) ylab = paste("% change since",nam[1],rownames(y)[ref])
  if(is.data.frame(y)) y = as.matrix(y)
  a = attributes(y)
  r = apply(y,2,function(x) pct.change(x,x[ref]))
  attributes(r) = a
  linechart(t(r),ylab=ylab,...)
}

period.grid <- function(start,p,lty="dotted",labels=T,cex=0.6,col="blue",...) {
  v = seq(start,par("usr")[2],by=p)
  abline(v=v,lty=lty,...)
  if(labels) {
    i = (v-start)/p + 1
    h = par("usr")[4] - yinch(0.03)
    n = length(v)
    # place text in center of each interval, except for last one
    x = c((v[1:(n-1)]+v[2:n])/2, v[n]+xinch(0.03))
    text(x,h,i,adj=c(0,1),cex=cex,...)
  }
}

format.pretty <- function(x,digits=2) {
  prettyNum(formatC(x,format="fg",digits=digits),big.mark=",")
}

cycle.matrix <- function(x,p,extended=F) {
  # x is a vector
  # p is cycle length
  n = length(x)
  inc = p/ceiling(p)
  if(ceiling(p) != p) {
    i = seq(1,n,by=inc)
    # linear interpolation
    x = approx(1:n,x,i)$y
    n = length(x)
    p = ceiling(p)
  }
  d = c(p,ceiling(n/p))
  m = array(x[1:prod(d)],d)
  rownames(m) = format.pretty(seq(1,nrow(m)*inc,len=nrow(m)))
  colnames(m) = 1:ncol(m)
  if(extended) {
    m2 = cbind(m[,2:ncol(m)],NA)
    m = rbind(m,m2)
  }
  names(dimnames(m)) = c("Time","Cycle")
  m = as.table(m)
  t(m)
}

plot.stack <- function(x,y,n=2,...) {
  xyplot(y ~ x | equal.count(x, n=n, overlap=0),
	        panel = function(x, y) {
	                panel.xyplot(x, y, type="l")
	        },
		strip = F,
	 	layout = c(1,n),
		scale = list(x = "free"))
}

plot.stack <- function(x,y,n=2,ylab="",asp=NULL,...) {
  real.asp = if(identical(asp,"auto")) auto.aspect(x,y) else asp
  opar <- par(mfrow=c(n,1),mar=auto.mar(xlab="",ylab=ylab))
  on.exit(par(opar))
  #b <- equal.count(x,num=n,overlap=0)$intervals
  b = break.equal(x,n)
  b = cbind(my.lag(b),b)[2:(n+1),]
  for(j in seq(nrow(b))) {
    i <- (x >= b[j,1]) & (x <= b[j,2])
    plot(x[i],y[i],ylab=ylab,asp=real.asp,...)
    grid(nx=0)
  }
}

spiral <- function(object, ...) UseMethod("spiral")
spiral.formula <- function(fmla, data=parent.frame(), ..., main) {
  x = model.frame.default(fmla,data,na.action=na.pass)
  resp = response.var(x)
  pred = predictor.vars(x)[1]
  if(missing(main)) main = resp
  spiral.default(x[,pred],x[,resp],...,main=main)
}
spiral.default <- function(x,y,p,offset=0,color.palette=YR.colors(64),
                           clockwise=T,
                           pch.palette=20,bg=gray(0.5),cex=par("cex"),
                           cex.label=0.8,label=T,interp=T,main=NULL,...) {
  # time flows clockwise by default, unlike other spiral plots
  dir = if(clockwise) -1 else 1
  theta = (x - min(x,na.rm=T) + offset)/p
  indent = floor(max(theta,na.rm=T)*0.25)
  theta = theta + indent

  opar = par(bg=bg,mar=auto.mar(main=main,xlab="",ylab="",xaxt="n",yaxt="n"))
  on.exit(par(opar))

  r = range(theta)
  # sample uniformly in circumference
  res = min(par("pin"))/inch.symbol(pch=pch.palette[1],cex=cex)*2.5
  theta.g = sqrt(seq(r[1]^2,r[2]^2,len=res*diff(r)))
  r = theta.g
  xr = r*cos(2*pi*theta.g*dir)
  yr = r*sin(2*pi*theta.g*dir)
  plot(xr,yr,asp=1,type="l",lwd=0.5,col="gray",xaxt="n",yaxt="n",main=main,...)

  if(label) {
    spiral.text(theta.g[1],format(x[1]),"inside",clockwise,cex=cex.label)
    n = length(xr)
    spiral.text(theta.g[n],format(x[length(x)]),"outside",clockwise,cex=cex.label)
  }

  if(interp) {
    y.g = approx(theta,y,theta.g)$y
    theta = theta.g
    y = y.g
  }

  nlev = length(color.palette)
  if(nlev > 1) {
    col = color.palette[as.numeric(cut.quantile(y,nlev))]
    pch = pch.palette[1]
  } else {
    col = color.palette[1]
    pch = pch.palette[as.numeric(cut.quantile(y,nlev))]
  }

  r = theta
  xr = r*cos(2*pi*theta*dir)
  yr = r*sin(2*pi*theta*dir)
  points(xr,yr,col=col,pch=pch,cex=cex,...)
}
spiral.text <- function(theta,label,type=c("inside","outside"),clockwise=T,...) {
  # place text on the spiral
  f <- function(theta) {
    theta = theta %% 1
    if(theta <= 1/4) theta*4
    else if(theta <= 1/2) 1
    else if(theta <= 3/4) (3/4 - theta)*4
    else 0
  }
  adj.inside <- function(theta) {
    c(1-f(theta),1-f(theta+1/4))
  }
  adj.outside <- function(theta) {
    c(f(theta), f(theta+1/4))
  }
  type = match.arg(type)
  adj = if(type == "inside") adj.inside(theta) else adj.outside(theta)
  if(!clockwise) adj[2] = 1-adj[2]
  dir = if(clockwise) -1 else 1
  spc = xinch(inch.symbol())
  spc = spc*(0.5-adj)
  r = theta
  xr = r*cos(2*pi*theta*dir)
  yr = r*sin(2*pi*theta*dir)
  text(xr+spc[1], yr+spc[2], label, adj=adj, ...)
}

# Functions for generating VRML scenes.
# The *.wrl functions are internal, not meant to be called from outside.

truename <- function(f) {
  # turn a path into an absolute path
  if((substring(f,2,2) == ":") || (substring(f,1,1) %in% c("\\","/")))
    return(f)
  paste(getwd(),f,sep="\\")
}

vrml.new.name <- function() {
  assign("vrml.counter",vrml.counter + 1,globalenv())
  paste("obj",vrml.counter,sep="")
}

vrml.color <- function(col) {
  t(col2rgb(col))/255
}
vrml.material <- list(diffuseColor=c(1,1,1),
                      emissiveColor=c(0,0,0),
                      specularColor=c(0,0,0),
                      transparency=0,shininess=0.2,
                      ambientIntensity=0.2)
vrml.material.wrl <- function(diffuseColor=vrml.material$diffuseColor,
                              emissiveColor=vrml.material$emissiveColor,
                              transparency=vrml.material$transparency,
                              specularColor=vrml.material$specularColor,
                              shininess=vrml.material$shininess,
                              ambientIntensity=vrml.material$ambientIntensity,
                              ...,file=vrml.file) {
  # writes an Appearance node
  if(is.null(diffuseColor)) diffuseColor = vrml.material$diffuseColor
  cat("Appearance { material Material {",file=file)
  if(any(emissiveColor != 0)) cat(" emissiveColor",emissiveColor,file=file)
  if(any(specularColor != 0)) cat(" specularColor",specularColor,file=file)
  if(any(diffuseColor != 0.8)) cat(" diffuseColor",diffuseColor,file=file)
  if(transparency != 0) cat(" transparency",transparency,file=file)
  if(ambientIntensity != 0.2) cat(" ambientIntensity",ambientIntensity,file=file)
  # shininess does nothing without specularColor
  if(shininess != 0.2) cat(" shininess",shininess,file=file)
  cat(" } } ",file=file)
}
vrml.clear.material.cache <- function() {
  assign("vrml.materials",list(),globalenv())
}
vrml.as.material <- function(diffuseColor=vrml.material$diffuseColor,
                             emissiveColor=vrml.material$emissiveColor,
                             transparency=vrml.material$transparency,
                             specularColor=vrml.material$specularColor,
                             shininess=vrml.material$shininess,
                             ambientIntensity=vrml.material$ambientIntensity,
                             ...) {
  if(is.null(diffuseColor)) diffuseColor = vrml.material$diffuseColor
  list(diffuseColor=diffuseColor,
       emissiveColor=emissiveColor,
       transparency=transparency,
       specularColor=specularColor,
       shininess=shininess,
       ambientIntensity=ambientIntensity)
}
vrml.material.key <- function(...) {
  color.spec = vrml.as.material(...)
  color.spec = lapply(color.spec, function(x) {attributes(x) = NULL; x})
  paste(deparse(color.spec),collapse=" ")
}
vrml.material.wrl.cache <- function(...,file=vrml.file) {
  # writes an appearance field, with USE or DEF as appropriate
  if(!exists("vrml.materials")) vrml.clear.material.cache()
  key = vrml.material.key(...)
  the.appearance = vrml.materials[[key]]
  if(is.null(the.appearance)) {
    the.appearance = vrml.new.name()
    cat("appearance DEF",the.appearance,"",file=file)
    vrml.material.wrl(...,file=file)
    vrml.materials[[key]] = the.appearance
    assign("vrml.materials",vrml.materials,globalenv())
  } else {
    cat("appearance USE",the.appearance,"",file=file)
  }
}
vrml.material.keys <- function(color,...) {
  key = c()
  for(i in 1:nrow(color)) {
    key[i] = vrml.material.key(color[i,],...)
  }
  key
}

vrml.set.scale <- function(xlim,ylim,zlim,scale=c(1,1,1)) {
  # sets the globals vrml.offset and vrml.scale,
  # each vectors of 3 numbers
  s = scale/c(diff(range(xlim)),diff(range(ylim)),diff(range(zlim)))
  offset = -c(xlim[1],ylim[1],zlim[1])
  assign("vrml.scale",s,globalenv())
  assign("vrml.offset",offset,globalenv())
  assign("vrml.symbol.radius",0.01,globalenv())
}
vrml.to.scale <- function(p) {
  # p is matrix of three columns
  if(ncol(p) != 3) stop("p must be a matrix of 3 columns")
  p = p + rep.row(vrml.offset,nrow(p))
  scale.cols(p, vrml.scale)
}

# http://www.cs.cmu.edu/afs/cs/academic/class/16741-s02/www/lecture6.pdf
# http://www.cs.cmu.edu/afs/cs/academic/class/16741-s02/www/lecture7.pdf
vrml.quaternion.rotation <- function(rot) {
  rot[1:3] = rot[1:3]/sqrt(sum(rot[1:3]^2))
  c(cos(rot[4]/2),sin(rot[4]/2)*rot[1:3])
}
vrml.rotation.quaternion <- function(q) {
  qv = q[2:4]
  if(all(qv == 0)) c(0,1,0,0)
  else c(qv/norm(qv), 2*atan2(norm(qv),q[1]))
}
vrml.compose.quaternion <- function(q1,q2) {
  q1v = q1[2:4]
  q2v = q2[2:4]
  c(q1[1]*q2[1] - sum(q1v*q2v),q1[1]*q2v+q2[1]*q1v+vec.cross(q1v,q2v))
}
vrml.compose.rotation <- function(rot1,rot2) {
  q1 = vrml.quaternion.rotation(rot1)
  q2 = vrml.quaternion.rotation(rot2)
  q = vrml.compose.quaternion(q1,q2)
  vrml.rotation.quaternion(q)
}
vec.cross <- function(p1,p2) {
  # cross-product of vectors
  # http://mathworld.wolfram.com/CrossProduct.html
  cbind(p1[2]*p2[3] - p1[3]*p2[2],
        p1[3]*p2[1] - p1[1]*p2[3],
        p1[1]*p2[2] - p1[2]*p2[1])
}
vec.angle <- function(p1,p2) {
  # angle between two vectors
  acos(sum(p1*p2)/sqrt(sum(p1*p1)*sum(p2*p2)))
}
vrml.rotate <- function(p,rot) {
  angle = rot[4]
  v = rot[1:3]
  vp = vec.cross(v,p)
  vvp = vec.cross(v,vp)
  p + sin(angle)*vp + (1-cos(angle))*vvp
}
vrml.rotation.between <- function(p1,p2) {
  # returns a four-element vector specifying an axis-angle rotation
  # that takes p1 to p2
  if(!is.null(dim(p1)) && nrow(p1) > 1) {
    rot = vrml.rotation.between(p1[1,,drop=F],p2[1,,drop=F])
    p1[2,] = vrml.rotate(p1[2,],rot)
    rot2 = vrml.rotation.between(p1[2,,drop=F],p2[2,,drop=F])
    vrml.compose.rotation(rot2,rot)
  } else {
    if(all(p1 == 0) || all(p2 == 0)) stop("impossible rotation (to/from the origin)")
    axis = vec.cross(p1,p2)
    if(sum(axis*axis) == 0) axis = cbind(0,1,0)
    else axis = axis/sqrt(sum(axis*axis))
    c(axis,vec.angle(p1,p2))
  }
}

#############################################################################

vrml.viewpoint.wrl <- function(pos,lookat,file=vrml.file) {
  default.orient = cbind(0,0,-1)
  rot = vrml.rotation.between(default.orient,lookat-pos)
  cat("Viewpoint { position",pos,"orientation",rot,"}\n",file=file)
}

vrml.file = ""

vrml.open <- function(xlim,ylim,zlim,scale=c(1,1,1),light=F,bg=c(1,1,1),...,
                      file.name=NULL) {
  if(is.null(file.name)) file.name = file.path(tempdir(),"Rplot")
  if(substring(file.name,nchar(file.name)-3,nchar(file.name)) != ".wrl") {
    file.name = paste(file.name,"wrl",sep=".")
  }
  if(nchar(file.name) > 4) {
    con = file(file.name,"w")
  } else {
    con = ""
  }
  assign("vrml.file",con,globalenv())
  assign("vrml.file.name",file.name,globalenv())
  assign("vrml.counter",0,globalenv())
  vrml.clear.material.cache()
  cat("#VRML V2.0 utf8\n",file=vrml.file)
  #cat("WorldInfo { title \"",match.call()[[1]],"\" }\n",sep="",file=vrml.file)
  cat("NavigationInfo { type \"EXAMINE\" speed 10 ",file=vrml.file)
  if(light) cat("headlight FALSE ",file=vrml.file)
  cat("}\n",file=vrml.file)
  # white background is easier on the eyes
  if(length(bg) != 3) bg = col2rgb(bg)[,]
  cat("Background { skyColor [",bg,"] }\n",file=vrml.file)
  cat("DEF black Appearance { material Material {} }\n",file=vrml.file)
  # camera looks at the center of the back wall
  p = cbind(0.5,-1,0.5)*scale - cbind(0,0.5,0)
  lookat = cbind(0.5,0,0.5)*scale
  vrml.viewpoint.wrl(p,lookat)
  vrml.set.scale(xlim,ylim,zlim,scale)
}

vrml.close <- function() {
  if(inherits(vrml.file,"connection")) {
    close(vrml.file)
    cat("wrote to",vrml.file.name,"\n")
    shell.exec = if(exists("shell.exec")) shell.exec else browseURL
    shell.exec(truename(vrml.file.name))
  }
}

# baseline specifies the plane in which text is drawn
vrml.baseline = rbind(c(1,0,0),c(0,1,0))
vrml.text.wrl <- function(p,labels,adj=c(0,0),color=NULL,
                          baseline=vrml.baseline,cex=1,billboard=F,...,
                          file=vrml.file) {
  # if only one vector given for baseline, use the default height vector
  if(nrow(baseline) == 1) baseline = rbind(baseline,vrml.baseline[2,])
  rot = vrml.rotation.between(vrml.baseline,baseline)
  s = cex*vrml.symbol.radius*6
  if(adj[1] > 0) {
    # BUG: this requires an existing plot
    w = strwidth(labels,units="i")*5
    dim(w) = c(length(w),1)
    p = p - (w %*% baseline[1,,drop=F]*adj[1]*s)
  }
  if(adj[2] > 0) {
    h = strheight(labels,units="i")*6
    dim(h) = c(length(h),1)
    p = p - (h %*% baseline[2,,drop=F]*adj[2]*s)
  }
  for(i in 1:nrow(p)) {
    cat("Transform { translation",p[i,],"\n",
        "rotation",rot,"\n",
        "scale",rep(s,3),"\n",
        "children ",file=file)
    if(billboard) {
      cat("Billboard { axisOfRotation 0 0 0 children ",file=file)
    }
    cat("Shape { ",file=file)
    color.i = if(is.null(color) || nrow(color) == 1) color else color[i,]
    vrml.material.wrl.cache(color.i,...,file=file)
    txt = paste("\"",labels[i],"\"",sep="")
    cat("geometry Text { string",txt,"} }",file=file)
    if(billboard) cat(" } ",file=file)
    cat("}\n",file=file)
  }
}
vrml.text <- function(x,y,z,col=1,...) {
  p = vrml.to.scale(cbind(x,y,z))
  vrml.text.wrl(p,color=vrml.color(col),...)
}

vrml.sphere <- function(x,y,z,cex=1,col=1,...) {
  p = vrml.to.scale(cbind(x,y,z))
  r = cex*vrml.symbol.radius
  vrml.sphere.wrl(p,r,color=vrml.color(col),...)
}
vrml.sphere.wrl <- function(p,r,color=NULL,...,
                            file=vrml.file) {
  if(length(r) == 1) r = rep(r,nrow(p))
  for(i in 1:nrow(p)) {
    cat("Transform { translation",p[i,],"\n",
        "scale",rep(r[i],3),"\n",
        "children Shape { ",file=file)
    color.i = if(is.null(color) || nrow(color) == 1) color else color[i,]
    vrml.material.wrl.cache(color.i,...,file=file)
    if(i == 1) {
      the.shape = vrml.new.name()
      cat("geometry DEF",the.shape,"Sphere {} ",file=file)
    } else {
      cat("geometry USE",the.shape,"",file=file)
    }
    cat("}\n}\n",file=file)
  }
}
vrml.cube <- function(x,y,z,s=c(1,1,1),cex=1,col=1,...) {
  p = vrml.to.scale(cbind(x,y,z))
  s = s*vrml.symbol.radius*cex
  vrml.cube.wrl(p,s,color=vrml.color(col),...)
}
vrml.cube.wrl <- function(p,s=c(1,1,1),color=NULL,...,
                          file=vrml.file) {
  # this would be more efficient with an IndexedFaceSet
  #if(is.null(dim(color))) dim(color) = c(1,3)
  if(length(s) == 1) s = rep(s,3)
  p = na.omit(p)
  omitted = attr(p,"na.action")
  if(is.null(dim(s)) || nrow(s)==1) s = rep.row(s,nrow(p))
  else if(!is.null(omitted)) s = s[omitted,]
  for(i in 1:nrow(p)) {
    cat("Transform { translation",p[i,],"\n",
        "scale",s[i,],"\n",
        "children Shape { ",file=file)
    color.i = if(is.null(color) || nrow(color) == 1) color else color[i,]
    vrml.material.wrl.cache(color.i,...,file=file)
    if(i == 1) {
      the.shape = vrml.new.name()
      cat("geometry DEF",the.shape,"Box {} ",file=file)
    } else {
      cat("geometry USE",the.shape,"",file=file)
    }
    cat("}\n}\n",file=file)
  }
}

vrml.cylinder.wrl <- function(p1,p2,r=vrml.symbol.radius,file=vrml.file) {
  midp = (p1+p2)/2
  dp = p2-p1
  h = sqrt(sum(dp*dp))
  rot = vrml.rotation.between(cbind(0,1,0),dp)
  cat("Transform { translation",midp[1,],"\n",
      "rotation",rot,"\n",
      "children Shape { geometry Cylinder { radius",r,"height",h,"} }\n",
      "}\n",
      file=file)
}

vrml.axes <- function(xlim,ylim,zlim,xlab=NULL,ylab=NULL,zlab=NULL,cex=2,...,
                      file=vrml.file) {
  # origin
  po = vrml.to.scale(cbind(xlim[1],ylim[1],zlim[1]))
  vrml.cube.wrl(po,vrml.symbol.radius,...,file=file)
  px = vrml.to.scale(cbind(xlim[2],ylim[1],zlim[1]))
  py = vrml.to.scale(cbind(xlim[1],ylim[2],zlim[1]))
  pz = vrml.to.scale(cbind(xlim[1],ylim[1],zlim[2]))
  if(F) {
    # cylinders for axes
    vrml.cylinder.wrl(po,px,file=file)
    vrml.cylinder.wrl(po,py,file=file)
    vrml.cylinder.wrl(po,pz,file=file)
  } else {
    vrml.line.wrl(po,px,file=file)
    vrml.line.wrl(po,py,file=file)
    vrml.line.wrl(po,pz,file=file)
  }
  # axis names
  mar = vrml.symbol.radius
  if(!is.null(xlab)) {
    vrml.text.wrl(px/2-pz*mar,xlab,cex=cex,adj=c(0.5,1),base=rbind(px,pz),file=file)
    #vrml.text.wrl(px/2-py*mar+pz,xlab,cex=cex,adj=c(0.5,1),base=px,file=file)
    #vrml.text.wrl(px/2-py*mar,xlab,cex=cex,adj=c(0.5,1),base=-px,file=file)
  }
  if(!is.null(ylab)) {
    vrml.text.wrl(py/2-pz*mar,ylab,cex=cex,adj=c(0.5,1),base=rbind(py,pz),file=file)
    #vrml.text.wrl(py/2-px*mar+pz,ylab,cex=cex,adj=c(0.5,1),base=py,file=file)
    #vrml.text.wrl(py/2-px*mar,ylab,cex=cex,adj=c(0.5,1),base=-py,file=file)
  }
  if(!is.null(zlab)) {
    vrml.text.wrl(pz/2-px*mar,zlab,cex=cex,adj=c(0.5,0),base=rbind(pz,-px),file=file)
    #vrml.text.wrl(pz/2-py*mar+px,zlab,cex=cex,adj=c(0.5,1),base=-pz,file=file)
    #vrml.text.wrl(pz/2-py*mar,zlab,cex=cex,adj=c(0.5,1),base=pz,file=file)
  }
}
vrml.line.wrl <- function(p1,p2,file=vrml.file) {
  cat("Shape { appearance USE black geometry IndexedLineSet {\n",
      "coord Coordinate {\n",
      "point [",t(rbind(p1,p2)),"]\n }\n",
      "coordIndex [",c(0,1,-1),"]\n}}\n",
      file=file)
}
vrml.points <- function(x,y,z,col=1,...) {
  p = vrml.to.scale(cbind(x,y,z))
  vrml.points.wrl(p,color=vrml.color(col),...)
}
vrml.points.wrl <- function(p,color=NULL,...,
                            file=vrml.file) {
  cat("Shape { ",file=file)
  color.i = if(is.null(color) || nrow(color) == 1) color else color[1,]
  vrml.material.wrl.cache(color.i,...,file=file)
  cat("geometry PointSet {\n",file=file)
  if(nrow(color) == nrow(p)) {
    cat("color Color { color [",t(color),"] }\n",file=file)
  }
  cat("coord Coordinate { point [",t(p),"] }\n",
      "}}\n",
      file=file)
}

vrml.box.wrl <- function(p,file=vrml.file) {
  points = array(0,c(8,3))
  g = data.matrix(expand.grid(1:2,1:2,1:2))
  for(i in 1:nrow(g)) {
    for(j in 1:3) points[i,j] = p[g[i,j],j]
  }
  index = c(0,1,5,5,4,0,-1,2,3,7,6,2,-1,
    0,2,-1,1,3,-1,5,7,-1,4,6,-1)
  cat("Shape { appearance USE black geometry IndexedLineSet {\n",
      "coord Coordinate {\n",
      "point [",t(points),"]\n }\n",
      "coordIndex [",index,"]\n}}\n",
      file=file)
}
vrml.box <- function(xlim,ylim,zlim,...) {
  vrml.box.wrl(vrml.to.scale(cbind(xlim,ylim,zlim)),...)
}


vrml.plot3d <- function(object, ...) UseMethod("vrml.plot3d")
vrml.plot3d.formula <- function(formula,data=parent.frame(),...) {
  x = model.frame.default(formula,data,na.action=na.pass)
  vrml.plot3d.data.frame(x,...)
}
vrml.plot3d.data.frame <- function(x,labels=NULL,...,xlab,ylab,zlab) {
  resp = response.var(x)
  pred = predictor.vars(x)
  if(missing(xlab)) xlab = pred[1]
  if(missing(ylab)) ylab = pred[2]
  if(missing(zlab)) zlab = resp
  if(identical(labels,TRUE)) labels = rownames(x)
  vrml.plot3d.default(x[,pred[1]],x[,pred[2]],x[,resp],labels=labels,
                  xlab=xlab,ylab=ylab,zlab=zlab,...)
}
vrml.plot3d.default <- function(x,y,z,data=parent.frame(),labels=NULL,
                                xlab,ylab,zlab,
                                pch=0,col=3,
                                scale=c(1,1,1),asp=NULL,
                                file.name=NULL,
                            cex=1,cex.axis=2,light=T,axes=T,type="p",...) {
  if(missing(xlab)) xlab <- deparse(substitute(x))
  if(missing(ylab)) ylab <- deparse(substitute(y))
  if(missing(zlab)) zlab <- deparse(substitute(z))
  x <- eval(substitute(x),data)
  y <- eval(substitute(y),data)
  z <- eval(substitute(z),data)

  xlim = range(x,na.rm=T); ylim = range(y,na.rm=T); zlim = range(z,na.rm=T)
  xlim = extend.pct(xlim); ylim = extend.pct(ylim); zlim = extend.pct(zlim)
  if(!is.null(asp)) {
    if(length(asp) == 1) asp = rep(asp,2)
    scale = c(diff(xlim),diff(ylim)*asp[1],diff(zlim)*asp[2])
    scale = scale/diff(xlim)
  }
  vrml.open(xlim,ylim,zlim,scale,light,file.name=file.name,...)
  if(light) {
    cat("PointLight { ambientIntensity 0.5 intensity 0.25 location 0 0 1 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.25 location 1 1 1 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.25 location 0 1 1 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.25 location 1 0 1 }\n",
        file=vrml.file)
  }
  if(axes) {
    cat("Group { children [\n",file=vrml.file)
    vrml.box(xlim,ylim,zlim)
    vrml.axes(xlim,ylim,zlim,xlab,ylab,zlab,cex=cex.axis,...)
    cat("]}\n",file=vrml.file)
  }
  cat("Group { children [\n",file=vrml.file)
  if(any(scale/cex >= 5) || length(x) > 5000) {
    pch = "."
    cat("Using pch=. to save time\n")
  }
  if(is.null(labels)) {
    if(pch == ".") {
      vrml.points(x,y,z,col=col,...)
    } else if(pch==1) {
      vrml.sphere(x,y,z,col=col,cex=cex,...)
    } else {
      # pch=0 is a cube
      vrml.cube(x,y,z,col=col,cex=cex,...)
    }
  } else {
    vrml.text(x,y,z,labels=labels,col=col,cex=cex,billboard=T,...)
  }
  if(type == "o") {
    p = vrml.to.scale(cbind(x,y,z))
    cat("Shape { appearance USE black geometry IndexedLineSet {\n",
        "coord Coordinate {\n",
        "point [",t(p),"]\n }\n",
        "coordIndex [",c(1:nrow(p),0)-1,"]\n}}\n",
        file=vrml.file)
  }
  cat("]}\n",file=vrml.file)
  vrml.close()
}
vrml.color.plot <- function(x,y,nlevels=4,color.palette=YlGnBu.colors,...) {
  if(missing(y)) {
    resp = response.var(x)
    y = x[[resp]]
    pred = predictor.terms(x)
    x = x[pred]
  }
  if(!is.factor(y)) y = cut.quantile(y,nlevels)
  else {
    color.palette = default.colors
    nlevels = length(levels(y))
  }
  if(is.function(color.palette)) color.palette = color.palette(nlevels)
  col = color.palette[y]
  vrml.plot3d(formula(x),x,col=col,light=F,...)
}

vrml.surface <- function(object,...) UseMethod("vrml.surface")
vrml.surface.loess <- function(object,res=20,
                        xlim,ylim,zlim,clip=T,xlab,ylab,zlab,...) {
  if(length(res) == 1) res <- rep(res,2)
  x <- model.frame(object)
  resp <- response.var(object)
  pred <- predictor.vars(object)
  if(missing(xlim)) xlim <- range(x[[pred[1]]],na.rm=T)
  x1 <- seq(xlim[1],xlim[2],length=res[1])
  if(missing(ylim)) ylim <- range(x[[pred[2]]],na.rm=T)
  x2 <- seq(ylim[1],ylim[2],length=res[2])
  xt <- expand.grid(x1,x2)
  names(xt) <- pred[1:2]
  z <- predict(object,xt)
  if(!missing(zlim)) {
    z[z < zlim[1]] <- zlim[1]
    z[z > zlim[2]] <- zlim[2]
  }
  if(clip != F) {
    if(clip == T) {
      i <- chull(x[pred])
      clip <- x[pred][i,]
    }
    z[!in.polygon(clip,xt)] <- NA
  }
  dim(z) <- c(length(x1),length(x2))
  if(missing(xlab)) xlab <- pred[1]
  if(missing(ylab)) ylab <- pred[2]
  if(missing(zlab)) zlab <- resp
  vrml.surface.default(x1,x2,z,xlab=xlab,ylab=ylab,zlab=zlab,...)
}
vrml.surface.default <- function(x,y,z,xlab=NULL,ylab=NULL,zlab=NULL,
                                 col="gray",scale=c(1,1,1),file.name=NULL,
                                 cex.axis=2,light=F,...) {
  if(is.null(xlab)) xlab <- deparse(substitute(x))
  if(is.null(ylab)) ylab <- deparse(substitute(y))
  if(is.null(zlab)) zlab <- deparse(substitute(z))
  #x <- eval(substitute(x),data)
  #y <- eval(substitute(y),data)
  #z <- eval(substitute(z),data)

  xlim = range(x,na.rm=T); ylim = range(y,na.rm=T); zlim = range(z,na.rm=T)
  xlim = extend.pct(xlim); ylim = extend.pct(ylim); zlim = extend.pct(zlim)
  if(identical(scale,F)) {
    scale = c(diff(range(xlim)),diff(range(ylim)),diff(range(zlim)))
    scale = scale/sum(scale)*3
  }
  vrml.open(xlim,ylim,zlim,scale,light,file.name=file.name,...)
  if(light) {
    # dir 1 1 -1
    cat("DirectionalLight { ambientIntensity 0 intensity 1",
        "direction 1 1 -0.1 }\n",file=vrml.file)
  }
  cat("Group { children [\n",file=vrml.file)
  vrml.box(xlim,ylim,zlim)
  vrml.axes(xlim,ylim,zlim,xlab,ylab,zlab,cex=cex.axis,...)
  cat("]}\n",file=vrml.file)
  cat("Group { children [\n",file=vrml.file)
  vrml.surface.wrl(x,y,z,color=vrml.color(col),...)
  cat("]}\n",file=vrml.file)
  vrml.close()
}
vrml.surface.wrl <- function(x,y,z,border=F,
                             creaseAngle=10,...,file=vrml.file) {
  # turn off LOD in the viewer or you will get weird effects
  # could also use an ElevationGrid
  points = expand.grid(x,y)
  points = cbind(points,as.vector(z))
  points = vrml.to.scale(points)
  face = 0
  nx = length(x)
  ny = length(y)
  index = array(0,c((nx-1)*(ny-1)*2,4))
  for(i in 0:(ny-2)) {
    for(j in 0:(nx-2)) {
      face = face + 1
      index[face,] = c(i*nx + j+1, i*nx + j, (i+1)*nx + j, -1)
      if(any(is.na(points[index[face,1:3]+1,3]))) face = face - 1
      face = face + 1
      index[face,] = c(i*nx + j+1, (i+1)*nx + j, (i+1)*nx + j+1, -1)
      if(any(is.na(points[index[face,1:3]+1,3]))) face = face - 1
    }
  }
  index = index[1:face,]
  cat("Shape { appearance ",file=file)
  vrml.material.wrl(...,file=file)
  cat("geometry IndexedFaceSet { solid FALSE creaseAngle",creaseAngle,"\n",
      "coord Coordinate {\n",
      "point [",t(na.dummy(points)),"]\n }\n",
      "coordIndex [",t(index),"]\n}}\n",file=file)

  if(border) {
    points[,3] = points[,3] + 4e-3
    index = array(0,dim=c((nx-1)*(ny-1),5))
    face = 0
    for(i in 0:(ny-2)) {
      for(j in 0:(nx-2)) {
        face = face + 1
        index[face,] = c(i*nx + j+1, i*nx + j, (i+1)*nx + j, (i+1)*nx + j+1, -1)
        if(any(is.na(points[index[face,1:4]+1,3]))) face = face - 1
      }
    }
    index = index[1:face,]
    cat("Shape { appearance ",file=file)
    # only emissiveColor is used
    vrml.material.wrl(emissiveColor=vrml.color(1),file=file)
    cat("geometry IndexedLineSet {",
        "coord Coordinate {\n",
        "point [",t(na.dummy(points)),"]\n }\n",
        "coordIndex [",t(index),"]\n}}\n",file=file)
  }
}
vrml.array.wrl <- function(m,file=vrml.file) {
  if(is.data.frame(m)) m = data.matrix(m)
  for(i in 1:nrow(m)) {
    cat(m[i,],", ",file=file)
  }
}
vrml.cone.grid <- function() {
  ntheta = 24
  points = array(0,dim=c(6*ntheta,3))
  col = c()
  index = c()
  k = 0
  for(z in seq(0,1,len=6)) {
    index = c(index,k+(1:ntheta)-1,-1)
    for(theta in seq(0,1,len=ntheta)) {
      x = z*cos(2*pi*theta)/2 + 0.5
      y = z*sin(2*pi*theta)/2 + 0.5
      k = k + 1
      points[k,] = c(x,y,z)
      col[k] = hsv(theta,1,z)
    }
  }
  cat("Shape { appearance ",file=vrml.file);
  vrml.material.wrl(vrml.color(grey(0.5)),file=vrml.file)
  cat("geometry IndexedLineSet {",
      "coord Coordinate {\n",
      "point [",t(points),"]\n }\n",
      "colorPerVertex TRUE\n",
      "color Color { color [",t(vrml.color(col)),"]\n}",
      "coordIndex [",t(index),"]\n}}\n",file=vrml.file)
}
color.cone <- function(col,cex=2,light=T,...) {
  # convert to HSY
  h = col2hsv(col)
  # replace Value by Lightness
  r = col2rgb(col)/255
  if(F) {
    r = r^3
    h[3,] = array(c(0.2126,0.7152,0.0722),c(1,3)) %*% r
    # lightness
    h[3,] = h[3,]^(1/3)
  } else {
    # approximation to lightness
    h[3,] = sqrt(colSums(r*r)/3)
  }
  h[3,] = h[3,]*255
  r <- hsv2cone(h)/255
  r = as.data.frame(t(r))
  x = r$x
  y = r$y
  z = r$z
  xlim = c(-1,1)
  ylim = c(-1,1)
  zlim = c(0,1)
  vrml.open(xlim,ylim,zlim,...)
  if(light) {
    cat("PointLight { ambientIntensity 0.5 intensity 0.5 location 0 0 0 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.5 location 1 1 0 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.5 location 0 1 0 }\n",
        file=vrml.file)
    cat("PointLight { intensity 0.5 location 1 0 0 }\n",
        file=vrml.file)
  }
  cat("Group { children [\n",file=vrml.file)
  vrml.axes(xlim,ylim,zlim,"x","y","z")
  if(T) {
    # outline of the cone
    cat("Transform { translation",cbind(0.5,0.5,0.5),"\n",
        "rotation",cbind(1,0,0,-pi/2),"\n",
        "scale",cbind(1,1,1)/2,"\n",
        "children Shape { ",file=vrml.file)
    cat("appearance ",file=vrml.file)
    vrml.material.wrl(diffuseColor=vrml.color(1),transparency=0.9,
                      specularColor=vrml.color(8),shininess=1,
                      file=vrml.file)
    cat("geometry Cone {}",file=vrml.file)
    cat("}\n}\n",file=vrml.file)

    # grey line
    vrml.line.wrl(cbind(0.5,0.5,0),cbind(0.5,0.5,1))
    vrml.cone.grid()
  }
  cat("]}\n",file=vrml.file)
  cat("Group { children [\n",file=vrml.file)
  vrml.sphere(x,y,z,col=col,cex=cex,...)
  p1 = vrml.to.scale(cbind(x[1],y[1],z[1]))
  n = length(x)
  for(i in 2:n) {
    p2 = vrml.to.scale(cbind(x[i],y[i],z[i]))
    vrml.line.wrl(p1,p2)
    p1 = p2
  }
  cat("]}\n",file=vrml.file)
  vrml.close()
}
test.color.cone <- function() {
  x = expand.grid(seq(0,1),seq(0,1),seq(0,1))
  col = rgb(x[,1],x[,2],x[,3])
  color.cone(col)
}
