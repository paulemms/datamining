# functions to convert colour representations

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


#' Sort levels of a factor
#'
#' The levels of a factor are sorted according to a summary statistic.
#' @param f a factor
#' @param x a vector, same length as \code{f}
#' @param fun function to use to summarize groups of \code{x}
#' @return
#'   A copy of \code{f} with re-ordered levels.
#' @author Tom Minka
#' @examples
#' data(OrchardSprays)
#' sort_levels(OrchardSprays$treatment, OrchardSprays$decrease)
#' sort_levels(OrchardSprays$treatment, OrchardSprays$decrease, fun=mean)
#' @export
sort_levels <- function(f,x,fun=median) {
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

#' @export
is.dimOrdered <- function(d) {
  # dimensions are ordered by default
  a <- attr(d,"ordered")
  is.null(a) || a
}
#' @export
dimOrdered <- function(x) {
  dn = dimnames(x)
  # never return NULL
  if(is.null(dn)) dn = dim(x)
  sapply(dn,is.dimOrdered)
}
#' @export
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
                     las=par("las"),xtick=NULL,ytick=NULL,cex.axis=1,
                     mar.scale=NULL,key=F,...) {
  # note: cex=1 actually means cex=par("cex")
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
  if(xaxt != "n" && (las == 3 || las == 2)) {
    if(!is.null(xtick)) {
      # this requires a figure to already be set up
      mar[1] = max(strwidth(xtick,units="inches",cex=cex.axis))
      # convert to "lines of text"
      #mar[1] = mar[1]/min(strheight(ytick,units="inches",cex=cex))
      mar[1] = mar[1]/0.1875 + 1
    } else {
      mar[1] = mar[1]*mar.scale
    }
  }
  if(yaxt != "n" && (las == 1 || las == 2)) {
    if(!is.null(ytick)) {
      # this requires a figure to already be set up
      mar[2] = max(strwidth(ytick,units="inches",cex=cex.axis))
      # convert to "lines of text"
      #mar[2] = mar[2]/min(strheight(ytick,units="inches",cex=cex))
      mar[2] = mar[2]/0.1875 + 1
    } else {
      mar[2] = mar[2]*mar.scale
    }
  }
  if(!is.null(sub) && sub != "") {
    mar[1] = mar[1] + 2
  }
  mar
}

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

