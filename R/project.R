

#' Project data into fewer dimensions
#'
#' Reduces the dimensionality of a data set.
#'
#' @param x a data frame
#' @param w a matrix with named rows and columns
#' @details
#'   Each column of \code{w} specifies a new variable, which is to be
#'   constructed by combining existing variables according to the given weights.
#' @return
#'   A data frame where the variables named in the rows of \code{w} are
#'   replaced by new variables named in the columns of \code{w}.
#'   Other variables are left unchanged.
#' @author Tom Minka
#' @seealso
#'   \code{\link{pca}}, \code{\link{projection}}, \code{\link{plot.axes}}
#' @examples
#' data(iris)
#' #w = projection(iris, k=2)
#' # w only involves the continuous attributes
#' # the new variables are h1 and h2
#' #x = project(iris, w)
#' #color.plot(x)
#' #plot.axes(w)
#' @export
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


#' Principal Component Analysis
#'
#' Computes a projection matrix that preserves spread.
#'
#' @param x a data frame.
#' @param k the number of dimensions to project down to.
#' @details
#'   The projection is chosen so that the projected data is as "spread out"
#'   as possible, measured according to the determinant of the covariance
#'   matrix.
#'   This turns out to be the top \code{k} eigenvectors of the data
#'   covariance matrix.
#'
#'   The projection is "stabilized" so that small changes in the data do
#'   not cause sign flips in the projection.
#' @return
#'   A matrix with named rows matching the numeric columns of \code{x} and
#'   columns named \code{h1}, ..., \code{hk}.
#'   Each column denotes a new dimension to be obtained as a linear combination
#'   of the numeric variables in \code{x}.
#'
#' @author Tom Minka
#' @seealso \code{\link{project}}, \code{\link{projection}},
#'   \code{pca} in the \code{multiv} package
#' @examples
#' data(Housing)
#' w <- pca(HousingT, k=2)
#' #plot(project(HousingT, w), asp=1)
#' #plot.axes(w)
#' @export
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


#' Discriminative projection
#'
#' Finds a projection matrix that separates data groups.
#'
#' @param x a data frame of variables to project.
#' @param y a factor or numeric vector which defines the groups to separate.
#'     If \code{y} is not given, it is
#'     taken to be the response variable of \code{x} (the last column if
#'                                                    \code{x} is not a \code{model.frame}).
#' @param k the number of dimensions to project \code{x} down to.
#' @param given A matrix specifying axes to avoid.  The projection matrix
#'     will be orthogonal to \code{given}.
#' @param type see below.
#' @param ... additional parameters depending on \code{type}.
#' @details
#'   This function only uses the second-order statistics of the data (means
#'                                                                    and covariances of the groups).
#'
#'   If \code{type="m"}, the within-group covariances are assumed equal
#'   and the projection will try to separate the projected means.
#'
#'   If \code{type="v"}, the within-group means are assumed equal
#'   and the projection will try to separate the projected covariances.
#'
#'   If \code{type="mv"}, the projection will try to separate the projected
#'   means and covariances, by maximizing the divergence between the
#'   projected classes (as Gaussians).
#'
#'   If \code{y} is a numeric vector, overlapping classes are defined by
#'   grouping data points with similar values of \code{y}.
#'   The optional argument \code{span} controls how big the classes will
#'   be (as a percentage of the dataset), and \code{res} controls the
#'   amount of overlap.  The total number of classes will be
#'   \code{res}/\code{span}.  The default values are usually acceptable.
#'
#'   The projection is "stabilized" so that small changes in the data do
#'   not cause sign flips in the projection.
#' @return
#'   A matrix suitable for input to \code{\link{project}},
#'   with named rows matching columns of \code{x} and
#'   columns named \code{h1}, ..., \code{hk}.
#'   Each column denotes a new dimension to be obtained as a linear combination
#'   of the variables in \code{x}.
#' @author Tom Minka
#' @references
#'   m-projection is Fisher's linear discriminant analysis.
#'   mv-projection is heteroscedastic discriminant analysis:
#'
#'   N. Kumar and A.G. Andreou.
#'   Heteroscedastic discriminant analysis and
#'   reduced rank HMMs for improved speech recognition.
#'   \emph{Speech Communication} 26: 283-297, 1998.
#'
#'   The case when \code{y} is numeric is sliced inverse
#'   regression:
#'
#'   K.-C. Li.  Sliced inverse regression for dimension reduction.
#'   \emph{Journal of the American Statistical Association} 86(414):
#'   316-327, 1991.
#' @seealso \code{\link{project}},\code{\link{pca}}
#' @examples
#' # illustrate difference between (m,v,mv)
#' library(MASS)
#' m1 <- c(6,6)
#' v1 <- array(c(2,1.9,1.9,2),c(2,2))
#' #v1 <- array(c(1,0,0,1),c(2,2))
#' x1 <- mvrnorm(100,m1,v1)
#' m2 <- c(0,0)
#' v2 <- array(c(20,0,0,10),c(2,2))
#' x2 <- mvrnorm(300,m2,v2)
#' x = as.data.frame(rbind(x1,x2))
#' y = factor(c(rep(1,nrow(x1)),rep(2,nrow(x2))))
#' plot(x[,1],x[,2],col=1,xlab="",ylab="",asp=1)
#' points(x2[,1],x2[,2],col=2)
#' w = projection(x,y,type="m")
#' abline(0,w[2]/w[1],col=3)
#' w = projection(x,y,type="v")
#' abline(0,w[2]/w[1],col=4)
#' w = projection(x,y,type="mv")
#' abline(0,w[2]/w[1],col=5)
#' my.legend(1,c("m","v","mv"),col=3:5,lty=1)
#'
#' # regression projection
#' x1 <- 2*runif(200)-1
#' x2 <- 2*runif(200)-1
#' y <- x1^2/2 + x2^2
#' x <- data.frame(x1,x2)
#' #color.plot(x[,1],x[,2],y)
#' w = projection(x,y)
#' abline(0,w[2]/w[1],col=4)
#' @export
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

#' Plot axes under projection
#'
#' Shows how unit vectors along the original dimensions appear
#'   under a projection.
#' @param w a numeric array with two columns and named rows
#' @param col color of arrows
#' @param origin If T, arrows emerge from (0,0).  If F, arrows emerge
#'     from the center of the figure.  If not given, arrows emerge from
#'     (0,0) if (0,0) is visible, otherwise from the figure center.
#' @param keep a length in inches, below which an arrow is not plotted.
#' @details
#'   Each row of the array specifies a location in the figure.  This
#'   location can be understood as the place where a unit vector
#'   along one of the original dimensions would appear.  An arrow is
#'   drawn from the origin pointing toward this location and labeled with
#'   the row name.  The relative length of the arrow is determined by the
#'   distance of the location from the origin, but all arrows are scaled as
#'   much as possible, to reduce crowding, while remaining inside the
#'   figure limits.  Arrows which are very short are not plotted at all.
#' @return Uses \code{\link{arrows}} to draw arrows on top of the current plot.
#' @author{Tom Minka}
#' @seealso
#'   \code{\link{project}}, \code{\link{pca}}, \code{\link{projection}}
#' @export
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
  s = eigen(A,symm=T)
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
    r <- rbind_extend(r,data.frame(ri))
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

