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


#' Internal routines
#'
#' Internal routines for \code{\link{merge.table}}
#'
#'
#' @aliases merge.table.cost merge.table.cost.inorder merge.table.cells
#' @author Tom Minka
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


#' Table merging
#'
#' Merges similar rows and columns of a contingency table.
#'
#' The desired table dimensions are achieved by successively merging the two
#' most similar slices. (`Slice' generalizes `row' and `column' to
#' higher-dimensional tables.) The distance between slices is measured
#' according to the chi-square statistic. Merging two slices means adding
#' together their counts, and concatenating their labels with a comma in
#' between. If a dimension is ordered (according to \code{\link{dim.ordered}}),
#' only adjacent slices are considered for merging, and their labels are
#' concatenated with a dash in between.
#'
#' @param x a \code{\link{table}}
#' @param bins the desired number of levels for each dimension being merged. a
#' numeric vector, the same length as \code{ds}.
#' @param ds a vector of dimensions to merge, either by name or number. default
#' is all of them.
#' @return A merged \code{\link{table}}.  The total count is the same as
#' \code{x}. A merging trace is plotted which shows, for each merge, the
#' chi-square distance of the slices which were merged. This is useful for
#' determining the appropriate dimensions.  An interesting number is one that
#' directly precedes a sudden jump in the chi-square distance.
#' @author Tom Minka
#' @seealso \code{\link{sort.table}}, \code{\link{mosaicplot}},
#' \code{\link{linechart}}
#' @examples
#'
#' i <- factor(c(1,1,2,2,3,3,4,4))
#' j <- factor(c(3,4,3,4,1,2,1,2))
#' x <- table(i,j)
#' merge.table(x,c(2,2))
#'
#' i <- factor(c(1,1,3,3,2,2,4,4))
#' j <- factor(c(2,4,2,4,1,3,1,3))
#' x <- table(i,j)
#' merge.table(x,c(2,2))
#'
#' # one ordered dimension
#' data(education)
#' merge.table(education,c(3,2))
#'
#' data(occupation)
#' merge.table(occupation,c(3,4))
#'
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
    diag(d) = 0
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

