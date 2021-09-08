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

#' @export
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

#' @export
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
  if(is.null(orth)) e <- eigen(iba,symmetric=F)
  else {
    if(is.vector(orth)) orth <- array(orth,c(length(orth),1))
    n <- nrow(orth)
    orth.k <- ncol(orth)
    orth.m <- eye(n) - (orth %*% t(orth))
    x <- orth.m %*% iba %*% orth.m
    e <- eigen(x,symmetric=F)
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

