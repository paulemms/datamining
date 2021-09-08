# bug fixes

replaceInNamespace <- function(name,value) {
  # if(exists("asNamespace")) {
  #   environment(value) = environment(get(name))
  #   assignInNamespace(name,value,environment(value))
  #   #assignInNamespace(name,value,env=which.environment(name))
  # } else {
  #   assign(name,value,env=.GlobalEnv)
  # }
}

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
#library(ts,warn=F)
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

# new functions
# also see cluster.r

# lag.default should really be lag.ts
lag.ts = getFromNamespace("lag.default","stats")
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

