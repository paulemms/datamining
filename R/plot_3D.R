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
draw.ellipse <- function(x,y,rx,ry,res=64,fill=F,...) {
  theta <- seq(0,2*pi,len=res)
  if(length(rx) == 1) rx <- rep(rx,length(x))
  if(length(ry) == 1) ry <- rep(ry,length(x))
  xt <- cos(theta)
  yt <- sin(theta)
  for(i in 1:length(x)) {
    if(fill)
      polygon(x[i]+xt*rx[i],y[i]+yt*ry[i],...)
    else
      lines(x[i]+xt*rx[i],y[i]+yt*ry[i],...)
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

if(!exists("persp.default"))
  persp.default = getFromNamespace("persp.default","graphics")
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
