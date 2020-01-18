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

#formula.nnet = formula.lm

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
# R 1.9.0
#coef.logistic = coef.glm
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


#' Linechart
#' 
#' Plot each data row as a curve.
#' 
#' If \code{xscale="linear"}, the columns are placed to make each curve as
#' straight as possible.  If \code{xscale="equal"}, the columns are placed
#' similar to \code{"linear"} but with the constraint that they must be equally
#' spaced.  If \code{xscale="none"}, the columns are placed in the order that
#' they appear in the matrix.  This is automatic if \code{y} has ordered
#' columns (see \code{\link{dim.ordered}}).  If \code{se != NULL}, error bars
#' are drawn around each point.
#' 
#' Linecharts are a replacement for dotcharts and mosaics.
#' 
#' @param y a named vector or matrix.
#' @param se a vector or matrix, same size as \code{y}, of error bounds.
#' Alternatively, \code{y} can be \code{list(y,se)}.
#' @param effects If \code{TRUE}, the columns are shifted to have zero mean.
#' @param med If \code{TRUE} and \code{effects=TRUE}, the columns are shifted
#' to have zero median.
#' @param xscale describes how the columns should be placed.
#' @param ... additional arguments to \code{\link{labeled.curves}}.
#' @author Tom Minka
#' @seealso \code{\link{dotchart}}, \code{\link{mosaicplot}}
#' @examples
#' 
#' # compare to a dotchart
#' data(VADeaths)
#' dotchart(VADeaths, main = "Death Rates in Virginia - 1940")
#' dimOrdered(VADeaths)[2] = F
#' linechart(VADeaths)
#' linechart(t(VADeaths))
#' 
#' # compare to a mosaicplot
#' data(HairEyeColor)
#' x <- margin.table(HairEyeColor,c(1,2))
#' dimOrdered(x) = F
#' mosaicplot(x)
#' x = t(x)
#' col = c("brown","blue","red","green")
#' linechart(row.probs(x),color.pal=col)
#' linechart(row.probs(x,se=T),color.pal=col)
#' linechart(row.probs(x,se=T),jitter=0.02,color.pal=col)
#' mosaicplot(x)
#' linechart(row.probs(t(x),se=T))
#' 
#' data(blood)
#' dimOrdered(blood) = F
#' linechart(row.probs(blood,se=T))
#' 
#' data(antacids)
#' dimOrdered(antacids) = F
#' linechart(row.probs(antacids,se=T))
#' mosaicplot(t(antacids))
#' 
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



#' Plot labeled curves
#' 
#' Draws multiple curves, each with a properly-placed label and unique
#' color/linestyle.
#' 
#' Point \code{j} in curve \code{i} is at \code{(x[j],y[j,i])}.  Thus all
#' curves must be the same length, and have points at the same horizontal
#' positions.  If \code{y[i,j]=NA}, then the curve has a break at the previous
#' point, and resumes at the next \code{non-NA} point.
#' 
#' If \code{se} is given, then an error bar will be placed around each point.
#' 
#' @param x a numeric vector giving the horizontal position of each point.
#' @param y a matrix or data frame giving the vertical position of each point.
#' Each column defines one curve.
#' @param se a matrix or data frame, the same size as \code{y}, giving an error
#' bound on each value.
#' @param labels a character vector of labels for the lines.  Defaults to the
#' column names of \code{y}.
#' @param xlab,ylab axis labels, by default the dimnames of \code{y}.
#' @param type indicates whether to draw points, lines, or both.  See
#' \code{\link{plot}}.
#' @param group a numeric vector specifying how to group cases when assigning
#' colors.
#' @param color.palette a vector of colors, arbitrary length, or a function
#' with integer argument which generates a vector of colors.  Used if
#' \code{col} is not specified.  If shorter than the number of curves, colors
#' will be recycled and the line style will change.
#' @param col a vector of colors, as in a call to \code{\link{plot}}.  Used to
#' explicitly set the color of each curve.
#' @param lty.palette a vector of line styles, arbitrary length.  The line
#' style will rotate through these when there aren't enough colors.
#' @param lty a vector of line styles, as in a call to \code{\link{plot}}.
#' Used to explicitly set the style of each curve.
#' @param lwd line width.
#' @param jitter the amount by which to jitter the error bars, to avoid
#' overplotting.  \code{jitter=0.02} is usually sufficient.
#' @param legend If \code{TRUE}, the labels are placed in a legend.  Otherwise
#' the labels are placed next to the lines.
#' @param cex character expansion factor.
#' @param horizontal If \code{TRUE}, the axes are laid out horizontally, so
#' that the curves run vertically.
#' @param mar a numerical vector giving the lines of margin on the four sides
#' of the plot (see \code{\link{par}}).  If missing, it is set according to
#' \code{\link{auto.mar}}.
#' @param bg background color.
#' @param ylim desired plotting limits.
#' @param main title for the plot.
#' @param ... extra arguments for low-level plotting commands.
#' @author Tom Minka
#' @seealso \code{\link{parallel.plot}}
labeled.curves <- function(x,y,se=NULL,labels,xtick=names(x),
                           xlab=NULL,ylab,type="o",
                           group,
                           color.palette=default.colors(6),col,
                           lty.palette=1:6,lty,lwd=2,jitter=0,
                           legend.=NULL,move=T,dryrun=F,
                           cex=1,horizontal=F,
                           mar,bg=par("bg"),ylim,main="",
                           xaxt=par("xaxt"),yaxt=par("yaxt"),...) {
  # must include (xaxt,yaxt) because they are swapped for horizontal
  # x is vector
  # y is a list or matrix of columns
  if(horizontal && missing(legend.)) legend. = 1
  if(horizontal) x = -x
  nam = names(dimnames(y))
  if(is.null(xtick)) xtick = rownames(y)
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
      auto.mar(main=main,xlab=ylab,ylab=xlab,xaxt=yaxt,yaxt=xaxt,ytick=xtick,...)
    else
      auto.mar(main=main,xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,xtick=xtick,...)
  }
  if(!missing(bg)) opar = par(mar=mar,bg=bg)
  else opar = par(mar=mar)
  on.exit(par(opar))
  # plot.new here so that yinch and strwidth work
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
  #if(dryrun) return(list(mar,xlim,ylim))
  if(horizontal) {
    plot.window(xlim=ylim,ylim=xlim)
    if(is.null(xtick)) {
      axis(2,yaxt=xaxt,...)
      grid(nx=0)
    } else {
      #cat("columns are",names(x),"\n")
      axis(2,x, labels=xtick,yaxt=xaxt,...)
      # grid
      abline(h=x,lty="dotted",col="lightgray")
    }
    axis(1,xaxt=yaxt,...)
  }
  else {
    plot.window(xlim=xlim,ylim=ylim)
    if(is.null(xtick)) {
      axis(1,xaxt=xaxt,...)
      grid(ny=0)
    } else {
      #cat("columns are",names(x),"\n")
      axis(1,x, labels=xtick,xaxt=xaxt,...)
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

errorbars <- function(x,se,las=1,color.palette=1,margin=1/3,...) {
  # a dotplot with error bars.
  # if x has more than one column, multiple charts will be drawn,
  # side-by-side via split.screen().
  #layout(matrix(c(1,1,2),1,3))
  #layout.show(2)
  if(F) {
    mxy = linechart(t(x[1]),t(se[1]),labels=NULL,ylab=colnames(x)[1],
                    horizontal=T,dryrun=T,
                    las=las,color.pal=color.palette,type="p",legend=NULL,
                    xaxt="s",...)
    print(mxy)
    mar = mxy[[1]]; xlim = mxy[[2]]; ylim = mxy[[3]]
    opar = par(mar=mar,mfrow=c(1,length(x)+1))
    on.exit(par(opar))
    plot.new()
    plot.window(xlim=c(1,1),ylim=xlim)
    print(par("usr"))
    axis(2,at=-(1:nrow(x)),labels=rownames(x),las=las)
    box()
  }
  #opar = par(mfrow=c(1,length(x)))
  #on.exit(par(opar))
  #split.screen(c(1,length(x)))
  n = length(x)
  if(n > 1) {
    m = matrix(0,n,4)
    x1 = 0; x2 = margin + (1-margin)/n
    for(i in 1:length(x)) {
      m[i,1] = x1
      m[i,2] = x2
      m[i,3] = 0
      m[i,4] = 1
      x1 = x2
      x2 = x1 + (1-margin)/n
    }
    split.screen(m)
    on.exit(close.screen(all=TRUE))
  }
  for(i in 1:n) {
    if(n > 1) screen(i)
    linechart(t(x[i]),t(se[i]),labels=NULL,ylab=colnames(x)[i],horizontal=T,
              las=las,color.pal=color.palette,type="p",legend=NULL,
              xaxt=if(i == 1) "s" else "n",...)
  }
}
# Extensions and uses of loess

smooth <- function(x,y) {
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

