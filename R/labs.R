# Functions for document and image retrieval
# By Tom Minka

# Functions for lab 1

strip_text <- function(txt) {
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

#' @export
read_doc <- function(fname, remove.header = TRUE) {
  txt <- readLines(fname)

  if(remove.header) {
    i <- which(txt == "")
    if(length(i) > 0) {
      txt <- txt[-(1:i[1])]
    }
  }
  strip_text(txt)
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
  doc.freq <- colMeans(div_by_sum(x))
  w <- sqrt(1/doc.freq)
  scale.cols(x,w)
}

#' @export
div_by_sum <- function(x) {
  scale.rows(x,1/(rowSums(x)+1e-16))
}

#' @export
div_by_euc_length <- function(x) {
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

#' @export
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
