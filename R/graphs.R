# Routines for manipulating directed graphs

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

