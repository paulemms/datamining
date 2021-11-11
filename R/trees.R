# extra routines for using regression trees
# requires: tree package, cluster.r

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
