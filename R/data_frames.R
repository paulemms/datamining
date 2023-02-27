# data.frame codes

# put data into the environment of formula
# from Brian Ripley
formula_with_data <- function(fmla, data) {
  if(identical(as.character(fmla[[3]]),".")) {
    # dot must be expanded
    resp = response.var(fmla)
    pred = setdiff(names(data),resp)
    fmla = formula(paste(resp,"~",paste(pred,collapse="+")))
  }
  env <- new.env()
  for(i in names(data)) assign(i, data[[i]], envir=env)
  environment(fmla) <- env
  fmla
}


#' Sort rows of a data frame
#'
#' Re-orders data frame rows so that a column is sorted.
#' @param x a data frame
#' @param f the name or number of a column of \code{x}, or a numeric vector
#'     to play the role of a column of \code{x}. Default: Last column.
#' @return
#'   A reordered version of \code{x} where the values in column \code{f} are in
#'   ascending order.
#' @author Tom Minka
#' @seealso
#'   \code{\link{sort_cells}}
#' @examples
#' data(mtcars)
#' sort(mtcars, f = "mpg")
#' @exportS3Method
sort.data.frame <- function(x, decreasing=FALSE, f=ncol(x), ...) {
  if(length(f) == 1) f <- x[[f]]
  else if(length(f) != nrow(x))
    stop("length of sort vector doesn't match nrows")
  x[order(f, decreasing = decreasing, ...),,drop=F]
}


#' Sort an array
#'
#' Sort the cells of an array.
#' @param x an array
#' @details Converts `x` to a data frame and calls
#'   \code{\link{sort.data.frame}}.
#' @return A data frame representation of `x`, sorted by cell value.
#' @author Tom Minka
#' @seealso [sort.data.frame()]
#' @examples
#' data(Titanic)
#' sort_cells(Titanic)
#' data(HairEyeColor)
#' sort_cells(HairEyeColor)
#' @export
sort_cells <- function(x) {
  sort.data.frame(as.data.frame(x))
}

# in R 4.1.1 can use list2DF - only 2 refs here but best leave
# e.g. list2DF(list(a = character(0), b = character(0)))
# use list2DF if you want the data.frame populated
empty_data_frame <- function(col_names) {
  stopifnot(is.character(col_names))
  if (missing(col_names)) y <- list()
  else {
    y <- sapply(col_names, function(x) NULL)
  }
  structure(y, class="data.frame")
}

rbind_extend2 <- function(df,df2) {
  # like rbind, but allows different numbers of columns (NAs inserted)
  # this version loops through all the column names and creates a new data.frame
  x <- list()
  col.names <- unique(c(names(df),names(df2)))
  for(j in 1:length(col.names)) {
    k <- col.names[j]
    # must check in names first to avoid abbreviation match
    if(k %in% names(df)) v1 <- df[[k]]
    else v1 <- factor(rep(NA,nrow(df)))
    if(k %in% names(df2)) v2 <- df2[[k]]
    else v2 <- factor(rep(NA,nrow(df2)))
    # requires cfactor
    # force dispatch on second argument too
    x[[j]] <- cfactor(v1,v2)
  }
  names(x) <- col.names
  row.names <- make.unique(c(rownames(df),rownames(df2)))
  attr(x,"row.names") <- row.names
  class(x) <- "data.frame"
  x
}

rbind_extend <- function(df,df2) {
  # like rbind, but allows different numbers of columns (NAs inserted)
  # this version modifies the data.frames with new NA columns if necessary
  if(nrow(df) == 0) return(df2)
  if(nrow(df2) == 0) return(df)
  not.1 <- setdiff(names(df2),names(df))
  not.2 <- setdiff(names(df),names(df2))
  if(length(not.1) > 0) df[not.1] <- rep(NA,nrow(df))
  if(length(not.2) > 0) df2[not.2] <- rep(NA,nrow(df2))
  col.names <- unique(c(names(df),names(df2)))
  x <- rbind(df[col.names],df2[col.names])
  row.names <- make.unique(c(rownames(df),rownames(df2)))
  attr(x,"row.names") <- row.names
  x
}

# not used
cbind_extend <- function(df,df2) {
  # like cbind.data.frame, but pads with NA until rownames match
  if(is.null(df) || nrow(df) == 0) return(df2)
  if(is.null(df2) || nrow(df2) == 0) return(df)
  not.1 = setdiff(rownames(df2),rownames(df))
  not.2 = setdiff(rownames(df),rownames(df2))
  if(length(not.1) > 0) df[not.1,] = NA
  if(length(not.2) > 0) df2[not.2,] = NA
  r = cbind(df,df2[rownames(df),])
  names(r) = c(names(df),names(df2))
  r
}


# Ripley suggests to get response
# https://stat.ethz.ch/pipermail/r-help/2007-November/145649.html


#' Get response variable
#'
#' @param obj formula, model frame, or fitted model
#'
#' @return character string
#'
#' @examples
#' response_var(a ~ b + c)
#' # response_var(~ b + c) # gives error
#' response_var(Height + Volume ~ Girth)
#' response_var(log(Volume) ~ Girth)
#' @export
response_var <- function(obj) {
  f <- formula(obj)
  if(length(f) < 3) stop("no response")
  resp <- f[[2]]
  deparse(resp)
}

#' TODO: remove this
#' Get response variable
#'
#' This works for
#' @param object formula, model frame, or fitted model
#'
#' @return the name of the response variable
#' @export
#'
#' @examples
#' response.var(a ~ b + c)
#' response.var(~ b + c) # gives '.'
#' response.var(Height + Volume ~ Girth) # incorrect
#' response.var(log(Volume) ~ Girth) # incorrect
#'
response.var <- function(object) {
  if(is.null(object) || is.array(object)) return(NULL)
  if(length(class(object)) == 1 && class(object) == "formula") return(all.vars(update(object, . ~ NULL))[1])
  if(inherits(object, "terms")) {
    a <- attributes(object)
    if(!a$response) return(character(0))
    return(as.character(a$variables[2]))
  }
  if(is.null(attr(object,"terms"))) {
    if(is.data.frame(object)) {
      # shortcut to avoid make.names
      return(names(object)[length(object)])
    }
    if(is.table(object)) return(response.var(terms(object)))
    #if(is.list(object) || is.vector(object) || is.array(object)) return(NULL)
  }
  response.var(terms(object))
}

#' Variables out of which the predictors are constructed
#'
#' If the object inherits from a `data.frame` then it is assumed the
#' response is in the first column.
#' @param object that inherits from a data.frame or terms object
#'
#' @return vector of predicts
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ .^2, data=mtcars)
#' predictor_vars(fit)
#' predictor_vars(mtcars)
predictor_vars <- function(object) {
  if(inherits(object, "terms")) {
    a <- attributes(object)
    pred = rownames(a$factors)[apply(a$factors, 1, any)]
    return(pred)
  }
  if(inherits(object, "data.frame") && is.null(attr(object, "terms"))) {
    # shortcut to avoid make.names
    return(names(object)[-1])
  }
  predictor.vars(terms(object))
}


# returns the variables out of which the predictors are constructed
predictor.vars <- function(object) {
  if(inherits(object, "terms")) {
    a <- attributes(object)
    pred = rownames(a$factors)[apply(a$factors,1,any)]
    # skip cross-products
    #pred = a$term.labels[a$order == 1]
    # drop I() terms
    #pred = pred[substring(pred,1,2) != "I("]
    return(pred)
  }
  if(inherits(object,"data.frame") && is.null(attr(object,"terms"))) {
    # shortcut to avoid make.names
    return(names(object)[-length(object)])
  }
  predictor.vars(terms(object))
}

# returns all terms on the rhs, including higher-order terms
predictor.terms <- function(object) {
  attr(terms(object),"term.labels")
}

as.data.frame.col <- function(x,n="V1") {
  frame = as.data.frame(as.matrix(x))
  names(frame) = n
  frame
}
as.data.frame.row <- function(x,row.name="") {
  extra.args <- list(check.names=F,row.names=row.name)
  do.call("data.frame",append(as.list(x),extra.args))
}

# nocheck.data.frame <- function(...,row.names=NULL) {
#   # must be done this way to prevent check.names (bug in data.frame)
#   # example: data.frame(list("a b"=3),check.names=F)
#   do.call("data.frame",append(as.list(...),list(row.names=row.names,check.names=F)))
# }

apply.df <- function(x, fun) {
  # apply a function to a data.frame and return a data.frame
  y <- lapply(x, fun)
  if(length(y[[1]]) == length(x[[1]]))
    data.frame(y, row.names=rownames(x), check.names=FALSE)
  else
    data.frame(y, check.names=FALSE)
}

my.model.frame <- function(...) {
  # bugfix for model.frame - puts response at end
  x = model.frame(...)
  a = attr(x,"terms")
  x = data.frame(c(x[-1],x[1]),check.names=F,row.names=rownames(x))
  attr(x,"terms") = a
  x
}

set.response.var <- function(fmla,v) {
  # v is character string
  vars = c(response.var(fmla),predictor.terms(fmla))
  i = grep(v,vars)
  if(length(i)>0) vars = vars[-i[1]]
  formula(paste(v,paste(vars,collapse="+"),sep="~"))
}

# concatenate factors (should be built in)
# unfortunately "sort" requires the original c(), which drops labels
cfactor <- function(f1,f2=NULL) {
  if(length(f1) == 0) return(f2)
  else if(length(f2) == 0) return(f1)
  if(!is.factor(f1) || !is.factor(f2)) return(c(f1,f2))
  all.levs <- unique(c(levels(f1),levels(f2)))
  factor(c(as.character(f1),as.character(f2)),levels=all.levs)
}

# should be built into R
sum.data.frame <- function(x,...) sapply(x,sum,...)

RVersionString <- function() {
  ver = R.Version()
  paste(ver$major,ver$minor,sep=".")
}

# bug fixes

which.environment <- function(x,...) {
  env = environment(...)
  while(!is.null(env) && !exists(x,envir=env,inherits=F)) {
    env = parent.env(env)
  }
  env
}

if(!exists("assignInNamespace")) {
  assignInNamespace <- function(subx,x,ns) {
    # if(!exists("asNamespace")) return(assign(subx,x,env=ns))
    # # from fixInNamespace
    # if (bindingIsLocked(subx, ns)) {
    #   unlockBinding(subx, ns)
    #   assign(subx, x, env = ns)
    #   w <- options("warn")
    #   on.exit(options(w))
    #   options(warn = -1)
    #   lockBinding(subx, ns)
    # }
    # else assign(subx, x, env = ns)
    # if (isNamespace(ns) && !isBaseNamespace(ns)) {
    #   S3 <- getNamespaceInfo(ns, "S3methods")
    #   if (!length(S3))
    #     return(invisible(NULL))
    #   S3names <- sapply(S3, function(x) x[[3]])
    #   if (subx %in% S3names) {
    #     i <- match(subx, S3names)
    #     genfun <- get(S3[[i]][[1]])
    #     defenv <- if (typeof(genfun) == "closure")
    #       environment(genfun)
    #     else .BaseNamespaceEnv
    #     S3Table <- get(".__S3MethodsTable__.", envir = defenv)
    #     if (exists(subx, envir = S3Table, inherits = FALSE))
    #       assign(subx, x, S3Table)
    #   }
    # }
  }
}

my.aggregate.data.frame <- function(x, by, FUN, ...) {
  if(!is.data.frame(x))
    x <- as.data.frame(x)
  if(!is.list(by))
    stop("`by' must be a list")
  if(is.null(names(by)))
    names(by) <- paste("Group", seq(along = by), sep = ".")
  else {
    nam <- names(by)
    ind <- which(nchar(nam) == 0)
    names(by)[ind] <- paste("Group", ind, sep = ".")
  }
  y <- lapply(x, tapply, by, FUN, ..., simplify = FALSE)
  if(any(sapply(unlist(y, recursive = FALSE), length) > 1))
    stop("`FUN' must always return a scalar")
  z <- y[[1]]
  d <- dim(z)
  w <- list()
  for (i in seq(along = d)) {
    j <- rep(rep(seq(1 : d[i]),
                 prod(d[seq(length = i - 1)]) * rep(1, d[i])),
             prod(d[seq(from = i + 1, length = length(d) - i)]))
    #w <- cbind(w, dimnames(z)[[i]][j])
    # minka: preserve original levels
    f <- factor(dimnames(z)[[i]][j],levels=levels(by[[i]]))
    w[[names(by)[i]]] = f
  }
  # minka: not sure what this is for
  #w <- w[which(!unlist(lapply(z, is.null))), ]
  y <- data.frame(w, lapply(y, unlist, use.names = FALSE))
  #names(y) <- c(names(by), names(x))
  y
}
replaceInNamespace("aggregate.data.frame",my.aggregate.data.frame)

#if(!exists("make.unique")) {
# desiderata for make.unique:
# 1. if A is unique, make.unique(c(A,B)) preserves A
# 2. make.unique(c(A,B)) = make.unique(c(make.unique(A),B))

# internal version
# does not satisfy desideratum #2
# make.unique(c("a","a","a")) != make.unique(c(make.unique(c("a","a")),"a"))


#' Make character strings unique
#'
#' Makes the elements of a character vector unique by appending sequence
#' numbers to duplicates.
#'
#' The algorithm used by \code{make.unique} has the property that
#' \code{make.unique(c(A,B)) = make.unique(c(make.unique(A),B))}.
#'
#' In other words, you can append one string at a time to a vector, making it
#' unique each time, and get the same result as applying \code{make.unique} to
#' all of the strings at once.
#'
#' If character vector \code{A} is already unique, then
#' \code{make.unique(c(A,B))} preserves \code{A}.
#'
#' @param names a character vector
#' @param sep a character string used to separate a duplicate name from its
#' sequence number.
#' @return a character vector of same length as \code{names} with duplicates
#' changed.
#' @author Thomas P Minka
#' @seealso \code{\link{make.names}}
#' @keywords character
#' @examples
#'
#' make.unique(c("a","a","a"))
#' make.unique(c(make.unique(c("a","a")),"a"))
#'
#' make.unique(c("a","a","a.2","a"))
#' make.unique(c(make.unique(c("a","a")),"a.2","a"))
#'
#' rbind(data.frame(x=1),data.frame(x=2),data.frame(x=3))
#' rbind(rbind(data.frame(x=1),data.frame(x=2)),data.frame(x=3))
#'
make.unique <- function(names) {
  # names is character vector
  if(!is.character(names)) stop("names must be a character vector")
  while(any(dups <- duplicated(names))) {
    names[dups] <- paste(names[dups], seq(length = sum(dups)), sep = "")
  }
  names
}
# satifies both desiderata
# make.unique(c("a","a","a.2","a")) ==
# make.unique(c(make.unique(c("a","a")),"a.2","a"))
make.unique <- function(names,sep=".",start=2) {
  # names is character vector
  if(!is.character(names)) stop("names must be a character vector")
  repeat {
    dups <- which(duplicated(names))
    if(length(dups) == 0) break
    # loop duplicates
    for(j in dups) {
      # for each duplicate, find the lowest value of cnt which makes it
      # different from previous names.
      i <- start
      repeat {
        newnam <- paste(names[j],i,sep=sep)
        # compare to previous elements only
        if(!any(is.element(newnam,names[1:(j-1)]))) break
        i <- i + 1
      }
      names[j] <- newnam
    }
    # repeat in case new duplicates have been created (see examples)
  }
  names
}



#' Make Syntactically Valid Names
#'
#' Make syntactically valid names out of character vectors.
#'
#' A syntactically valid name consists of letters, numbers, and the dot
#' character and starts with a letter or the dot.
#'
#' All invalid characters are translated to \code{"."}.  A missing value is
#' translated to \code{"NA"}.  Names which match R keywords have a dot appended
#' to them.
#'
#' If \code{unique = TRUE} then \code{\link{make.unique}} is used to append
#' sequence numbers to duplicates (after coercion).
#'
#' @param names character (vector) to be coerced to syntactically valid names.
#' @param unique logical; if \code{TRUE}, the resulting elements are unique.
#' This may be desired for, e.g., column names.
#' @return A character vector of same length as \code{names} with each changed
#' to a syntactically valid name.
#' @seealso \code{\link{names}}, \code{\link{character}},
#' \code{\link{data.frame}}.
#' @keywords character
#' @examples
#'
#' make.names(c("a and b", "a_and_b"), unique=TRUE)
#' # "a.and.b"  "a.and.b1"
#'
#' data(state)
#' state.name[make.names(state.name) != state.name]# those 10 with a space
#'
make.names <- function(names, unique=F) {
  names <- .Internal(make.names(as.character(names), unique))
  # minka: change keywords
  i <- is.element(names, c("for","while","repeat","if","else","function"))
  if(any(i)) names[i] <- paste(names[i],".",sep="")
  if(unique) names <- make.unique(names)
  names
}
replaceInNamespace("make.names",make.names)
#}


#
#  just rearrange data.frame and use stats functions?
#
#' # does NOT throw an error if no terms
#' #' @export
#' terms <- function(...) UseMethod("terms")
#'
#' #' @export
#' terms.default <- function(x,...) x$terms
#'
#' @exportS3Method
terms.data.frame <- function(x,env=parent.frame(),...) {
  fmla <- attr(x,"terms")
  if(is.null(fmla)) {
    # minka: assume the last column is the response
    #nm <- make.names(names(x))
    nm = names(x)
    if(length(nm) > 1) {
      lhs <- nm[length(nm)]
      rhs <- nm[-length(nm)]
    }
    else {
      lhs <- NULL
      rhs <- nm
    }
    fmla <- terms(formula(paste(lhs,"~",paste(rhs,collapse="+")),env=env,...))
  }
  fmla
}

#' #' @export
#' terms.table <- function(x,...) {
#'   terms.data.frame(as.data.frame(x),...)
#' }
#'
#' #' @export
#' formula <- function(...) UseMethod("formula")
#'
#' #' @export
#' formula.default <- function(x,env=parent.frame(), ...)
#' {
#'   if (!is.null(x$formula))		eval(x$formula)
#'   else if (!is.null(x$call$formula))	eval(x$call$formula)
#'   # minka: always return formula, not terms
#'   else if (!is.null(x$terms))		formula.terms(x$terms)
#'   else if (!is.null(attr(x, "formula"))) attr(x, "formula")
#'   else {form<-switch(mode(x),
#'                      NULL = structure(NULL, class = "formula"),
#'                      character = formula(eval(parse(text = x)[[1]])),
#'                      call = eval(x), stop("invalid formula"))
#'   environment(form)<-env
#'   form
#'   }
#' }
#'
#' #' @export
#' formula.data.frame <- function(x,env=parent.frame(),...) {
#'   formula(terms(x,env=env),...)
#' }

# Ripley says the original update is not broken (just non-modular)
# related problem:
# http://maths.newcastle.edu.au/~rking/R/help/02b/1318.html
my.update <-
  function (object, formula., ..., evaluate = TRUE)
  {
    call <- object$call
    if (is.null(call))
      stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
      call$formula <- update.formula(formula(object), formula.)
    if(length(extras) > 0) {
      existing <- !is.na(match(names(extras), names(call)))
      ## do these individually to allow NULL to remove entries.
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if(any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }
    if(evaluate) {
      # minka: use environment of formula instead of parent.frame
      # see the man page for formula
      env<-environment(call$formula)
      if (is.null(env)) env<-parent.frame()
      eval(call,env)
    }
    else call
  }

plot.default <- function(x, y=NULL, type="p", xlim=NULL, ylim=NULL,
                         log="", main=NULL, sub=NULL, xlab=NULL, ylab=NULL,
                         ann=par("ann"), axes=TRUE, frame.plot=axes,
                         panel.first=NULL, panel.last=NULL,
                         col=par("col"), bg=NA, pch=par("pch"),
                         cex = 1, lty=par("lty"), lab=par("lab"),
                         lwd=par("lwd"), asp=NA, ...)
{
  xlabel <- if (!missing(x)) deparse(substitute(x))
  ylabel <- if (!missing(y)) deparse(substitute(y))
  # minka
  # treat character as factor
  if(mode(y) == "character") y <- factor(y)
  if(mode(x) == "character") x <- factor(x)
  # this works even if x and y are factors
  xy <- xy.coords(x, y, xlabel, ylabel, log)
  if(all(is.na(xy$x))) xy$x <- 1:length(x)
  xlab <- if (is.null(xlab)) xy$xlab else xlab
  ylab <- if (is.null(ylab)) xy$ylab else ylab
  xlim <- if (is.null(xlim)) range(xy$x[is.finite(xy$x)]) else xlim
  ylim <- if (is.null(ylim)) range(xy$y[is.finite(xy$y)]) else ylim
  plot.new()
  plot.window(xlim, ylim, log, asp, ...)
  panel.first
  plot.xy(xy, type, col=col, pch=pch, cex=cex, bg=bg, lty=lty, lwd=lwd, ...)
  panel.last
  if (axes) {
    # minka
    if(is.factor(x)) {
      lev <- levels(x)
      axis(1, at=1:length(lev), labels=lev, ...)
    } else axis(1, ...)
    # minka
    if(is.factor(y)) {
      lev <- levels(y)
      axis(2, at=1:length(lev), labels=lev, ...)
    } else axis(2, ...)
  }
  if (frame.plot)
    box(...)
  if (ann)
    title(main=main, sub=sub, xlab=xlab, ylab=ylab, ...)
  invisible()
}

model.frame.multinom <- function(object) {
  oc <- object$call
  oc[[1]] <- NULL
  do.call("model.frame",as.list(oc))
}

model.frame.glm <-
  function (object, ...)
  {
    dots = list(...)
    nargs = dots[match(c("data", "na.action", "subset"), names(dots), 0)]
    if (any(nargs > 0) || is.null(object$model)) {
      fcall <- object$call
      fcall$method <- "model.frame"
      #fcall[[1]] <- as.name("glm")
      # minka
      if (!is.null(object$terms))
        env <- environment(object$terms)
      else
        env <- environment(fcall$formula)
      if (is.null(env))
        env <- parent.frame()
      eval(fcall, env)
    }
    else object$model
  }
model.frame.nnet <- model.frame.glm
# missing from modreg
model.frame.loess <- model.frame.glm

# bugfix

if(compareVersion(RVersionString(),"1.9.0") < 0) {
  my.loess <-
    function(formula, data=NULL, weights, subset, na.action, model = FALSE,
             span = 0.75, enp.target, degree = 2, parametric = FALSE,
             drop.square = FALSE, normalize = TRUE,
             family = c("gaussian", "symmetric"),
             method = c("loess", "model.frame"),
             control = loess.control(...), ...)
    {
      family <- match.arg(family)
      method <- match.arg(method)
      mt <- terms(formula, data = data)
      mf <- match.call(expand.dots=FALSE)
      mf$model <- mf$span <- mf$enp.target <- mf$degree <-
        mf$parametric <- mf$drop.square <- mf$normalize <- mf$family <-
        mf$control <- mf$... <- NULL
      # minka: bugfix
      mf$method <- NULL
      mf[[1]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
      if (match.arg(method) == "model.frame") return(mf)
      na.act <- attr(mf, "na.action")
      y <- model.response(mf, "numeric")
      w <- model.weights(mf)
      if(is.null(w)) w <- rep(1, length(y))
      nmx <- as.character(attr(mt, "variables"))[-(1:2)]
      x <- mf[, nmx, drop=FALSE]
      if(any(sapply(x, is.factor))) stop("predictors must all be numeric")
      x <- as.matrix(x)
      D <- ncol(x)
      nmx <- colnames(x)
      names(nmx) <- nmx
      drop.square <- match(nmx, nmx[drop.square], 0) > 0
      parametric <- match(nmx, nmx[parametric], 0) > 0
      if(!match(degree, 0:2, 0)) stop("degree must be 0, 1 or 2")
      iterations <- if(family=="gaussian") 1 else control$iterations
      if(!missing(enp.target))
        if(!missing(span))
          warning("both span and enp.target specified: span will be used")
      else {				# White book p.321
        tau <- switch(degree+1, 1, D+1, (D+1)*(D+2)/2) - sum(drop.square)
        span <- 1.2 * tau/enp.target
      }
      fit <- simpleLoess(y, x, w, span, degree, parametric, drop.square,
                         normalize, control$statistics, control$surface,
                         control$cell, iterations, control$trace.hat)
      fit$call <- match.call()
      fit$terms <- mt
      fit$xnames <- nmx
      fit$x <- x
      fit$y <- y
      fit$weights <- w
      if(model) fit$model <- mf
      if(!is.null(na.act)) fit$na.action <- na.act
      fit
    }
} else {
  my.loess = loess
}
# minka: robust by default
formals(my.loess)$family = "symmetric"

my.predict.loess <- function(object, newdata = NULL, se = FALSE, ...)
{
  if(!inherits(object, "loess"))
    stop("First argument must be a loess object")
  if(is.null(newdata) & (se == FALSE)) return(fitted(object))

  if(is.null(newdata)) newx <- object$x
  else {
    vars <- as.character(attr(delete.response(terms(object)),
                              "variables"))[-1]
    if(length(vars) > 1 || NCOL(newdata) > 1) {
      if(any(!match(vars, colnames(newdata), FALSE)))
        # minka
        newdata = model.frame.default(delete.response(terms(object)),newdata)
      else
        newdata = newdata[,vars,drop=F]
    }
    newx = as.matrix(newdata)
  }
  res <- predLoess(object$y, object$x, newx, object$s, object$weights,
                   object$pars$robust, object$pars$span, object$pars$degree,
                   object$pars$normalize, object$pars$parametric,
                   object$pars$drop.square, object$pars$surface,
                   object$pars$cell, object$pars$family,
                   object$kd, object$divisor, se=se)
  if(se)
    res$df <- object$one.delta^2/object$two.delta
  res
}

replaceInNamespace("loess",my.loess)
environment(my.predict.loess) = environment(loess)
#assignInNamespace("predict.loess",my.predict.loess,environment(loess))

# bug fix

# minka: added `key' and `main' arguments
my.filled.contour <-
  function (x = seq(0, 1, len = nrow(z)),
            y = seq(0, 1, len = ncol(z)),
            z,
            xlim = range(x, finite=TRUE),
            ylim = range(y, finite=TRUE),
            zlim = range(z, finite=TRUE),
            levels = pretty(zlim, nlevels), nlevels = 20,
            color.palette = cm.colors,
            col = color.palette(length(levels) - 1),
            plot.title, plot_axes, key.title, key.axes, key=T,
            asp = NA, xaxs="i", yaxs="i", las = 1, axes = TRUE,
            frame.plot = axes, main="", ...)
  {
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq(0, 1, len = nrow(z))
        }
      }
      else stop("no `z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing x and y values expected")

    if(key) {
      mar.orig <- (par.orig <- par(c("mar","las","mfrow")))$mar
      on.exit(par(par.orig))
      par(las = las)

      w <- (3 + mar.orig[2]) * par('csi') * 2.54
      layout(matrix(c(2, 1), nc=2), widths=c(1, lcm(w)))
      par(las = las)
      ## Plot the `plot key' (scale):
      mar <- mar.orig
      mar[4] <- mar[2]
      mar[2] <- 1
      par(mar = mar)
      plot.new()
      plot.window(xlim=c(0,1), ylim=range(levels), xaxs="i", yaxs="i")
      rect(0, levels[-length(levels)], 1, levels[-1], col = col)
      if (missing(key.axes)) {
        if (axes)
          axis(4)
      }
      else key.axes
      box()
      if (!missing(key.title))
        key.title
      mar <- mar.orig
      mar[4] <- 1
      par(mar=mar)
    }
    plot.new()
    plot.window(xlim, ylim, "", xaxs=xaxs, yaxs=yaxs, asp=asp)

    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
      stop("no proper `z' matrix specified")
    if (!is.double(z))
      storage.mode(z) <- "double"
    .Internal(filledcontour(as.double(x),
                            as.double(y),
                            z,
                            as.double(levels),
                            col = col))
    if (axes) {
      if (missing(plot_axes)) {
        # minka: no title here
        axis(1)
        axis(2)
      } else plot_axes
      if (missing(plot.title))
        title(main=main,...)
      else
        plot.title
    }
    if(frame.plot) box()
    invisible()
  }
replaceInNamespace("filled.contour",my.filled.contour)


# Routines for contingency tables
# Tom Minka

barplot.table <- function(m,col=default.colors.w(nrow(m)),
                          xlab,main,
                          ylab=response.var(m),legend.text=T,...) {
  nam <- names(dimnames(m))
  opar <- par(mar=auto.mar(main=nam[1]))
  on.exit(par(opar))
  if(missing(xlab)) xlab = nam[2]
  if(missing(main)) main = nam[1]
  barplot.default(m,legend.text=legend.text,col=col,...,
                  xlab=xlab,ylab=ylab,main=main)
}

flipud <- function(m) {
  m[nrow(m):1,]
}
fliplr <- function(m) {
  m[,ncol(m):1]
}

image.table <- function(m,col=YR.colors(64),bg=gray(0.5),mar,mar.scale=NULL,
                        ...) {
  m <- flipud(m)
  b = break.quantile(as.vector(m),length(col))
  #b = break.equal(as.vector(m),length(col))
  if(is.null(colnames(m))) colnames(m) = 1:ncol(m)
  if(is.null(rownames(m))) rownames(m) = 1:nrow(m)
  nam <- names(dimnames(m))
  if(missing(mar)) mar = auto.mar(xlab=nam[2],ylab=nam[1],mar.scale=mar.scale,...)
  opar = par(bg=bg,mar=mar)
  on.exit(par(opar))
  image.default(t(m),axes=F,col=col,breaks=b,...)
  box()
  # axis 1 is horizontal
  axis(1,at=seq(0,1,len=ncol(m)),labels=colnames(m),...)
  # axis 2 is vertical
  axis(2,at=seq(0,1,len=nrow(m)),labels=rownames(m),...)
  title(xlab=nam[2],ylab=nam[1])
}

if(!exists('pie.default')) {
  pie.default = pie
  pie = function(x,...) UseMethod("pie")
}

# draws a pie for each column
pie.table <- function(m,layout,col=default.colors.w(nrow(m)),...) {
  if(length(dim(m)) == 1) return(pie.default(m,col=col,...))
  if(missing(layout)) layout <- auto.layout(ncol(m))
  opar <- par(mfrow=layout,mar=c(1,1,1,0)*2)
  on.exit(par(opar))
  nam <- names(dimnames(m))
  for(i in 1:ncol(m)) {
    pie(m[,i],col=col,...)
    title(paste(nam[1],"\n",nam[2],"=",colnames(m)[i]))
  }
}


#' Sort rows and columns of a contingency table
#'
#' Score rows and columns via correspondence analysis, then sort
#' the scores.
#'
#' Only for 2D tables right now.
#' @param x a contingency table
#' @seealso \code{\link{mosaicplot}}
#' @examples
#' data(blood)
#' sort_table(blood)
#'
#' @export
sort_table <- function(x,k=NULL) {
  if(!is.null(k)) {
    if(k == 1) {
      # sort rows only
      n <- margin.table(x,k)+1e-15
      v2 <- 1:ncol(x)
      v1 <- drop(x %*% v2)/n
      return(x[order(v1),])
    }
    if(k == 2) {
      # sort cols only
      n <- margin.table(x,k)+1e-15
      v1 <- 1:nrow(x)
      v2 <- drop(v1 %*% x)/n
      return(x[,order(v2)])
    }
  }

  n1 <- margin.table(x,1)
  n2 <- margin.table(x,2)
  n <- sum(n1)
  if(n == 0) return(x)

  # iterate
  v2 <- rnorm(ncol(x))
  iter <- 1
  repeat {
    old.v2 <- v2
    # update rows
    # add eps in case there are zeros
    v1 <- drop(x %*% v2)/(n1+eps)
    # shift to zero mean
    v1 <- v1 - sum(v1*n1)/n
    # scale to unit sd
    v1 <- v1/sqrt(sum(v1^2*n1)/n)
    # update cols
    v2 <- drop(v1 %*% x)/(n2+eps)
    if(max(abs(v2 - old.v2)) < 1e-8) break
    iter <- iter + 1
  }
  # remove sign ambiguity
  nz <- which(n1 > 0)[1]
  s <- sign(v1[nz])
  v1 <- v1*s
  v2 <- v2*s

  # empty rows/cols get placed at end
  z1 <- which(n1 == 0)
  z2 <- which(n2 == 0)
  v1[z1] <- 1e10
  v2[z2] <- 1e10
  i1 <- order(v1)
  i2 <- order(v2)
  y <- x[i1,i2]
  class(y) <- class(x)
  y
  #list(v1=v1,v2=v2)
}

indep.fit <- function(x) {
  ## Compute the expected cell counts under the independence model.
  loglin(x, as.list(1:length(dim(x))), fit = TRUE, print = FALSE)$fit
}


#' Associations in a contingency table
#'
#' Find unusually frequent variable combinations.
#'
#' Enumerates all two-way marginal tables, sorts the cells by lift, and returns
#' the top.  The formula for lift is \deqn{\frac{p(i,j)}{p(i)p(j)}}
#'
#' @param x a contingency table
#' @param top the number of top associations to return
#' @param targets a character vector of table variables, one of which must be
#' included in every returned association.
#' @return A data frame where the first column is \code{Lift} and the remaining
#' columns correspond to the variables of \code{x}. Each row describes a cell
#' of \code{x}, where marginalized variables have value \code{NA}. The lift
#' value describes the ratio of the actual count in the cell versus the
#' expected count under independence.
#' @author Tom Minka
#' @export
#' @examples
#'
#' data(Titanic)
#' mine.associations(Titanic)
#' # Females are twice as likely to survive as the average person.
#' # Members of 3rd class are twice as likely to be children as the average person.
#' # Etc.
#'
#' # focus on associations with survival
#' mine.associations(Titanic,target="Survived")
#'
#' # focus on children
#' mine.associations(Titanic[,,1,],target="Survived")
#'
mine.associations <- function(x,top=10,targets=NULL,z=1) {
  if(!inherits(x,"class")) x <- as.table(x)
  dm <- dim(x)
  nd <- length(dm)
  dn <- attr(dimnames(x),"names")
  if(is.null(targets)) targets <- 1:nd
  else if(is.character(targets)) {
    targets <- pmatch(targets,dn)
    if(any(is.na(targets))) stop("unrecognized target")
  }
  not.targets <- setdiff(1:nd, targets)

  # create a data frame with 0 rows
  res <- empty_data_frame(c("Lift",dn))

  # loop rows
  for(i in targets) {
    predictors <- c(not.targets, targets[targets > i])
    for(j in predictors) {
      y <- margin.table(x, c(i,j))
      df <- sort_cells((y - z*sqrt(y))/(indep.fit(y)+eps))
      names(df)[3] <- "Lift"
      # take top k lifts
      if(nrow(df) > top) df <- df[nrow(df)+1-(top:1),]
      # extend cols
      for(v in 1:nd) {
        if(v != i && v != j) df[[dn[v]]] <- factor(rep(NA,nrow(df)))
      }
      res <- rbind_extend(res,df)
      # ensure unique rownames
      #rownames(res) <- nrow(res):1
    }
  }
  res <- sort.data.frame(res, f = 1)
  if(nrow(res) > top) res <- res[nrow(res)+1-(top:1),]
  rownames(res) <- nrow(res):1
  res
}
