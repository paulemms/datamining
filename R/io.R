# functions for input and output

#' Read and write arrays
#'
#' Read data into an array from the console or file
#'
#' The array can have one of several formats.  The preferred format,
#'   produced by \code{\link{write.array}}, looks like \code{
#'     cvar.nam
#'     rvar.1.nam   ... rvar.k.nam    cvar.lev.1 ... cvar.lev.l
#'     rvar.1.lev.1 ... rvar.k.lev.1  ...        ... ...
#'   }
#'
#' Modified from rw1030/library/base/R/base/read.ftable
#' @param file either a character string naming a file or a connection which the data are to be read from or written to. "" indicates input from the console for reading and output to the console for writing.
#' @param sep the field separator string. Values on each line of the file are separated by this string.
#' @param quote the set of quoting characters as a single character string.
#' @param skip the number of lines of the input file to skip before beginning to read data.
#' @return An array.
#' @aliases write_array
#' @author Tom Minka
#' @examples
#' data(HairEyeColor)
#' fn <- tempfile()
#' write_array(HairEyeColor, fn)
#' read_array(fn)
#' @export
read_array <- function(file, sep = "", quote = "\"", skip = 0)
{
  z <- count.fields(file, sep, quote, skip)
  if(z[2] != max(z)) stop("unknown array representation")
  if(z[1] == 1) {
    ## Variable names and levels.  File looks like
    ##                                cvar.nam
    ## rvar.1.nam   ... rvar.k.nam    cvar.lev.1 ... cvar.lev.l
    ## rvar.1.lev.1 ... rvar.k.lev.1  ...        ... ...
    ##
    # read the header
    n.col.vars <- 1
    col.vars <- vector("list", length = n.col.vars)
    s <- scan(file, what = "", sep = sep, quote = quote,
              nlines = 2, skip = skip, quiet = TRUE)
    names(col.vars) <- s[1]
    s <- s[-1]
    n.row.vars <- z[max(which(z == max(z)))] - z[length(z)] + 1
    row.vars <- vector("list", length = n.row.vars)
    i <- 1:n.row.vars
    names(row.vars) <- s[i]
    col.vars[[1]] <- s[-i]
    z <- z[3:length(z)]
    skip <- skip + 2
  } else if(z[1] == min(z[2:length(z)])-1) {
    ## Variable levels, no names.  File looks like
    ##                                cvar.lev.1 ... cvar.lev.l
    ## rvar.1.lev.1 ... rvar.k.lev.1  ...        ... ...
    n.col.vars <- 1
    col.vars <- vector("list", length = n.col.vars)
    s <- scan(file, what = "", sep = sep, quote = quote,
              nlines = 1, skip = skip, quiet = TRUE)
    col.vars[[1]] <- s
    n.row.vars <- z[max(which(z == max(z)))] - z[length(z)] + 1
    row.vars <- vector("list", length = n.row.vars)
    z <- z[2:length(z)]
    skip <- skip + 1
  }
  else if(all(z == max(z))) {
    ## Just a block of numbers.
    values <- scan(file, sep = sep, quote = quote, quiet = TRUE, skip = skip)
    dim(values) <- c(z[1],length(z))
    return(t(values))
  } else {
    stop("unknown array representation")
  }
  # read the rest of the file
  s <- scan(file, what = "", sep = sep, quote = quote, quiet = TRUE,
            skip = skip)
  # select the data portion
  is.row.lab <- rep(rep(c(TRUE, FALSE), length(z)),
                    c(rbind(z - min(z) + 1, min(z) - 1)))
  values <- as.numeric(s[!is.row.lab])
  # label portion
  tmp <- s[is.row.lab]
  len <- length(tmp)
  p <- 1
  n <- integer(n.row.vars)
  for(k in seq(from = 1, to = n.row.vars)) {
    # n[k] is the cardinality of variable k
    n[k] <- sum(z >= max(z) - k + 1) / p
    # minka: added this line
    p <- n[k]
    # i indexes the lines
    i <- seq(from = 1, to = len, by = len / n[k])
    row.vars[[k]] <- unique(tmp[i])
    tmp <- tmp[seq(from = 2, to = len / n[k])]
    len <- length(tmp)
  }
  # minka: transposed from R version
  dim(values) <- c(sapply(col.vars, length),sapply(rev(row.vars), length))
  dimnames(values) <- c(col.vars,rev(row.vars))
  aperm(values)
}

# from rw1030/library/base/R/base/write.ftable
# modified to accept arrays
#' @export
write_array <- function(x, file = "", quote = TRUE,
                        digits = getOption("digits"))
{
  ox <- x
  charQuote <- function(s) {
    if(length(s) > 1) sapply(s,charQuote)
    else {
      if(!is.null(s) && length(grep(" ",s))>0)
        paste("\"", s, "\"", sep = "")
      else s
    }
  }
  makeLabels <- function(lst) {
    lens <- sapply(lst, length)
    cplensU <- c(1, cumprod(lens))
    cplensD <- rev(c(1, cumprod(rev(lens))))
    y <- NULL
    for (i in rev(seq(along = lst))) {
      ind <- 1 + seq(from = 0, to = lens[i] - 1) * cplensD[i + 1]
      tmp <- character(length = cplensD[i])
      tmp[ind] <- charQuote(lst[[i]])
      y <- cbind(rep(tmp, times = cplensU[i]), y)
    }
    y
  }
  dn <- dimnames(x)
  if(is.null(dn)) {
    x <- format(unclass(x),digits=digits)
    cat(t(x), file = file, sep = c(rep(" ", ncol(x) - 1), "\n"))
    return(invisible(ox))
  }
  xrv <- dn[1:(length(dn)-1)]
  xcv <- dn[length(dn)]
  LABS <- rbind(matrix("", nr = length(xcv), nc = length(xrv)),
                charQuote(names(xrv)),
                makeLabels(xrv))
  x <- aperm(x)
  dim(x) <- c(prod(sapply(xcv,length)),prod(sapply(xrv,length)))
  x <- t(x)
  DATA <- rbind(t(makeLabels(xcv)),
                format(unclass(x), digits = digits))
  if(!is.null(names(xcv))) {
    DATA <- rbind(c(charQuote(names(xcv)),rep("",ncol(x)-1)),DATA)
  }
  x <- cbind(apply(LABS, 2, format, justify = "left"),
             apply(DATA, 2, format, justify = "right"))
  cat(t(x), file = file, sep = c(rep(" ", ncol(x) - 1), "\n"))
  invisible(ox)
}
