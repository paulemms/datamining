

#' Cars, Housing, States, and VA Deaths
#' 
#' Cleaned-up versions of some popular R datasets.
#' 
#' The variables are given more informative names, some irrelevant variables
#' are dropped, erroneous values are fixed, and some variables are transformed
#' (in the standardized case).
#' 
#' @aliases Cars Housing States va.deaths
#' @return The specified dataset is loaded as a \code{\link{data.frame}}.  The
#' first three also define a transformed dataset (ending in `T') where the
#' variables are standardized.
#' @author Tom Minka
#' @seealso
#' \code{\link{Cars93}},\code{\link{Boston}},\code{\link{state}},\code{\link{VADeaths}}
#' @examples
#' 
#' data(States)
#' hist(States)
#' hist(StatesT)
#' w = pca(StatesT,2)
#' text.plot(project(StatesT,w),asp=1,cex=.6)
#' plot.axes(w)
#' 
NULL





#' Is a dimension ordered?
#' 
#' Query whether each array dimension is ordered or not.
#' 
#' Intuitively, a dimension is ordered if its values are in a natural
#' progression, such as (1,2,3) or ("one","two","three").  Programmatically, a
#' dimension is considered ordered if and only if the corresponding
#' \code{\link{dimnames}} vector has an \code{ordered} attribute set to
#' \code{T}.
#' 
#' @param x an array
#' @return A named logical vector, with an entry for each dimension, indicating
#' whether it is ordered or not.
#' @author Tom Minka
NULL



