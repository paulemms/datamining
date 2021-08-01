#' Cars, Housing, States, and VA Deaths
#'
#' Cleaned-up versions of some popular R datasets.
#'
#' The variables are given more informative names, some irrelevant variables
#' are dropped, erroneous values are fixed, and some variables are transformed
#' (in the standardized case).
#'
#' @name Cars
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







