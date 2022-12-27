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
#' dev.new()
#' hist(States)
#' hist(StatesT)
#' w <- pca(StatesT, 2)
#' text_plot(project(StatesT, w), asp=1, cex=.6)
#' plot_axes(w)
#'
NULL


#' U.S. state area 1990
#' @name state.area
#' @details
#'   \code{state.area} is a named vector of the 1990 U.S. state land areas,
#'   including the District of Columbia, in square miles.
#'   The numbers are more accurate and up-to-date than the ones in
#'   \code{data(state)}.
#' @seealso
#'   \code{\link{state.pop}}, \code{\link{state}}
#' @references
#'   The data was obtained from
#'   \url{http://www.census.gov/population/censusdata/90den_stco.txt}
#'   and edited into R format.
#'   County areas can also be obtained from the above link.
NULL


#' U.S. state population 1790-1990
#' @name state.pop
#' @details
#'   \code{state.pop} is a matrix containing the population of
#'   each U.S. state, including the District of Columbia,
#'   from 1790 to 1990 in ten-year increments.
#'   The rows are the years "1790" to "1990" and the columns are state names.
#' @seealso
#'   \code{\link{state.area}}, \code{\link{state.name}}, \code{\link{state}}
#' @references
#'   The data was obtained from
#'   \url{http://merrill.olm.net/pdocs/feas/pop/pop1790_1990/pii.txt}
#'   and edited into R format.
#' @examples
#'   dotchart(sort(log(state.pop["1990",])))
#'
#'   # distribution of growth rates over time
#'   rates <- diff(log(state.pop))
#'   boxplot(as.data.frame(t(rates)))
NULL




