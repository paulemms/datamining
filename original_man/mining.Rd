\author{Tom Minka}
\name{Cars}
\alias{Cars}
\alias{Housing}
\alias{States}
\alias{va.deaths}
\title{Cars, Housing, States, and VA Deaths}
\description{Cleaned-up versions of some popular R datasets.}
\usage{
data(Cars)
data(Housing)
data(States)
data(va.deaths)
}
\details{
  The variables are given more informative names, some irrelevant
  variables are dropped, erroneous values are fixed, and some variables
  are transformed (in the standardized case).
}
\value{
  The specified dataset is loaded as a \code{\link{data.frame}}.
  The first three also define a
  transformed dataset (ending in `T') where the variables are standardized.
}
\seealso{\code{\link{Cars93}},\code{\link{Boston}},\code{\link{state}},\code{\link{VADeaths}}}
\examples{
data(States)
hist(States)
hist(StatesT)
w = pca(StatesT,2)
text.plot(project(StatesT,w),asp=1,cex=.6)
plot.axes(w)
}
\eof
\author{Tom Minka}
\name{auto.aspect}
\alias{auto.aspect}
\title{Choose aspect ratio}
\description{Uses Cleveland's "banking to 45" rule to determine an optimal
  aspect ratio for a plot.}
\usage{
auto.aspect(x,y)
}
\arguments{
  \item{x,y}{numeric vectors of coordinates defining a continuous curve,
    or multiple curves delimited by \code{NA}s.  Alternatively, \code{x}
    can be a list of two vectors.}
}
\value{
  An aspect ratio, suitable for the \code{asp} parameter to
  \code{\link{plot}} or \code{\link{plot.window}}.
}
\examples{
data(lynx)
plot(lynx)
plot(lynx,asp=auto.aspect(time(lynx),lynx))
}
\eof
\author{Tom Minka}
\name{bhist}
\alias{bhist}
\title{Histogram with confidence intervals}
\description{Same as \code{\link{hist}} except a confidence interval is drawn
around each bin height.}
\details{The width of the interval for height p is
sqrt(p*(1-p)/n)*exp(-1/6/p/n).}
\usage{
bhist(x,b=NULL,z=1)
}
\arguments{
\item{x}{a numerical vector}
\item{b}{the number of bins, or a vector of break locations.  If NULL, chosen automatically by \code{\link{hist}}.}
}
\value{None.}
\examples{
x <- c(rnorm(100,-2,0.5),rnorm(100,2,0.5))
b <- seq(-4,4,by=0.25)
bhist(x,b)
}
\eof
\author{Tom Minka}
\name{boxplot.hclust}
\alias{boxplot.hclust}
\title{Boxplot with hierarchical cluster breaks}
\description{A representation of a hierarchical clustering of predefined groups}
\usage{
boxplot.hclust(h,x,k=2:5,...)
}
\arguments{
\item{h}{an \code{hclust} object}
\item{x}{the list of vectors that was clustered to produce \code{h}
(typically via \code{\link{ward}})}
\item{k}{a vector of the cluster cardinalities to plot}
\item{...}{arguments passed to \code{boxplot}}
}
\value{
A boxplot of \code{x} is shown with blue lines cutting the x-axis.
The tallest lines correspond to divisions made at the top of the hierarchy.
By reading top to bottom, you can see how each cluster is subdivided.
This can be far more illuminating than a plot of the hierarchy as a tree.
}
\seealso{
\code{\link{hist.hclust}}, \code{\link{ward}}, \code{\link{break.ward}}
}
\eof
\author{Tom Minka}
\name{break.quantile}
\alias{break.quantile}
\title{Equal-count breaks of a dataset}
\description{Computes a set of breaks using the equal-count algorithm.}
\usage{

}
\eof
\author{Tom Minka}
\name{break.ts}
\alias{break.ts}
\title{Change-point analysis by clustering}
\description{Divide a time-series into homogenous segments.}
\usage{
break.ts(x,n=2,trace=T,same.var=T)
}
\arguments{
\item{x}{a numerical vector or \code{ts} object}
\item{n}{the desired number of segments}
\item{trace}{If TRUE, shows a merging trace via
\code{\link{plot.hclust.trace}}}
\item{same.var}{argument passed to \code{\link{ward}}}
}
\details{
Calls \code{\link{ward}} with \code{sortx=F} to cluster the series into
segments.  Only the marginal distribution of data is used;
temporal smoothness, for example, is ignored.
}
\value{A vector of time breaks.  The breaks are also plotted visually via
\code{\link{plot.segments.ts}}.}
\seealso{
\code{\link{plot.segments.ts}}, \code{\link{plot.breaks}}
}
\examples{
library(ts)
data(LakeHuron)
# single major change
break.ts(LakeHuron,2)
# merging trace suggests n=6 is also interesting:
break.ts(LakeHuron,6)
# interesting oscillation

data(treering)
break.ts(treering[1:500],9,same=T)
break.ts(treering[1:100],7,same=T)
# interesting multiscale structure

x <- c(rnorm(100),rnorm(300)*3,rnorm(200)*2)
b <- break.ts(x,3,same=F)
plot(x,type="l")
plot.breaks(b)
}
\eof
\author{Tom Minka}
\name{break.ward}
\alias{break.ward}
\alias{break.kmeans}
\alias{break.hclust}
\title{Quantize by clustering}
\description{Quantize a one-dimensional variable by calling a clustering
routine.}
\usage{
break.ward(x,n=2,plot=T)
break.hclust(x,n=2,method="ward",plot=T)
break.kmeans(x,n=2,plot=T)
}
\arguments{
\item{x}{a numerical vector}
\item{n}{the desired number of bins}
\item{method}{argument given to \code{\link{hclust}}}
\item{plot}{If TRUE, a histogram with break lines is plotted
(\code{\link{hist.hclust}} or \code{\link{plot.breaks}}).
For \code{break.ward} and \code{break.hclust}, also shows a merging trace
(\code{\link{plot.hclust.trace}}).}
}
\value{A break vector, suitable for use in \code{\link{cut}}
or \code{\link{hist}}.}
\details{These are convenience routines which simply call the appropriate
clustering routine (\code{\link{ward}}, \code{\link{hclust}}, or
\code{\link{kmeans}}), convert
the output to a break vector, and make plots.}
\examples{
x <- c(rnorm(700,-2,1.5),rnorm(300,3,0.5))
break.ward(x,2)
break.hclust(x,2,method="complete")
break.kmeans(x,2)

x <- c(rnorm(700,-2,0.5),rnorm(1000,2.5,1.5),rnorm(500,7,0.1))
break.ward(x,3)
break.hclust(x,3,method="complete")
break.kmeans(x,3)
}
\eof
\author{Tom Minka}
\name{color.cone}
\alias{color.cone}
\title{Depict colors geometrically}
\description{Plots named colors as points in a three-dimensional cone.}
\usage{
color.cone(col,cex=2,light=T,...)
}
\arguments{
  \item{col}{a vector of strings, naming colors.}
  \item{cex}{a number controlling the size of the points.}
  \item{light}{If \code{TRUE}, the cone will be surrounded by point
    light sources.}
}
\details{
  The colors are mapped into Hue-Saturation-Lightness (HSL).
  Hue gives the angle around the cone, Saturation the radial distance
  from the center line, and Lightness the height from the base of the
  cone.
}
\note{
  The high-lightness, high-saturation part of the cone will always be
  empty, because these colors are not achievable on computer displays
  (due to the RGB representation).
}
\value{
  Produces a VRML file which is opened by a VRML viewer.
}
\references{
  Rich Franzen's Wheel of Saturation, Intensity, and Hue.
  \url{http://home.att.net/~rocq/SIHwheel.html}

  Charles Poynton's Color FAQ.
  \url{http://www.poynton.com/notes/colour_and_gamma/ColorFAQ.html}
}
\examples{
color.cone(YlGnBu.colors(8))
color.cone(YR.colors(16))
color.cone(RYB.colors(7))
color.cone(topo.colors(20))
# reveals how topo.colors is not sequential in lightness
}
\eof
\author{Tom Minka}
\name{color.key}
\alias{color.key}
\title{Add a thin key to the top of a plot}.
\description{Similar to \code{\link{legend}} but takes up less space.}
\usage{
color.key(col=NULL,pch=NULL,labels=NULL,breaks=NULL,digits=2,cex=0.75)
}
\arguments{
  \item{col}{a vector of colors, one for each group.  Assumed black if
    none given.}
  \item{pch}{a vector of symbols, one for each group.  If less
    than two distinct symbols are given, no symbols are put in the key.}
  \item{labels}{a character vector giving the label of each group.}
  \item{breaks}{a numeric vector of boundaries between groups.}
  \item{digits}{number of digits to use for \code{breaks}.}
  \item{cex}{character expansion factor.}
}
\value{
  A thin bar is placed at the top of the plot, just inside the plotting
  area.  The bar is divided into equal-sized lengths, colored according
  to \code{col}.  If \code{pch} is given, that symbol is plotted in the
  center of each segment.  If \code{labels} is given, each label is
  placed at the center of the corresponding segment, just above the
  plotting area (using \code{\link{mtext}}).  If \code{breaks} is given,
  the boundary between segments is also labeled with a number.
}
\seealso{\code{\link{color.plot}}}
\examples{
data(iris)
y = as.numeric(iris$Species)
plot(Sepal.Width ~ Sepal.Length, iris,col=y,pch=y)
color.key(1:3,1:3,levels(iris$Species))
}
\eof
\author{Tom Minka}
\name{color.plot.data.frame}
\alias{color.plot.data.frame}
\alias{color.plot.formula}
\title{Plot cases as colored points}
\description{Like \code{\link{plot}} and \code{\link{text.plot}}
  but colors according to the response variable.}
\usage{
color.plot.data.frame(data,labels=F,...)
color.plot.formula(formula, data=parent.frame(), ...)
}
\arguments{
  \item{data}{a data frame.}
  \item{formula}{a formula specifying a response and two predictors from \code{data}}
  \item{labels}{If NULL, cases are plotted as points.  If T, cases are
    plotted as labels, according to \code{rownames}.}
  \item{...}{Extra arguments passed to \code{\link{color.plot.default}}.}
}
\details{
  Calls \code{\link{color.plot.default}} with \code{x} as the first
  predictor in \code{data},
  \code{y} as the second predictor, and \code{z} as the response.
  To get a different predictor/response breakdown than the default, use
  \code{color.plot(formula, x, ...)}, which is shorthand for
  \code{color.plot(model.frame(formula, x), ...)}.

  Each case is plotted with a color determined by the
  response.  If the response is a factor, each factor level is in a
  different color.  If the response is numeric, then it is color-coded
  by assigning each quantile a different color.
}
\seealso{\code{\link{color.plot.default}}}
\examples{
data(iris)
color.plot(iris)
color.plot(Species ~ Petal.Length + Petal.Width, iris)
color.plot(Species ~ Petal.Length, iris)
color.plot(Species ~ Petal.Length, iris,jitter=T)
color.plot(iris, col=1)
color.plot(iris, col=c(1,2))

data(state)
x <- data.frame(state.x77)
color.plot(Murder ~ Frost + Illiteracy, x, labels=T, cex=0.5)
}
\eof
\author{Tom Minka}
\name{color.plot.loess}
\alias{color.plot.loess}
\title{Contour plot of a regression surface}
\usage{
color.plot.loess(object,data,res=50,fill=F,add=fill,clip=T,...)
}
\arguments{
  \item{object}{a \code{\link{loess}} object.}
  \item{data}{data to use instead of \code{model.frame(object)}.}
  \item{res}{resolution of the sampling grid in each direction.}
  \item{fill}{passed to \code{\link{contour.plot}}.}
  \item{add}{If \code{TRUE}, the contours are added to an existing
    plot.  Otherwise a new plot is created.}
  \item{clip}{a polygon over which the surface is to be defined.
    Possible values are \code{FALSE} (no clipping), \code{TRUE}
    (clip to the convex hull of the data), or a
    matrix with two columns specifying (x,y) coordinates.}
  \item{...}{extra arguments to \code{\link{color.plot.data.frame}}}.
}
\details{
  The regression surface is evaluated at all points on a grid,
  clipped values are set to \code{NA}, and \code{contour.plot} is
  used to plot the contours.

  If \code{add=FALSE}, the data is plotted on top using
  \code{color.plot.data.frame}.
}
\value{
  A plot is produced.
}
\seealso{\code{\link{contour.plot}}}
\examples{
data(Housing)
fit = loess(Price ~ Rooms + Low.Status, Housing)
color.plot(fit)
}
\eof
\author{Tom Minka}
\name{color.plot}
\alias{color.plot}
\alias{color.plot.default}
\title{Plot subgroups as colored points}
\description{Like \code{\link{plot}} and \code{\link{text.plot}}
  but colors according to a third variable.}
\usage{
color.plot.default(x,y,z,labels=NULL,data=parent.frame(),
                   xlab,ylab,zlab,
                   xlim=NULL,ylim=NULL,
                   axes=T,key=T,add=F,nlevels=4,
                   color.palette=default.colors,col,
                   pch.palette=c(1,4,2,22,16,17,15),pch,
                   jitter=F,digits=2,mar,bg,...)
}
\arguments{
  \item{x,y}{numeric vectors of the same length.}
  \item{z}{a numeric vector or factor, same length as
    \code{x} and \code{y}.}
  \item{labels}{a character vector of labels, same length as \code{z}.
    If NULL, cases are plotted as points (default).}
  \item{axes}{If FALSE, no axes are plotted (but there may still be a key).}
  \item{key}{If FALSE, no key is plotted.}
  \item{add}{If FALSE, a new plot is made.  If TRUE, points or labels
    are added to the existing plot.}
  \item{nlevels}{an integer.  If \code{z} is numeric, it is
    color-coded using this many levels.  (If \code{z} is a factor,
    color-coding follows the factor levels.)}
  \item{color.palette}{a vector of colors, arbitrary length, or a
    function with integer
    argument which generates a vector of colors (e.g.
    \code{\link{YlGnBu.colors}}).
    Used if \code{col}
    is not specified.  If shorter than the number of
    levels, colors will be recycled and the plotting symbol will change.}
  \item{col}{a vector of colors, as in a call to \code{\link{plot}}.  Used
    to explicitly set the color of each point.}
  \item{pch.palette}{a vector of plotting symbols, arbitrary length.  If
    \code{labels=NULL}, the plot symbol will
    rotate through these when there aren't enough colors.}
  \item{pch}{a vector of plotting symbols, as in a call to \code{\link{plot}}.
    Used to explicitly set the symbol of each point.}
  \item{digits}{the number of digits to use in the color key when \code{z}
    is numeric.}
  \item{bg}{the background color}
  \item{mar}{figure margins (see \code{par}).  If not specified,
    appropriate margins will be chosen automatically, which will not
    necessarily match the current value of \code{par("mar")}.}
}
\details{
  Each (x,y) point is plotted with a color determined by \code{z}.
  If \code{z} is a factor, each factor level is in a
  different color.  If \code{z} is numeric, then it is color-coded
  by assigning each quantile a different color.

  Despite the name, this function can also make plots using different
  symbols instead of colors.  For example, if \code{color.palette=1} then
  all points will be black but use different symbols.

  When \code{labels != NULL}, the result is equivalent to
  \code{\link{text.plot}} with colors.

  If \code{key=TRUE} and there is more than one color or symbol on the
  plot, a key is displayed at the top of the figure.
}
\note{
  This function sets the figure margins permanently, so that you can
  draw on the color plot.  Unfortunately, this also means future plots
  will use the same margins, until you change them with \code{par("mar")}.
}
\value{
  A plot is produced.
}
\seealso{\code{\link{color.plot.data.frame}},\code{\link{color.plot.loess}},\clode{\link{color.plot.glm}},\code{\link{color.plot.knn}},\code{\link{color.plot.tree}},\code{\link{YlGnBu.colors}}}
\examples{
# See the examples for color.plot.data.frame
}
\eof
\author{Tom Minka}
\name{colors}
\alias{default.colors}
\alias{default.colors.w}
\alias{YR.colors}
\alias{YlGnBu.colors}
\alias{OrRd.colors}
\alias{gray.colors}
\alias{RYB.colors}
\alias{BrBg.colors}
\alias{RC.colors}
\alias{GM.colors}
\title{Color schemes}
\description{Compute a color scheme with a specified number of levels.}
\usage{
default.colors(n)
default.colors.w(n)

YR.colors(n)
YlGnBu.colors(n)
OrRd.colors(n)
gray.colors(n)

RYB.colors(n)
BrBg.colors(n)
RC.colors(n)
GM.colors(n)
}
\arguments{
  \item{n}{the number of colors desired}
}
\details{
  The first two schemes are \emph{categorical}, providing maximum
  separation in hue, intended for depicting unordered categories.
  \code{default.colors} has only dark colors, good for coloring points,
  while \code{default.colors.w} includes light colors, good for filling
  regions.

  The next four schemes are \emph{sequential}, from light to dark,
  with variation in hue to increase discrimination.
  They are intended for depicting ordered levels.
  The sequential order is more easily perceived with these schemes than
  with the built-in palettes \code{heat.colors},
  \code{terrain.colors},
  and \code{topo.colors}.
  The ordering can also be seen by the color-blind and when printed in black
  and white.

  The main difference between the sequential schemes is the variation in
  hue, with \code{YR.colors} having the most variation and
  \code{gray.colors} having the least.
  Generally you should choose the amount of variation according to the
  number of levels.
  I recommend \code{OrRd.colors} for three levels, \code{YlGnBu.colors}
  for four to eight levels, and \code{YR.colors} beyond eight levels.

  The last four schemes are \emph{double-ended} or \emph{diverging}
  schemes,
  which progress from one hue to a second hue, passing through white
  in the middle.  They are intended for representing signed ordered
  levels, such as residuals.  The main difference between them is the
  amount of separation between colors, so generally you use
  \code{GM.colors} when you want a few levels and
  \code{RYB.colors} when you want many.

  These functions can be used as the \code{color.palette} parameter to
  \code{\link{filled.contour}} and \code{\link{color.plot}}, for example.
}
\value{
  A vector of strings, naming colors.
}
\seealso{\code{\link{colors}},\code{\link{rainbow}},\code{\link{color.cone}}}
\references{
  The schemes in
  \code{YlGnBu.colors}, \code{OrRd.colors}, \code{RYB.colors}, and
  \code{BrBg.colors} are from ColorBrewer.
  The scheme in \code{YR.colors} is from Howard Seltman.

  Mark A. Harrower and Cynthia A. Brewer.  ColorBrewer: An
  Online Tool for Selecting Color Schemes for Maps, \emph{The Cartographic
  Journal}, in press.
  \url{http://www.colorbrewer.org/},
  \url{http://www.personal.psu.edu/faculty/c/a/cab38/ColorBrewerBeta.html}

  Generalized color schemes for Mapping and Visualization.
  From Cynthia Brewer, Color Use Guidelines for Mapping and
  Visualization.
  Reprinted at the Gallery of Data Visualization by Michael Friendly.
  \url{http://www.math.yorku.ca/SCS/Gallery/images/S12-fullstructureClean.gif}

  Dan Carr.
  Color perception, the importance of gray and residuals, on a choropleth map.
  \emph{Statistical Computing & Graphics Newsletter} 5(1):16-20, 1994
  \url{http://cm.bell-labs.com/cm/ms/who/cocteau/newsletter/issues/back/v51.pdf}
}
\examples{
data(Housing)
color.plot(Price ~ Rooms + Low.Status, Housing, bg=gray(0.5),
           color.palette=YlGnBu.colors)
color.plot(Price ~ Rooms + Low.Status, Housing, bg=gray(0.5),
           color.palette=YR.colors)
color.plot(Price ~ Rooms + Low.Status, Housing, bg=gray(0.5),
           color.palette=RYB.colors,nlevels=5)

# also see examples for color.cone
}
\eof
\author{Tom Minka}
\name{contour.plot}
\alias{contour.plot}
\title{Display contours}
\description{Create a filled or unfilled contour plot.}
\usage{
contour.plot(x,y,z,fill=F,
        levels=NULL,nlevels=if(is.null(levels)) 4 else length(levels),
                         level.data=z,
                         zlim,equal=F,pretty=T,lwd=2,drawlabels=F,
                         color.palette=YlGnBu.colors,bg=grey(0.5),
                         key=T,main="",...)
}
\arguments{
  \item{x,y}{locations at which the values in \code{z} are measured.}
  \item{z}{a matrix containing the values to be plotted (NAs are allowed).}
  \item{fill}{If \code{TRUE}, makes a filled contour plot.}
  \item{levels}{numeric vector of levels at which to draw contour
    lines.  The next few arguments only apply to the case when
    \code{levels==NULL}.}
  \item{nlevels}{The number of contour lines to draw.}
  \item{level.data}{numeric vector to use for computing \code{levels}.}
  \item{zlim}{The range that the levels should cover.  Defaults to the
    range of \code{level.data}.}
  \item{equal}{If \code{TRUE}, the levels are equally spaced over
    \code{zlim}.  Otherwise they match the quantiles of
    \code{level.data}.}
  \item{pretty}{If \code{TRUE}, the levels are chosen to be round numbers.}
  \item{lwd}{width of the contour lines.}
  \item{drawlabels}{If \code{TRUE}, the contour lines are labeled with
    their level.}
  \item{color.palette}{a vector of colors, or a function which takes a
    number and returns a vector of that length.}
  \item{bg}{the background color.}
  \item{key}{If \code{TRUE}, a color key is drawn at the top of the
    plot.  If \code{key==2}, the key of \code{filled.contour}
    is used.}
  \item{main}{the title of the plot.}
  \item{...}{extra arguments passed to \code{contour} or
    \code{filled.contour}.}
}
\details{
  This is a wrapper function for \code{\link{contour}} and
  \code{\link{filled.contour}} whose main purpose is to provide a
  uniform interface and provide a decent automatic choice of levels.

  If \code{equal=FALSE}, the levels are chosen according to the
  equal-count algorithm of \code{\link{break.quantile}}.
}
\seealso{\code{\link{contour}},\code{\link{filled.contour}},
  \code{\link{color.plot.loess}}}
\examples{
# see the examples for color.plot.loess
}
\eof
\author{Tom Minka}
\name{data.image}
\alias{data.image}
\title{Data Image}
\description{Each value in a data matrix is represented by a colored pixel.}
\usage{
data.image(x,reorder=T,scale=T,col=YR.colors(64),...)
}
\arguments{
  \item{x}{a matrix or data frame.}
  \item{reorder}{If \code{TRUE}, the rows and columns are reordered so
    that the data along each row and each column follows a linear trend.
    Can also be a vector of two logicals, indicating separately
    whether to reorder rows and/or columns.
  }
  \item{scale}{If \code{TRUE}, each column is scaled to have minimum
    0 and maximum 1.}
  \item{col}{a vector of colors to use for representing values.}
  \item{...}{extra arguments to \code{\link{image.default}}.}
}
\seealso{\code{\link{YR.colors}}}
\references{
M. C. Minnotte and R. W. West.
"The Data Image: A Tool For Exploring High Dimensional Data Sets."
Proceedings of the ASA Section on Statistical Graphics, 1998.
\url{http://citeseer.nj.nec.com/72085.html}
\url{http://math.usu.edu/~minnotte/research/pubs.html}
}
\examples{
data(USJudgeRatings)
data.image(USJudgeRatings,col=RC.colors(32))
}
\eof
\author{Tom Minka}
\name{dim.ordered}
\alias{dim.ordered}
\title{Is a dimension ordered?}
\description{Query whether each array dimension is ordered or not.}
\usage{
dim.ordered(x)
}
\arguments{
\item{x}{an array}
}
\details{
  Intuitively, a dimension is ordered if its values are in a natural
  progression, such as (1,2,3) or ("one","two","three").
  Programmatically, a dimension is considered ordered if and only if the
  corresponding \code{\link{dimnames}} vector has an \code{ordered}
  attribute set to \code{T}.
}
\value{
  A named logical vector, with an entry for each dimension, indicating
  whether it is ordered or not.
}
\eof
\author{Tom Minka}
\name{hist.hclust}
\alias{hist.hclust}
\title{Histogram with hierarchical cluster breaks}
\description{A representation of a one-dimensional hierarchical
clustering}
\usage{
hist.hclust(h,x,k=2:5)
}
\arguments{
\item{h}{an \code{hclust} object}
\item{x}{the numerical vector that was clustered to produce \code{h}}
\item{k}{a vector of the cluster cardinalities to plot}
}
\value{
A histogram of \code{x} is shown with blue lines cutting the x-axis.
The tallest lines correspond to divisions made at the top of the hierarchy.
By reading top to bottom, you can see how each cluster is subdivided.
This can be far more illuminating than a plot of the hierarchy as a tree.
}
\seealso{
\code{\link{boxplot.hclust}}, \code{\link{ward}}, \code{\link{break.ward}}
}
\eof
\author{Tom Minka}
\name{labeled.curves}
\alias{labeled.curves}
\title{Plot labeled curves}
\description{Draws multiple curves, each with a properly-placed label
  and unique color/linestyle.}
\usage{
labeled.curves(x,y,se=NULL,labels,xlab,ylab,type="o",group,
               color.palette=default.colors(6),col,
               lty.palette=1:6,lty,lwd=2,jitter=0,
               legend=F,
               cex=par("cex"),horizontal=F,
               mar,bg=par("bg"),ylim,main="",...)
}
\arguments{
  \item{x}{a numeric vector giving the horizontal position of each point.}
  \item{y}{a matrix or data frame giving the vertical position of each point.
    Each column defines one curve.}
  \item{se}{a matrix or data frame, the same size as \code{y}, giving an
    error bound on each value.}
  \item{labels}{a character vector of labels for the lines.  Defaults to
    the column names of \code{y}.}
  \item{xlab,ylab}{axis labels, by default the dimnames of \code{y}.}
  \item{type}{indicates whether to draw points, lines, or both.
    See \code{\link{plot}}.}
  \item{group}{a numeric vector specifying how to group cases when
    assigning colors.}
  \item{color.palette}{a vector of colors, arbitrary length, or a
    function with integer
    argument which generates a vector of colors.  Used if \code{col}
    is not specified.  If shorter than the number of
    curves, colors will be recycled and the line style will change.}
  \item{col}{a vector of colors, as in a call to \code{\link{plot}}.  Used
    to explicitly set the color of each curve.}
  \item{lty.palette}{a vector of line styles, arbitrary length.
    The line style will rotate through these when there aren't enough colors.}
  \item{lty}{a vector of line styles, as in a call to \code{\link{plot}}.
    Used to explicitly set the style of each curve.}
  \item{lwd}{line width.}
  \item{jitter}{the amount by which to jitter the error bars, to avoid
    overplotting.  \code{jitter=0.02} is usually sufficient.}
  \item{legend}{If \code{TRUE}, the labels are placed in a legend.
    Otherwise the labels are placed next to the lines.}
  \item{cex}{character expansion factor.}
  \item{horizontal}{If \code{TRUE}, the axes are laid out horizontally,
    so that the curves run vertically.}
  \item{mar}{a numerical vector giving the lines of margin on the four
    sides of the plot (see \code{\link{par}}).
    If missing, it is set according to \code{\link{auto.mar}}.}
  \item{bg}{background color.}
  \item{ylim}{desired plotting limits.}
  \item{main}{title for the plot.}
  \item{...}{extra arguments for low-level plotting commands.}
}
\details{
  Point \code{j} in curve \code{i} is at \code{(x[j],y[j,i])}.
  Thus all curves must be the same length, and have points at the same
  horizontal positions.
  If \code{y[i,j]=NA}, then the curve has a break at the previous point, and
  resumes at the next \code{non-NA} point.

  If \code{se} is given, then an error bar will be placed around each point.
}
\seealso{
  \code{\link{parallel.plot}}}, \code{\link{linechart}}
}

\eof
\author{Tom Minka}
\name{linechart}
\alias{linechart}
\title{Linechart}
\description{Plot each data row as a curve.}
\usage{
linechart(y,se=NULL,effects=F,med=F,
          xscale=c("equal","linear","none"),...)
}
\arguments{
  \item{y}{a named vector or matrix.}
  \item{se}{a vector or matrix, same size as \code{y}, of error bounds.
    Alternatively, \code{y} can be \code{list(y,se)}.}
  \item{effects}{If \code{TRUE}, the columns are shifted to have zero
    mean.}
  \item{med}{If \code{TRUE} and \code{effects=TRUE}, the columns are
    shifted to have zero median.}
  \item{xscale}{describes how the columns should be placed.}
  \item{...}{additional arguments to \code{\link{labeled.curves}}.}
}
\details{
  If \code{xscale="linear"}, the columns are placed to make each curve
  as straight as possible.
  If \code{xscale="equal"}, the columns are placed similar to
  \code{"linear"} but with the constraint that they must be equally
  spaced.
  If \code{xscale="none"}, the columns are placed in the order that
  they appear in the matrix.  This is automatic if \code{y} has ordered
  columns (see \code{\link{dim.ordered}}).
  If \code{se != NULL}, error bars are drawn around each point.

  Linecharts are a replacement for dotcharts and mosaics.
}
\seealso{
  \code{\link{dotchart}},
  \code{\link{mosaicplot}}
}
\examples{
# compare to a dotchart
data(VADeaths)
dotchart(VADeaths, main = "Death Rates in Virginia - 1940")
dimOrdered(VADeaths)[2] = F
linechart(VADeaths)
linechart(t(VADeaths))

# compare to a mosaicplot
data(HairEyeColor)
x <- margin.table(HairEyeColor,c(1,2))
dimOrdered(x) = F
mosaicplot(x)
x = t(x)
col = c("brown","blue","red","green")
linechart(row.probs(x),color.pal=col)
linechart(row.probs(x,se=T),color.pal=col)
linechart(row.probs(x,se=T),jitter=0.02,color.pal=col)
mosaicplot(x)
linechart(row.probs(t(x),se=T))

data(blood)
dimOrdered(blood) = F
linechart(row.probs(blood,se=T))

data(antacids)
dimOrdered(antacids) = F
linechart(row.probs(antacids,se=T))
mosaicplot(t(antacids))
}
\eof
\name{make.names}
\alias{make.names}
\title{Make Syntactically Valid Names}
\description{
  Make syntactically valid names out of character vectors.
}
\usage{
make.names(names, unique = FALSE)
}
\arguments{
 \item{names}{character (vector) to be coerced to syntactically valid
   names.}
 \item{unique}{logical; if \code{TRUE}, the resulting elements are
     unique.  This may be desired for, e.g., column names.}
}
\value{
  A character vector of same length as \code{names} with each changed to
  a syntactically valid name.
}
\details{
  A syntactically valid name consists of letters, numbers, and the dot
  character and starts with a letter or the dot.

  All invalid characters are translated to \code{"."}.  A missing value
  is translated to \code{"NA"}.  Names which match R keywords have a dot
  appended to them.

  If \code{unique = TRUE} then \code{\link{make.unique}} is used to
  append sequence numbers to duplicates (after coercion).
}
\seealso{
  \code{\link{names}},
  \code{\link{character}},
  \code{\link{data.frame}}.
}
\examples{
make.names(c("a and b", "a_and_b"), unique=TRUE)
# "a.and.b"  "a.and.b1"

data(state)
state.name[make.names(state.name) != state.name]# those 10 with a space
}
\keyword{character}
\eof
\author{Thomas P Minka}
\name{make.unique}
\alias{make.unique}
\title{Make character strings unique}
\description{Makes the elements of a character vector unique by
  appending sequence numbers to duplicates.}
\usage{
make.unique(names,sep=".")
}
\arguments{
  \item{names}{a character vector}
  \item{sep}{a character string used to separate a duplicate name from
    its sequence number.}
}
\value{
  a character vector of same length as \code{names} with duplicates changed.
}
\details{
  The algorithm used by \code{make.unique} has the
  property that
  \code{make.unique(c(A,B)) = make.unique(c(make.unique(A),B))}.

  In other words, you can append one string at a time to a vector,
  making it unique each time, and get the same result as applying
  \code{make.unique} to all of the strings at once.

  If character vector \code{A} is already unique, then
  \code{make.unique(c(A,B))} preserves \code{A}.
}
\seealso{
  \code{\link{make.names}}
}
\examples{
make.unique(c("a","a","a"))
make.unique(c(make.unique(c("a","a")),"a"))

make.unique(c("a","a","a.2","a"))
make.unique(c(make.unique(c("a","a")),"a.2","a"))

rbind(data.frame(x=1),data.frame(x=2),data.frame(x=3))
rbind(rbind(data.frame(x=1),data.frame(x=2)),data.frame(x=3))
}
\keyword{character}
\eof
\author{Tom Minka}
\name{merge.factor}
\alias{merge.factor}
\title{Merge factor levels}
\description{Merges factor levels with similar response distributions
(assumed normal).}
\usage{
merge.factor(f,x,n,same.var=T,trace=T,xlab=NA,ylab=NA)
}
\arguments{
\item{f}{a factor}
\item{x}{a numerical vector, same length as \code{f}}
\item{n}{the desired number of factor levels}
\item{same.var}{argument passed to \code{\link{ward}}}
\item{trace}{If TRUE, a merging trace is plotted
(\code{\link{plot.hclust.trace}})}
\item{xlab,ylab}{axis labels.  If NA, taken from f and x arguments.}
}
\details{Calls \code{ward(split(x,f))} to get a tree, cuts the tree, and
constructs a new factor.  The tree is shown via \code{\link{boxplot.hclust}}.}
\value{A new factor, same length as \code{f}, but with \code{n} levels.}
\examples{
n <- 20
x <- c(rnorm(n)+1, rnorm(n)+2, rnorm(n)*4+2)
f <- gl(3,n)
levels(f) <- c("a","b","c")
merge.factor(f,x,2,same.var=T)
merge.factor(f,x,2,same.var=F)

# an ordered factor
data(va.deaths)
merge.factor(va.deaths$Age,va.deaths$Rate,2)
}
\eof
\author{Tom Minka}
\name{merge.hist}
\alias{merge.hist}
\title{Merge histogram bins}
\description{Quantize a variable by merging similar histogram bins.}
\details{
The desired number of bins is achieved by successively merging the
two most similar histogram bins.  The distance between bins of height
(f1,f2) and width (w1,w2) is measured according to the
chi-square statistic
\deqn{w1*(f1-f)^2/f + w2*(f2-f)^2/f}
where f is the height of
the merged bin:
\deqn{f = (f1*w1 + f2*w2)/(w1 + w2)}
}
\usage{merge.hist(x,b=NULL,n=b)}
\arguments{
\item{x}{a numerical vector}
\item{b}{the starting number of bins, or a vector of starting break locations.
If NULL, chosen automatically by \code{\link{hist}}.}
\item{n}{the desired number of bins.}
}
\value{A vector of bin breaks, suitable for use in \code{\link{hist}},
\code{\link{bhist}}, or \code{\link{cut}}.
Two plots are shown: a \code{\link{bhist}} using the returned bin breaks,
and a merging trace.  The trace shows, for each merge, the chi-square
distance of the bins which were merged.  This is useful for determining
the appropriate number of bins.  An interesting number of bins is one
that directly precedes a sudden jump in the chi-square distance.
}
\examples{
x <- c(rnorm(100,-2,0.5),rnorm(100,2,0.5))
b <- seq(-4,4,by=0.25)
merge.hist(x,b,10)
# according to the merging trace, n=5 and n=11 are most interesting.

x <- runif(1000)
b <- seq(0,1,by=0.05)
merge.hist(x,b,10)
# according to the merging trace, n=6 and n=9 are most interesting.
# because the data is uniform, there should only be one bin,
# but chance deviations in density prevent this.
# a multiple comparisons correction in merge.hist may fix this.
}
\eof
\author{Tom Minka}
\name{merge.table.cost}
\alias{merge.table.cost}
\alias{merge.table.cost.inorder}
\alias{merge.table.cells}
\usage{
merge.table.cost(x,v,cost)
merge.table.cost.inorder(x,v,cost)
merge.table.cells(x,v,i,j)
}
\title{Internal routines}
\description{Internal routines for \code{\link{merge.table}}}
\eof
\author{Tom Minka}
\name{merge.table}
\alias{merge.table}
\title{Table merging}
\description{Merges similar rows and columns of a contingency table.}
\usage{
merge.table(x, bins=rep(2,length(ds)), ds=1:length(dim(x)))
}
\arguments{
\item{x}{a \code{\link{table}}}
\item{bins}{the desired number of levels for each dimension being merged.
a numeric vector, the same length as \code{ds}.}
\item{ds}{a vector of dimensions to merge, either by name or number.
default is all of them.}
}
\value{A merged \code{\link{table}}.  The total count is the same as \code{x}.
A merging trace is plotted which shows, for each merge, the chi-square
distance of the slices which were merged.
This is useful for determining
the appropriate dimensions.  An interesting number is one
that directly precedes a sudden jump in the chi-square distance.
}
\details{
The desired table dimensions are achieved by successively merging the
two most similar slices.
(`Slice' generalizes `row' and `column' to higher-dimensional tables.)
The distance between slices
is measured according to the chi-square statistic.
Merging two slices means adding together their counts, and
concatenating their labels with a comma in between.
If a dimension is ordered (according to \code{\link{dim.ordered}}), only
adjacent slices are considered for merging, and their labels are concatenated
with a dash in between.
}
\seealso{
\code{\link{sort.table}},
\code{\link{mosaicplot}},
\code{\link{linechart}}
}
\examples{
i <- factor(c(1,1,2,2,3,3,4,4))
j <- factor(c(3,4,3,4,1,2,1,2))
x <- table(i,j)
merge.table(x,c(2,2))

i <- factor(c(1,1,3,3,2,2,4,4))
j <- factor(c(2,4,2,4,1,3,1,3))
x <- table(i,j)
merge.table(x,c(2,2))

# one ordered dimension
data(education)
merge.table(education,c(3,2))

data(occupation)
merge.table(occupation,c(3,4))
}
\eof
\author{Tom Minka}
\name{mine.associations}
\alias{mine.associations}
\title{Associations in a contingency table}
\description{Find unusually frequent variable combinations.}
\usage{
mine.associations(x,top=10,targets=NULL)
}
\arguments{
\item{x}{a contingency table}
\item{top}{the number of top associations to return}
\item{targets}{a character vector of table variables, one of which must be
included in every returned association.}
}
\details{
Enumerates all two-way marginal tables, sorts the cells by lift, and returns
the top.  The formula for lift is \deqn{\frac{p(i,j)}{p(i)p(j)}}.
}
\value{
A data frame where the first column is \code{Lift} and the remaining columns
correspond to the variables of \code{x}.
Each row describes a cell of \code{x}, where marginalized
variables have value \code{NA}.
The lift value describes the ratio of the actual count in the cell versus
the expected count under independence.
}
\examples{
data(Titanic)
mine.associations(Titanic)
# Females are twice as likely to survive as the average person.
# Members of 3rd class are twice as likely to be children as the average person.
# Etc.

# focus on associations with survival
mine.associations(Titanic,target="Survived")

# focus on children
mine.associations(Titanic[,,1,],target="Survived")
}
\eof
\author{Tom Minka}
\name{parallel.plot}
\alias{parallel.plot}
\title{Parallel-coordinate plot}
\description{Plots a data matrix using parallel coordinates.}
\usage{
parallel.plot(x,yscale=c("linear","range","none"),
              xscale=c("equal","linear","none"),proto=T,...)
}
\arguments{
  \item{x}{a data matrix or frame.  Only numeric columns will be
    plotted.}
  \item{yscale}{describes how the columns should be scaled.}
  \item{xscale}{describes how the columns should be placed.}
  \item{proto}{If \code{TRUE} and the response variable is a factor,
    then only the mean value for each response is plotted.
    See \code{\link{prototypes}}.}
  \item{...}{Extra arguments to \code{\link{labeled.curves}}.}
}
\details{
  Let the rows of the matrix be \emph{cases} and the columns be
  \emph{dimensions}.
  Each dimension is laid out as a vertical axis, and the value for each
  case is plotted as a point on that axis.  Points corresponding to the
  same case are connected by lines.  Thus each case is represented by a
  curve or \emph{profile}, moving horizontally from dimension to dimension.

  To enhance readability, the dimensions should usually be shifted to have
  a common center and scaled to have
  comparable units.  This is controlled by \code{yscale}.
  If \code{yscale="linear"}, the dimensions are automatically shifted
  and scaled
  by the linear profiles method.  That is, the profiles are
  made to be as straight as possible.
  This may result in negative scaling for a dimension,
  i.e. the axis is reversed so that moving up means a lower value.
  When this happens, the dimension name is prepended with a minus sign.
  If \code{yscale="range"}, the dimensions are shifted and scaled (by a
  positive number) so that the data ranges from 0 to 1 on each axis.
  If \code{yscale="none"}, the dimensions are plotted as-is.

  Another important choice is the ordering and spacing of the axes
  on the plot.
  If \code{xscale="linear"}, the axes are placed according to the
  linear profiles method.
  If \code{xscale="equal"}, the axes are placed similarly to
  \code{"linear"} but with the constraint that they must be equally
  spaced.
  If \code{xscale="none"}, the axes are placed in the order
  that they appear in the matrix.

  If the data frame has a categorical response, the profiles are colored
  according to the response.  The color scheme can be adjusted using the
  arguments to \code{\link{labeled.curves}}.
}
\references{
A. Inselberg and B. Dimsdale.
"Parallel coordinates: A tool for visualizing multidimensional geometry."
Proc. of Visualization '90, p. 361-78, 1990.

E. J. Wegman.
"Hyperdimensional Data Analysis Using Parallel Coordinates."
JASA 85:664-675, 1990.
\url{http://links.jstor.org/sici?sici=0162-1459%28199009%2985%3A411%3C664%3AHDAUPC%3E2.0.CO%3B2-7}
}
\seealso{\code{\link{parallel.cases}},\code{\link{star.plot}}}
\examples{
data(iris)
parallel.plot(iris)
parallel.plot(iris,proto=F,labels=NULL)
}
\eof

\eof
\author{Tom Minka}
\name{plot.axes}
\alias{plot.axes}
\title{Plot axes under projection}
\description{Shows how unit vectors along the original dimensions appear
  under a projection.}
\usage{
plot.axes(w,col=2,origin,keep=0.2)
}
\arguments{
  \item{w}{a numeric array with two columns and named rows}
  \item{col}{color of arrows}
  \item{origin}{If T, arrows emerge from (0,0).  If F, arrows emerge
    from the center of the figure.  If not given, arrows emerge from
    (0,0) if (0,0) is visible, otherwise from the figure center.}
  \item{keep}{a length in inches, below which an arrow is not plotted.}
}
\details{
  Each row of the array specifies a location in the figure.  This
  location can be understood as the place where a unit vector
  along one of the original dimensions would appear.  An arrow is
  drawn from the origin pointing toward this location and labeled with
  the row name.  The relative length of the arrow is determined by the
  distance of the location from the origin, but all arrows are scaled as
  much as possible, to reduce crowding, while remaining inside the
  figure limits.  Arrows which are very short are not plotted at all.
}
\value{Uses \code{\link{arrows}} to draw arrows on top of the current plot.}
\seealso{
  \code{\link{project}}, \code{\link{pca}}, \code{\link{projection}}
}
\eof
\author{Tom Minka}
\name{plot.breaks}
\alias{plot.breaks}
\title{Draw break locations}
\description{Draws lines to indicate break locations.}
\usage{
plot.breaks(b)
}
\arguments{
\item{b}{a numeric vector}
}
\value{Uses \code{\link{segments}} to draw thick blue vertical lines on top
of the current plot, at the given locations.}
\seealso{
\code{\link{plot.hclust.breaks}},
\code{\link{plot.segments.ts}},
\code{\link{break.kmeans}}
}
\eof
\author{Tom Minka}
\name{plot.hclust.trace}
\alias{plot.hclust.trace}
\title{Plot a merging trace}
\description{Plots the cost of successive merges in a hierarchical clustering}
\usage{
plot.hclust.trace(h,k=1:10)
}
\arguments{
\item{h}{an \code{hclust} object}
\item{k}{a vector of the cluster cardinalities to plot}
}
\value{
The trace shows, for each merge, the
distance of the clusters which were merged.  This is useful for determining
the appropriate number of clusters.  An interesting number of clusters is one
that directly precedes a sudden jump in distance.
}
\seealso{
\code{\link{ward}}, \code{\link{break.ward}}
}
\eof
\author{Tom Minka}
\name{plot.segments.ts}
\alias{plot.segments.ts}
\title{Plot time-series segments}
\description{Shows how the mean level of a time-series varies between
segments.}
\usage{
plot.segments.ts(x,b)
}
\arguments{
\item{x}{a numeric vector or \code{ts} object}
\item{b}{a numeric vector of break times}
}
\value{Plots the time series as a line, then draws a blue horizontal line
through each segment, at the mean value of the segment.}
\seealso{
\code{\link{break.ts}},
\code{\link{plot.breaks}}
}
\eof
\author{Tom Minka}
\name{predict.plot}
\alias{predict.plot}
\alias{predict.plot.data.frame}
\alias{predict.plot.formula}
\title{Plot predictors versus response.}
\description{Makes a matrix of pairwise scatterplots with
  lowess-type trend lines.}
\usage{
predict.plot(...)
predict.plot.formula(formula,data=parent.frame(),...)
predict.plot.data.frame(x,layout=NULL,partial=NULL,
                           rtype=c("prob","logit","probit"),
                           highlight,se=T,
                           identify.pred=F,
                           mcol=0,mlwd=2,
                           scol="red",slwd=2,
                           span=2/3,degree=2,family="symmetric",
                           mar=NA,bg=par("bg"),
                           xaxt=par("xaxt"),yaxt=par("yaxt"),
                           col=1,asp=NULL,
                           given=NULL,given.lab=NULL,nlevels=2,
                           pretty=T,key=!is.null(given),
                           color.palette=default.colors,
                           pch.palette=default.pch,
                           main=NULL,xlab,ylab,...)
}
\arguments{
  \item{formula}{a formula specifying the response and predictor variables}
  \item{data,x}{a data frame with at least two columns}
  \item{partial}{a model from which to compute partial residuals (used
    by \code{\link{predict.plot.lm}}).}
  \item{mcol,mlwd}{If plotting partial residuals of an \code{lm},
    the color and width of the model predictions.}
  \item{layout}{a vector \code{c(rows,cols)} specifying the desired
    layout of panels.  Otherwise chosen automatically based on the size
    of the plotting window.}
  \item{highlight}{a logical vector specifying which predictors to highlight.}
  \item{se}{If \code{TRUE}, show standard errors in linecharts.}
  \item{scol,slwd}{color and width of trend lines.}
  \item{span,degree,family}{parameters for the trend line (see \code{loess}).}
  \item{rtype}{how a factor response should be handled when drawing a
    trend line.}
  \item{identify.pred}{A character vector of predictor names for which
    to interactively \code{\link{identify}} points.  If \code{TRUE},
    done for all predictors.}
  \item{mar}{margins within each panel}
  \item{xaxt,yaxt}{arguments to \code{par}}
  \item{col}{plotting color for symbols}
  \item{asp}{Aspect ratio for each panel.  If \code{"auto"}, the aspect
    ratio is chosen automatically based on the trend line and
    \code{\link{auto.aspect}}.}
  \item{given,given.lab,nlevels,pretty,key,bg,color.palette,pch.palette}{used for conditioning plots.}
  \item{main,xlab,ylab}{axis labels.}
  \item{...}{extra arguments passed to \code{predict.plot.data.frame} or
    \code{\link{plot}}.}
}
\value{
  If the predictor is numeric, makes a scatterplot with loess line on top.
  If the predictor is a factor, makes a \code{\link{linechart}}.
}
\seealso{\code{\link{loess}},\code{\link{model.plot}}}
\examples{
data(Cars)
predict.plot(Price~.,CarsT)
fit = lm(Price~.,CarsT)
predict.plot(Price~.,CarsT,partial=fit)
# same thing using predict.plot.lm
predict.plot(fit,partial=T)
}
\eof
\author{Tom Minka}
\name{predict.plot.lm}
\alias{predict.plot.lm}
\title{Plot predictors versus residuals.}
\description{Makes a matrix of pairwise scatterplots with
  lowess-type trend lines.}
\usage{
predict.plot.lm(object,data,partial=F,ylab,...)
}
\arguments{
  \item{object}{the output of \code{lm}.}
  \item{data}{a data frame to use instead of \code{model.frame(object)}.}
  \item{partial}{If \code{TRUE}, plot partial residuals instead of residuals.}
  \item{ylab}{axis label.}
  \item{...}{extra arguments for \code{\link{predict.plot.data.frame}}.}
}
\details{
  Partial residuals are computed by fitting a new model with the
  predictor removed, which is different from \code{residuals(type="partial")}.
}
\value{
  A plot similar to \code{\link{predict.plot}}, but where the vertical
  axis is residuals.  These plots can be used to judge which predictors
  should be added to the model.
}
\seealso{\code{\link{predict.plot}}}
\examples{
# see the examples for predict.plot
}
\eof
\author{Tom Minka}
\name{project}
\alias{project}
\title{Project data into fewer dimensions}
\description{Reduces the dimensionality of a data set.}
\usage{
project(x,w)
}
\arguments{
  \item{x}{a data frame}
  \item{w}{a matrix with named rows and columns}
}
\details{
  Each column of \code{w} specifies a new variable, which is to be
  constructed by combining existing variables according to the given weights.
}
\value{
  A data frame where the variables named in the rows of \code{w} are
  replaced by new variables named in the columns of \code{w}.
  Other variables are left unchanged.
}
\seealso{
  \code{\link{pca}}, \code{\link{projection}}, \code{\link{plot.axes}}
}
\examples{
data(iris)
w = projection(iris,k=2)
# w only involves the continuous attributes
# the new variables are h1 and h2
x = project(iris,w)
color.plot(x)
plot.axes(w)
}
\eof
\author{Tom Minka}
\name{projection}
\alias{projection}
\title{Discriminative projection}
\description{Finds a projection matrix that separates data groups.}
\usage{
projection(x,y,k=1,given=NULL,type=c("mv","m","v"),...)
}
\arguments{
  \item{x}{a data frame of variables to project.}
  \item{y}{a factor or numeric vector which defines the groups to
    separate.
    If \code{y} is not given, it is
    taken to be the response variable of \code{x} (the last column if
    \code{x} is not a \code{model.frame}).}
  \item{k}{the number of dimensions to project \code{x} down to.}
  \item{given}{A matrix specifying axes to avoid.  The projection matrix
    will be orthogonal to \code{given}.}
  \item{type}{see below.}
  \item{...}{additional parameters depending on \code{type}.}
}
\details{
  This function only uses the second-order statistics of the data (means
  and covariances of the groups).

  If \code{type="m"}, the within-group covariances are assumed equal
  and the projection will try to separate the projected means.

  If \code{type="v"}, the within-group means are assumed equal
  and the projection will try to separate the projected covariances.

  If \code{type="mv"}, the projection will try to separate the projected
  means and covariances, by maximizing the divergence between the
  projected classes (as Gaussians).

  If \code{y} is a numeric vector, overlapping classes are defined by
  grouping data points with similar values of \code{y}.
  The optional argument \code{span} controls how big the classes will
  be (as a percentage of the dataset), and \code{res} controls the
  amount of overlap.  The total number of classes will be
  \code{res}/\code{span}.  The default values are usually acceptable.

  The projection is "stabilized" so that small changes in the data do
  not cause sign flips in the projection.
}
\value{
  A matrix suitable for input to \code{\link{project}},
  with named rows matching columns of \code{x} and
  columns named \code{h1}, ..., \code{hk}.
  Each column denotes a new dimension to be obtained as a linear combination
  of the variables in \code{x}.
}
\references{
  m-projection is Fisher's linear discriminant analysis.
  mv-projection is heteroscedastic discriminant analysis:

  N. Kumar and A.G. Andreou.
  Heteroscedastic discriminant analysis and
  reduced rank HMMs for improved speech recognition.
  \emph{Speech Communication} 26: 283-297, 1998.

  The case when \code{y} is numeric is sliced inverse
  regression:

  K.-C. Li.  Sliced inverse regression for dimension reduction.
  \emph{Journal of the American Statistical Association} 86(414):
  316-327, 1991.
}
\seealso{\code{\link{project}},\code{\link{pca}}}
\examples{
# illustrate difference between (m,v,mv)
library(MASS)
m1 <- c(6,6)
v1 <- array(c(2,1.9,1.9,2),c(2,2))
#v1 <- array(c(1,0,0,1),c(2,2))
x1 <- mvrnorm(100,m1,v1)
m2 <- c(0,0)
v2 <- array(c(20,0,0,10),c(2,2))
x2 <- mvrnorm(300,m2,v2)
x = as.data.frame(rbind(x1,x2))
y = factor(c(rep(1,nrow(x1)),rep(2,nrow(x2))))
plot(x[,1],x[,2],col=1,xlab="",ylab="",asp=1)
points(x2[,1],x2[,2],col=2)
w = projection(x,y,type="m")
abline(0,w[2]/w[1],col=3)
w = projection(x,y,type="v")
abline(0,w[2]/w[1],col=4)
w = projection(x,y,type="mv")
abline(0,w[2]/w[1],col=5)
my.legend(1,c("m","v","mv"),col=3:5,lty=1)

# regression projection
x1 <- 2*runif(200)-1
x2 <- 2*runif(200)-1
y <- x1^2/2 + x2^2
x <- data.frame(x1,x2)
color.plot(x[,1],x[,2],y)
w = projection(x,y)
abline(0,w[2]/w[1],col=4)
}
\eof
\author{Tom Minka}
\name{read.array}
\alias{read.array}
\alias{write.array}
\title{Read and write arrays}
\description{Read data into an array from the console or file}
\details{The array can have one of several formats.  The preferred format,
produced by \code{\link{write.array}}, looks like \code{
                               cvar.nam
rvar.1.nam   ... rvar.k.nam    cvar.lev.1 ... cvar.lev.l
rvar.1.lev.1 ... rvar.k.lev.1  ...        ... ...
}
}
\usage{
read.array(file, sep = "", quote = "\"", skip = 0)
}
\argument{
\item{file}{either a character string naming a file or a connection which the data are to be read from or written to. "" indicates input from the console for reading and output to the console for writing.}
\item{sep}{the field separator string. Values on each line of the file are separated by this string.}
\item{quote}{the set of quoting characters as a single character string.}
\item{skip}{the number of lines of the input file to skip before beginning to read data.}
}
\value{An array.}
\examples{
#read.array("class/data/ex9-21.dat")
#x <- read.array("class/data/mm2-14.dat")
#x <- read.array("class/data/arizona.txt")
}
\eof
\author{Tom Minka}
\name{sort.cells}
\title{Sort an array}
\description{Sort the cells of an array}
\usage{
sort.cells(x)
}
\arguments{
\item{x}{an array}
}
\details{Converts \code{x} to a data frame and calls
\code{\link{sort.data.frame}}.}
\value{A data frame representation of \code{x}, sorted by cell value.}
\seealso{
\code{\link{sort.data.frame}}
}
\examples{
data(Titanic)
sort.cells(Titanic)

data(HairEyeColor)
sort.cells(HairEyeColor)
}
\eof
\author{Tom Minka}
\name{sort.data.frame}
\title{Sort rows of a data frame}
\description{Re-orders data frame rows so that a column is sorted.}
\usage{
sort.data.frame(df,f=ncol(df))
}
\arguments{
\item{df}{a data frame}
\item{f}{the name or number of a column of \code{df}, or a numeric vector
to play the role of a column of \code{df}.}
}
\value{
A reordered version of \code{df} where the values in column \code{f} are in
ascending order.
}
\seealso{
\code{\link{sort.cells}}
}
\examples{
data(mtcars)
sort.data.frame(mtcars,"mpg")
}
\eof
\author{Tom Minka}
\name{sort.levels}
\alias{sort.levels}
\title{Sort levels of a factor}
\description{The levels of a factor are sorted according to a summary
statistic.}
\usage{
sort.levels(f,x,fun=median)
}
\arguments{
\item{f}{a factor}
\item{x}{a vector, same length as \code{f}}
\item{fun}{function to use to summarize groups of \code{x}}
}
\value{
A copy of \code{f} with re-ordered levels.
}
\examples{
data(OrchardSprays)
sort.levels(OrchardSprays$treatment, OrchardSprays$decrease)
sort.levels(OrchardSprays$treatment, OrchardSprays$decrease, fun=mean)
}
\eof
\name{sort.table}
\title{Sort rows and columns of a contingency table}
\description{Score rows and columns via correspondence analysis, then sort
the scores.}
\usage{
sort.table(x)
}
\arguments{
\item{x}{a contingency table}
}
\seealso{
\code{\link{mosaicplot}}
}
\examples{
data(blood)
sort.table(blood)
}
\eof
\author{Tom Minka}
\name{star.plot}
\alias{star.plot}
\title{Reordered and scaled star plot}
\description{Automatically reorders and scales the variables to make a
  readable star plot.}
\usage{
star.plot(x,proto=T,reorder=T,scale=T,...)
}
\arguments{
  \item{x}{a matrix or data frame.  Each row makes one star.}
  \item{proto}{If \code{TRUE} and the response variable is a factor,
    then only the mean value for each response is plotted.
    See \code{\link{prototypes}}.}
  \item{reorder}{If \code{TRUE}, the rows and columns are reordered so
    that the data along each row and each column follows a linear trend.
    Can also be a vector of two logicals, indicating separately
    whether to reorder rows and/or columns.
  }
  \item{scale}{If \code{TRUE}, each column is scaled to have minimum
    0 and maximum 1.}
  \item{...}{additional arguments for \code{\link{stars}}.}
}
\seealso{\code{\link{stars}}}
\examples{
data(iris)
star.plot(iris)
}
\eof
\name{state.area}
\title{U.S. state area 1990}
\usage{
data(state.area)
}
\details{
\code{state.area} is a named vector of the 1990 U.S. state land areas,
including the District of Columbia, in square miles.
The numbers are more accurate and up-to-date than the ones in
\code{data(state)}.
}
\seealso{
\code{\link{state.pop}}, \code{\link{state}}
}
\references{
The data was obtained from
\url{http://www.census.gov/population/censusdata/90den_stco.txt}
and edited into R format.
County areas can also be obtained from the above link.
}
\eof
\name{state.pop}
\title{U.S. state population 1790-1990}
\usage{
data(state.pop)
}
\details{
\code{state.pop} is a matrix containing the population of
each U.S. state, including the District of Columbia,
from 1790 to 1990 in ten-year increments.
The rows are the years "1790" to "1990" and the columns are state names.
}
\seealso{
\code{\link{state.area}}, \code{\link{state.name}}, \code{\link{state}}
}
\references{
The data was obtained from
\url{http://merrill.olm.net/pdocs/feas/pop/pop1790_1990/pii.txt}
and edited into R format.
}
\examples{
dotchart(sort(log(state.pop["1990",])))

# distribution of growth rates over time
rates <- diff(log(state.pop))
boxplot(as.data.frame(t(rates)))
}
\eof
\author{Tom Minka}
\name{text.plot}
\alias{text.plot}
\alias{text.plot.default}
\alias{text.plot.data.frame}
\alias{text.plot.formula}
\title{Make a plot of text labels}
\description{Like \code{\link{text}} except it creates a new plot with
  limits chosen to make all labels visible.}
\usage{
text.plot.default(x,y,labels,xlab,ylab,xlim=NULL,ylim=NULL,
                  cex=par("cex"),adj=NULL,srt=0,axes=T,...)
text.plot.data.frame(data,labels=rownames(data),...)
text.plot.formula(formula,data=parent.frame(),...)
}
\arguments{
  \item{x,y}{numeric vectors, same length}
  \item{data}{a data frame with at least two columns}
  \item{formula}{a formula specifying a response and predictor variable}
  \item{labels}{character vector, same length as \code{x}}
  \item{srt}{rotation angle as in \code{\link{text}}.}
  \item{adj}{text justification as in \code{\link{text}}.}
}
\details{
  The main job of \code{\link{text.plot}} is finding the right plot
  limits so that all labels are visible.  It does this by computing the
  bounding box of each label (taking into account text rotation and
  justification) and solving a system of inequalities which ensure that
  each label fits entirely into the plot window.
}
\examples{
data(state)
x <- data.frame(state.x77)
text.plot(x$Frost, x$HS.Grad, rownames(x))
# same thing, using text.plot.formula
text.plot(HS.Grad ~ Frost, x)
# notice how the limits change
text.plot(HS.Grad ~ Frost, x, srt=45)
text.plot(HS.Grad ~ Frost, x, srt=90)
}
\eof
\author{Tom Minka}
\name{vrml.plot3d}
\alias{vrml.plot3d}
\alias{vrml.plot3d.default}
\alias{vrml.plot3d.data.frame}
\alias{vrml.plot3d.formula}
\title{View points in 3D}
\description{Creates a VRML scene with the given points.}
\usage{
vrml.plot3d.default(x,y,z,data=parent.frame(),xlab,ylab,zlab,
                    pch=0,col=3,scale=c(1,1,1),file.name=NULL,
                    cex=1,cex.axis=12,light=T,axes=T,...)
vrml.plot3d.data.frame(frame,...)
vrml.plot3d.formula(formula,data,...)
}
\arguments{
  \item{x,y,z}{numeric vectors giving the coordinates of the points.}
  \item{frame}{a data frame with three columns.}
  \item{data}{an environment in which to lookup x,y,z.}
  \item{xlab,ylab,zlab}{axis labels.}
  \item{pch}{a numeric code for the graphic object to place at
    each location.  See details.}
  \item{col}{a single color or vector of colors, used to color the points.}
  \item{scale}{a vector of 3 numbers, defining the size of the box into
    which the points are placed.  You can control the aspect ratio this way.}
  \item{file.name}{a filename (with or without the .wrl extension) to
    receive the VRML.  If \code{NULL}, a temporary name is chosen.}
  \item{cex}{a scale factor for the objects.}
  \item{cex.axis}{a scale factor for the axis labels.}
  \item{light}{If \code{TRUE}, the scene will contain a light at each
    corner of the cube, with no headlight.
    Otherwise, the scene will contain no lights and rely on the
    headlight only.}
  \item{axes}{If \code{TRUE}, a box with axis labels is drawn around the
    points.}
  \item{...}{additional arguments to pass to the internal drawing routines,
    or to pass through to \code{vrml.plot3d.default}
    (for the \code{data.frame} and \code{formula} methods).}
}
\details{
  A VRML scene description file is created, and opened with a browser if
  one is available.

  VRML is a standard language for describing 3D scenes, and the file
  produced by this function should be portable across all machines
  with a VRML browser.

  If \code{pch=0}, a cube is placed at each location.
  If \code{pch=1}, a sphere is placed at each location.
  If \code{pch="."}, a dot is placed at each location.
  The latter option is the fastest, and will be chosen automatically if
  the scale is small or the dataset is large.
}
\references{
  \url{http://web3d.vapourtech.com/tutorials/vrml97/}

  \url{www.wed3d.org/resources/vrml_ref_manual/Book.html}
}
\examples{
data(Housing)
w = pca(HousingT,k=3)
x = project(HousingT,w)
plot.new()
vrml.plot3d(x)
}
\eof
\author{Tom Minka}
\name{vrml.surface}
\alias{vrml.surface}
\alias{vrml.surface.default}
\alias{vrml.surface.loess}
\title{View a surface in 3D}
\description{Creates a VRML scene with a shaded surface.}
\usage{
vrml.surface.default(x,y,z,xlab=NULL,ylab=NULL,zlab=NULL,
                     col="gray",scale=c(1,1,1),file.name=NULL,
                     cex.axis=12,light=F,...)
}
\arguments{
  \item{x,y}{locations of grid lines at which the values in \code{z} are
    measured.}
  \item{z}{a matrix containing the values to be plotted (\code{NA}s are
    allowed).}
  \item{xlab,ylab,zlab}{axis labels.  If \code{NULL}, taken from the
    deparsed expressions for \code{x,y,z}.}
  \item{col}{the color of the surface.}
  \item{scale}{a vector of 3 numbers, defining the size of the box into
    which the surface is placed.  You can control the aspect ratio this way.}
  \item{file.name}{a filename (with or without the .wrl extension) to
    receive the VRML.  If \code{NULL}, a temporary name is chosen.}
  \item{cex.axis}{a scale factor for the axis labels.}
  \item{light}{If \code{TRUE}, the scene will contain its own light
    and no headlight.
    Otherwise, the scene will contain no lights and rely on the
    headlight only.}
  \item{border}{If \code{TRUE}, a grid will be drawn on top of the surface.}
  \item{creaseAngle}{a parameter controlling the smoothness of the
    surface.  When it is small, faces at large angles to each other will
    create visible "creases".}
  \item{...}{additional arguments for internal drawing routines.}
}
\details{
  A VRML scene description file is created, and opened with a browser if
  one is available.

  VRML is a standard language for describing 3D scenes, and the file
  produced by this function should be portable across all machines
  with a VRML browser.

  This function is similar to \code{\link{persp}} except it creates a 3D
  scene which can be manipulated with a VRML viewer.
}
\seealso{\code{\link{persp}},\code{\link{vrml.plot3d}}}
\examples{
#  The Obligatory Mathematical surface.
#  Rotated sinc function.
x <- seq(-10, 10, length= 30)
y <- x
f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
z <- outer(x, y, f)
z[is.na(z)] <- 1
vrml.surface(x,y,z)
vrml.surface(x,y,z,border=T)
}
\eof
\author{Tom Minka}
\name{ward}
\alias{ward}
\title{Create a hierarchy by Ward's method}
\description{Produces a hierarchical clustering of one-dimensional data via Ward's method.}
\usage{
ward(x, n=rep(1,length(x)), s=rep(1,length(x)), sortx=T, same.var=T)
}
\arguments{
\item{x}{a numerical vector, or a list of vectors.}
\item{n}{if x is a vector of cluster means, n is the size of each cluster.}
\item{s}{if x is a vector of cluster means, s is the sum of squares in each
cluster.  only needed if \code{same.var=F}.}
\item{sortx}{if \code{sortx=F}, only clusters which are adjacent in \code{x}
can be merged.  Used by \code{\link{break.ts}}.}
}
\item{same.var}{if \code{same.var=T}, clusters are assumed to have the same
true variance, otherwise not.  This affects the cost function for merging.}
}
\details{Repeatedly merges clusters in order to minimize the clustering cost.
By default, it is the same as \code{hclust(method="ward")}.
If \code{same.var=T}, the cost is the sum of squares:
\deqn{sum_c sum_{i in c} (x_i - m_c)^2}
The incremental cost of merging clusters ca and cb is
\deqn{(n_a*n_b)/(n_a+n_b)*(m_a - m_b)^2}
It prefers to merge clusters which are small and have similar means.

If \code{same.var=F}, the cost is the sum of log-variances:
\deqn{sum_c n_c*log(1/n_c*sum_{i in c} (x_i - m_c)^2)}
It prefers to merge clusters which are small, have similar means,
and have similar variances.

If \code{x} is a list of vectors, each vector is assumed to be a
cluster.  \code{n} and \code{s} are computed for each cluster and
\code{x} is replaced by the cluster means.
Thus you can say \code{ward(split(x,f))} to cluster the data for different
factors.
}
\value{The same type of object returned by \code{\link{hclust}}.}
\section{Bugs}{Because of the adjacency constraint used in implementation,
the clustering that results
from \code{sortx=T} and \code{same.var=F} may occasionally be suboptimal.}
\seealso{
\code{\link{hclust}},
\code{\link{plot.hclust.trace}},
\code{\link{hist.hclust}},
\code{\link{boxplot.hclust}},
\code{\link{break.ward}},
\code{\link{break.ts}},
\code{\link{merge.factor}}
}
\examples{
x <- c(rnorm(700,-2,1.5),rnorm(300,3,0.5))
hc <- ward(x)
opar <- par(mfrow=c(2,1))
plot.hclust.trace(hc)
hist.hclust(hc,x)
par(opar)

x <- c(rnorm(700,-2,0.5),rnorm(1000,2.5,1.5),rnorm(500,7,0.1))
hc <- ward(x)
opar <- par(mfrow=c(2,1))
plot.hclust.trace(hc)
hist.hclust(hc,x)
par(opar)

data(OrchardSprays)
x <- OrchardSprays$decrease
f <- factor(OrchardSprays$treatment)
# shuffle levels
#lev <- levels(OrchardSprays$treatment)
#f <- factor(OrchardSprays$treatment,levels=sample(lev))
hc <- ward(split(x,f))
# is equivalent to:
#n <- tapply(x,f,length)
#m <- tapply(x,f,mean)
#s <- tapply(x,f,var)*n
#hc <- ward(m,n,s)
boxplot.hclust(hc,split(x,f))
}
\eof
