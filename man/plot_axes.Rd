% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project.R
\name{plot_axes}
\alias{plot_axes}
\title{Plot axes under projection}
\usage{
plot_axes(
  w,
  col = 2,
  origin = NULL,
  keep = NULL,
  top = NULL,
  cex = par("cex"),
  labels,
  ...
)
}
\arguments{
\item{w}{a numeric array with two columns and named rows}

\item{col}{color of arrows}

\item{origin}{If T, arrows emerge from (0,0).  If F, arrows emerge
from the center of the figure.  If not given, arrows emerge from
(0,0) if (0,0) is visible, otherwise from the figure center.}

\item{keep}{a length in inches, below which an arrow is not plotted.}
}
\value{
Uses \code{\link{arrows}} to draw arrows on top of the current plot.
}
\description{
Shows how unit vectors along the original dimensions appear
under a projection.
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
\seealso{
\code{\link{project}}, \code{\link{pca}}, \code{\link{projection}}
}
\author{
{Tom Minka}
}
