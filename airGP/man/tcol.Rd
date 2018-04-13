\name{tcol}
\alias{tcol}
\title{Transparent color}
\description{Creates transparent variations of named \code{R} colors}
\usage{
tcol(col, alpha = 1)
}
\arguments{
\item{col}{A vector of R colors. Poissble choices for specifying the colors include a character string. e.g. "red", an integer value, e.g., 2 for red, or a "#RRGGBB" code form such as "#FF0000" for red. See \code{\link{par}} for more details and other choices.}
\item{alpha}{A transparency value between 0 (transparent) and 1 (opaque).}
}
\value{
A color code of the form "#RRGGBBTT", with the last two places encoding transparency.
}

\seealso{"Color Specification" under \code{\link{par}}.}
\examples{
tcol("red", 0.5)
tcol(2, 0.5)
tcol("#FF0000", 0.5)
## all give "#FF000080"

tcol(c("red", "blue", "darkgreen"), 0.5)
## returns
##         red        blue   darkgreen
## "#FF000080" "#0000FF80" "#00640080"
##
}
\keyword{programming}
