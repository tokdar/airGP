\name{summary.airGP}
\alias{summary.airGP}
\title{Summary Method for airGP Model Fit}
\description{Summarize model fit, including variable importance and interaction, for \code{\link{airGP}}}
\usage{
\method{summary}{airGP}(object, burn = 0.1, nmc = 200, ngraph = 10,
        use.names = FALSE, pars.plot = TRUE, graph.plot = TRUE,
        inclusion.plot = TRUE, keep.layout = FALSE, reorder = TRUE,
        ...)
}
\arguments{
\item{object}{a fitted model of the class 'airGP'.}
\item{burn}{fraction of MCMC draws to be discurded as burn-in}
\item{nmc}{number of samples to be used for Monte Carlo estimation from the pos-burn segment of the Markov chain}
\item{ngraph}{number of primary nodes to be shown in the variable interaction graphs}
\item{use.names}{logical indicating if variable names to be used in the graph. Default is \code{FALSE} for which variables are represented by their position number.}
\item{pars.plot}{logical indicating if a component-by-component posterior summary of the R-square and proximity parameters should be presented. Defaults to \code{TRUE}}
\item{graph.plot}{logical indicating if the interaction graph should be plotted. Defaults to \code{TRUE}}
\item{inclusion.plot}{logical indicating if the posterior variable inclusions probabilities should be plotted. Defaults to \code{TRUE}}
\item{keep.layout}{logical indicating if the existing plot layout should be preserved. Defaults to \code{FALSE}, for which a 2x2 plot grid is used}
\item{reorder}{Logical indicating whether components should be reordered from the most important to the least important for better visualization.}
\item{...}{Additional arguments. None currently implemented.}
}
\value{
Displays a visual summary of the posterior, including posterior spreads of the covariance parameters \code{Rsq} and \code{prox}, variable inclusion probability, and, a variable interaction graph.

The following quantities are returned invisibly as list items:
\item{graph}{the connectivity matrix for of the interaction graph}
\item{inclusion}{vector of inclusion probabilities for predictors}
}

\seealso{\code{\link{airGP}} and \code{\link{predict.airGP}}.}
\examples{
\dontrun{
require(ppls)
data(cookie)

train.set <- 1:40
freqs <- seq(1100, 2498, 2)

x.train <- as.matrix(cookie[train.set, 1:700])
y.train <- cookie[train.set, 701]
x.test <- as.matrix(cookie[-train.set, 1:700])
y.test <- cookie[-train.set, 701]

dimnames(x.train)[[2]] <- paste(freqs, "nm", sep = "")
dimnames(x.test)[[2]] <- paste(freqs, "nm", sep = "")

n <- nrow(x.train)
p <- ncol(x.train)

nsweep <- 20 ## for a real run try nsweep <- 1e3
ss <- 11:20 ## to be used for posterior summary

## fit model and summarize
fit <- airGP(x.train, y.train, nsamp = 500, thin = 10, lowrank = list(max.rank = 40))
summary(fit, ngraph = 10, use.names = TRUE)
}
}
\keyword{programming}
