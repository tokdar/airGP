\name{airGP1v}
\alias{airGP1v}
\title{
Single Variable Scoring in airGP
}
\description{
Likelihood scoring of each predictor
}
\usage{
airGP1v(x.train, y.train,
       lowrank = list(max.rank = nrow(x.train), tol = 1e-6, lpenalty = NULL),
       hyper = list(ERsq = 0.95, lam = c(9.9, 0.1), rho = c(1,1), sig = c(0,0),
           delta = 0.1, prox = c(0.99, 0.93, 0.75, 0.45, 0.15)))
}
\arguments{
\item{x.train}{numeric design matrix of predictors for training data. A column of ones should NOT be included.}

\item{y.train}{numeric vector of response values for training data.}

\item{lowrank}{ lowrank approximation parameters supplied as a list containing the following items: \code{max.rank}: maximum low rank to use, defaults to number of observations in the training data; \code{tol}: error tolerance for truncation of incomplete Cholesky factorization; \code{penalty}: penalty scores used in rank ordering observations for selection as knots.}

\item{hyper}{ parameters for Gaussian process (GP) priors used on each additive component, supplied as a list consisting of the items: \code{ERsq}: a priori expected \emph{R-squared} value used to define the prior mean on the signal-to-noise-ratio (SNR) parameter of the GP; \code{lam}: hyper-parameters for the prior on GP length-scale; \code{rho}: hyper-parameters for the prior on GP SNR; \code{sig}: hyper-parameters for noise variance; \code{delta}: mesh-size for the discretization of SNR to determine a grid of \emph{rho-squared} values; \code{prox}: vector of real positive numbers between 0 and 1. The length-scale or correlation range parameters of the component GPs are reparametrized in terms of "proximity" -- measured by the correlation they induce between the function values at a distance 0.1 apart. Defaults to c(0.99, 0.93, 0.75, 0.45, 0.15).}
}
\value{
Returns a list containing the following components:
\item{ls0}{A single number giving the log-likelihood score of for the null model (no predictors).}
\item{ls1}{A \code{p}-vector containing the log-likelihood score of each predictor.}
\item{runtime}{total runtime in seconds}
}
\references{
Qamar, S. and Tokdar, S.T. (2015). Additive Gaussian Process Regression. arXiv:1411.7009 [stat.ME]
}
\author{
Surya T Tokdar
}

\seealso{
\code{\link{airGP}} for full airGP regression.
}

\examples{
\dontrun{
require(ppls)
data(cookie)

train.set <- 1:40
freqs <- seq(1100, 2498, 2)

x.train <- as.matrix(cookie[train.set, 1:700])
y.train <- cookie[train.set, 701]

scores.gp <- airGP1v(x.train, y.train)
plot(exp(scores.gp$ls1 - scores.gp$ls0), ty = "h", ylab = "Single predictor scores", xlab = "Predictor index")
}
}
