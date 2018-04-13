\name{predict.airGP}
\alias{predict.airGP}
\title{Prediction Method for airGP Model Fit}
\description{Predict response values for test cases under a \code{\link{airGP}} model fit.}
\usage{ 
\method{predict}{airGP}(object, x.test, burn = 0.1, nmc = 200, ...)
} 
\arguments{
\item{object}{a fitted model of the class 'airGP'}
\item{x.test}{design matrix of test case}
\item{burn}{fraction of MCMC draws to be discurded as burn-in}
\item{nmc}{number of samples to be used in Monte Carlo etimation from the pos-burn segment of the Markov chain}
\item{...}{Additional arguments, non currently implemented.}
}
\value{
Return the following quantities:
\item{f.hat}{a \code{n + nnew}-by-\code{nsamp} matrix of the posterior mean values of the regression function at the training and test data. First \code{n} rows are training data fitted values, the last \code{nnew} rows are the predicted test cases}
\item{f.samp}{a \code{n + nnew}-by-\code{nsamp} matrix of the posterior sampled values of the regression function at the training and test data. First \code{n} rows are training data fitted values, the last \code{nnew} rows are the predicted test cases}
\item{sigma}{vector of length \code{nsamp} of the posterior draws of \code{sigma}, the noise standard deviation}
\item{subsamp}{The indices of the MCMC draws that were used in Monte Carlo averaging.}
}

\seealso{\code{\link{airGP}} and \code{\link{summary.airGP}}.}
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


## predict at x.test
## NB. returns both fitted values for training cases, followed by predicted values for test cases
pp <- predict(fit, x.test, burn = 0, nsamp = nsweep)

## for prediction on test cases, discard the first n rows which correspond to training cases
f.test.mean <- apply(pp$f.hat[,ss], 1, median)[-(1:n)]
f.test.lowr <- apply(pp$f.hat[,ss], 1, quantile, p = .025)[-(1:n)]
f.test.uppr <- apply(pp$f.hat[,ss], 1, quantile, p = .975)[-(1:n)]

## calculate RMSE
rmse <- function(a, b) return(sqrt(mean((a - b)^2)))
cat(agp.rmse <- round(rmse(y.test, f.test.mean), 4), "\n")

## plot observed vs predicted
plot(f.test.mean, y.test, ylim = range(y.test), xlim = range(y.test), ann = FALSE, ty = "n")
segments(f.test.mean, f.test.lowr, f.test.mean, f.test.uppr, col = "gray")
points(f.test.mean, y.test, pch = 19, cex = 0.5)
abline(0, 1, col = 2)
title(paste("airGP (", round(agp.rmse, 2), ")", sep = ""))
title(ylab = "Test y", xlab = "Predicted")
}
}
\keyword{programming}



