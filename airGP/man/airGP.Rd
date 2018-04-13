\name{airGP}
\alias{airGP}
\title{
Additive Interactive Regression with Gaussian Processes
}
\description{
Estimates a nonparametric regression function in high dimensional setting assuming a sparse additive-interactive structure. The regression function is taken to decompose into additive component functions, each of which depends on a small number of interactive predictors.
}
\usage{
airGP(x.train, y.train, nsamp = 100, thin = 10,
      sparsity = list(ncomp = NULL, dmax = 4, Esize = 1, pactive = c(1,1),
          shrink = c("shape", "rate"), decay.type = c("poly", "geo")),
      lowrank = list(max.rank = nrow(x.train), tol = 1e-6, lpenalty = NULL),
      hyper = list(ERsq = 0.95, lam = c(9.9, 0.1), rho = c(1, 5.0), sig = c(0,0),
          prox = c(0.99, 0.93, 0.75, 0.45, 0.15), nrho = 10),
      adapt = list(var.p = NULL, init.meth = c("corr", "lasso", "gp1"),
          taper.imp = 0.95, simMat = NULL, simFn = cor),
      state = NULL, verbose = TRUE, repeat.call = FALSE)
}
\arguments{
\item{x.train}{ numeric design matrix of predictors for training data. A column of ones should NOT be included.}

\item{y.train}{ numeric vector of response values for training data.}

\item{nsamp}{ number of Markov chain samples to store. Defaults to 100.}

\item{thin}{ thinning frequency of the Markov chain, defaults to 10. The Markov chain runs for \code{nsamp*thin} many iterations and samples are stored every \code{thin}-th iteration.}

\item{sparsity}{ sparsity parameters supplied as a list containing the following items: \code{ncomp}: maximum number of additive components (see Details for default); \code{dmax}: maximum interaction order, i.e., maximum number of predictors to be included in any of the additive components, defaults to 4; \code{Esize}: a priori average size of each (active) component, defaults to 1; \code{pactive}: A vetor of length 2 giving the 50-th and 99.9-th percentiles for a beta prior on the probability that a component is active; \code{shrink}: whether and how to shrink the signal strengths of the components by employing different gamma priors on the signal strength parameters. When set to "shape" or "rate"", the shape or rate parameters of the gamma distributions vary across components. No shrinkage happens when set to anything else. Defaults to shrinkage by shape; \code{decay.type}: if shrinkage is applied, then one can decide whether to decay the shape or inverse-rate parameters geoemtrically or polynomially across components, so that first component has, a priori, stochastically higher signal than the second component and so on. If set to anything other than "poly" or "geo" then all components are shrunk uniformly. Defaults to polynomial decay.}

\item{lowrank}{ lowrank approximation parameters supplied as a list containing the following items: \code{max.rank}: maximum low rank to use, defaults to number of observations in the training data; \code{tol}: error tolerance for truncation of incomplete Cholesky factorization; \code{penalty}: penalty scores used in rank ordering observations for selection as knots.}

\item{hyper}{ parameters for Gaussian process (GP) priors used on each additive component, supplied as a list consisting of the items: \code{ERsq}: a priori expected \emph{R-squared} value used to define the prior mean on the signal-to-noise-ratio (SNR) parameter of the GP; \code{lam}: hyper-parameters for the prior on GP length-scale; \code{rho}: hyper-parameters for the prior on GP SNR; \code{sig}: hyper-parameters for noise variance; \code{prox}: vector of real positive numbers between 0 and 1. The length-scale or correlation range parameters of the component GPs are reparametrized in terms of "proximity" -- measured by the correlation they induce between the function values at a distance 0.1 apart. Defaults to c(0.99, 0.93, 0.75, 0.45, 0.15); \code{nrho}: Number of draws to use in parameter augmentation step for \emph{rho}. Defaults to 10.}

\item{adapt}{ parameters controlling adaptation of the Markov chain, supplied as a list containing the following items: \code{var.p}: initial importance scores of the predictors; \code{init.meth}: method to use to assign initial importance scores; \code{taper.imp}: a tapering parameter that controls how many 'similar' predictors should get importance score boost when boosting the score of a given predictor; \code{simMat}: a similarity matrix could be provided to judge similarity between predictors; \code{simFn}: alternatively one can provide the function to measure similarity between predictors, defaults to absolute value of correlation coefficient.}

\item{state}{ a list giving the initial state of the Markov chain. If \code{NULL} (default) then this state is initialized from the prior. If supplied, one must provide a list containing all of the following items: \code{p.empty}: probability that a component is empty; \code{pincl.A}: odds of the expected probability of not including a predictor; \code{p.incl}: a vector of length \code{sparsity$ncomp} giving the predictor inclusion probabilities of each component; \code{comp.size}: an integer vector of length \code{sparsity$ncomp} giving the size of each component, size means the number of predictors included in the component; \code{comp.ix}: a character vector of length \code{sparsity$ncomp} giving the predictors included in each component -- each element of the vector is character string of predictor indices separate by "." (i.e., "2.5.19" means predictors 2, 5, and, 19 are included); \code{gp.par}: a vector of length \code{2*sparsity$ncomp} giving the length-scale and SNR parameters of the components; \code{fcomp}: a matrix of dim \code{n}-by-\code{ncomp} (or a long vector of same size) giving the GP function evaluation of each component at all training data points.}

\item{verbose}{ logical indicating whether to print progress report during the Markov chain run. Defaults to \code{TRUE}. Recommended to be set to \code{FALSE} when running in batch or in parallel.}

\item{repeat.call}{ logical indicating if the exact function call should be printed. This may help when trying out different possible settings. Defaults to \code{FALSE}.}

}

\value{
Returns a list of class "airGP" containing the following components:
\item{ix.store}{A \code{ncomp}-by-\code{nsamp} matrix of character strings. Each column lists the groups of active variables in all components.}
\item{pars.store}{A \code{2*ncomp}-by-\code{nsamp} matrix. Each consecuitive pair of columns gives saved draws of \code{(rho^2, lam^2)} of a component.}
\item{sigma.store}{A vector of length \code{nsamp} giving the draws of residual standard deviation for the scaled response data. To get draws of \code{sigma} one should multiply \code{sy}.}
\item{ftot.store}{A vector of length \code{nsamp} giving the draws of regression function at the observed predictors for the scaled response. To get draws of \code{f} one should multiply by \code{sy} and then add to it \code{my}.}
\item{pactive.store}{A vector of length \code{nsamp} giving the draws of \code{pactive} -- the probability that a component is active.}
\item{rhob.store}{A vector of length \code{nsapm} giving the draws of \code{rhob} -- the global signal-to-noise ratio parameter.}
\item{state}{A list giving the state of the Markov chain at the completion of the MCMC. See \code{state} under input variables.}
\item{nprop}{A vector of length 4 giving the counts of the different proposal types: "add", "remove", "swap", "do nothing". The last one is vacuous.}
\item{nacpt}{A vector of length 4 giving the acceptance counts for each proposal type}
\item{var.p}{a vector of length \code{p} giving the final importance score attached to each predictors.}
\item{my}{scalar giving the mean response}
\item{sy}{scalar giving the response standard deviation}
\item{sx}{an \code{n}-by-\code{p} matrix giving the scaled version of the design matrix used in model fitting.}
\item{x}{the original design matrix}
\item{obs}{the response vector scaled to have mean zero and variance one}
\item{hyper}{the list of hyper parameters}
\item{lamsq.grid}{the grid on \code{lam^2} used in model fitting, as deduced from the specified grid on \code{prox}}
\item{runtime}{total runtime in seconds}
\item{lowrank}{A list giving the lowrank approximatino details. See \code{lowrank} under input variables.}
\item{details}{A list giving details of the problem setup that are useful for making addiitonal runs.}
}
\references{
Qamar, S. and Tokdar, S.T. (2015). Additive Gaussian Process Regression. arXiv:1411.7009 [stat.ME]
}
\author{
Surya T Tokdar
}

\seealso{
\code{\link{summary.airGP}} and \code{\link{predict.airGP}} for visual summarization and prediction to new cases.
}
\examples{
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
fit <- airGP(x.train, y.train, nsamp = 20, thin = 1, lowrank = list(max.rank = 20))
## WARNINGT: the above run is too short for MCMC mixing
\dontrun{
fit <- airGP(x.train, y.train, nsamp = 500, thin = 10)
}
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
