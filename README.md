# airGP
Additive-Interactive regression with Gaussian processes

Contains R package for additive-interactive regression with Gaussian processes. This is primarily aimed at nonparametric smoothing based regression in the high-dimensional setting, where one may have many more predictors than observations. The additive-interactive regression framework (Yang and Tokdar, 2015) assumes that the response function f(x), with x being high dimensional, can be decomposed as an additive sum of component functions f_k(x), where each component function depends on a handful of the predictors and their interactions. The implementation done in this package is loosely based on Qamar and Tokdar (2015), with important technical differences that will be documented in a forthcoming paper.

## References
Yang, Y. and Tokdar, S. T. (2015). Minimax-Optimal Nonparametric Regression in High Dimensions. The Annals of Statistics, 43(2), 652-674.

Qamar, S. and Tokdar, S. T. (2016). Bayesian Additive Gaussian Process Regression. arXiv:1411.7009.


