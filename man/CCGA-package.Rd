\name{CCGA-package}
\alias{CCGA-package}
\docType{package}
\title{
Case-control genetic association analysis with adjustment of non-confounding covariates. 
}
\description{
In CCGA, the power of disease-SNP association analysis in case-control studies can be potentially improved through efficiently adjustmenting non-confounding covariates, compared with the standard logistic regression method with or without covariate adjustment.

Let Y, S, G, and X denote the case-control status, stratum variable, SNP genotype coded by 0, 1, or 2, and covariate(s), respectively. The model for relating the response Y and (S,G,X) is

g(pr(Y=1)) = alpha + beta_S S + beta_G G + beta_X X,

where g(.) is a given function, which can be the logit function or the probit function. The disease prevalence is required for each stratum.

The disease prevalence(s), Hardy-Weinberg equilibrium, and independence between  non-confounding covariates and SNP genotype are efficiently incorporated in the retrospective likelihood function, and the nonparametric distribution of the covariate(s) is profiled out through the application of Lagrange's multipler method. The multipliers can be directly estimated using the available data, which yields the so called profile maximum likelihood estimates of the regression parameter (pMLE). Alternatively, the data-dependent multiplier(s) can be replaced with the limiting value(s) to yield a modified profile likeihood function, which results in modified profile maximum likelihood estimates (mpMLE). Theoretically, mpMLE and pMLE are asymptotically equivalent in terms of efficiency, but mpMLE could be computationally much simplier and faster.

There are two candidate link functions ("logit" and "probit"), and the p

This package includes the main function "CCGA" and a simulated dataset "data" for illustration.
}
\details{
\tabular{ll}{
Package: \tab CCGA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-10-27\cr
License: \tab Artistic License 2.0\cr
}
}
\author{
Hong Zhang

Maintainer: Hong Zhang <zhanghfd@fudan.edu.cn>
}
\references{
Zhang H, Chatterjee N, Rader D, Chen J. (2016) Adjustment of Non-confounding Covariates in Case-control Genetic Association Studies. Annals of Applied Statistics (revised).
}


\keyword{package}
\keyword{case-control study}
\keyword{Hardy-Weinberg equilibrium}
\keyword{Lagrange multiplier}
\keyword{maximum likelihood estimate}
\keyword{non-confounding covariate}
\keyword{prevalence}
