\name{SingleSNP}
\alias{SingleSNP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of regression parameters and their standard errors for a single SNP.
}
\description{
This function returns parameter estimates and their standard errors for the parameters in the model relating disease status and (1) SNP genotype, (2) stratum variable, and (3) covariate(s). The estimated parameters include (1) regression parameters and (2) minor allele frequency. If the original profile likelihood is maximized, If the outputs also include (3) Lagrange's multipliers.

The imputs of this function include (1) disease status, (2) SNP genotype, (3) covariate(s), and (4) stratum indicator for each subject. Furthermore, the disease prevalence for each stratum is also required. If a covariate is categorical, dummy variables should be constructed before using the function.
}
\usage{
SingleSNP(Y, G, X, S, fs, par = NULL, link = "logit", modified = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Response variable (vector), which should take a value of 1 (case) or 0 (control).
}
  \item{G}{SNP genotype coded by 0, 1, or 2 according to the number of minor alleles.
}
  \item{X}{Covariate variable(s) (numerical vector for a single covariate or matrix for multiple covariates). Categorical variable should be coded as dummy variables.
}
  \item{S}{Stratum variable (vector), which should take a value of 1, 2, ..., or K, where K is the number of strata.
}
  \item{fs}{Prevalence(s) (vector). The ith entry is the prevalence for stratum i.
}
  \item{par}{Initial regression parameters (list) including alpha, betaG, betaX, and betaS, with the default value being NULL.
}
  \item{link}{Link function, which should be either 'logit' (default value) or 'probit'.
}
  \item{modified}{An indicator for modifiying the profile likelihood or not. If it is TRUE (default value), then the profile likelihood will be modified; otherwise the original profile likelihood function will be used.
}
}

\value{
  \item{LOGIT}{Standard logistic regression results including the MLEs of regression parameters and their standard errors.}
  \item{mpMLE}{Modified profile MLEs and their standard errors.}
  \item{pMLE}{The original prifile MLEs and their standard errors, available when "modified" is FALSE.}
}


\examples{

  data(data);

  attach(data);

  fs = c(0.010869, 0.000867, 0.001707);

  covariate = cbind(covariate.1,covariate.2,covariate.3);

  fit = SingleSNP(Y=status,G=genotype,X=covariate,S=stratum,fs=fs);

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

