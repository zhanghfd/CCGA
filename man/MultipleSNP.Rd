\name{MultipleSNP}
\alias{MultipleSNP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of log-ORs (SEs) and significance test p-values for multiple SNPs.
}
\description{
This function returns (1) log-OR estimates and the corresponding (2) standard errors and the significance test (3) p-values for multiple SNPs. The modified profile likelihood is used. Parallel computation can be used to speed up the analysis through specification of more than one CPU cores. Because the function "mcapply" in the R package "parallel" is used, the core number can be greater than 1 in UNIX-like OS but it can be only 1 in Windows OS.

If a covariate is categorical, dummy variables should be constructed before using the function.
}
\usage{
MultipleSNP(Gs, Y, Z, S, fs, par = NULL, link = "logit", modified=TRUE, cl.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Gs}{SNP genotype matrix coded by 0, 1, or 2 according to the number of minor alleles, one colume per one SNP.
}
  \item{Y}{Response variable (vector), which should take a value of 1 (case) or 0 (control).
}

  \item{Z}{Covariate variable(s) (numerical vector for a single covariate or matrix for multiple covariates). Categorical variable should be coded as dummy variables.
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
  \item{cl.cores}{CPU cores used, with a defaul value of 1.
}
}


\examples{

  data(data);

  status = data[,1];

  stratum = data[,2];

  covariate = data[,3:5];

  Gs = data[,-(1:5)];

  fs = c(0.010869, 0.000867, 0.001707);

  ## not run
  # res = MultipleSNP(Gs,Y=status,Z=covariate,S=stratum,fs=fs,cl.cores=1)

}
