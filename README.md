# CCGA
CCGA: Case-control genetic association analysis incorporating non-confounding covariates

Adjusting for non-confounding covariates in the standard case-control genetic association analyses may lead to decreased power when the phenotype is rare. On the other hand, the existing methods not adjusting for non-confounding covariates may be less powerful than those with adjustment. In the R package CCGA (short for Case-Control Genetic Association), a novel method is implemented, which exploits external phenotype prevalence data, Hardy-Weinberg equilibrium, and independence between non-confounding covariates and interested SNP. CCGA can fit two commonly used penetrance models, i.e., the logistic and probit regression models. It guarantees increased power through covariate adjustment regardless of whether the phenotype is rare or common.

The manual file is "CCGA-manual.pdf". 

Installation of CCGA in R:

> library(‘devtools’);

> install_github(‘zhanghfd/CCGA’);
