R statistical language function package mainly provides bivariate CR estimation of logOR and logHR from Weibull PH mixture cure model 
under Firth type penalized likelihood, by estimating the edge points for and plotting the CR. The primary bivariate CR is 
profile likelihood-based cR with the area of CR also being calculated, while the Wald-type CR can also be calculated. This package also provide function to calculate 
and plot confidence intervals under profile likelihood and Wald. For parameter estimation, both Firth-type penalized 
likelihood and standard maximum likelihood are included.
Created by Changchang Xu

Contact:changchang.xu@alumni.utoronto.ca

This package can be installed via the following R code:

devtools::install_github("ChangchangXu-LTRI2025/mixcuref.cr", build = TRUE, build_opts = c())

library(mixcuref.cr)
