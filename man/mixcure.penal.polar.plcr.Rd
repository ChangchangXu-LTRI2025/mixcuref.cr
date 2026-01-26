\name{mixcure.penal.polar.plcr}
\alias{mixcure.penal.polar.plcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Profile likelihood confidence region estimation by radial-angular method
}
\description{
Radial-angular method estimated profile likelihood confidence region edge points for logOR and logHR of the same covariate under the Weibull PH MC model, with the option of Firth-type penalized likelihood estimation or standard maximum likelihood estimation.
}
\usage{
mixcure.penal.polar.plcr(formula, data, apct = 0.05,init, pl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
MC model specification in the form of a conventional survival model, i.e., Surv()~. Note the whatever x variables specified in the latency part of the model are assumed to be also regressed on in the incidence part, and vice versa.
}
  \item{data}{
  a data object in the form of dataframe.
  }
  \item{init}{
  a vector of initial values input for the optimization of the vector of parameters specified in the likelihood. Usually, specified in the order of 1) regression parameters for the incidence part, 2) regression parameters for the latency part, and 3) the shape parameter of Weibull hazard function for the latency part.
  }
  \item{pl}{
  if TRUE, uses FT-PL; otherwise uses the usual likelihood for the model.
  }
  \item{apct}{
  specifies the alpha level of the bivariate confidence region, default is 0.05.
  }
}
\details{
For categorical variables of 3 or more categories, pre-treatment for transforming it to binary dummy variables are needed, when being specified in the formula.
}
\value{
The mixcure.penal.est function returns a list of three sub-lists containing:
  \item{coefficients}{contains 3 tables of parameter estimates with standard errors, as well as Wald-type test statistics, p values and interval estimations: i) 'CURE' table for regression parameters of the incidence part; ii) 'SURVIVAL' table for regression parameters of the latency part; iii) 'ALPHA' table for shape parameter of the Weibull baseline hazard function in the latency part.}
  \item{plcr}{containg a dataframe of radial specific diameters of all covariates as well as the varying values of radius from calculation of radial-angular PLCR calculation}
  \item{cov}{variance-covariance matrix of all the parameters in the order of incidence part, latency part and shape parameter.}

}
\references{
20.	Jaeger A. Computation of Two- and Three-Dimensional Confidence Regions with the Likelihood Ratio. The American Statistician (Statistical Computing and Graphics). 2016; 70: 395-398.
}
\author{
Changchang Xu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Low event rate, 5 variable
data(ANNbcBMdat5)
# Calculate the PLCR of the MC model parameters using maximum likelihoodlizards.
mc.plcr.mle5 <- mixcure.penal.polar.plcr(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,data=ANNbcBMdat5,init=c(5,rep(0,5),-10,rep(0,5),0.1), pl=F)
# Calculate the PLCR of the MC model parameters using Firth-type penalized likelihoodlizards.
mc.plcr.ple5 <- mixcure.penal.polar.plcr(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,data=ANNbcBMdat5,init=c(5,rep(0,5),-10,rep(0,5),0.1), pl=T)

mc.plcr.mle5$coefficients
mc.plcr.mle5$plcr
mc.plcr.ple55$coefficients
mc.plcr.ple5$plcr

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
