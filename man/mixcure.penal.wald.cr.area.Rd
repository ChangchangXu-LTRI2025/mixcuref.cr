\name{mixcure.penal.wald.cr.area}
\alias{mixcure.penal.wald.cr.area}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wald-type 2d test statistic and confidence region area
}
\description{
A function that calculates Wald-type 2d test statistic, p values and confidence region areas for each of the covariates under MC model
}
\usage{
mixcure.penal.wald.cr.area(formula, data, apct = 0.05, init, true.val, pl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  MC model specification in the form of a conventional survival model, i.e., Surv()~. Note the whatever x variables specified in the latency part of the model are assumed to be also regressed on in the incidence part, and vice versa.
  }
  \item{data}{
  a data object in the form of dataframe.
  }
  \item{apct}{
  identify the significance level for Wald test statistic.
  }
  \item{init}{
  a vector of initial values input for the optimization of the vector of parameters specified in the likelihood. Usually, specified in the order of 1) regression parameters for the incidence part, 2) regression parameters for the latency part, and 3) the shape parameter of Weibull hazard function for the latency part.
  }
  \item{true.val}{
  identify the null values for constructing the Wald test statistic.
  }
  \item{pl}{
  if TRUE, uses FT-PL; otherwise uses the usual likelihood for the model.
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a table of parameter estimates, se, covariance, 2 df chisq value, 2df p value, and CR area.
}
\references{
%% ~put references to the literature/web site here ~
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
# Calculate the Wald test statistic under the null values, 2df p value and CR area of the MC model parameters using maximum likelihoodlizards.
mc.plcr.mle5 <- mixcure.penal.wald.cr.area(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,
data=ANNbcBMdat5,
true.val = c(rep(0,5),rep(0,5)), pl=F)

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
