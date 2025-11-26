\name{plcr.surface}
\alias{plcr.surface}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D interactive Profile likelihood confidence surface for a covariate under MC
}
\description{
A function to generate the interactive 3D Profile likelihood confidence surface of the joint confidence region (PLCR) estimation of logOR and logHR of the same covariate under MC model
}
\usage{
plcr.surface(formula, covar,indata, init, pl, alpha_lo=-2, alpha_hi=2, beta_lo=-2, beta_hi=2, a_step=0.1, b_step=0.1, title_name=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
MC model specification in the form of a conventional survival model, i.e., Surv()~. Note the whatever x variables specified in the latency part of the model are assumed to be also regressed on in the incidence part, and vice versa.
}
  \item{data}{
  a data object in the form of dataframe.
  }
  \item{covar}{
  identify the order number of single covariate in the specified MC model.
  }
  \item{init}{
  a vector of initial values input for the optimization of the vector of parameters specified in the likelihood. Usually, specified in the order of 1) regression parameters for the incidence part, 2) regression parameters for the latency part, and 3) the shape parameter of Weibull hazard function for the latency part.
  }
  \item{pl}{
  if TRUE, uses FT-PL; otherwise uses the usual likelihood for the model.
  }
  \item{alpha_lo}{
  defined lower bound for logOR.
  }
  \item{alpha_hi}{
  defined upper bound for logOR.
  }
  \item{beta_lo}{
  defined lower bound for logHR.
  }
  \item{beta_hi}{
  defined upper bound for logHR.
  }
  \item{a_step}{
  defined incremental step size for logOR.
  }
  \item{b_step}{
  defined incremental step size for logHR.
  }
  \item{title_name}{
  Title of the plot. Note default is NULL.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns a 3D PLCR interactive surface plot in a html file for visualization.
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
data(ANNbcBMdat5)
# Plot the 3D PLCR of single covariate Her2 under MC model parameters using maximum likelihoodlizards.
mixcure.penal.polar.plcr(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,
covar = 1, data = ANNbcBMdat5,init = c(5,rep(0,5),-10,rep(0,5),0.1), pl=F)

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
