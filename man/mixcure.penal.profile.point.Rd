\name{mixcure.penal.profile.point}
\alias{mixcure.penal.profile.point}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D data points generating for interactive Profile likelihood confidence surface under MC
}
\description{
A function to generate 3D data points for interactive Profile likelihood confidence surface under MC. Called under plcr.surface() function.
}
\usage{
mixcure.penal.profile.point(formula, data, init, k, kval, pl, iterlim = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
MC model specification in the form of a conventional survival model, i.e., Surv()~. Note the whatever x variables specified in the latency part of the model are assumed to be also regressed on in the incidence part, and vice versa.
}
  \item{data}{
  a data object in the form of dataframe.
  }
  \item{k}{
  identify the order number of single covariate in the specified MC model.
  }
    \item{kval}{
  a two-element vector (i,j) that points to the value of point (logOR, logHR) to be considered for estimating the profile likelihood function.
  }
  \item{init}{
  a vector of initial values input for the optimization of the vector of parameters specified in the likelihood. Usually, specified in the order of 1) regression parameters for the incidence part, 2) regression parameters for the latency part, and 3) the shape parameter of Weibull hazard function for the latency part.
  }
  \item{pl}{
  if TRUE, uses FT-PL; otherwise uses the usual likelihood for the model.
  }
  \item{iterlim}{
  identify uplimt of iteration numbers for likelihood optimization under nlm().
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns a list object, in which the first item is a table containing the matrix of covariate number, logOR, logHR and loglikelihood value.
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
## For the first covariate (second counting the intercept) of MC model
## To find out the points for PLCR surface distribution between (-0.5, 0.5) for its logOR and (-1, 1) for its logHR

for (i in seq(-0.5, 0.5, by = 0.05)) {
  for (j in seq(-1, 1, by = 0.1)) {

    single.prof.point <- mixcure.cr::mixcure.penal.profile.point(formula=formula,
                                                     data = indata, k = covar+1,
                                                     kval = c(i,j), init=init, pl=pl)
    outdat <- rbind.data.frame(outdat, single.prof.point$coefficients[[1]][1:4])

  }
}
colnames(outdat) <- c("var", "kval1","kval2","loglik")
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
