\name{mixcure.penal.plcr.area}
\alias{mixcure.penal.plcr.area}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Area calculation of profile likelihood confidence region under MC
}
\description{
A function calculates the area of PLCR for each of all the covariates under the specified MC model
}
\usage{
mixcure.penal.plcr.area(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
is specified as the output of mixcure.penal.polar.plcr() function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns a table listing the covariate names and corresponding area of PLCR for each of the covariate
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

mc.plcr.mle5 <- mixcure.penal.polar.plcr(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,data=ANNbcBMdat5,init=c(5,rep(0,5),-10,rep(0,5),0.1), pl=F)
mixcure.penal.plcr.area(mc.plcr.mle5)

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
