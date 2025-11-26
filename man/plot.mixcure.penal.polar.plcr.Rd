\name{plot.mixcure.penal.polar.plcr.R}
\alias{plot.mixcure.penal.polar.plcr.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting profile likelihood confidence region (PLCR) estimation by radial-angular method
}
\description{
A function that uses the output from mixcure.penal.polar.plcr() to plot PLCR by Radial-angular method with the option of overlaying it with marginal PLCIs for the same covariates. The default is plotting PLCR, and Wald-CR can also be plotted as an alternative.
}
\usage{
plot.mixcure.penal.polar.plcr.R(obj_cr, cr_type="PLCR", title_name =NULL, plci=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj_cr}{
The output object from  mixcure.penal.polar.plcr().
}
  \item{cr_type}{
If plotting PLCRs, input "PLCR", else specify "Wald" for plotting Wald-CRs.
}
  \item{title_name}{
If title's name needs to be specified, can be filled here.
}
  \item{plci}{
If NULL is specified, no overlay of PLCIs will be plotted, otherwise PLCIs will be plotted.
}}
\details{
For categorical variables of 3 or more categories, pre-treatment for transforming it to binary dummy variables are needed, when being specified in the formula.
}
\value{
This ggplot based plotting function returns the plots of confidence region of logOR and logHR of the same covariates.
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

# Calculate the PLCR of the MC model parameters using maximum likelihoodlizards.
mc.plcr.mle5 <- mixcure.penal.polar.plcr(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,data=ANNbcBMdat5,init=c(5,rep(0,5),-10,rep(0,5),0.1), pl=F)

plot.mixcure.penal.polar.plcr.R(mc.plcr.mle5, cr_type="PLCR", title_name =NULL, plci=NULL)
plot.mixcure.penal.polar.plcr.R(mc.plcr.mle5, cr_type="Wald", title_name =NULL, plci=NULL)

# Calculate the PLCR of the MC model parameters using Firth-type penalized likelihoodlizards.
mc.plcr.ple5 <- mixcure.penal.polar.plcr(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,data=ANNbcBMdat5,init=c(5,rep(0,5),-10,rep(0,5),0.1), pl=T)

plot.mixcure.penal.polar.plcr.R(mc.plcr.ple5, cr_type="PLCR", title_name =NULL, plci=NULL)
plot.mixcure.penal.polar.plcr.R(mc.plcr.ple5, cr_type="Wald", title_name =NULL, plci=NULL)


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
