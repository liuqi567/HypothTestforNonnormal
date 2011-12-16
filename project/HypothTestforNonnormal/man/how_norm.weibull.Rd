\name{how_norm.weibull}
\alias{how_norm.weibull}
\title{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Weibull distribution parameters.}
\usage{how_norm.weibull(n, shape=1, scale=1, n_2=0,
    shape_2=1, scale_2=1, sim=10^4)}
\description{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Weibull distribution parameters.}
\value{QQ plot of one sample mean or difference of two sample mean vs Normal}
\note{the shape and scale parameters of Weibull distribution must be big than 0; The Weibull distribution with shape parameter a and scale parameter b has density given by f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a) for x>0}
\arguments{\item{n}{the sample size for sample one}
\item{shape}{the location parameter for sample one with default value of 0}
\item{scale}{the scale parameter for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{shape_2}{the location parameter for sample two with default value of 0}
\item{scale_2}{the scale parameter for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^4}}
\examples{how_norm.weibull(n=5, shape=2, scale=1)
how_norm.weibull(n=5, shape=1,scale=2, n_2=7,shape_2=1, scale_2=1, sim=10^4)}
