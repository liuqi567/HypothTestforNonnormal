\name{how_norm.gamma}
\alias{how_norm.gamma}
\title{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Gamma distribution parameters.}
\usage{how_norm.gamma(n, shape=1, scale=1, n_2=0,
    shape_2=1, scale_2=1, sim=10^4)}
\description{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Gamma distribution parameters.}
\value{QQ plot of one sample mean or difference of two sample mean vs Normal}
\note{the shape and scale parameters of Gamma distribution must be big than 0; The Gamma distribution with parameters shape = a and scale = s has density f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s); E(X) = a*s and Var(X) = a*s^2}
\arguments{\item{n}{the sample size for sample one}
\item{shape}{the shape parameter for sample one with default value of 1}
\item{scale}{the scale parameter for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{shape_2}{the shape parameter for sample two with default value of 1}
\item{scale_2}{the scale parameter for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^4}}
\examples{how_norm.gamma(n=5, shape=2, scale=1)
how_norm.gamma(n=5, shape=1,scale=2, n_2=7,shape_2=1, scale_2=1, sim=10^4)}
