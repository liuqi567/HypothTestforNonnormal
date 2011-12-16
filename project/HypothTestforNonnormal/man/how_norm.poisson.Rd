\name{how_norm.poisson}
\alias{how_norm.poisson}
\title{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Poisson distribution parameter.}
\usage{how_norm.poisson(n, lambda=1, n_2=0, lambda_2=1,
    sim=10^4)}
\description{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Poisson distribution parameter.}
\value{QQ plot of one sample mean or difference of two sample mean vs Normal}
\note{the lambda parameter of Poisson distribution must be big than 0}
\arguments{\item{n}{the sample size for sample one}
\item{lambda}{the mean for sample one with default value of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{lambda_2}{the mean for sample two with default value of 1}
\item{sim}{the number of simulations with default value of 10^4}}
\examples{how_norm.poisson(n=5, lambda=1/2)
how_norm.poisson(n=5, lambda=1,n_2=7, lambda_2=0.2, sim=10^4)}
