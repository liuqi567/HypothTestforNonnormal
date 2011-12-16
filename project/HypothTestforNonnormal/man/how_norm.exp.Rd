\name{how_norm.exp}
\alias{how_norm.exp}
\title{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and exponential distribution rate parameter.}
\usage{how_norm.exp(n, rate=1, n_2=0, rate_2=1,
    sim=10^4)}
\description{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and exponential distribution rate parameter.}
\value{QQ plot of one sample mean or difference of two sample mean vs Normal}
\note{rate of exponential distribution must be big than 0}
\arguments{\item{n}{the sample size for sample one with default value of 1}
\item{rate}{the rate parameter for sample one}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{rate_2}{the rate parameter for sample two with default value of 1}
\item{sim}{the number of simulations with default value of 10^4}}
\examples{how_norm.exp(n=5, rate=1/2)
how_norm.exp(n=5, rate=1, n_2=7, rate_2=1/2, sim=10^3)}
