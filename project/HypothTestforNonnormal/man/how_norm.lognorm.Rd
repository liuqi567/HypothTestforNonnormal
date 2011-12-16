\name{how_norm.lognorm}
\alias{how_norm.lognorm}
\title{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and lognormal distribution parameters.}
\usage{how_norm.lognorm(n, meanlog=0, sdlog=1, n_2=0,
    meanlog_2=0, sdlog_2=1, sim=10^4)}
\description{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and lognormal distribution parameters.}
\value{QQ plot of one sample mean or difference of two sample mean vs Normal}
\note{sdlog parameter of lognormal distribution must be big than 0}
\arguments{\item{n}{the sample size for sample one}
\item{meanlog}{the mean of the distribution on the log scale for sample one with default value of 0}
\item{sdlog}{the standard deviation of the distribution on the log scale for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{meanlog_2}{the mean of the distribution on the log scale for sample two with default value of 0}
\item{sdlog_2}{the standard deviation of the distribution on the log scale for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^4}}
\examples{how_norm.lognorm(n=5, meanlog=0, sdlog=1)
how_norm.lognorm(n=5, meanlog=-1,sdlog=2, n_2=7,meanlog_2=0, sdlog_2=1, sim=10^3)}
