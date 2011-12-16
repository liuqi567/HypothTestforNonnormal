\name{how_norm.cauchy}
\alias{how_norm.cauchy}
\title{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Cauchy distribution parameters.}
\usage{how_norm.cauchy(n, location=0, scale=1, n_2=0,
    location_2=0, scale_2=1, sim=10^4)}
\description{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Cauchy distribution parameters.}
\value{QQ plot of one sample mean or difference of two sample mean vs Normal}
\note{scale parameter of Cauchy distribution must be big than 0}
\arguments{\item{n}{the sample size for sample one}
\item{location}{the location parameter for sample one with default value of 0}
\item{scale}{the scale parameter for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{location_2}{the location parameter for sample two with default value of 0}
\item{scale_2}{the scale parameter for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^4}}
\examples{how_norm.cauchy(n=5, location=-1, scale=1)
how_norm.cauchy(n=5, location=-1,scale=2, n_2=7,location_2=0, scale_2=1, sim=10^3)}
