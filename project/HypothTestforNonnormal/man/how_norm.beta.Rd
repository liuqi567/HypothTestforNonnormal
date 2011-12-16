\name{how_norm.beta}
\alias{how_norm.beta}
\title{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Beta distribution parameters.}
\usage{how_norm.beta(n, shape1=1, shape2=1, n_2=0,
    shape1_2=1, shape2_2=1, sim=10^4)}
\description{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Beta distribution parameters.}
\value{QQ plot of one sample mean or difference of two sample mean vs Normal}
\note{the shape parameters of Beta distribution must be big than 0;}
\arguments{\item{n}{the sample size for sample one}
\item{shape1}{the first shape parameter for sample one with default value of 1}
\item{shape2}{the second shape parameter for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{shape1_2}{the first shape parameter for sample two with default value of 1}
\item{shape2_2}{the second shape parameter for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^4}}
\examples{how_norm.beta(n=5, shape1=2, shape2=1)}
