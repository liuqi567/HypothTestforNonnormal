\name{how_norm.bernoulli}
\alias{how_norm.bernoulli}
\title{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Bernoulli distribution parameter.}
\usage{how_norm.bernoulli(n, prob=1/2, n_2=0,
    prob_2=1/2, sim=10^4)}
\description{Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Bernoulli distribution parameter.}
\value{QQ plot of one sample mean or difference of two sample mean vs Normal}
\note{the prob parameter of Bernoulli distribution must be between 0 and 1}
\arguments{\item{n}{the sample size for sample one}
\item{prob}{the probability of success for sample one with default value of 1/2}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{prob_2}{the probability of success for sample two with default value of 1/2}
\item{sim}{the number of simulations with default value of 10^4}}
