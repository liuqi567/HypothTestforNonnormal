\name{rank_t.test}
\alias{rank_t.test}
\title{Transforms data to rank, then performs one or two sample t_tests on rank.}
\usage{rank_t.test(x, y, ...)}
\description{Transforms data to rank, then performs one or two sample t_tests on rank.}
\value{A list with class "htest" containing the value of t-statisc, parameter, p-value, confidence interval, estimate, null hypothesis, alternative hypothesis, method, and data name.}
\arguments{\item{x}{a (non-empty) numeric vector of data values}
\item{y}{an optional (non-empty) numeric vector of data values}
\item{...}{further arguments to be passed to t-test}}
\examples{rank_t.test(1:10, y=c(7:20))}
