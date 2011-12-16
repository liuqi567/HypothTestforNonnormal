\name{t_vs_wilcox.poisson}
\alias{t_vs_wilcox.poisson}
\title{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Poisson samples.}
\usage{t_vs_wilcox.poisson(n, lambda=1, n_2=0,
    lambda_2=1, sim=10^3, mu=0, conf.level=0.95)}
\description{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Poisson samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.}
\value{simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis}
\note{lambda>0 It may produce warnings since Poisson sample may have ties.}
\arguments{\item{n}{the sample size for sample one}
\item{lambda}{the mean for sample one with default value of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{lambda_2}{the mean for sample two with default value of 1}
\item{sim}{the number of simulations with default value of 10^3}
\item{conf.level}{confidence level with default value of 0.95}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0}}
\examples{t_vs_wilcox.poisson(n=10, lambda=1, n_2=10, lambda_2=1, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.poisson(n=30, lambda=1,  n_2=30, lambda_2=0.2, sim=10^3,mu=0, conf.level=0.95)}
