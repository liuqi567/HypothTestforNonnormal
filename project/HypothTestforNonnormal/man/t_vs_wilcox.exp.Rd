\name{t_vs_wilcox.exp}
\alias{t_vs_wilcox.exp}
\title{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Exponential samples.}
\usage{t_vs_wilcox.exp(n, rate, n_2=0, rate_2=1, mu=0,
    sim=10^3, conf.level=0.95)}
\description{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Exponential samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.}
\value{simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis}
\note{rate of exponential distribution must be big than 0}
\arguments{\item{n}{the sample size for sample one with default value of 1}
\item{rate}{the rate parameter for sample one}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{rate_2}{the rate parameter for sample two with default value of 1}
\item{sim}{the number of simulations with default value of 10^3}
\item{conf.level}{confidence level with default value of 0.95}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0}}
\examples{t_vs_wilcox.exp(n=30,rate=1,n_2=30,rate_2=1, sim=10^3, conf.level=0.95)
t_vs_wilcox.exp(n=30,rate=0.1,n_2=30,rate_2=1, sim=10^3)}
