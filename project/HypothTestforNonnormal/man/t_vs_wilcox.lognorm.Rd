\name{t_vs_wilcox.lognorm}
\alias{t_vs_wilcox.lognorm}
\title{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Lognormal samples.}
\usage{t_vs_wilcox.lognorm(n, meanlog=0, sdlog=1,
    n_2=0, meanlog_2=0, sdlog_2=1, sim=10^3,
    mu=0, conf.level=0.95)}
\description{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Lognormal samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.}
\value{simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis}
\note{sdlog parameter of lognormal distribution must be big than 0}
\arguments{\item{n}{the sample size for sample one}
\item{meanlog}{the mean of the distribution on the log scale for sample one with default value of 0}
\item{sdlog}{the standard deviation of the distribution on the log scale for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{meanlog_2}{the mean of the distribution on the log scale for sample two with default value of 0}
\item{sdlog_2}{the standard deviation of the distribution on the log scale for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^3}
\item{conf.level}{confidence level with default value of 0.95}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0}}
\examples{t_vs_wilcox.lognorm(n=10, meanlog=0, sdlog=1, n_2=10, meanlog_2=0, sdlog_2=1,sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.lognorm(n=30, meanlog=0, sdlog=1, n_2=30, meanlog_2=4, sdlog_2=1,sim=10^3,mu=0, conf.level=0.95)}
