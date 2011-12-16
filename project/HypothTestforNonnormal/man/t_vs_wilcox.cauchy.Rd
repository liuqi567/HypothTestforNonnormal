\name{t_vs_wilcox.cauchy}
\alias{t_vs_wilcox.cauchy}
\title{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Cauchy samples.}
\usage{t_vs_wilcox.cauchy(n, location=0, scale=1,
    n_2=0, location_2=0, scale_2=1, sim=10^3,
    mu=0, conf.level=0.95)}
\description{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Cauchy samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.}
\value{simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis}
\note{t-test might perform badly on Cauchy samples}
\arguments{\item{n}{the sample size for sample one}
\item{location}{the location parameter for sample one with default value of 0}
\item{scale}{the scale parameter for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{location_2}{the location parameter for sample two with default value of 0}
\item{scale_2}{the scale parameter for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^3}
\item{conf.level}{confidence level with default value of 0.95}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0}}
\examples{t_vs_wilcox.cauchy(n=10, location=0, scale=1, n_2=10, location_2=0, scale_2=1,sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.cauchy(n=30, location=0, scale=1, n_2=30, location_2=4, scale_2=1,sim=10^3,mu=0, conf.level=0.95)}
