\name{t_vs_wilcox.weibull}
\alias{t_vs_wilcox.weibull}
\title{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Weibull samples.}
\usage{t_vs_wilcox.weibull(n, shape=1, scale=1, n_2=0,
    shape_2=1, scale_2=1, sim=10^3, mu=0,
    conf.level=0.95)}
\description{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Weibull samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.}
\value{simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis}
\note{the shape and scale parameters of Weibull distribution must be big than 0; The Weibull distribution with shape parameter a and scale parameter b has density given by f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a) for x>0}
\arguments{\item{n}{the sample size for sample one}
\item{shape}{the shape parameter for sample one with default value of 1}
\item{scale}{the scale parameter for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{shape_2}{the shape parameter for sample two with default value of 1}
\item{scale_2}{the scale parameter for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^3}
\item{conf.level}{confidence level with default value of 0.95}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0}}
\examples{t_vs_wilcox.weibull(n=10, shape=1, scale=1, n_2=10, shape_2=1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.weibull(n=50, shape=10, scale=5, n_2=50, shape_2=0.1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)}
