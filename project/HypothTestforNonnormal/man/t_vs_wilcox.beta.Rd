\name{t_vs_wilcox.beta}
\alias{t_vs_wilcox.beta}
\title{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Beta samples.}
\usage{t_vs_wilcox.beta(n, shape1=1, shape2=1, n_2=0,
    shape1_2=1, shape2_2=1, sim=10^3, mu=0,
    conf.level=0.95)}
\description{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Beta samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.}
\value{simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis}
\note{the shape parameters of Beta distribution must be big than 0;}
\arguments{\item{n}{the sample size for sample one}
\item{shape1}{the first shape parameter for sample one with default value of 1}
\item{shape2}{the second shape parameter for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{shape1_2}{the first shape parameter for sample two with default value of 1}
\item{shape2_2}{the second shape parameter for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^3}
\item{conf.level}{confidence level with default value of 0.95}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0}}
\examples{t_vs_wilcox.beta(n=10, shape1=1, shape2=1, n_2=10, shape1_2=1, shape2_2=1, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.beta(n=30, shape1=1, shape2=1, n_2=30, shape1_2=3, shape2_2=1, sim=10^3,mu=0, conf.level=0.95)}
