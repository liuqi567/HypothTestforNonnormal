\name{t_vs_wilcox.gamma}
\alias{t_vs_wilcox.gamma}
\title{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Gamma samples.}
\usage{t_vs_wilcox.gamma(n, shape=1, scale=1, n_2=0,
    shape_2=1, scale_2=1, sim=10^3, mu=0,
    conf.level=0.95)}
\description{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Gamma samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.}
\value{simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis}
\note{shape>0, scale>0 The Gamma distribution with parameters shape = a and scale = s has density f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s); E(X) = a*s and Var(X) = a*s^2}
\arguments{\item{n}{the sample size for sample one}
\item{shape}{the shape parameter for sample one with default value of 1}
\item{scale}{the scale parameter for sample one with default values of 1}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{shape_2}{the shape parameter for sample two with default value of 1}
\item{scale_2}{the scale parameter for sample two with default values of 1}
\item{sim}{the number of simulations with default value of 10^3}
\item{conf.level}{confidence level with default value of 0.95}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0}}
\examples{t_vs_wilcox.gamma(n=10, shape=1, scale=1, n_2=10, shape_2=1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.gamma(n=30, shape=1, scale=1, n_2=30, shape_2=5, scale_2=1, sim=10^3,mu=0, conf.level=0.95)}
