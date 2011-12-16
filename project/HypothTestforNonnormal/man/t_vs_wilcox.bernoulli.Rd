\name{t_vs_wilcox.bernoulli}
\alias{t_vs_wilcox.bernoulli}
\title{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Bernoulli samples.}
\usage{t_vs_wilcox.bernoulli(n, prob=1/2, n_2=0,
    prob_2=1/2, sim=10^3, mu=0, conf.level=0.95)}
\description{Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Bernoulli samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.}
\value{simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis}
\note{0< prob <1 It may produce warnings since Bernoulli sample has large portion of ties.}
\arguments{\item{n}{the sample size for sample one}
\item{prob}{the probability of success for sample one with default value of 1/2}
\item{n_2}{the sample size for sample two with default value of 0 (for one sample test case)}
\item{prob_2}{the probability of success for sample two with default value of 1/2}
\item{sim}{the number of simulations with default value of 10^3}
\item{conf.level}{confidence level with default value of 0.95}
\item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0}}
\examples{t_vs_wilcox.bernoulli(n=10, prob=1/2, n_2=10, prob_2=1/2, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.bernoulli(n=50, prob=0.25,  n_2=50, prob_2=0.75, sim=10^3,mu=0, conf.level=0.95)}
