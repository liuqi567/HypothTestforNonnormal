library(roxygen)
#' A package to assess the normality of test statistcs when sample from common nonnormal disttributions
#'
#' \tabular{11}{
#' Package:\tab HypothTestforNonnormal\cr
#' Type:\tab Package\cr
#' Version: \tab 0.1\cr
#' Date:\tab 2011-11-15\cr
#' License: \tab GPL (>=2)\cr
#' Lazyload: \tab yes\cr
#' }
#'
#' how_norm functions simulate sample mean from given distribution and give QQ plot vs Normal
#' 
#' @name  HypothTestforNonnormal
#' @aliases HypothTestforNonnormal
#' @docType package
#' @title Assess the normality of of test statistcs when sample from nonnormal disttributions
#' @author Qi Liu \email{qi.liu.1@@vanderbilt.edu}
#' @keywords package
#' @examples
#' how_norm.exp(n=5, rate=1, n_2=9, rate_2=1/2,sim=10^4)
#' how_norm.cauchy(n=1000,location=0, scale=1)
roxygen()

#' Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and exponential distribution rate parameter.
#' 
#' @param n the sample size for sample one with default value of 1
#' @param rate the rate parameter for sample one
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param rate_2 the rate parameter for sample two with default value of 1
#' @param sim the number of simulations with default value of 10^4
#' @return QQ plot of one sample mean or difference of two sample mean vs Normal
#' @note rate of exponential distribution must be big than 0
#' @examples
#' how_norm.exp(n=5, rate=1/2)
#' how_norm.exp(n=5, rate=1, n_2=7, rate_2=1/2, sim=10^3)

how_norm.exp<- function(n, rate=1, n_2=0, rate_2=1, sim=10^4){
  TS <- rep(NA, sim)
  x_bar <- rep(NA, sim)
  for(i in 1:sim){
    x <- rexp(n, rate)
    x_bar[i] <- mean(x)
    }
  
  if(n_2) {
    y_bar <- rep(NA, sim)
    for(i in 1:sim) {
      y <- rexp(n_2, rate_2)
      y_bar[i] <- mean(y)
    }
    x_bar <- x_bar-y_bar
  }
  TS <- x_bar
  qqnorm(TS, main="Q-Q Plot: Sampling Distribution from Exponential vs. Normal", xlab="Normal Quantiles", 
         ylab="Sample Quantiles of Mean")
}

#' Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and lognormal distribution parameters.
#' 
#' @param n the sample size for sample one
#' @param meanlog the mean of the distribution on the log scale for sample one with default value of 0
#' @param sdlog the standard deviation of the distribution on the log scale for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param meanlog_2 the mean of the distribution on the log scale for sample two with default value of 0
#' @param sdlog_2 the standard deviation of the distribution on the log scale for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^4
#' @return QQ plot of one sample mean or difference of two sample mean vs Normal
#' @note sdlog parameter of lognormal distribution must be big than 0
#' @examples
#' how_norm.lognorm(n=5, meanlog=0, sdlog=1)
#' how_norm.lognorm(n=5, meanlog=-1,sdlog=2, n_2=7,meanlog_2=0, sdlog_2=1, sim=10^3)
how_norm.lognorm<- function(n, meanlog=0, sdlog=1, n_2=0, meanlog_2=0, sdlog_2=1,sim=10^4){
  TS <- rep(NA, sim)
  x_bar <- rep(NA, sim)
  for(i in 1:sim){
    x <- rlnorm(n, meanlog, sdlog)
    x_bar[i] <- mean(x)
    }
  
  if(n_2) {
    y_bar <- rep(NA, sim)
    for(i in 1:sim) {
      y <- rlnorm(n_2, meanlog_2, sdlog_2)
      y_bar[i] <- mean(y)
    }
    x_bar <- x_bar-y_bar
  }
  
  TS <- x_bar
  qqnorm(TS, main="Q-Q Plot: Sampling Distribution of Mean from Lognorm vs. Normal", xlab="Normal Quantiles", 
         ylab="Sample Quantiles of Mean")
}


#' Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Cauchy distribution parameters.
#' 
#' @param n the sample size for sample one
#' @param location the location parameter for sample one with default value of 0
#' @param scale the scale parameter for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param location_2 the location parameter for sample two with default value of 0
#' @param scale_2 the scale parameter for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^4
#' @return QQ plot of one sample mean or difference of two sample mean vs Normal
#' @note scale parameter of Cauchy distribution must be big than 0
#' @examples
#' how_norm.cauchy(n=5, location=-1, scale=1)
#' how_norm.cauchy(n=5, location=-1,scale=2, n_2=7,location_2=0, scale_2=1, sim=10^3)
how_norm.cauchy<- function(n, location=0, scale=1, n_2=0, location_2=0, scale_2=1,sim=10^4){
  TS <- rep(NA, sim)
  x_bar <- rep(NA, sim)
  for(i in 1:sim){
    x <- rcauchy(n, location, scale)
    x_bar[i] <- mean(x)
    }
  
  if(n_2) {
    y_bar <- rep(NA, sim)
    for(i in 1:sim) {
      y <- rcauchy(n_2, location_2, scale_2)
      y_bar[i] <- mean(y)
    }
    x_bar <- x_bar-y_bar
  }
  
  TS <- x_bar
  qqnorm(TS, main="Q-Q Plot: Sampling Distribution of Mean from Cauchy vs. Normal", xlab="Normal Quantiles", 
         ylab="Sample Quantiles of Mean")
}

#' Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Bernoulli distribution parameter.
#' 
#' @param n the sample size for sample one
#' @param prob the probability of success for sample one with default value of 1/2
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param prob_2 the probability of success for sample two with default value of 1/2
#' @param sim the number of simulations with default value of 10^4
#' @return QQ plot of one sample mean or difference of two sample mean vs Normal
#' @note the prob parameter of Bernoulli distribution must be between 0 and 1

how_norm.bernoulli<- function(n, prob=1/2, n_2=0, prob_2=1/2,sim=10^4){
  TS <- rep(NA, sim)
  x_bar <- rep(NA, sim)
  for(i in 1:sim){
    x <- rbinom(n, size=1, prob)
    x_bar[i] <- mean(x)
    }
  
  if(n_2) {
    y_bar <- rep(NA, sim)
    for(i in 1:sim) {
      y <- rbinom(n_2, size=1, prob_2)
      y_bar[i] <- mean(y)
    }
    x_bar <- x_bar-y_bar
  }
  
  TS <- x_bar
  qqnorm(TS, main="Q-Q Plot: Sampling Distribution of Mean from Bernoulli vs. Normal", xlab="Normal Quantiles", 
         ylab="Sample Quantiles of Mean")
}

#' Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Poisson distribution parameter.
#' 
#' @param n the sample size for sample one
#' @param lambda the mean for sample one with default value of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param lambda_2 the mean for sample two with default value of 1
#' @param sim the number of simulations with default value of 10^4
#' @return QQ plot of one sample mean or difference of two sample mean vs Normal
#' @note the lambda parameter of Poisson distribution must be big than 0 
#' @examples
#' how_norm.poisson(n=5, lambda=1/2)
#' how_norm.poisson(n=5, lambda=1,n_2=7, lambda_2=0.2, sim=10^4)
how_norm.poisson<- function(n, lambda=1, n_2=0, lambda_2=1,sim=10^4){
  TS <- rep(NA, sim)
  x_bar <- rep(NA, sim)
  for(i in 1:sim){
    x <- rpois(n, lambda)
    x_bar[i] <- mean(x)
    }
  
  if(n_2) {
    y_bar <- rep(NA, sim)
    for(i in 1:sim) {
      y <- rpois(n_2, lambda_2)
      y_bar[i] <- mean(y)
    }
    x_bar <- x_bar-y_bar
  }
  
  TS <- x_bar
  qqnorm(TS, main="Q-Q Plot: Sampling Distribution of Mean from Poisson vs. Normal", xlab="Normal Quantiles", 
         ylab="Sample Quantiles of Mean")
}

#' Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Gamma distribution parameters.
#' 
#' @param n the sample size for sample one
#' @param shape the shape parameter for sample one with default value of 1
#' @param scale the scale parameter for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param shape_2 the shape parameter for sample two with default value of 1
#' @param scale_2 the scale parameter for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^4
#' @return QQ plot of one sample mean or difference of two sample mean vs Normal
#' @note the shape and scale parameters of Gamma distribution must be big than 0; The Gamma distribution with parameters shape = a and scale = s has density f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s); E(X) = a*s and Var(X) = a*s^2
#' @examples
#' how_norm.gamma(n=5, shape=2, scale=1)
#' how_norm.gamma(n=5, shape=1,scale=2, n_2=7,shape_2=1, scale_2=1, sim=10^4)
how_norm.gamma<- function(n, shape=1, scale=1, n_2=0, shape_2=1, scale_2=1,sim=10^4){
  TS <- rep(NA, sim)
  x_bar <- rep(NA, sim)
  for(i in 1:sim){
    x <- rgamma(n, shape, scale)
    x_bar[i] <- mean(x)
    }
  
  if(n_2) {
    y_bar <- rep(NA, sim)
    for(i in 1:sim) {
      y <- rgamma(n_2, shape_2, scale_2)
      y_bar[i] <- mean(y)
    }
    x_bar <- x_bar-y_bar
  }
  
  TS <- x_bar
  qqnorm(TS, main="Q-Q Plot: Sampling Distribution of Mean from Gamma vs. Normal", xlab="Normal Quantiles", 
         ylab="Sample Quantiles of Mean")
}

#' Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Beta distribution parameters.
#' 
#' @param n the sample size for sample one
#' @param shape1 the first shape parameter for sample one with default value of 1
#' @param shape2 the second shape parameter for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param shape1_2 the first shape parameter for sample two with default value of 1
#' @param shape2_2 the second shape parameter for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^4
#' @return QQ plot of one sample mean or difference of two sample mean vs Normal
#' @note the shape parameters of Beta distribution must be big than 0; 
#' @examples
#' how_norm.beta(n=5, shape1=2, shape2=1)

how_norm.beta<- function(n, shape1=1, shape2=1, n_2=0, shape1_2=1, shape2_2=1,sim=10^4){
  TS <- rep(NA, sim)
  x_bar <- rep(NA, sim)
  for(i in 1:sim){
    x <- rbeta(n, shape1, shape2)
    x_bar[i] <- mean(x)
    }
  
  if(n_2) {
    y_bar <- rep(NA, sim)
    for(i in 1:sim) {
      y <- rbeta(n_2, shape1_2, shape2_2)
      y_bar[i] <- mean(y)
    }
    x_bar <- x_bar-y_bar
  }
  
  TS <- x_bar
  qqnorm(TS, main="Q-Q Plot:Sampling Distribution of Mean from Beta vs. Normal", xlab="Normal Quantiles", 
         ylab="Sample Quantiles of Mean")
}

#' Outputs QQ plot of one sample mean or difference of two sample mean for given sample size and Weibull distribution parameters.
#' 
#' @param n the sample size for sample one
#' @param shape the location parameter for sample one with default value of 0
#' @param scale the scale parameter for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param shape_2 the location parameter for sample two with default value of 0
#' @param scale_2 the scale parameter for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^4
#' @return QQ plot of one sample mean or difference of two sample mean vs Normal
#' @note the shape and scale parameters of Weibull distribution must be big than 0; The Weibull distribution with shape parameter a and scale parameter b has density given by f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a) for x>0
#' @examples
#' how_norm.weibull(n=5, shape=2, scale=1)
#' how_norm.weibull(n=5, shape=1,scale=2, n_2=7,shape_2=1, scale_2=1, sim=10^4)
how_norm.weibull<- function(n, shape=1, scale=1, n_2=0, shape_2=1, scale_2=1,sim=10^4){
  TS <- rep(NA, sim)
  x_bar <- rep(NA, sim)
  for(i in 1:sim){
    x <- rweibull(n, shape, scale)
    x_bar[i] <- mean(x)
    }
  
  if(n_2) {
    y_bar <- rep(NA, sim)
    for(i in 1:sim) {
      y <- rweibull(n_2, shape_2, scale_2)
      y_bar[i] <- mean(y)
    }
    x_bar <- x_bar-y_bar
  }
  
  TS <- x_bar
  qqnorm(TS, main="Sampling Distribution of Mean from Weibull vs. Normal", xlab="Normal Quantiles", 
         ylab="Sample Quantiles of Mean")
}

#' Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Exponential samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.
#' 
#' @param n the sample size for sample one with default value of 1
#' @param rate the rate parameter for sample one
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param rate_2 the rate parameter for sample two with default value of 1
#' @param sim the number of simulations with default value of 10^3
#' @param conf.level confidence level with default value of 0.95
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0
#' @return simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis
#' @note rate of exponential distribution must be big than 0
#' @examples
#' t_vs_wilcox.exp(n=30,rate=1,n_2=30,rate_2=1, sim=10^3, conf.level=0.95)
#' t_vs_wilcox.exp(n=30,rate=0.1,n_2=30,rate_2=1, sim=10^3)

t_vs_wilcox.exp <- function(n, rate, n_2=0, rate_2=1,mu=0, sim=10^3, conf.level=0.95){
  p.t <- rep(NA, sim)
  p.wilcox <- rep(NA, sim)
  for(i in 1: sim){
    x <- rexp(n, rate)
    if(n_2){
      y <- rexp(n_2, rate_2)
      } else y <- NULL
    T1 <- t.test(x,y, mu=mu)
    p.t[i] <- ifelse(T1$p.value<(1-conf.level), 1, 0)
    T2 <- wilcox.test(x,y,mu=mu)
    p.wilcox[i] <- ifelse( T2$p.value<(1-conf.level), 1, 0)
  }
  compare<- list(t.test=mean(p.t), wilcox.test=mean(p.wilcox))
  return(compare)
}

#' Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Lognormal samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.
#' 
#' @param n the sample size for sample one
#' @param meanlog the mean of the distribution on the log scale for sample one with default value of 0
#' @param sdlog the standard deviation of the distribution on the log scale for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param meanlog_2 the mean of the distribution on the log scale for sample two with default value of 0
#' @param sdlog_2 the standard deviation of the distribution on the log scale for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^3
#' @param conf.level confidence level with default value of 0.95
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0
#' @return simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis
#' @note sdlog parameter of lognormal distribution must be big than 0
#' @examples
#' t_vs_wilcox.lognorm(n=10, meanlog=0, sdlog=1, n_2=10, meanlog_2=0, sdlog_2=1,sim=10^3,mu=0, conf.level=0.95)
#' t_vs_wilcox.lognorm(n=30, meanlog=0, sdlog=1, n_2=30, meanlog_2=4, sdlog_2=1,sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.lognorm <- function(n, meanlog=0, sdlog=1, n_2=0, meanlog_2=0, sdlog_2=1,sim=10^3,mu=0, conf.level=0.95){
  p.t <- rep(NA, sim)
  p.wilcox <- rep(NA, sim)
  for(i in 1: sim){
    x <- rlnorm(n, meanlog, sdlog)
    if(n_2){
      y <- rlnorm(n_2, meanlog_2, sdlog_2)
      } else y <- NULL
    T1 <- t.test(x,y, mu=mu)
    p.t[i] <- ifelse(T1$p.value<(1-conf.level), 1, 0)
    T2 <- wilcox.test(x,y,mu=mu)
    p.wilcox[i] <- ifelse( T2$p.value<(1-conf.level), 1, 0)
  }
  compare<- list(t.test=mean(p.t), wilcox.test=mean(p.wilcox))
  return(compare)
}


#' Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Cauchy samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.
#' 
#' @param n the sample size for sample one
#' @param location the location parameter for sample one with default value of 0
#' @param scale the scale parameter for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param location_2 the location parameter for sample two with default value of 0
#' @param scale_2 the scale parameter for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^3
#' @param conf.level confidence level with default value of 0.95
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0
#' @return simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis
#' @note t-test might perform badly on Cauchy samples
#' @examples
#' t_vs_wilcox.cauchy(n=10, location=0, scale=1, n_2=10, location_2=0, scale_2=1,sim=10^3,mu=0, conf.level=0.95)
#' t_vs_wilcox.cauchy(n=30, location=0, scale=1, n_2=30, location_2=4, scale_2=1,sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.cauchy <- function(n, location=0, scale=1, n_2=0, location_2=0, scale_2=1,sim=10^3,mu=0, conf.level=0.95){
  p.t <- rep(NA, sim)
  p.wilcox <- rep(NA, sim)
  for(i in 1: sim){
    x <- rcauchy(n, location, scale)
    if(n_2){
      y <- rcauchy(n_2, location_2, scale_2)
      } else y <- NULL
    T1 <- t.test(x,y, mu=mu)
    p.t[i] <- ifelse(T1$p.value<(1-conf.level), 1, 0)
    T2 <- wilcox.test(x,y,mu=mu)
    p.wilcox[i] <- ifelse( T2$p.value<(1-conf.level), 1, 0)
  }
  compare<- list(t.test=mean(p.t), wilcox.test=mean(p.wilcox))
  return(compare)
}

#' Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Bernoulli samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.
#' 
#' @param n the sample size for sample one
#' @param prob the probability of success for sample one with default value of 1/2
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param prob_2 the probability of success for sample two with default value of 1/2
#' @param sim the number of simulations with default value of 10^3
#' @param conf.level confidence level with default value of 0.95
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0
#' @return simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis
#' @note 0< prob <1 It may produce warnings since Bernoulli sample has large portion of ties.
#' @examples
#' t_vs_wilcox.bernoulli(n=10, prob=1/2, n_2=10, prob_2=1/2, sim=10^3,mu=0, conf.level=0.95)
#' t_vs_wilcox.bernoulli(n=50, prob=0.25,  n_2=50, prob_2=0.75, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.bernoulli <- function(n, prob=1/2, n_2=0, prob_2=1/2,sim=10^3,mu=0, conf.level=0.95){
  p.t <- rep(NA, sim)
  p.wilcox <- rep(NA, sim)
  for(i in 1: sim){
    x <- rbinom(n, size=1, prob)
    if(n_2){
      y <- rbinom(n_2, size=1, prob_2)
      } else y <- NULL
    T1 <- t.test(x,y, mu=mu)
    p.t[i] <- ifelse(T1$p.value<(1-conf.level), 1, 0)
    T2 <- wilcox.test(x,y,mu=mu)
    p.wilcox[i] <- ifelse( T2$p.value<(1-conf.level), 1, 0)
  }
  compare<- list(t.test=mean(p.t), wilcox.test=mean(p.wilcox))
  return(compare)
}

#' Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Poisson samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.
#' 
#' @param n the sample size for sample one
#' @param lambda the mean for sample one with default value of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param lambda_2 the mean for sample two with default value of 1
#' @param sim the number of simulations with default value of 10^3
#' @param conf.level confidence level with default value of 0.95
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0
#' @return simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis
#' @note lambda>0 It may produce warnings since Poisson sample may have ties.
#' @examples
#' t_vs_wilcox.poisson(n=10, lambda=1, n_2=10, lambda_2=1, sim=10^3,mu=0, conf.level=0.95)
#' t_vs_wilcox.poisson(n=30, lambda=1,  n_2=30, lambda_2=0.2, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.poisson <- function(n, lambda=1, n_2=0, lambda_2=1,sim=10^3,mu=0, conf.level=0.95){
  p.t <- rep(NA, sim)
  p.wilcox <- rep(NA, sim)
  for(i in 1: sim){
    x <- rpois(n, lambda)
    if(n_2){
      y <- rpois(n_2, lambda_2)
      } else y <- NULL
    T1 <- t.test(x,y, mu=mu)
    p.t[i] <- ifelse(T1$p.value<(1-conf.level), 1, 0)
    T2 <- wilcox.test(x,y,mu=mu)
    p.wilcox[i] <- ifelse( T2$p.value<(1-conf.level), 1, 0)
  }
  compare<- list(t.test=mean(p.t), wilcox.test=mean(p.wilcox))
  return(compare)
}

#' Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Gamma samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.
#' 
#' @param n the sample size for sample one
#' @param shape the shape parameter for sample one with default value of 1
#' @param scale the scale parameter for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param shape_2 the shape parameter for sample two with default value of 1
#' @param scale_2 the scale parameter for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^3
#' @param conf.level confidence level with default value of 0.95
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0
#' @return simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis
#' @note shape>0, scale>0 The Gamma distribution with parameters shape = a and scale = s has density f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s); E(X) = a*s and Var(X) = a*s^2
#' @examples
#' t_vs_wilcox.gamma(n=10, shape=1, scale=1, n_2=10, shape_2=1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
#' t_vs_wilcox.gamma(n=30, shape=1, scale=1, n_2=30, shape_2=5, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.gamma <- function(n, shape=1, scale=1, n_2=0, shape_2=1, scale_2=1,sim=10^3,mu=0, conf.level=0.95){
  p.t <- rep(NA, sim)
  p.wilcox <- rep(NA, sim)
  for(i in 1: sim){
    x <- rgamma(n, shape, scale)
    if(n_2){
      y <- rgamma(n_2, shape_2, scale_2)
      } else y <- NULL
    T1 <- t.test(x,y, mu=mu)
    p.t[i] <- ifelse(T1$p.value<(1-conf.level), 1, 0)
    T2 <- wilcox.test(x,y,mu=mu)
    p.wilcox[i] <- ifelse( T2$p.value<(1-conf.level), 1, 0)
  }
  compare<- list(t.test=mean(p.t), wilcox.test=mean(p.wilcox))
  return(compare)
}

#' Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Beta samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.
#' 
#' @param n the sample size for sample one
#' @param shape1 the first shape parameter for sample one with default value of 1
#' @param shape2 the second shape parameter for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param shape1_2 the first shape parameter for sample two with default value of 1
#' @param shape2_2 the second shape parameter for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^3
#' @param conf.level confidence level with default value of 0.95
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0
#' @return simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis
#' @note the shape parameters of Beta distribution must be big than 0; 
#' @examples
#' t_vs_wilcox.beta(n=10, shape1=1, shape2=1, n_2=10, shape1_2=1, shape2_2=1, sim=10^3,mu=0, conf.level=0.95)
#' t_vs_wilcox.beta(n=30, shape1=1, shape2=1, n_2=30, shape1_2=3, shape2_2=1, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.beta <- function(n, shape1=1, shape2=1, n_2=0, shape1_2=1, shape2_2=1,sim=10^3,mu=0, conf.level=0.95){
  p.t <- rep(NA, sim)
  p.wilcox <- rep(NA, sim)
  for(i in 1: sim){
    x <- rbeta(n, shape1, shape2)
    if(n_2){
      y <- rbeta(n_2, shape1_2, shape2_2)
      } else y <- NULL
    T1 <- t.test(x,y, mu=mu)
    p.t[i] <- ifelse(T1$p.value<(1-conf.level), 1, 0)
    T2 <- wilcox.test(x,y,mu=mu)
    p.wilcox[i] <- ifelse( T2$p.value<(1-conf.level), 1, 0)
  }
  compare<- list(t.test=mean(p.t), wilcox.test=mean(p.wilcox))
  return(compare)
}

#' Compare the power or type I error of t-test or Wilcoxon rank sum test for simulated Weibull samples. If input null hypothesis parameters, it campares type I error. If input alternative parameters, it campares power.
#' 
#' @param n the sample size for sample one
#' @param shape the shape parameter for sample one with default value of 1
#' @param scale the scale parameter for sample one with default values of 1
#' @param n_2 the sample size for sample two with default value of 0 (for one sample test case)
#' @param shape_2 the shape parameter for sample two with default value of 1
#' @param scale_2 the scale parameter for sample two with default values of 1
#' @param sim the number of simulations with default value of 10^3
#' @param conf.level confidence level with default value of 0.95
#' @param mu a number indicating the true value of the mean (or difference in means if you are performing a two sample test) with default value of 0
#' @return simulated type I error of t-test and Wilcoxon rank sum test with parameters under null hypothesis; simulated power of t-test and Wilcoxon rank sum test with parameters under alternative hypothesis
#' @note the shape and scale parameters of Weibull distribution must be big than 0; The Weibull distribution with shape parameter a and scale parameter b has density given by f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a) for x>0
#' @examples
#' t_vs_wilcox.weibull(n=10, shape=1, scale=1, n_2=10, shape_2=1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
#' t_vs_wilcox.weibull(n=50, shape=10, scale=5, n_2=50, shape_2=0.1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
t_vs_wilcox.weibull <- function(n, shape=1, scale=1, n_2=0, shape_2=1, scale_2=1,sim=10^3,mu=0, conf.level=0.95){
  p.t <- rep(NA, sim)
  p.wilcox <- rep(NA, sim)
  for(i in 1: sim){
    x <- rweibull(n, shape, scale)
    if(n_2){
      y <- rweibull(n_2, shape_2, scale_2)
      } else y <- NULL
    T1 <- t.test(x,y, mu=mu)
    p.t[i] <- ifelse(T1$p.value<(1-conf.level), 1, 0)
    T2 <- wilcox.test(x,y,mu=mu)
    p.wilcox[i] <- ifelse( T2$p.value<(1-conf.level), 1, 0)
  }
  compare<- list(t.test=mean(p.t), wilcox.test=mean(p.wilcox))
  return(compare)
}

#' Transforms data to rank, then performs one or two sample t_tests on rank.

#' @param x a (non-empty) numeric vector of data values
#' @param y an optional (non-empty) numeric vector of data values
#' @param ... further arguments to be passed to t-test
#' @return A list with class "htest" containing the value of t-statisc, parameter, p-value, confidence interval, estimate, null hypothesis, alternative hypothesis, method, and data name.
#' @examples
#' rank_t.test(1:10, y=c(7:20))
rank_t.test <- function(x, y = NULL, ...){
  data <- c(x,y)
  rank.data <- rank(data, ties.method="min")
  index.x <- 1:length(x)
  rank.x <- rank.data[index.x]
  rank.y <- NULL
  if(length(y)!=0){
    index.y <- (length(x)+1) : length(data)
    rank.y <- rank.data[index.y]
  }
  t.test(rank.x, rank.y, ...)
}

