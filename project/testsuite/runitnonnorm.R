#Test function test.how_norm.-- test whether sampling distribution converge to Normal when large N
#test.how_norm.-- test for both one sample mean and two sample mean difference. plot for one sample mean; plot_2 for two sample mean difference
#The sampling distribution of mean from Cauchy does not converge to Normal even large N

#The test function test.t_vs_wilcox.-- test for reasonal type I error and power with moderate sample size or approprate effect size.
#t-test performs badly on Cauchy or Weibull distribution. test for power waivered for these two distributions
#Bernoulli and Poisson distribution provide discrete sample, and may produce warning due to ties during test.
library(RUnit)

test.how_norm.bernoulli <- function(){
  plot <- how_norm.bernoulli(n=1000,prob=0.5, sim=10^4)
  checkEquals(cor(unlist(plot[1]), unlist(plot[2])),1, tolerance=0.001)
  plot_2 <- how_norm.bernoulli(n=200,prob=0.5, n_2=200,prob_2=0.5, sim=10^4)
  checkEquals(cor(unlist(plot_2[1]), unlist(plot_2[2])),1, tolerance=0.001)
}

test.how_norm.beta <- function(){
  plot <- how_norm.beta(n=1000, shape1=1, shape2=1, sim=10^4)
  checkEquals(cor(unlist(plot[1]), unlist(plot[2])),1, tolerance=0.001)
  plot_2 <- how_norm.beta(n=200, shape1=1, shape2=1, n_2=200, shape1_2=1,shape2_2=1,sim=10^4)
  checkEquals(cor(unlist(plot_2[1]), unlist(plot_2[2])),1, tolerance=0.001)
}

test.how_norm.cauchy <- function(){
  plot <- how_norm.cauchy(n=1000,location=0, scale=1, sim=10^4)
  checkTrue(cor(unlist(plot[1]), unlist(plot[2])) < 0.5)
  plot_2 <- how_norm.cauchy(n=200,location=0, scale=1, n_2=200, location_2=0, scale_2=1,sim=10^4)
  checkTrue(cor(unlist(plot_2[1]), unlist(plot_2[2])) < 0.5)

}

test.how_norm.exp <- function(){
  plot <- how_norm.exp(n=1000, rate=1,sim=10^4)
  checkEquals(cor(unlist(plot[1]), unlist(plot[2])),1, tolerance=0.001)
  plot_2 <- how_norm.exp(n=200, rate=1, n_2=200, sim=10^4)
  checkEquals(cor(unlist(plot_2[1]), unlist(plot_2[2])),1, tolerance=0.001)
}

test.how_norm.gamma <- function(){
  plot <- how_norm.gamma(n=1000, shape=1, scale=3, sim=10^4)
  checkEquals(cor(unlist(plot[1]), unlist(plot[2])),1, tolerance=0.001)
  plot_2 <- how_norm.gamma(n=200, shape=1, scale=3, n_2=200, shape_2=1,scale_2=3,sim=10^4)
  checkEquals(cor(unlist(plot_2[1]), unlist(plot_2[2])),1, tolerance=0.001)
}

test.how_norm.lognorm <- function(){
  plot <- how_norm.lognorm(n=1000, meanlog=0, sdlog=1, sim=10^4)
  checkEquals(cor(unlist(plot[1]), unlist(plot[2])),1, tolerance=0.005)
  plot_2 <- how_norm.lognorm(n=200, meanlog=0, sdlog=1,n_2=200,meanlog_2=0, sdlog_2=1,sim=10^4)
  checkEquals(cor(unlist(plot_2[1]), unlist(plot_2[2])),1, tolerance=0.001)
}

test.how_norm.poisson <- function(){
  plot <- how_norm.poisson(n=1000,lambda=1, sim=10^4)
  checkEquals(cor(unlist(plot[1]), unlist(plot[2])),1, tolerance=0.001)
  plot_2 <- how_norm.poisson(n=200,lambda=1, n_2=200, lambda_2=1,sim=10^4)
  checkEquals(cor(unlist(plot_2[1]), unlist(plot_2[2])),1, tolerance=0.001)
}

test.how_norm.weibull <- function(){
  plot <- how_norm.weibull(n=1000, shape=1, scale=1, sim=10^4)
  checkEquals(cor(unlist(plot[1]), unlist(plot[2])),1, tolerance=0.001)
  plot_2 <- how_norm.weibull(n=200, shape=1, scale=1, n_2=200, shape_2=1,scale_2=1,sim=10^4)
  checkEquals(cor(unlist(plot_2[1]), unlist(plot_2[2])),1, tolerance=0.001)
}

test.t_vs_wilcox.exp <- function(){
  alpha <- t_vs_wilcox.exp(n=30,rate=1,n_2=30,rate_2=1, sim=1000, conf.level=0.95)
  checkTrue(alpha$t.test<0.1)
  checkTrue(alpha$wilcox.test<0.1)
  power <- t_vs_wilcox.exp(n=30,rate=0.1,n_2=30,rate_2=1, sim=1000)
  checkTrue(power$t.test>0.8) 
  checkTrue(power$wilcox.test>0.8)
}

test.t_vs_wilcox.lognorm <- function(){
  alpha <- t_vs_wilcox.lognorm(n=10, meanlog=0, sdlog=1, n_2=10, meanlog_2=0, sdlog_2=1,sim=10^3,mu=0, conf.level=0.95)
  checkTrue(alpha$t.test<0.1)
  checkTrue(alpha$wilcox.test<0.1)
  power<- t_vs_wilcox.lognorm(n=30, meanlog=0, sdlog=1, n_2=30, meanlog_2=4, sdlog_2=1,sim=10^3,mu=0, conf.level=0.95)
  checkTrue(power$t.test>0.8) 
  checkTrue(power$wilcox.test>0.8)
}

test.t_vs_wilcox.cauchy <- function(){
  alpha <- t_vs_wilcox.cauchy(n=10, location=0, scale=1, n_2=10, location_2=0, scale_2=1,sim=10^3,mu=0, conf.level=0.95)
  checkTrue(alpha$t.test<0.1)
  checkTrue(alpha$wilcox.test<0.1)
  power<- t_vs_wilcox.cauchy(n=30, location=0, scale=1, n_2=30, location_2=4, scale_2=1,sim=10^3,mu=0, conf.level=0.95)
  checkTrue(power$wilcox.test>0.8)
}

test.t_vs_wilcox.bernoulli <- function(){
  alpha <- t_vs_wilcox.bernoulli(n=10, prob=1/2, n_2=10, prob_2=1/2, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(alpha$t.test<0.1)
  checkTrue(alpha$wilcox.test<0.1)
  power<- t_vs_wilcox.bernoulli(n=50, prob=0.25,  n_2=50, prob_2=0.75, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(power$t.test>0.8)
  checkTrue(power$wilcox.test>0.8)
}

test.t_vs_wilcox.poisson <- function(){
  alpha <- t_vs_wilcox.poisson(n=10, lambda=1, n_2=10, lambda_2=1, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(alpha$t.test<0.1)
  checkTrue(alpha$wilcox.test<0.1)
  power<- t_vs_wilcox.poisson(n=30, lambda=1,  n_2=30, lambda_2=0.2, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(power$t.test>0.8)
  checkTrue(power$wilcox.test>0.8)
}

test.t_vs_wilcox.gamma <- function(){
  alpha <- t_vs_wilcox.gamma(n=10, shape=1, scale=1, n_2=10, shape_2=1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(alpha$t.test<0.1)
  checkTrue(alpha$wilcox.test<0.1)
  power<- t_vs_wilcox.gamma(n=30, shape=1, scale=1, n_2=30, shape_2=5, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(power$t.test>0.8)
  checkTrue(power$wilcox.test>0.8)
}

test.t_vs_wilcox.beta <- function(){
  alpha <- t_vs_wilcox.beta(n=10, shape1=1, shape2=1, n_2=10, shape1_2=1, shape2_2=1, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(alpha$t.test<0.1)
  checkTrue(alpha$wilcox.test<0.1)
  power<- t_vs_wilcox.beta(n=30, shape1=1, shape2=1, n_2=30, shape1_2=3, shape2_2=1, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(power$t.test>0.8)
  checkTrue(power$wilcox.test>0.8)
}

test.t_vs_wilcox.weibull<- function(){
  alpha <- t_vs_wilcox.weibull(n=10, shape=1, scale=1, n_2=10, shape_2=1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(alpha$t.test<0.1)
  checkTrue(alpha$wilcox.test<0.1)
  power<- t_vs_wilcox.weibull(n=50, shape=10, scale=5, n_2=50, shape_2=0.1, scale_2=1, sim=10^3,mu=0, conf.level=0.95)
  checkTrue(power$wilcox.test>0.8)
}