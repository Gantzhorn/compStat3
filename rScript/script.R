library(tidyverse)
theme_set(theme_bw())
library(CSwR)
library(microbenchmark)
library(Rcpp)
library(profvis)
library(numDeriv)
#---------------------------------Simulation-----------------------------------#
#Set parameters
mu <- 2
sigma_sq <- 3
nu <- 1

initpar <- c(2, 4)

#Simulate using hint
N <- 100000 #number of simulations
W <- rchisq(N, df = nu)
X <- rnorm(N, mean = mu, sd = sqrt(nu*sigma_sq/W))

#Simulate raw X (non-standard t-distribution)
# raw_X <- mu + rt(N, df = nu)*sqrt(sigma_sq)

 ggplot2::ggplot(data = tibble(x = X), aes(x = x)) +
   geom_histogram(fill = "Firebrick", col = "Black")

#Basic EM
source("~/Desktop/Skole/Compstat/Assignment3/compStat3/rScript/EM_1.R")

# Test of basic EM algorithm
EM_tracer <- CSwR::tracer(c("par", "par0"))
EM_1(data = X, par = initpar, nu = nu, trace = EM_tracer$tracer, maxiter = 50)
summary(EM_tracer) %>% as_tibble() %>% mutate(mu = abs(par.1-par0.1), sigma_sq = abs(par.2- par0.2), time = .time) %>% 
  select(time, mu, sigma_sq) %>% pivot_longer(-time, names_to = "Parameter", values_to = "stepsize") %>% 
  ggplot(aes(x = time, y = log(stepsize), col = Parameter)) + geom_point(size = 2) + xlab("Time") + 
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12))


EM_tracer <- CSwR::tracer(c("par", "par0"))
EM_1(data = X, par = initpar2, nu = nu, trace = EM_tracer$tracer)
summary(EM_tracer) %>% as_tibble() %>% mutate(mu = abs(par.1-par0.1), sigma_sq = abs(par.2- par0.2), time = .time) %>% 
  select(time, mu, sigma_sq) %>% pivot_longer(-time, names_to = "Parameter", values_to = "stepsize") %>% 
  ggplot(aes(x = time, y = log(stepsize), col = Parameter)) + geom_point(size = 2) + xlab("Time") + 
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12))


#Alternative implementation - Does not work yet
source("~/Desktop/Skole/Compstat/Assignment3/compStat3/rScript/gradient_descent.R")
par_grad <- gradient_descent(initpar, X, nu = nu)
#Rcpp implementation
Rcpp::sourceCpp("EM_CPP.cpp")
est_cpp <- EM_CPP(X,initpar2, nu, 1e-5, maxiter = 50)
#profiling
profvis::profvis(source("~/Desktop/Skole/Compstat/Assignment3/compStat3/rScript/EM_1.R"))

#Benchmarking
init_bench <-microbenchmark::microbenchmark(EM_1(data = X, par = initpar, nu = nu),
                                   EM_CPP(X,initpar, nu, 1e-5, maxiter = 50), times = 100)

as_tibble(init_bench) %>% mutate(Implementation = ifelse(str_detect(expr, "EM_1"), "R", "Rcpp"), time = time/1000000000) %>% 
  filter(time< 4) %>% ggplot(aes(x = time, fill = Implementation)) + geom_density(alpha = .7)

# Benchmarking as a function of N
# Simulate datasets
N1 <- 100
W1 <- rchisq(N1, df = nu)
X1 <- rnorm(N1, mean = mu, sd = sqrt(nu*sigma_sq/W1))

N2 <- 1000
W2 <- rchisq(N2, df = nu)
X2 <- rnorm(N2, mean = mu, sd = sqrt(nu*sigma_sq/W2))

N3 <- 5000
W3 <- rchisq(N3, df = nu)
X3 <- rnorm(N3, mean = mu, sd = sqrt(nu*sigma_sq/W3))

N4 <- 10000
W4 <- rchisq(N4, df = nu)
X4 <- rnorm(N4, mean = mu, sd = sqrt(nu*sigma_sq/W4))

N5 <- 25000
W5 <- rchisq(N5, df = nu)
X5 <- rnorm(N5, mean = mu, sd = sqrt(nu*sigma_sq/W5))

N6 <- 50000
W6 <- rchisq(N6, df = nu)
X6 <- rnorm(N6, mean = mu, sd = sqrt(nu*sigma_sq/W6))


more_bench <- microbenchmark::microbenchmark(EM_1(data = X1, par = initpar, nu = nu),
                                             EM_CPP(X1,initpar, nu, 1e-5, maxiter = 50),
                                             EM_1(data = X2, par = initpar, nu = nu),
                                             EM_CPP(X2,initpar, nu, 1e-5, maxiter = 50),
                                             EM_1(data = X3, par = initpar, nu = nu),
                                             EM_CPP(X3,initpar, nu, 1e-5, maxiter = 50),
                                             EM_1(data = X4, par = initpar, nu = nu),
                                             EM_CPP(X4,initpar, nu, 1e-5, maxiter = 50),
                                             EM_1(data = X5, par = initpar, nu = nu),
                                             EM_CPP(X5,initpar, nu, 1e-5, maxiter = 50),
                                             EM_1(data = X6, par = initpar, nu = nu),
                                             EM_CPP(X6,initpar, nu, 1e-5, maxiter = 50))

more_bench_sum <- summary(more_bench)

p1 <- more_bench_sum %>% mutate(N = c(100, 100, 1000, 1000, 5000, 5000, 10000, 10000, 25000, 25000, 50000, 50000),
                          Implementation = ifelse(str_detect(expr, "CPP"), "Rcpp", "R")) %>% 
  select(lq, median, uq, N, Implementation) %>% 
  ggplot(aes(x = N, y = median, col = Implementation)) + geom_line(size = 1)  + 
  ggtitle("Well-behaved") +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        plot.title = element_text(face = "bold", size = 16)) + ylab("Median")

well <- more_bench_sum %>% mutate(N = c(100, 100, 1000, 1000, 5000, 5000, 10000, 10000, 25000, 25000, 50000, 50000),
                                  Implementation = ifelse(str_detect(expr, "CPP"), "Rcpp", "R"),
                                  Type = "Well") %>% 
  select(lq, median, uq, N, Implementation, Type)


# Benchmarking as a function of with bad true values
# Simulate datasets
nu1 <- 0.3

N1 <- 100
W1 <- rchisq(N1, df = nu1)
X1 <- rnorm(N1, mean = mu, sd = sqrt(nu1*sigma_sq/W1))

N2 <- 1000
W2 <- rchisq(N2, df = nu1)
X2 <- rnorm(N2, mean = mu, sd = sqrt(nu1*sigma_sq/W2))

N3 <- 5000
W3 <- rchisq(N3, df = nu1)
X3 <- rnorm(N3, mean = mu, sd = sqrt(nu1*sigma_sq/W3))

N4 <- 10000
W4 <- rchisq(N4, df = nu1)
X4 <- rnorm(N4, mean = mu, sd = sqrt(nu1*sigma_sq/W4))

N5 <- 25000
W5 <- rchisq(N5, df = nu1)
X5 <- rnorm(N5, mean = mu, sd = sqrt(nu1*sigma_sq/W5))

N6 <- 50000
W6 <- rchisq(N6, df = nu1)
X6 <- rnorm(N6, mean = mu, sd = sqrt(nu1*sigma_sq/W6))


more_bench_ill <- microbenchmark::microbenchmark(EM_1(data = X1, par = initpar, nu = nu1, maxiter = 1000),
                                             EM_CPP(X1,initpar, nu1, 1e-5, maxiter = 1000),
                                             EM_1(data = X2, par = initpar, nu = nu1, maxiter = 1000),
                                             EM_CPP(X2,initpar, nu1, 1e-5, maxiter = 1000),
                                             EM_1(data = X3, par = initpar, nu = nu1, maxiter = 1000),
                                             EM_CPP(X3,initpar, nu1, 1e-5, maxiter = 1000),
                                             EM_1(data = X4, par = initpar, nu = nu1, maxiter = 1000),
                                             EM_CPP(X4,initpar, nu1, 1e-5, maxiter = 1000),
                                             EM_1(data = X5, par = initpar, nu = nu1, maxiter = 1000),
                                             EM_CPP(X5,initpar, nu1, 1e-5, maxiter = 1000),
                                             EM_1(data = X6, par = initpar, nu = nu1, maxiter = 1000),
                                             EM_CPP(X6,initpar, nu1, 1e-5, maxiter = 1000))

more_bench_ill_sum <- summary(more_bench_ill)

p2 <- more_bench_ill_sum %>% mutate(N = c(100, 100, 1000, 1000, 5000, 5000, 10000, 10000, 25000, 25000, 50000, 50000),
                          Implementation = ifelse(str_detect(expr, "CPP"), "Rcpp", "R")) %>% 
  select(lq, median, uq, N, Implementation) %>% 
  ggplot(aes(x = N, y = median, col = Implementation)) + geom_line(size = 1) + 
  ggtitle("Bad true values") +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        plot.title = element_text(face = "bold", size = 16)) + ylab("Median") 

ill <- more_bench_ill_sum %>% mutate(N = c(100, 100, 1000, 1000, 5000, 5000, 10000, 10000, 25000, 25000, 50000, 50000),
                                     Implementation = ifelse(str_detect(expr, "CPP"), "Rcpp", "R"),
                                     Type = "badTrueVal") %>% 
  select(lq, median, uq, N, Implementation, Type)

# Benchmarking as a function of with bad true values
# Simulate datasets

initpar2 <- c(-1000, 400)
N1 <- 100
W1 <- rchisq(N1, df = nu)
X1 <- rnorm(N1, mean = mu, sd = sqrt(nu*sigma_sq/W1))

N2 <- 1000
W2 <- rchisq(N2, df = nu)
X2 <- rnorm(N2, mean = mu, sd = sqrt(nu*sigma_sq/W2))

N3 <- 5000
W3 <- rchisq(N3, df = nu)
X3 <- rnorm(N3, mean = mu, sd = sqrt(nu*sigma_sq/W3))

N4 <- 10000
W4 <- rchisq(N4, df = nu)
X4 <- rnorm(N4, mean = mu, sd = sqrt(nu*sigma_sq/W4))

N5 <- 25000
W5 <- rchisq(N5, df = nu)
X5 <- rnorm(N5, mean = mu, sd = sqrt(nu*sigma_sq/W5))

N6 <- 50000
W6 <- rchisq(N6, df = nu)
X6 <- rnorm(N6, mean = mu, sd = sqrt(nu*sigma_sq/W6))


more_bench_badstart <- microbenchmark::microbenchmark(EM_1(data = X1, par = initpar2, nu = nu, maxiter = 50),
                                                 EM_CPP(X1,initpar2, nu, 1e-5, maxiter = 50),
                                                 EM_1(data = X2, par = initpar2, nu = nu, maxiter = 50),
                                                 EM_CPP(X2,initpar2, nu, 1e-5, maxiter = 50),
                                                 EM_1(data = X3, par = initpar2, nu = nu, maxiter = 50),
                                                 EM_CPP(X3,initpar2, nu, 1e-5, maxiter = 50),
                                                 EM_1(data = X4, par = initpar2, nu = nu, maxiter = 50),
                                                 EM_CPP(X4,initpar2, nu, 1e-5, maxiter = 50),
                                                 EM_1(data = X5, par = initpar2, nu = nu, maxiter = 50),
                                                 EM_CPP(X5,initpar2, nu, 1e-5, maxiter = 50),
                                                 EM_1(data = X6, par = initpar2, nu = nu, maxiter = 50),
                                                 EM_CPP(X6,initpar2, nu, 1e-5, maxiter = 50))

more_bench_badstart_sum <- summary(more_bench_badstart)

p3 <- more_bench_badstart_sum %>% mutate(N = c(100, 100, 1000, 1000, 5000, 5000, 10000, 10000, 25000, 25000, 50000, 50000),
                              Implementation = ifelse(str_detect(expr, "CPP"), "Rcpp", "R")) %>% 
  select(lq, median, uq, N, Implementation)
  ggplot(aes(x = N, y = median, col = Implementation)) + geom_line(size = 1)  +
  ggtitle("Poor starting values") +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12),
        plot.title = element_text(face = "bold", size = 16)) + ylab("Median")

  
badstart <- more_bench_badstart_sum %>% mutate(N = c(100, 100, 1000, 1000, 5000, 5000, 10000, 10000, 25000, 25000, 50000, 50000),
                                    Implementation = ifelse(str_detect(expr, "CPP"), "Rcpp", "R"),
                                    Type = "badStart") %>% 
  select(lq, median, uq, N, Implementation, Type)

bind_rows(well, ill, badstart) %>% mutate(Type = factor(Type, levels = c("Well", "badTrueVal", "badStart"))) %>%
  ggplot(aes(x = log(N), y = log(median), col = Implementation)) + geom_line(size = 1.2) +
  facet_grid(~Type) + theme(axis.title = element_text(face = "bold", size = 14),
                            axis.text = element_text(face = "bold", size = 12),
                            legend.title = element_text(face = "bold", size = 16),
                            legend.text = element_text(face = "bold", size =12),
                            strip.text.x = element_text(face = "bold", size = 12))
