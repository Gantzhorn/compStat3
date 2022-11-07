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
sigma_sq <- 4
nu <- 5

initpar <- c(1, 5)

#Simulate using hint
N <- 10000 #number of simulations
W <- rchisq(N, df = nu)
X <- rnorm(N, mean = mu, sd = sqrt(nu*sigma_sq/W))

# ggplot2::ggplot(data = tibble(x = X), aes(x = x)) +
#   geom_histogram(fill = "Firebrick", col = "Black")

#Basic EM
source("~/Desktop/Skole/Compstat/Assignment3/compStat3/rScript/EM_1.R")

# Test of basic EM algorithm
EM_tracer <- CSwR::tracer(c("par", "par0"))
EM_1(data = X, par = initpar, nu = nu, trace = EM_tracer$tracer)
summary(EM_tracer) %>% as_tibble() %>% mutate(mu = abs(par.1-par0.1), sigma_sq = abs(par.2- par0.2), time = .time) %>% 
  select(time, mu, sigma_sq) %>% pivot_longer(-time, names_to = "Parameter", values_to = "stepsize") %>% 
  ggplot(aes(x = time, y = log(stepsize), col = Parameter)) + geom_point(size = 2) + xlab("Time") + 
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12))

initpar2 <- c(-100000, 400)
EM_tracer <- CSwR::tracer(c("par", "par0"))
EM_1(data = X, par = initpar2, nu = nu, trace = EM_tracer$tracer)
summary(EM_tracer) %>% as_tibble() %>% mutate(mu = abs(par.1-par0.1), sigma_sq = abs(par.2- par0.2), time = .time) %>% 
  select(time, mu, sigma_sq) %>% pivot_longer(-time, names_to = "Parameter", values_to = "stepsize") %>% 
  ggplot(aes(x = time, y = log(stepsize), col = Parameter)) + geom_point(size = 2) + xlab("Time") + 
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size =12))


#Alternative implementation
source("~/Desktop/Skole/Compstat/Assignment3/compStat3/rScript/gradient_descent.R")
#Rcpp implementation
Rcpp::sourceCpp("EM_CPP.cpp")
est_cpp <- EM_CPP(X,initpar, nu, 1e-5, maxiter = 50)
#profiling
profvis::profvis(source("~/Desktop/Skole/Compstat/Assignment3/compStat3/rScript/EM_1.R"))

#Benchmarking
init_bench <-microbenchmark::microbenchmark(EM_1(data = X, par = initpar, nu = nu),
                                   EM_CPP(X,initpar, nu, 1e-5, maxiter = 50), times = 100)

as_tibble(init_bench) %>% mutate(Implementation = ifelse(str_detect(expr, "EM_1"), "R", "Rcpp"), time = time/1000000) %>% 
  filter(time< 4) %>% ggplot(aes(x = time, fill = Implementation)) + geom_density(alpha = .7)
