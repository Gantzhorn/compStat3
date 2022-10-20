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


#Simulate using hint
N <- 10000 #number of simulations
W <- rchisq(N, df = nu)
X <- rnorm(N, mean = mu, sd = sqrt(nu*sigma_sq/W))

ggplot2::ggplot(data = tibble(x = X), aes(x = x)) +
  geom_histogram(fill = "Firebrick", col = "Black")

source("~/Desktop/Skole/Compstat/Assignment3/rScript/assignment3/EM_1.R")
source("~/Desktop/Skole/Compstat/Assignment3/rScript/assignment3/gradient_descent.R")
#profiling
profvis::profvis(source("~/Desktop/Skole/Compstat/Assignment3/rScript/assignment3/EM_1.R"))

#Benchmarking
a <-microbenchmark::microbenchmark(EM_1(data = X, par = initpar, nu = nu), times = 10000)

ggplot2::autoplot(a)