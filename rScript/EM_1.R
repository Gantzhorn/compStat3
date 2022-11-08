#E-step
Estep0 <- function(x, mu = 0, sig_sq = 1, nu = 5){
  Q1 <- (nu*sig_sq*(nu+1))/(nu*sig_sq+(x-mu)^2)
  Q2 <- -log(1/2*(1+(x-mu)^2/(nu*sig_sq))) + digamma((nu+1)/2)
  out <- list(Q1 = Q1, Q2 = Q2)
  out
}

Qfunc <- function(par, data, nu){
  N <- length(data)
  E <- Estep0(data, par[1], par[2], nu)
  -N/2*log(pi*nu*par[2]) + sum((nu-1)/2*E$Q2) + sum(-1/2*E$Q1*(1+(data-par[1])^2/(nu*par[2])))
}

Mstep0 <- function(x, E, nu = 5){
  mu_est <- sum(x*E)/sum(E)
  sigma_sq_est <- 1/(length(x)*nu)*sum(E*(x-mu_est)^2)
  out <- c(mu_est, sigma_sq_est)
  out
}

EM_1 <- function(data,
                 par,
                 nu,
                 epsilon = 1e-6,
                 maxiter = 50,
                 trace = NULL,
                 prints = FALSE){
  for(i in 1:maxiter){
    par0 <- par

    par <- Mstep0(data,
                  Estep0(
                    data,
                    par0[1],
                    par0[2],
                    nu)$Q1,
                  nu)
    if(!is.null(trace)) trace()
    if(sum((par-par0)^2) <= epsilon * (sum(par^2)+epsilon)) break
  }
  if(i == maxiter) {base::warning("Maximum number of iterations reached")}
  if(prints == TRUE && i < maxiter){
    print(glue::glue("Convergence reached after ", i, " steps."))
  }
  par
}

est <- EM_1(data = X, par = initpar, nu = nu, prints = TRUE)

grad_Qfunc <- function(par,
                       data,
                       nu){
  N <- length(data)
  E <- Estep0(data, par[1], par[2], nu)
  wrt_mu <- sum(E$Q1*(data-par[1])/(nu* par[2]))
  wrt_sigma_sq <- -N/(2* par[2])+1/2*sum(E$Q1*(data-par[1])^2/(par[2]^2*nu))
  out <- c(wrt_mu,wrt_sigma_sq)
  out
}

#Calculate Fisher-information

Phi <- function(par0){
  Mstep0(X,Estep0(X,par0[1],par0[2],nu = nu)$Q1,nu = nu)
}

DPhi <- numDeriv::jacobian(Phi, est)

fisher1 <- function(grad, estimates, Data, ...){
  iY <- -numDeriv::jacobian(grad, estimates, data = Data, ...)
  iX <- (diag(1,2) - t(DPhi)) %*% iY
  solve(iX)
}

fisher2 <- function(grad, estimates, Data, ...){
  iY <- -numDeriv::jacobian(grad, estimates, data = Data, ...)
  iYinv <- solve(iY)
  iYinv + iYinv %*% t(solve(diag(1, 2) - DPhi, DPhi))
}

iXinv1 <- fisher1(grad_Qfunc, est, X, nu = nu)

iXinv2 <- fisher2(grad_Qfunc, est, X, nu = nu)


confint_est <- matrix(c(est[1]+c(-1,1)*1.96*sqrt(iXinv1[1,1]), est[2] + c(-1,1)*1.96*sqrt(iXinv1[2,2])),
                      nrow = 2, byrow = T)

# #Benchmarking
# fish_speed <- microbenchmark::microbenchmark(fisher1(grad_Qfunc, est10, X10, nu = nu),
#                                fisher2(grad_Qfunc, est10, X10, nu = nu),
#                                fisher1(grad_Qfunc, est100, X100, nu = nu),
#                                fisher2(grad_Qfunc, est100, X100, nu = nu),
#                                fisher1(grad_Qfunc, est1000, X1000, nu = nu),
#                                fisher2(grad_Qfunc, est1000, X1000, nu = nu),
#                                fisher1(grad_Qfunc, est10000, X10000, nu = nu),
#                                fisher2(grad_Qfunc, est10000, X10000, nu = nu),
#                                fisher1(grad_Qfunc, est100000, X100000, nu = nu),
#                                fisher2(grad_Qfunc, est100000, X100000, nu = nu),
#                                fisher1(grad_Qfunc, est1000000, X1000000, nu = nu),
#                                fisher2(grad_Qfunc, est1000000, X1000000, nu = nu))
# 
# fisher_tib <- summary(fish_speed) %>% as_tibble()
# 
# fisher_tib %>% mutate(N = c(10, 10 , 100, 100, 1000, 1000, 10000, 10000, 100000, 100000, 1000000, 1000000),
#                       Implementation = ifelse(str_detect(expr, "fisher1"), "Fisher1", "Fisher2")) %>% 
#   ggplot(aes(x = log(N), y = log(median), col = Implementation)) + geom_line(size = 1, alpha = 1)
