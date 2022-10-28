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
    # if(!is.null(trace)) CSwR::tracer()
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

#Jacobian
iY <- -numDeriv::jacobian(grad_Qfunc, est, data = X, nu = nu)

-optimHess(est, Qfunc, grad_Qfunc, data = X, nu = nu)

iX <- (diag(1,2) - t(DPhi)) %*% iY

iXinv1 <- solve(iX) #One way

iYinv <- solve(iY)
iXinv2  <-  iYinv + iYinv %*% t(solve(diag(1, 2) - DPhi, DPhi))

confint_est <- matrix(c(est[1]+c(-1,1)*1.96*sqrt(iXinv1[1,1]), est[2] + c(-1,1)*1.96*sqrt(iXinv1[2,2])),
                      nrow = 2, byrow = T)