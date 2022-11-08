# Can reuse E step

loss <- function(theta, data){sum(
    c(sum(data*Estep0(data, theta[1], theta[2], theta[3])$Q1)/sum(Estep0(data, theta[1], theta[2], theta[3])$Q1) - theta[1],
      
    1/theta[3]*mean(Estep0(data, theta[1], theta[2], theta[3])$Q1*(data-theta[1])^2) - theta[2],
    
    1/(2*theta[3]) + 1/2*log(2) + digamma(theta[3]/2) -
      1/2*mean(Estep0(data, theta[1], theta[2], theta[3])$Q2) -
      1/2*mean(Estep0(data, theta[1], theta[2], theta[3])$Q1*((data-theta[1])^2)/(theta[3]^2*theta[2]))
    )^2)
}
loss(parinit3, X)

numDeriv::grad(loss, x = parinit3, data =X)

nlm(loss, c(5.020901, 4.064585, 2.968535), data = X)

parinit3 <- c(2,3,1)

EM_mu <- function(data,
                 par,
                 epsilon = 1e-6,
                 maxiter = 50,
                 trace = NULL,
                 prints = FALSE){
  for(i in 1:maxiter){
    par0 <- par
    par <- nlm(loss, par, data = data)$estimate
    print(i)
    if(!is.null(trace)) trace()
    if(sum((par-par0)^2) <= epsilon * (sum(par^2)+epsilon)) break
  }
  if(i == maxiter) {base::warning("Maximum number of iterations reached")}
  if(prints == TRUE && i < maxiter){
    print(glue::glue("Convergence reached after ", i, " steps."))
  }
  par
}

EM_mu(X, parinit3)

