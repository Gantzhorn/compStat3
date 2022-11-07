gradient_descent <- function(parinit,
                             #grad1,
                             data,
                             nu,
                             learning_rate = 0.001,
                             maxiter = 50,
                             epsilon = 1e-06){
  par <- parinit
  for (i in 1:50){
    print(par)
    par0 <- par
    par <- par0 - learning_rate*grad_Qfunc(par0, data, nu)
    print(par)
    if(sum(grad_Qfunc(par0, data, nu)) <= epsilon) break
  }
  if(i == maxiter)(base::warning("Maximum number of iterations reached"))
  else{print(glue::glue("Convergence reached after ", i, " steps."))}
  par
}

gradient_descent(parinit = initpar, data = raw_X, nu = nu)
