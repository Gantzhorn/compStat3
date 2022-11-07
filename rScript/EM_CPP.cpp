#include <Rcpp.h>
// [[Rcpp::export]]
Rcpp::NumericVector eStepcpp(Rcpp::NumericVector x, double mu, double sig_sq, double nu) {
  Rcpp::NumericVector Q = (nu*sig_sq*(nu+1))/(nu*sig_sq+pow(x-mu,2));
  return Q;
}

Rcpp::NumericVector mStepcpp(Rcpp::NumericVector x, Rcpp::NumericVector E, double nu) {
  int N = x.size();
  double mu_est = sum(E*x)/(sum(E));
  double sigma_sq_est = 1/(N*nu)*sum(E*pow(x-mu_est,2));
  Rcpp::NumericVector out(2);
  out[0] = mu_est; out[1] = sigma_sq_est;
  return out;
}
// [[Rcpp::export]]
Rcpp::NumericVector EM_CPP(Rcpp::NumericVector x,
                           Rcpp::NumericVector par,
                           double nu,
                           double epsilon,
                           int maxiter
                           ){
  for(int i = 0; i<maxiter+1; ++i){
    Rcpp::NumericVector par0 = par;
    par = mStepcpp(x,
                   eStepcpp(x,
                      par0[0],
                      par0[1],
                          nu),
                    nu
                   );
    if(sum(pow(par-par0,2)) <= epsilon * (sum(pow(par,2))+epsilon)) break;
  }
  return par;
}


