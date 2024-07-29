#include <RcppArmadillo.h>
#include <R.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace R;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
mat Stein(const mat& coords, 
          const double& alpha, const double& nu, const double& eps, 
          const vec& Z){
  
  int n = coords.n_rows;
  mat out(n, n, arma::fill::zeros);
  double cnst = M_PI*pow(alpha, 2);
  double sigsq = cnst*pow(1e-8, nu)*R::bessel_k(1e-8, nu, 1)/pow(2, nu-1)/R::gammafn(nu+1);
  
  for (int i=0; i < n; i++) {
    double x1 = coords(i,0);
    double y1 = coords(i,1);
    double t1 = coords(i,2);
    
    for (int j=i+1; j < n; j++) {
      double x2 = coords(j,0);
      double y2 = coords(j,1);
      double t2 = coords(j,2);
      
      double xstar = x1-x2;
      double ystar = y1-y2;
      double tstar = t1-t2;
      double abs_tstar = abs(tstar);
      
      double r = alpha*sqrt(pow(xstar - eps*tstar*Z(0), 2) + 
                            pow(ystar - eps*tstar*Z(1), 2));
      
      if (r > 0) {
        out(i, j) = cnst*pow(r, nu+abs_tstar)*R::bessel_k(r, nu+abs_tstar, 1)/pow(2, nu+abs_tstar-1)/R::gammafn(nu+abs_tstar+1);
      } else {
        out(i, j) = sigsq;
      }
    }
  }
  out += out.t();
  out.diag().fill(sigsq);
  
  return out;
}



