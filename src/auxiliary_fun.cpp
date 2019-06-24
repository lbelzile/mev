#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]
using namespace Rcpp;

//' Distance matrix with geometric anisotropy
//'
//' The function computes the distance between locations, with geometric anisotropy.
//' The parametrization assumes there is a scale parameter, so that \code{scale}
//' is the distortion for the second component only. The angle \code{rho} must lie in
//' \eqn{[-\pi/2, \pi/2]}.
//'
//' @param loc a \code{d} by 2 matrix of locations giving the coordinates of a site per row.
//' @param scale numeric vector of length 1, greater than 1.
//' @param rho angle for the anisotropy, must be larger than \eqn{\pi/2} in modulus.
//' @return a \code{d} by \code{d} square matrix of pairwise distance
//' @export
// [[Rcpp::export]]
arma::mat distg(arma::mat loc, NumericVector scale, NumericVector rho){
  int d = loc.n_rows;
  int m = loc.n_cols;
  arma::mat aniso(2, 2);
  if(loc.n_cols > 2){
    stop("Invalid location matrix; only geometric anisotropy for bivariate is supported");
  }
  if(rho.size() > 1 || scale.size() > 1){
    stop("Invalid length for `scale` or `rho`");
  }
  if(std::abs(rho[0]) > M_PI){
    stop("Invalid `rho` argument: angle must be in [-pi/2, pi/2]");
  }
  if(scale[0] < 1){
    stop("Scale parameter should be larger than 1 for identifiability");
  }
    aniso(0,0) = cos(rho)[0];
      aniso(0,1) = sin(rho)[0];
      aniso(1,0) = -scale[0]*aniso(0,1);
      aniso(1,1) = scale[0]*aniso(0,0);
      arma::mat distmat(d, d);
      distmat.zeros();
      arma::vec dist_ij(loc.n_cols);
      for(int i=0; i<d-1; i++){
        for(int j=i+1; j < d; j++){
          for(int k=0; k < m; k++){
            dist_ij(k) = loc(i,k) - loc(j,k);
          }
          distmat(i,j) = arma::norm(aniso * dist_ij);
          distmat(j,i) = distmat(i,j);
        }
      }
      return distmat;
}



