// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]
# include <RcppArmadillo.h>
using namespace Rcpp;
//using namespace arma;

// [[Rcpp::export(.EuclideanWeights)]]
arma::vec EuclideanWeights(arma::mat x, arma::rowvec mu){
  if(mu.n_cols!=x.n_cols){
    Rcpp::stop("Invalid argument for Euclidean likelihood estimator: nonconformal mean and sample");
  }
  arma::rowvec xbar = sum(x,0) / (double)x.n_rows;
  arma::mat xcent = x-repmat(xbar, x.n_rows,1);
  arma::mat S = cov(x,1); // Equivalent to xcent.t() * xcent/(double)x.n_rows;
  arma::vec weights = (1 - xcent * S.i() * trans(xbar-mu)) / (double)x.n_rows;
  return(weights);
}


// Degree 4 polynomial approx to log
arma::vec h(arma::vec y, arma::vec cvals){
  arma::vec ans(y.size());
  ans.fill(cvals(0));
  for(int j=1; j<5; j++){
    ans = ans + pow(y,j)*cvals(j);
  }
  return ans;
}

// Derivative of hp, from approximation at pt
arma::vec hp(arma::vec y, double pt){
  arma::vec ans(y.size(), arma::fill::zeros);
  for(int j=0; j<4; j++){
    ans = ans + pow(-y/pt, j);
  }
  ans = ans / (-pt);
  return ans;
}


// Second derivative of h at y, from approx at pt
arma::vec hpp(arma::vec y, double pt){
  arma::vec ans(y.size(), arma::fill::zeros);
  for(int j=0; j<3; j++){
    ans = ans + (j+1.0)*pow(-y/pt, j);
  }
  ans = ans / (pt * pt);
  return ans;
}

// //' 	Self-concordant Taylor approximation logstar function
// //'
// //' 	Minus log and its first \code{der} derivatives, on \code{eps < x < M}, with
// //' 	fourth order Taylor series approximation to the left of \code{eps} and to the right of \code{M}
// //'
// //' 	@param x column vector of observations
// //' 	@param eps lower tolerance
// //' 	@param M maximum value
// //' 	@param der derivative, either 0, 1 or 2.
// //'
// //' 	@return a matrix with \eqn{der+1} columns containing values
arma::mat mllog(arma::colvec x, double eps, double M, int der=0){
  if(!(der==0 || der==1 || der==2)){
    Rcpp::stop("Only first two derivatives and log functions implemented");
  }
  if(eps > M){
    Rcpp::stop("Thresholds out of order");
  }
  //Coefficients for 4th order Taylor approx below eps
  arma::vec coefs(5);
  arma::vec Coefs(5);
  coefs(0) = -log(eps);
  Coefs(0) = -log(M);
  for(int i=1; i<5; i++){
    coefs(i) = R_pow(-eps,-i)/(double)i;
    Coefs(i) = R_pow(-M,-i)/(double)i;
  }

  arma::uvec lo = find(x < eps);
  arma::uvec hi = find(x > M);

  //Function value
  arma::vec f(x.size());
  f = -log(x);
  f(lo) = h(x(lo)-eps, coefs);
  f(hi) = h(x(hi)-M,   Coefs);

  if(der < 1){
    return f;
  }
  //First derivative

  arma::vec fp(x.size());
  fp = -1.0/x;
  fp(lo) = hp(x(lo)-eps, eps);
  fp(hi) = hp(x(hi)-M,   M);

  if(der < 2){
    return join_rows(f, fp);
  }
  //Second derivative
  arma::vec fpp(x.size());
  fpp = pow(x,-2);
  fpp(lo) = hpp(x(lo)-eps, eps);
  fpp(hi) = hpp(x(hi)-M,   M);

  return join_rows(join_rows(f,fp),fpp);
}

// Singular value decomposition least squares solution for linear regression
arma::vec svdlm (arma::mat X, arma::colvec y){
  // Linear model regression coefficient via SVD

  // Tolerances for generalized inverse via SVD
  double RELTOL = 1e-9;
  double ABSTOL = 1e-100;
  //SVD decomposition and handling of failures
  arma::vec d;
  arma::mat U;
  arma::mat V;
  arma::svd(U, d, V, X);
  U.resize(y.n_elem,d.n_elem);
  double maxd = max(d);
  for(int i=0; i<d.n_elem; i++){
    if(d(i) < (RELTOL * max(d) + ABSTOL)){
      d(i) = 0;
    } else{
      d(i)=1 / d(i);
    }
  }
  return V * diagmat(d) * U.t()  * y;
}
//[[Rcpp::export(.emplik_intern)]]
List emplik_intern(arma::mat z, arma::colvec mu, arma::vec lam, double eps,
	double M = 1e30, double thresh = 1e-12, int itermax = 1000){
//(arma::mat z, arma::vec mu  = vec(z.n_cols,fill::zeros), double eps = 1/z.nrows, double M = datum::inf);
// # Backtracking line search parameters [Tweak only with extreme caution.]
// # See Boyd and Vandenberghe, pp 464-466.
double ALPHA = 0.3; // seems better than 0.01 on some 2d test data (sometimes fewer iters)
double BETA  = 0.8;
// # We need  0 < ALPHA < 1/2   and 0 < BETA < 1
// # Backtrack threshold: you can miss by this much.
double BACKEPS = 0;
// # Consider replacing 0 by 1e-10 if backtracking seems to be
// # failing due to round off.
int n = z.n_rows;
int d = z.n_cols;
 z.each_row() -= trans(mu);
//  Use lam = 0 or initial lam, whichever is best
arma::vec onen = arma::vec(n,arma::fill::ones);
arma::mat init0 = mllog(onen, eps, M, 2);
arma::mat init(init0.n_rows, init0.n_cols);
if(!any(lam)){
  init = init0;
} else {
  init = mllog(onen+z*lam, eps, M, 2);
	if(sum(init0.col(0) < sum(init.col(0)))){
    lam = arma::rowvec(z.n_cols,arma::fill::zeros);
    init = init0;
  }
}
// # Initial f, g
double fold = sum(init.col(0));
arma::rowvec gold = sum(z.each_col() % init.col(1),0);
bool converged = false;
int iter = 0;
arma::mat oldvals = init;
arma::mat newvals(oldvals.n_rows, oldvals.n_cols);
arma::vec rootllpp;
arma::mat zt(z.n_rows, z.n_cols);
arma::vec yt(z.n_rows);
arma::vec wts;
double fnew; double targ;
double logelr;
arma::vec step(z.n_cols);
int s; double ndec; double gradnorm;
bool backtrack = false;
while(!converged){
  iter += 1;
  //   # Get Newton Step
  rootllpp = sqrt(oldvals.col(2));  //# sqrt 2nd deriv of -llog lik
  zt = z;
     for(int j=0; j<d; j++){
      zt.col(j) %= rootllpp;
    }
    yt   = oldvals.col(1) / rootllpp;
    step = -svdlm(zt,yt);
    backtrack = false;
    s = 1;
    while(!backtrack){
         newvals = mllog(onen + z * (lam+s*step), eps, M, 2);
        fnew = sum(newvals.col(0));
        targ = fold + ALPHA * s * sum(trans(gold) % step) + BACKEPS; // (BACKEPS for roundoff, should not be needed)
         if(fnew <= targ){ // backtracking has converged
          backtrack = true;
          oldvals = newvals;
          fold = fnew;
     			gold = sum(z.each_col() % oldvals.col(1),0);
          lam = lam + s*step;
         } else{
          s = s * BETA;
         }
       }
    //   # Newton decrement and gradient norm
      ndec     = sqrt(sum(square(step % trans(gold))));
      gradnorm = sqrt(sum(square(gold)));
      converged = (ndec * ndec <= thresh);
       if( iter > itermax )
       	break;
     }
     wts = pow(1 + z * lam, -1)/n;
		logelr = sum(mllog(onen + z * lam, eps, M, 0).col(0));
    return Rcpp::List::create(Named("logelr")=logelr, Named("lam") = lam, Named("wts") = wts,
        Named("conv") = converged, Named("niter") = iter,Named("ndec") = ndec, Named("gradnorm") = gradnorm);
     }

// [[Rcpp::export(.Pickands_emp)]]
NumericVector Pickands_emp(NumericVector s, NumericVector ang, NumericVector wts){
  if(wts.size()!=ang.size()){
    warning("Only implemented in the bivariate case");
    stop("Non-conformal arguments; size of angles does not match weights.");
  }
  NumericVector pick(s.size());
  for(int i=0;i<s.size();i++){
    pick[i] = 2*sum(pmax((1-s[i])*ang,s[i]*(1-ang))*wts);
  }
  return pick;
}

// Log of Dirichlet function
// [[Rcpp::export(ldirfn)]]
double ldirfn(NumericVector param){
  double res = lgamma(sum(param))-sum(lgamma(param));
  return res;
}



// [[Rcpp::export(.loocvdens)]]
NumericVector loocvdens(double nu, NumericMatrix ang, NumericVector wts, NumericMatrix loowts) {
  NumericVector result(1);
  int n = loowts.ncol();
  NumericVector sumbeta(1);
  for(int i = 0; i < n; i++){
    sumbeta[0] = 0;
    for(int j = 0; j < n; j++){
      if(i == j){
        continue;
      }
      sumbeta[0] = sumbeta[0] + exp(log(loowts(j, i)) + ldirfn( nu * ang(j,_) ) + sum((nu * ang(j,_) - 1.0) * log(ang(i,_))));
    }
    result[0] = result[0] - log(sumbeta[0]);
  }
  return result;
}
