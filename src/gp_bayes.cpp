#include <Rcpp.h>
using namespace Rcpp;
// Zhang and Stephens (2009) and Zhang (2010) routines for function gp.fit

// Profile log-likelihood for GPD, parametrized in function of (theta, k)
NumericVector pll(NumericVector x, NumericVector theta){
  NumericVector q(2);
  q[1] = sum(log(1-theta[0]*x))/x.size(); //xi
  q[0] = x.size()*(log(-theta[0]/q[1])-q[1]-1); //theta
  return q;
}
// GPD prior of Zhang and Stephens
double prZS(NumericVector theta, double bound, double scale0, double xi0){
  return -log(scale0)-(1.0/xi0 + 1.0)*log(1 + xi0*(bound-theta[0])/scale0);
}

// Adaptive Metropolis-Hastings algorithm for posterior estimates of GPD profile likelihood model
// [[Rcpp::export]]
List Zhang_Stephens(NumericVector x, NumericVector init, NumericVector adapt_sd=0.1, bool adapt=true, int burnin=1000,
                             int niter = 10000, int thin=1, int method=1) { // int verbose = 200 not implemented
  // http://www.johndcook.com/blog/standard_deviation/
  //Sort data to retrieve first quartile and sample maximum
  std::sort(x.begin(), x.end());
  //Initialize containers
  // NumericMatrix chain = NumericMatrix(niter,2);
  // NumericVector loglikvals = NumericVector(niter);
  NumericMatrix summary = NumericMatrix(2,2); //Mean first row, var second row
  NumericVector cur = clone(init);
  NumericVector prop = clone(init);
  NumericVector accept = NumericVector(init.size());
  NumericVector attempt = NumericVector(init.size());
  NumericVector cur_ll(1);
  NumericVector prop_ll(1);
  NumericVector cur_kt(2); // To reduce calculations, save intermediate k values
  NumericVector prop_kt(2); // idem
  NumericVector temp(2); //temporary container for running var
  int counter = 0;
  double bound; double scale0; double xi0;
  //Quantities used in the empirical Bayes prior
  if(method==1){
    bound = 1.0/x(x.size()-1);
    //First quartile of the data Xstar
    scale0 = 1.0/6.0*x(floor((double)x.size()/4.0+0.5)-1);
    xi0 = 0.5;
  } else if(method==2){
    bound = (x.size()-1.0)/(x.size()+1.0)*1.0/x(x.size()-1);
    double p; double kp;
    NumericVector sigmap(7);
    for(int j=0; j<sigmap.size();j++){
    p = 0.3 + j*0.1;
    kp = log(x(floor((1-p*p)*x.size()+0.5)-1)/x(floor((1-p)*x.size()+0.5)-1)-1.0)/log(p);
    sigmap[j] = kp*x(floor((1-p)*x.size()+0.5)-1)/(1-pow(p,kp));
    }
    std::sort(sigmap.begin(), sigmap.end());
    scale0 = 1.0/(2*sigmap[3]);
    xi0 = 1.0;
  }

  //Initialization of the algorithm
  cur_kt = pll(x, init);
  cur_ll[0] = cur_kt[0] + prZS(init, bound, scale0, xi0);

  for(int i=0; i<(burnin+niter*thin); i++){
    //Metropolis-Hastings step
    //Random walk Normal proposals
      prop[0] = rnorm(1,cur[0], adapt_sd[0])[0];
      prop_kt = pll(x, prop);
      prop_ll[0] = prop_kt[0] + prZS(prop, bound, scale0, xi0);
    if(log(runif(1)[0])<(prop_ll[0]-cur_ll[0])){ // - dnorm(cur,0.0,adapt_sd[0],1)[0]+ dnorm(prop,0.0,adapt_sd[0],1)[0]
      //Increase acceptance counter and save current proposal
        accept[0]=accept[0]+1.0;
        cur[0] = prop[0];
        cur_ll[0] = prop_ll[0];
        cur_kt = prop_kt;
    }
    attempt = attempt + 1.0;
    // Reporting the results

     if((i + 1)%5000 == 0){
 						checkUserInterrupt();
//       Rcpp::Rcout << "Iteration " << i+1<< ": prop_ll "<< prop_ll[0] << " versus cur_ll " << cur_ll[0]  << std::endl;
//       Rcpp::Rcout << " Parameter values: " << cur[0] << std::endl;
    }

    // Adapting the proposals during burnin to get proposals toward 0.44 acceptance rate
    // Note: in general d-dimensional problems, 0.234 is optimal
       if(i<burnin*0.95){
          if(adapt==true){
          if(attempt[0]>50.0){
            if(accept[0]/attempt[0]>0.9){
              adapt_sd[0] = adapt_sd[0]*1.5;
              accept[0] = 0.0; attempt[0] = 0.0;
            } else if(accept[0]/attempt[0]<0.05){
              adapt_sd[0] = adapt_sd[0]*0.5;
              accept[0] = 0.0; attempt[0] = 0.0;
            } else if(accept[0]/attempt[0]>0.64){
              adapt_sd[0] = adapt_sd[0]*1.25;
              accept[0] = 0.0; attempt[0] = 0.0;
            } else if(accept[0]/attempt[0]<0.24){
              adapt_sd[0] = adapt_sd[0]*0.75;
              accept[0] = 0.0; attempt[0] = 0.0;
            } else if(accept[0]/attempt[0]>0.49){
              adapt_sd[0] = adapt_sd[0]*1.1;
              accept[0] = 0.0; attempt[0] = 0.0;
            } else if(accept[0]/attempt[0]<0.39){
              adapt_sd[0] = adapt_sd[0]*0.95;
              accept[0] = 0.0; attempt[0] = 0.0;
            } else if(accept[0]/attempt[0]>0.47){
            	adapt_sd[0] = adapt_sd[0]*1.05;
            	accept[0] = 0.0; attempt[0] = 0.0;
            }
          }
        }
       }

    // Reset the acceptance rate after burnin
    if(i==(burnin-1)){
    	// Increase burnin period if too short
    	if((accept[0]/attempt[0] > 0.48) || (accept[0]/attempt[0] < 0.38)){
    	 		burnin = burnin + 1000;
    		//Put an upper bound for adaptation
    		if(burnin>10000){
    			burnin = 10000;
    		}
    	} else{
        accept[0] = 0.0; attempt[0] = 0.0;
      }
    }
    if(i>=burnin && i%thin==0){
      // Saving the output
    if(counter==0){
      summary(0,0) = -cur_kt[1]/cur[0];  //mean of sigma
      summary(0,1) = cur_kt[1]; //mean of xi
      summary(1,0) = 0.0; //var of sigma
      summary(1,1) = 0.0; //var of xi
    } else{
    // Running mean and variances: Welford (1962)
    temp[0] = summary(0,0); temp[1] = summary(0,1);
    summary(0,0) = summary(0,0) - (cur_kt[1]/cur[0]+summary(0,0))/(counter+1.0);
    summary(0,1) = summary(0,1) + (cur_kt[1]-summary(0,1))/(counter+1.0);
    summary(1,0) = summary(1,0) + (-cur_kt[1]/cur[0]-temp[0])*(-cur_kt[1]/cur[0]-summary(0,0));
    summary(1,1) = summary(1,1) + (cur_kt[1]-temp[1])*(cur_kt[1]-summary(0,1));
    }

//       chain(counter,1) = cur_kt[1];//Shape parameter
//       chain(counter,0) = -cur_kt[1]/cur[0];//Scale parameter
//       loglikvals[counter] = cur_ll[0];

      counter++;
    }

  }
  //Rescale components to get variance
  summary(1,0) = summary(1,0)/(niter-1);
  summary(1,1) = summary(1,1)/(niter-1);
  //NumericVector mean
  return List::create(Named("rate")=accept/attempt, Named("adapted.sd")=adapt_sd,Named("summary")=summary,
                            Named("thin")=thin,Named("burnin")=burnin,Named("niter")=niter);
  // out.attr("class") = "gpdbayes";
  // out;
//   return List::create(Named("posterior")=chain, Named("rate")=accept/attempt, Named("adapted.sd")=adapt_sd, Named("thin")=thin,
//                       Named("param.init")=init, Named("ll.vals")=loglikvals);
 }




