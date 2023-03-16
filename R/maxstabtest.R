
#' P-P plot for testing max stability
#'
#' The diagnostic, proposed by Gabda, Towe, Wadsworth and Tawn,
#' relies on the fact that, for max-stable vectors on the unit Gumbel scale,
#' the distribution of the maxima is Gumbel distribution with a location parameter equal to the exponent measure.
#' One can thus consider tuples of size \code{m} and estimate the location parameter via maximum likelihood
#' and transforming observations to the standard Gumbel scale. Replicates are then pooled and empirical quantiles are defined.
#' The number of combinations of \code{m} vectors can be prohibitively large, hence only \code{nmax} randomly selected
#' tuples are selected from all possible combinations. The confidence intervals are obtained by a
#' nonparametric bootstrap, by resampling observations with replacement observations for the selected tuples and re-estimating the
#' location parameter. The procedure can be computationally intensive as a result.
#' @references Gabda, D.; Towe, R. Wadsworth, J. and J. Tawn, Discussion of ``Statistical Modeling of Spatial Extremes'' by A. C. Davison, S. A. Padoan and M. Ribatet. \emph{Statist. Sci.} \bold{27} (2012), no. 2, 189--192.
#' @param dat matrix or array of max-stable observations, typically block maxima. The first dimension should consist of replicates
#' @param m integer indicating how many tuples should be aggregated.
#' @param nmax maximum number of pairs. Default to 500L.
#' @param B number of nonparametric bootstrap replications. Default to 1000L.
#' @param ties.method string indicating the method for \code{\link[base]{rank}}. Default to \code{"random"}.
#' @param plot logical indicating whether a graph should be produced (default to \code{TRUE}).
#' @return a Tukey probability-probability plot with 95% confidence intervals obtained using a nonparametric bootstrap
#' @export
#' @examples
#' \dontrun{
#' dat <- mev::rmev(n = 250, d = 100, param = 0.5, model = "log")
#' maxstabtest(dat, m = 100)
#' maxstabtest(dat, m = 2, nmax = 100)
#' dat <- mev::mvrnorm(n = 250, Sigma = diag(0.5, 10) + matrix(0.5, 10, 10), mu = rep(0, 10))
#' maxstabtest(dat, m = 2, nmax = 100)
#' maxstabtest(dat, m = ncol(dat))
#' }
maxstabtest <- function(dat, m = prod(dim(dat)[-1]), nmax = 500L, B = 1000L,
                        ties.method = "random", plot = TRUE){
  if(is.null(dim(dat))){
    stop("\"dat\" should be an array or a matrix, with replicates")
  }
  d <- length(dim(dat))
  if(!d %in% 2:3){
    stop("\"dat\" must be a n by D matrix or else an (n by D by D) array if observations are gridded")
  }
  D <- prod(dim(dat)[2:d])
  n <- dim(dat)[1]
  m <- as.integer(m)
  if(B < 20){
   stop("Not enough bootstrap replications to compute 95% confidence intervals.")
  }
  if(m < 1 && m > D){
    stop("Invalid number of tuples for the test; should be an integer")
  }


  #Transform observations to Gumbel margins
  #Use of apply -> observations are now first
  gdat <- matrix(apply(dat, 2:d,
                       function(x){-log(-log(rank(x, na.last = NA, ties.method = ties.method)/(n+1)# + runif(length(x), min = -0.5/n, max = 0.5/n)
                                        )) }), nrow = n, ncol = D)
  #discreteness of ranks means that for strongly dependent processes

  # Selection of m-tuples
  ntuples <- min(nmax, choose(D, m))
  if(ntuples * n > 10000L){
    warning("The number of statistics is higher than 10000L: consider reducing \"nmax\".")
  }
  #TODO Choose pairs that are not far apart - using some form of preferential sampling based on distance
  if(ntuples < nmax){
    indm <- as.matrix(combn(x = 1:D, m = m)) #replicates per column
  } else{
    # Select subsample of combn
    # https://stackoverflow.com/questions/18292144/select-a-subset-of-combinations
    # @author Alex
    .combn_sub <- function (x, m, nset = 5000, simplify = TRUE, ...) {
      stopifnot(length(m) == 1L)
      if (m < 0)
        stop("m < 0", domain = NA)
      if (is.numeric(x) && length(x) == 1L && x > 0 && trunc(x) ==
          x)
        x <- seq_len(x)
      n <- length(x)
      if (n < m)
        stop("n < m", domain = NA)
      m <- as.integer(m)
      e <- 0
      h <- m
      a <- seq_len(m)
      len.r <- length(r <-  x[a] )
      count <- as.integer(round(choose(n, m)))
      if( count < nset ) nset <- count
      dim.use <- c(m, nset)

      ##-----MOD 1: Change the output matrix size--------------
      out <- matrix(r, nrow = len.r, ncol = nset)

      if (m > 0) {
        i <- 2L
        nmmp1 <- n - m + 1L

        ##----MOD 2: Select a subset of indices
        #set.seed(seed)
        samp <- sort(c(1, sample( 2:count, nset - 1 )))

        ##----MOD 3: Start a counter.
        counter <- 2L

        while (a[1L] != nmmp1 ) {
          if (e < n - h) {
            h <- 1L
            e <- a[m]
            j <- 1L
          }
          else {
            e <- a[m - h]
            h <- h + 1L
            j <- 1L:h
          }
          a[m - h + j] <- e + j

          #-----MOD 4: Whenever the counter matches an index in samp,
          #a combination of row indices is produced and stored in the matrix `out`
          if(samp[i] == counter){
            out[, i] <- x[a]
            if( i == nset ) break
            i <- i + 1L
          }
          #-----Increase the counter by 1 for each iteration of the while-loop
          counter <- counter + 1L
        }
      }
      array(out, dim.use)
    }
    indm <- .combn_sub(x = 1:D, m = m, nset = nmax)
  }
  if(ntuples * n < 1000){
    p <- 1:(ntuples*n)/(ntuples*n+1)
  } else{
    p <- seq(0.01, 0.99, by=0.005)
  }
  bsamp <- matrix(0, nrow = B, ncol = length(p))
  csamp <- numeric(length = ntuples * n)
  #MLE for mu has closed form log(n/sum(exp(-z)))
  for(i in 1:ntuples){
    csamp[(n*(i-1)+1):(n*i)] <- apply(gdat[,indm[,i]], 1, max)
    csamp[(n*(i-1)+1):(n*i)] <- csamp[(n*(i-1)+1):(n*i)] -
      min(max(log(n)-log(sum(exp(-csamp[(n*(i-1)+1):(n*i)]))), 0), log(m))
  }
  if(ntuples * n < 1000){
    bsamp[1,] <- sort(mev::pgev(csamp, loc = 0, scale = 1, shape = 0))
  } else{
    bsamp[1,] <- quantile(mev::pgev(csamp, loc = 0, scale = 1, shape = 0), p, type = 3)
  }

  for(b in 2:B){
    # if(ntuples == 1 && m == D){ #case m = D, taking maxima of all spatial replications
    #   csamp[1:n] <- mev::rgev(n = n, loc = 0, shape = 0) #transform to ranks also
    #   csamp[1:n] <- csamp[1:n] - (log(n)-log(sum(exp(-csamp[1:n]))))
    # } else{
      #Nonparametric bootstrap procedure
    for(i in 1:ntuples){
      csamp[(n*(i-1)+1):(n*i)] <- apply(gdat[sample.int(n, n, replace = TRUE),indm[,i]], 1, max)
      csamp[(n*(i-1)+1):(n*i)] <- csamp[(n*(i-1)+1):(n*i)] -
        min(max(log(n)-log(sum(exp(-csamp[(n*(i-1)+1):(n*i)]))), 0), log(m))
    #}
    }
    if(ntuples * n < 1000){
      bsamp[b,] <- sort(mev::pgev(csamp, loc = 0, scale = 1, shape = 0))
    } else{
      bsamp[b,] <- quantile(mev::pgev(csamp, loc = 0, scale = 1, shape = 0), p, type = 3)
    }
  }
  #Problem with size of matrix allocation, since we have ntuples * B * n elements in the matrix.
  quants <- t(apply(bsamp, 2, quantile, c(0.025,0.975)))
  estimate <- cbind(percentiles = bsamp[1,], lower = quants[,1], upper = quants[,2])
  matplot(p, estimate - p, type = "l", lty = c(1,2,2),
          col = c(1, "grey80","grey80"), ylim = c(max(-1, 1.2*min(estimate - p)), min(1.2*max(estimate - p),1)),
          panel.first = {abline(h=0)}, bty = "l",  xaxs = "i", ylab = "Probability difference")
  return(invisible(list(p = p, estimate = estimate)))
}

