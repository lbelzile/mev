## -----------------------------------------------------------------------------
library(mev)
#Sample of size 1000 from a 5-dimensional logistic model
x <- rmev(n=1000, d=5, param=0.5, model="log")
#Marginal parameters are all standard Frechet, meaning GEV(1,1,1)
apply(x, 2, function(col){ismev::gev.fit(col, show=FALSE)$mle})


#Sample from the corresponding spectral density
w <- rmevspec(n=1000, d=5, param=0.5, model="log")
#All rows sum to 1 by construction
head(rowSums(w))
#The marginal mean is 1/d
round(colMeans(w),2)

