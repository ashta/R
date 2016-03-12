# Local-EM estimator

# density_estimator - density estimator
# datapts - locations of the mismeasured data
# tgtpts - locations at which to evaluate the smoothed fit
# bw - bandwidth
density_estimator <- function(datapts, tgtpts, bw) {
  sapply( tgtpts, function(x) mean(dnorm((datapts - x) / bw))/bw )
}
# lEM_estimator - local-EM estimator for Gaussian mismeasured error
# datapts - locations of the mismeasured data
# niter - number of iterations
# bw - bandwidth
# mu - known mean parameter of the mismeasurement error
# sigma - known standard deviation parameter of the mismeasurement error
# tgtpts - locations at which to evaluate the smoothed fit
lEM_estimator <- function(datapts, niter = 1, bw, mu, sigma, tgtpts = NULL) {
  
  if(is.null(tgtpts)) {
    errpts <- 0.25 * (max(datapts) - min(datapts))
    tgtpts <- seq(min(datapts) - errpts, max(datapts) + errpts, length = max(100, length(datapts)))
  }
  # 1st iteration of estimator, mathematically the same as density estimator with different arguments
  tgtlem <- density_estimator(datapts, tgtpts + mu, sqrt(bw^2 + sigma^2))
  result <- tgtlem
  
  # nth iteration of estimator
  if(niter > 1) {
    require(statmod)
    # locations to evaluate integral
    sdpts <- max(sd(datapts),sd(tgtpts))
    limitpts <- c(min(datapts,tgtpts) - 25 * sdpts, max(datapts,tgtpts) + 25 * sdpts)
    # Gauss Quadrature approximations
    evalpts <- 0.5 * (diff(limitpts) * gauss.quad(500)$nodes + sum(limitpts))
    evalwgts <- gauss.quad(500)$weights 
    # 1st iteration of locations of integral
    newlem <- density_estimator(datapts, evalpts + mu, sqrt(bw^2 + sigma^2))
    
    # exp_estimator - estimator for local-EM expectation
    # pts - locations at which to evaluate the integral
    # wgts - weights associated with locations
    # prev - previous ith local-EM estimator for evaluating locations
    # datapt - data location
    # tgtpt - target location
    # bw - bandwidth
    # mu - known mean parameter of the mismeasurement error
    # sigma - known standard deviation parameter of the mismeasurement error
    exp_estimator <- function(datapt, tgtpt,pts, wgts, prev,  bw, mu, sigma) {
      numer <- sum(dnorm((pts - tgtpt)/bw)/bw * dnorm(datapt - pts, mu, sigma) * prev * wgts)
      denom <- sum(dnorm(datapt - evalpts, mu, sigma) * prev * wgts)
      return(numer/denom)
    }
    iter <- 1
    while(iter < niter) {
      result <- colMeans(outer(datapts, tgtpts, Vectorize(function(a, b) 
        exp_estimator(a, b, evalpts, evalwgts, newlem, bw, mu, sigma))))
      # update estimator for locations of integral
      iter <- iter + 1
      if(iter < niter) {
        prevlem <- newlem
        newlem <- colMeans(outer(datapts, evalpts, Vectorize(function(a, b)
          exp_estimator(a,b, evalpts, evalwgts, prevlem, bw, mu, sigma))))  
      }
    }
  }
  return(result)
}