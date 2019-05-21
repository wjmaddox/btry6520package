soft_threshold = function(theta, lambda = 1) sign(theta) * pmax(0, abs(theta) - lambda)

#' generates the simulated data for a sparse system
#'
#' @param n number of data points
#' @param p dimensionality of system
#' @param noise noise in system
#' @param seed enables setting seed for reproducibility
#'
#' @examples
#' output = gen_fs_system()
#' output$y
#'
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
gen_fs_system = function(n = 20, p = 20, noise = 1, seed = 2019){
  set.seed(seed)

  epsilon = rnorm(n, mean = 0, sd = noise)

  design = diag(p)
  design[1,2] = 0.8
  design[5,6] = 0.8

  x = mvrnorm(n, rep(0, p), design)

  beta = c(2, 0, -2, 0, 1, 0, -1, rep(0, p - 7))

  y = x%*% beta + epsilon
  return(list(
    y = y,
    x = x,
    beta = beta
  ))
}

#' implementation of elastic net coordinate descent as described in Statistical Learning
#' with Sparsity, Hastie, Tibshirani, and Wainwright.
#' eg. hat{beta} = argmin_{beta} ||y - x beta||_2^2 + lambda(alpha ||beta||_1 + (1-alpha) ||beta||_2)
#'
#' returns the parameters estimated
#'
#' @param y response vector
#' @param x input feature matrix
#' @param alpha sparsity parameter for elastic net (default: 0, so lasso)
#' @param maxSteps maximum steps for iterative coordinate descent (default: 100)
#' @param minTol minimum convergence tolerance (default: 1e-8)
#' @param lambda shrinkage parameter
#'
#' @export
#'
#' @examples
#' data = gen_fs_system()
#' elnet_coord(data$y, data$x, lambda = 0.5, alpha = 0.5)
elnet_coord = function(y, x, lambda = 1, alpha = 0, maxSteps = 100, minTol = 1e-8){
  N = length(y)
  p = ncol(x)
  beta = rep(0, p)

  partial_residual = function(j, beta){
    adj_beta = beta
    adj_beta[j] = 0

    # compute residual of beta[-j]
    y - x %*% adj_beta
  }

  x_vars = apply(x, 2, function(i) mean(i^2))

  nsteps = 0
  beta_old = rep(1e5, p)
  tol = 1e5

  while(nsteps < maxSteps & tol > minTol){
    tol = sum((beta_old - beta)^2)
    #print(paste("Step: ", nsteps, "Change: ", tol))

    beta_old = beta

    # iterate through all
    for(j in 1:p){
      r_j = partial_residual(j, beta)
      # 2/N r^j x[,j]
      mean_residual = 2 * (t(r_j) %*% x[,j]) / N

      #\beta[j] = S_{\alpha \lambda}(residual) / (2 \sum x^2 + \lambda (1 - \alpha))
      beta[j] = soft_threshold(mean_residual, lambda = alpha * lambda) / (2 * x_vars[j] + lambda * (1 - alpha))
    }

    nsteps = nsteps + 1
  }
  if(nsteps > maxSteps - 2){print(paste("Exited at step: ", nsteps, "Tolerance: ", tol))
  }

  return(beta)
}
