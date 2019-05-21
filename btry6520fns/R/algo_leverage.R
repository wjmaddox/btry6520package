wl_ols = function(y, X, probability_weights, size = 10, return_indices = FALSE){
  print(size)
  N = nrow(X)

  # draw random sample
  indices_to_keep = sample(1:N, size = size, replace = TRUE, prob = probability_weights)

  # compute model on kept indices with weights that are given
  model = lm(y[indices_to_keep] ~ X[indices_to_keep,] - 1,
             weights = probability_weights[indices_to_keep])

  # return either indices or coefficients
  if(return_indices){
    return(list(coef = coef(model), indices = indices_to_keep))
  }else{
    return(coef(model))
  }

}


unif_leveraging = function(y, X, size = 10, return_indices = FALSE){
  #weights = 1/N
  weights = rep(1/nrow(X), nrow(X))

  return(
    wl_ols(y, X, weights, size = size, return_indices = return_indices)
  )
}

hat_leveraging = function(y, X, size = 10, return_indices = FALSE){
  #weights \propto sum(U_{i.}^2)
  #eq. 14 of ma and sun
  u = svd(X)$u
  weights = apply(svd(X)$u, 1, function(x) sum(x^2))
  return(
    wl_ols(y, X, weights, size = size, return_indices = return_indices)
  )
}

#' Leveraging techniques with two different methods (uniform and hat).
#' Returns the coefficients from a linear model call.
#'
#' @param y response (vector or matrix)
#' @param X predictors (matrix)
#' @param size size option for uniform sub-sampling
#' @param method either uniform (random subsampling) or hat (SVD-based) leveraging.
#' See Algorithm 1 of Ma and Sun, "Leveraging for big data regression," WIREs Comp. Stat. 2014.
#'
#' @examples
#' X = matrix(rnorm(100), ncol = 2)
#' y = rnorm(50)
#' algo_leverage(y, X, method = "hat")
#'
#' @export
#' @import stats
algo_leverage = function(y, X, size = 10, method = "hat"){
  if(method == "hat"){
    print("running method")
    return(
      hat_leveraging(y, X, size = size)
    )
  }else if(method == "uniform"){
    return(
      unif_leveraging(y, X, size = size)
    )
  }else{
    stop("Only hat or uniform leveraging implemented")
  }
}
