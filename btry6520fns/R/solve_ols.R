#library(foreach)

#' generates a tri-diagonal matrix of size n
#'
#' @param alpha diagonal of matrix
#' @param n size of matrix
#'
#' @examples
#' A = gen_tridiag_mat(4, 10)
#' A
#'
#' @export
gen_tridiag_mat = function(alpha = 1, n = 100){
  a = matrix(0, nrow = n, ncol = n)
  diag(a) = alpha
  for(i in 1:(n-1)){
    a[i,i+1] = -1
    a[i+1,i] = -1
  }
  return(a)
}

#' @import foreach
jacobi_iterative_parallel = function(A, b, nit = 10000){
  # start time at initialization
  ptm = proc.time()

  # convergence condition is \rho( -D^{-1}(L + U) ) < 1
  dinv = diag(1 / diag(A))
  nj = A
  diag(nj) = 0
  condition = max(eigen(dinv %*% nj, only.values = TRUE)$values)

  #elapsed_time = rep(0, nit)
  #relative_err = rep(0, nit)

  n = nrow(A)
  x = rep(0, n)

  for(k in 1:nit){
    xcurr = x
    x <- foreach(i=1:n, .combine=c, .export = NULL) %dopar%{
      sgn = rep(1, n)
      sgn[i] = 0
      ( b[i] - A[i,] %*% (sgn * xcurr) ) / A[i,i]
    }

  }
  elapsed_time = (proc.time() - ptm)[3]
  soln = x
  attr(soln, which = "elapsed_time") <- elapsed_time
  return(soln)
  #relative_err = relative_err / norm2(soln)
  # return(list(error = relative_err,
  #             time = elapsed_time,
  #             condition = condition))
}


jacobi_vector = function(A, b, soln, nit = 10000){
  # start time at initialization
  ptm = proc.time()

  #elapsed_time = rep(0, nit)
  #relative_err = rep(0, nit)

  n = nrow(A)
  x = rep(0, n)

  l = A
  u = A
  l[upper.tri(A, diag = TRUE)] = 0
  u[lower.tri(A, diag = TRUE)] = 0

  nj = -(l + u)

  dinv = diag(1 / diag(A))

  # convergence condition is \rho( -D^{-1}(L + U) ) < 1
  condition = max(eigen(dinv %*% nj, only.values = TRUE)$values)

  for(k in 1:nit){
    # from pg 511 of golub and van loan, matrix computations
    # x^{k} = D^{-1} (-(L + U) x^{k-1} + b)
    x = dinv %*% (nj %*% x + b)

    #elapsed_time[k] = (proc.time() - ptm)[3]
    #relative_err[k] = norm2(x - soln)
  }
  elapsed_time = (proc.time() - ptm)[3]
  soln = x
  attr(soln, which = "elapsed_time") <- elapsed_time
  return(soln)
}

gauss_seidel = function(A, b, nit = 10000){
  # start time at initialization
  ptm = proc.time()

  #elapsed_time = rep(0, nit)
  #relative_err = rep(0, nit)

  n = nrow(A)
  x = rep(0, n)

  l = A
  u = A
  l[upper.tri(A, diag = FALSE)] = 0
  u[lower.tri(A, diag = TRUE)] = 0

  ng = -u

  # l^{-1}
  mg = solve(l)

  # convergence condition is \rho(-L^{-1} U) < 1
  condition = max(Re(
    eigen(mg %*% ng, only.values = TRUE)$values)
  )
  #print(max(Re(eigen(mg %*% ng, only.values = TRUE)$values)))

  for(k in 1:nit){
    #from pg 511 of golub and van loan, matrix computations
    #x^k = L^{-1} (-u x^{k-1} + b)
    x = mg %*% (ng %*% x + b)
    #elapsed_time[k] = (proc.time() - ptm)[3]
    #relative_err[k] = norm2(x - soln)
  }

  elapsed_time = (proc.time() - ptm)[3]
  soln = x
  attr(soln, which = "elapsed_time") <- elapsed_time
  return(soln)
}

#' solves ordinary least squares problems (Ax = b) with eitehr Jacobi or Gauss-Seidel method.
#' Returns the soln with an attribute of time
#'
#' @param A lhs of system (must be n x n matrix type)
#' @param b rhs of system; should be n x 1 matrix/vector
#' @param nIt number of iterations
#' @param method either "Jacobi" or "GS" (Gauss-Seidel)
#' @param parallel for Jacobi: whether or not to construct a parallel backend with doParallel
#' @param nCores number of cores for parallel backend
#'
#' @examples
#' A = gen_tridiag_mat(4, n = 100) #generate a tridiagonal matrix
#' v = rep(c(1, 0), nrow(A) / 2)
#' b = A %*% v
#' solve_ols(A, b, nIt = 100) #should match v
#'
#' @export
#'
#' @importFrom foreach %dopar%
#' @importFrom parallel detectCores makeCluster stopCluster
#'
solve_ols = function(A, b, nIt = 10000, method = "Jacobi", parallel = FALSE, nCores = 4){
  if(method == "Jacobi" & parallel == FALSE){

    return(jacobi_vector(A, b, nIt))

  }else if(method == "GS"){
    if(parallel==TRUE){
      warning("Parallel unused for Gauss-Seidel")
    }
    return(gauss_seidel(A, b, nIt))

  }else if(method == "Jacobi" & parallel == TRUE){
    #setup parallel backend to use many processors
    cores = parallel::detectCores()
    cl <- parallel::makeCluster(nCores) #not to overload your computer
    doParallel::registerDoParallel(cl)

    soln = jacobi_iterative_parallel(A, b, nIt)

    parallel::stopCluster(cl)
    return(soln)
  }
}

