% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_ols.R
\name{solve_ols}
\alias{solve_ols}
\title{solves ordinary least squares problems (Ax = b) with eitehr Jacobi or Gauss-Seidel method.
Returns the soln with an attribute of time}
\usage{
solve_ols(A, b, nIt = 10000, method = "Jacobi", parallel = FALSE,
  nCores = 4)
}
\arguments{
\item{A}{lhs of system (must be n x n matrix type)}

\item{b}{rhs of system; should be n x 1 matrix/vector}

\item{nIt}{number of iterations}

\item{method}{either "Jacobi" or "GS" (Gauss-Seidel)}

\item{parallel}{for Jacobi: whether or not to construct a parallel backend with doParallel}

\item{nCores}{number of cores for parallel backend}
}
\description{
solves ordinary least squares problems (Ax = b) with eitehr Jacobi or Gauss-Seidel method.
Returns the soln with an attribute of time
}
\examples{
A = gen_tridiag_mat(4, n = 100) #generate a tridiagonal matrix
v = rep(c(1, 0), nrow(A) / 2)
b = A \%*\% v
solve_ols(A, b, nIt = 100) #should match v

}
