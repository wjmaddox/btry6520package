HW2 Functions
================

Installatation
--------------

``` r
devtools::install_github('wjmaddox/btry6520package/btry6520fns')
```

Then load the package

``` r
library(btry6520fns)
```

Leveraging Examples
-------------------

Here's a plot of the two different types of algorithmic leveraging. Pretty clearly, this is random data so there should be no effects, but the uniform leveraging procedure seems to pick up a higher variance of effects than hat based leveraging.

``` r
X = matrix(rnorm(10000), ncol = 20)
y = rnorm(nrow(X))

# algorithmic leveraging
hat_lev = algo_leverage(y, X, size = 100, method = "hat")
#> [1] "running method"
#> [1] 100

# uniform leveraging
unif_lev = algo_leverage(y, X, size = 100, method = "uniform")
#> [1] 100

plot(hat_lev, unif_lev)
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

Elastic net coordinate descent
------------------------------

Here's an example of medium n, large p (it's not ImageNet but it's okay...)

``` r
data = gen_fs_system(n = 100, p = 500)
beta_pred = elnet_coord(data$y, data$x, lambda = 0.5, alpha = 0.5)

# we're able to pick up the true components in even a large p system
print(data$beta[1:7] - beta_pred[1:7])
#> [1]  0.6042922  0.0000000 -0.6115097  0.0000000  0.3793045  0.0000000
#> [7] -0.5853123

# number of nonzero elements
print(sum(beta_pred > 0))
#> [1] 22
```

solving OLS systems
-------------------

``` r
A = gen_tridiag_mat(4, n = 100) #generate a tridiagonal matrix
v = rep(c(1, 0), nrow(A) / 2)
b = A %*% v

jac_soln = solve_ols(A, b, nIt = 100, method = "Jacobi", parallel = FALSE)
jac_error = sum(jac_soln - v)
print(jac_error)
#> [1] 0

gs_soln = solve_ols(A, b, nIt = 100, method = "GS", parallel = FALSE)
gs_error = sum(gs_soln - v)
print(gs_error)
#> [1] 0

print("Times")
#> [1] "Times"
print(attr(jac_soln,"elapsed_time"))
#> elapsed 
#>   0.257
print(attr(gs_soln,"elapsed_time"))
#> elapsed 
#>   0.012
```

``` r
A = gen_tridiag_mat(4, n = 1000) #generate a tridiagonal matrix
v = rep(c(1, 0), nrow(A) / 2)
b = A %*% v

jac_soln = solve_ols(A, b, nIt = 100, method = "Jacobi", parallel = TRUE)
jac_error = sum(jac_soln - v)
print(jac_error)
#> [1] 0

print(attr(jac_soln,"elapsed_time"))
#> elapsed 
#>  47.113
```
