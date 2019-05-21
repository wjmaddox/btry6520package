# btry6520package
test R package for BTRY6520 HW assignment

to install use 
```{r}
devtools::install_github('wjmaddox/btry6520package/btry6520fns')
```

Then, check out the wonderful vignette in `btry6520fns/vignettes/description.Rmd`.

I've included examples for:
- algorithmic leveraging (and uniform)
- elastic net coordinate descent
- iterative linear systems solving with Jacobi (and parallel Jacobi) and Gauss-Seidel.

Unfortunately, my system seemed to freeze for the system solver for parallel Jacobi after about N > 5000. 
