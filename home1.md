Homework1
================

Q1
--

Create an R-function that takes a data matrix *X* ∈ ℝ<sup>*n* × *p*</sup>, *n* &gt; *p* as an argument and returns the unbiased estimator of the covariance matrix. Do not use the built-in functions cov(), cor() or any additional R-packages.

Code of Q1
----------

``` r
D1 <- as.matrix(read.table("data/Data1.txt", sep = "\t", header = F))

covariance <- function(x){
  ave <- apply(x, 2, mean)
  cent <- sweep(x, 2, ave, "-")
  cova <- t(cent) %*% cent / (dim(x)[1]-1)
  return(cova)
}

D1.cov <- covariance(D1)
cov(D1)
```

    ##              V1          V2           V3           V4
    ## V1  0.085600783 0.007107367 -0.008906472 -0.004708558
    ## V2  0.007107367 0.069159452  0.015104514  0.013638556
    ## V3 -0.008906472 0.015104514  0.104474247 -0.003473715
    ## V4 -0.004708558 0.013638556 -0.003473715  0.080600878

``` r
D1.cov
```

    ##              V1          V2           V3           V4
    ## V1  0.085600783 0.007107367 -0.008906472 -0.004708558
    ## V2  0.007107367 0.069159452  0.015104514  0.013638556
    ## V3 -0.008906472 0.015104514  0.104474247 -0.003473715
    ## V4 -0.004708558 0.013638556 -0.003473715  0.080600878

Q2
--

Create an R-function that takes a full-rank covariance matrix *A* ∈ ℝ<sup>*p* × *p*</sup>as an argument and returns the square root of the inverse matrix such that $A^{-\\frac{1}{2}}A^{-\\frac{1}{2}} = A^{-1}$

``` r
trans_square_root <- function(x){
  lam <- eigen(x)$values
  gam <- eigen(x)$vectors
  
  eigen_value <- sqrt(solve(diag(lam)))
  
  return(gam %*% eigen_value %*% t(gam))
}

trans_square_root(D1.cov)
```

    ##            [,1]       [,2]       [,3]       [,4]
    ## [1,]  3.4524853 -0.2186977  0.1846031  0.1347682
    ## [2,] -0.2186977  3.9304827 -0.3442340 -0.3735435
    ## [3,]  0.1846031 -0.3442340  3.1484650  0.1164939
    ## [4,]  0.1347682 -0.3735435  0.1164939  3.5810484

Q3
--

Create an R-function that takes a full-rank covariance matrix A as an argument and returns the corresponding correlation matrix.

``` r
correlation <- function(x){
  x.diag <- diag(diag(x))
  
  corr <- sqrt(solve(x.diag)) %*% x %*% sqrt(solve(x.diag))
  
  return(corr)
}

correlation(D1.cov)
```

    ##             [,1]      [,2]        [,3]        [,4]
    ## [1,]  1.00000000 0.0923728 -0.09418076 -0.05668644
    ## [2,]  0.09237280 1.0000000  0.17769547  0.18267232
    ## [3,] -0.09418076 0.1776955  1.00000000 -0.03785468
    ## [4,] -0.05668644 0.1826723 -0.03785468  1.00000000

``` r
cor(D1)
```

    ##             V1        V2          V3          V4
    ## V1  1.00000000 0.0923728 -0.09418076 -0.05668644
    ## V2  0.09237280 1.0000000  0.17769547  0.18267232
    ## V3 -0.09418076 0.1776955  1.00000000 -0.03785468
    ## V4 -0.05668644 0.1826723 -0.03785468  1.00000000
