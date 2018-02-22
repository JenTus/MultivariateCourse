Bivariate Correspondence Analysis
================

Bivariate Correspondence Analysis
---------------------------------

The data SMOKING1.txt contains a two-dimensional frequency table, where the employees of a company have been categorized according to their position (5 cat- egories: SM = Senior Managers, JM= Junior Managers, SE = Senior Employees, JE = Junior Employees, SC = Secretaries). The smoking of the employees have 4 categories (None, Light, Medium, Heavy).

Q1
--

Form the row and column profiles.

``` r
smok <- read.table("data/SMOKING.txt", header = TRUE, row.names = 1)
smok <- smok[(1:5),(1:4)]
row.sum <- apply(smok, 1, sum)
row.profile <- sweep(smok, 1, row.sum, "/")

col.sum <- apply(smok, 2, sum)
col.profile <- sweep(smok, 2, col.sum, "/")
```

Q2
--

How much of the variation is explained by the combination of components 1 and 3. Give the answer in percentages relative to the total variation.

``` r
library(ca)
smok.ca <- ca(smok)
```

Q3
--

Produce the BCA graph with respect to the first two components.

``` r
plot(smok.ca, arrows=c(TRUE, T))
```

![](home6_files/figure-markdown_github/unnamed-chunk-3-1.png)

Q4
--

How much of the variation of the modality Heavy is explained by the first two components? How much of the variation of the modality Medium is explained by the first two components? Give the answers in percentages relative to the total variation. Hint: use the quality of representation.

``` r
n <- sum(smok)
ca.row <- smok.ca$rowcoord
f_row <- row.sum / n

ca.col <- smok.ca$colcoord
f_col <- col.sum / n 

#the modality Heavy explained by the first two components.
sum(ca.col["Heavy",c(1,2)]*ca.col["Heavy",c(1,2)])/sum(ca.col["Heavy",]*ca.col["Heavy",])
```

    ## [1] 0.752805

``` r
#the modality Medium explained by the first two components.
sum(ca.col["Medium",c(1,2)]*ca.col["Medium",c(1,2)])/sum(ca.col["Medium",]*ca.col["Medium",])
```

    ## [1] 0.2465586