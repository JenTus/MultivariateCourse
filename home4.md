home4
================

Q1
--

Plot the empirical influence function of the median.

``` r
library(ggplot2)

x <- sort(runif(100, 0, 100))
y <- seq(from = 0, to = 100, length.out = 100)
IF <- seq(from = 0, to = 0, length.out = 100)

for (i in 1:100){
  IF[i] <- 101*(median(c(x, y[i])) - median(x))
}

df <- data.frame(cbind(y, IF))

ggplot(df, aes(df$y, df$IF)) +
  geom_point() +
  labs(x = "y", y = "IF")
```

![](home4_files/figure-markdown_github/cars-1.png)
