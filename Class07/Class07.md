Class 7 R function and packages
================
Pin-Chung (Tony) Cheng
10/23/2019

``` r
source("http://tinyurl.com/rescale-R")
```

Try rescle function
    again

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale(c(3,10,NA,7))
```

    ## [1] 0.0000000 1.0000000        NA 0.5714286

``` r
rescale2(c(3,10,NA,7))
```

    ## [1] 0.0000000 1.0000000        NA 0.5714286

``` r
both_na4 = function(x, y) {
  if (length(x) != length(y)) {
    stop("The two vectors are not the same length!!!")
    }
  sum( is.na(x) & is.na(y) )
}
```

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA)

both_na4(x, y2)
```

    ## [1] 2

``` r
a = "tony cheng12345 hello"
sub("cheng)", "\\1", a)
```

    ## [1] "tony cheng12345 hello"

``` r
a = c(100, 100, 100, 100, 100, 90)
r = which.min(a)

a[-r]
```

    ## [1] 100 100 100 100 100

``` r
grade = function(x) {
  
  x_minus_one = x[-which.min(x)]
  sum(x_minus_one, na.rm= TRUE) / length(x_minus_one)
  
}
```

``` r
# student 1
s1 = c(100, 100, 100, 100, 100, 100, 100, 90)
# student 2
s2 = c(100, NA, 90, 90, 90, 90, 97, 80)
```

``` r
grade(s1)
```

    ## [1] 100

``` r
grade(s2)
```

    ## [1] 79.57143
