---
title: "Class 6 are functions!"
author: "Pin-Chung (Tony) Cheng"
date: "2019-10-18"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# This is a H1
This is my class 6 work. **hello**

## hello hello
```{r}
plot(1:10)
```

```{r}
test1 = read.csv("test1.txt")
test2 = read.table("test2.txt", sep= "$", header = TRUE)
test3 = read.table("test3.txt")
# The default sep="" means separate by one or more spaces
```

our function!!
```{r}
add <- function(x, y=1) {
 # Sum the input x and y
 x + y
}
```

```{r}
add(4)
add(c(5,4,4), c(2,1))
```
```{r}
rescale = function(x){ 
  xmin = min(x)
  (x - xmin) / (max(x) - xmin)
  }
```

```{r}
rescale2 = function(y){
  rng = range(y, na.rm = TRUE)
  (y - rng[1]) / (rng[2] - rng[1])
}
```

```{r}
rescale(1:10)
rescale2(1:10)
rescale2(c(1,2,NA,3,10))

```
```{r}
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
 if(na.rm) {
 rng <-range(x, na.rm=na.rm)
 } else {
 rng <-range(x)
 }
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
 return(answer)
}

```

```{r}
rescale3(1:10)
```

