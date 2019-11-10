---
title: "class08_RMD"
author: "Pin-Chung (Tony) Cheng"
date: "10/25/2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


# K mean example
```{r}
tmp <- c(rnorm(30,-3), rnorm(30,3))
x <- cbind(x=tmp, y=rev(tmp))
plot(x)
```

Use kmeans function
```{r}
k = kmeans(x, centers=2, nstart = 20)
```

print it out
```{r}
print(k)
```

check ur cluster
```{r}
table(k$cluster)
```


plot with k cluster, add cluster centers as blue points
```{r}
plot(x, col = k$cluster)
# points function add features to existing plot
points(k$centers, col="blue", pch=15)
```

Hierarchical clustering in R
```{r}
hc = hclust(dist(x))
hc
```

plot hc
```{r}
plot(hc)

# draw a line to cut tree! and here u get two big branch, 2 clusters!
# height is like the distance difference
abline(h=6, col="red")
abline(h=3, col= "blue")
# return a vector of groups
grp1 = cutree(hc, h=6)
grp2 = cutree(hc, h=4)
grp3 = cutree(hc, h=3)

# define how many k group
grp4 = cutree(hc, k=2)
```
check ut group with table()
```{r}
table(grp4)
```

plot
```{r}
plot(x, col=grp4)

```

Slide sample
```{r}
# Step 1. Generate some example data for clustering
x <- rbind(
 matrix(rnorm(100, mean=0, sd = 0.3), ncol = 2), # c1
 matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2), # c2
 matrix(c(rnorm(50, mean = 1, sd = 0.3), # c3
 rnorm(50, mean = 0, sd = 0.3)), ncol = 2))
colnames(x) <- c("x", "y")
# Step 2. Plot the data without clustering
plot(x)
# Step 3. Generate colors for known clusters
# (just so we can compare to hclust results)
col <- as.factor( rep(c("c1","c2","c3"), each=50) )
plot(x, col=col)
```

Make it into hc
```{r}
hc = hclust(dist(x))
plot(hc)

grp5 = cutree(hc, k=3)
table(grp5, col)
plot(x, col=grp5)
```

```{r}
mydata <- read.csv("https://tinyurl.com/expression-CSV", row.names = 1) 
head(mydata)
str(mydata)
```

prcomp function 
```{r}
pca = prcomp(t(mydata), scale = TRUE)
# check whats in pca
attributes(pca)
```

```{r}
# get pca components in x , and column means pca1, pca2 ....
plot(pca$x[,1], pca$x[,2])
```

```{r}
summary(pca)
```

```{r}
## Variance captured per PC
pca.var <- pca$sdev^2 
# try plot pca
plot(pca)

```

```{r}
plot(pca$x[,1:2], col=c("red","red","red","red","red","blue","blue","blue","blue","blue"))
```

```{r}
col = c(rep("green",5), rep("blue",5))
plot(pca$x[,1:2], col=col)       
```

Hands on exercise
```{r}
x <- read.csv("UK_foods.csv", row.names = 1)
dim(x)
```

```
```{r}
pairs(x, col=rainbow(10), pch=16)
```

pca
```{r}
pca = prcomp(t(x))
summary(pca)
```

```{r}
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2", xlim=c(-270,500))
text(pca$x[,1], pca$x[,2], colnames(x), col=c("red", "blue", "green","orange"))
```

```{r}
par(mar=c(10, 3, 0.35, 0))
barplot( pca$rotation[,1], las=2 )
```

```{r}
## Lets focus on PC2 as it accounts for > 90% of variance 
par(mar=c(10, 3, 0.35, 0))
barplot( pca$rotation[,2], las=2 )
```

biplot
```{r}
biplot(pca)
```

