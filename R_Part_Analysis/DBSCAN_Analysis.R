##### Essai DBSCAN

library(dbscan)

## Example 1: use dbscan on the iris data set
data(iris)
iris <- as.matrix(iris[, 1:4])

## Find suitable DBSCAN parameters:
## 1. We use minPts = dim + 1 = 5 for iris. A larger value can also be used.
## 2. We inspect the k-NN distance plot for k = minPts - 1 = 4
kNNdistplot(iris, minPts = 5)

## Noise seems to start around a 4-NN distance of .7
abline(h=.7, col = "red", lty = 2)

## Cluster with the chosen parameters
res <- dbscan(iris, eps = .7, minPts = 5)
res

pairs(iris, col = res$cluster + 1L)

## Use a precomputed frNN object
fr <- frNN(iris, eps = .7)
dbscan(fr, minPts = 5)

## Example 2: use data from fpc
set.seed(665544)
n <- 100
x <- cbind(
  x = runif(10, 0, 10) + rnorm(n, sd = 0.2),
  y = runif(10, 0, 10) + rnorm(n, sd = 0.2)
)

res <- dbscan(x, eps = .3, minPts = 3)
res

## plot clusters and add noise (cluster 0) as crosses.
plot(x, col = res$cluster)
points(x[res$cluster == 0, ], pch = 3, col = "grey")
hullplot(x, res)

## Predict cluster membership for new data points
## (Note: 0 means it is predicted as noise)
newdata <- x[1:5,] + rnorm(10, 0, .3)
hullplot(x, res)
points(newdata, pch = 3 , col = "red", lwd = 3)
text(newdata, pos = 1)

pred_label <- predict(res, newdata, data = x)
pred_label
points(newdata, col = pred_label + 1L,  cex = 2, lwd = 2)

## Compare speed against fpc version (if microbenchmark is installed)
## Note: we use dbscan::dbscan to make sure that we do now run the
## implementation in fpc.
if (FALSE) {
  if (requireNamespace("fpc", quietly = TRUE) &&
      requireNamespace("microbenchmark", quietly = TRUE)) {
    t_dbscan <- microbenchmark::microbenchmark(
      dbscan::dbscan(x, .3, 3), times = 10, unit = "ms")
    t_dbscan_linear <- microbenchmark::microbenchmark(
      dbscan::dbscan(x, .3, 3, search = "linear"), times = 10, unit = "ms")
    t_dbscan_dist <- microbenchmark::microbenchmark(
      dbscan::dbscan(x, .3, 3, search = "dist"), times = 10, unit = "ms")
    t_fpc <- microbenchmark::microbenchmark(
      fpc::dbscan(x, .3, 3), times = 10, unit = "ms")
    
    r <- rbind(t_fpc, t_dbscan_dist, t_dbscan_linear, t_dbscan)
    r
    
    boxplot(r,
            names = c('fpc', 'dbscan (dist)', 'dbscan (linear)', 'dbscan (kdtree)'),
            main = "Runtime comparison in ms")
    
    ## speedup of the kd-tree-based version compared to the fpc implementation
    median(t_fpc$time) / median(t_dbscan$time)
  }}

## Example 3: manually create a frNN object for dbscan (dbscan only needs ids and eps)
nn <- structure(list(ids = list(c(2,3), c(1,3), c(1,2,3), c(3,5), c(4,5)), eps = 1),
                class =  c("NN", "frNN"))
nn
dbscan(nn, minPts = 2)


