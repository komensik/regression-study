library(rdrobust)

set.seed(12231)

# simulate the data
x <- rnorm(2000, 50, 20)
x <- x[x >= 0 & x < 100]

# c_0 at 50
D <- ifelse(x > 50, 1, 0)

# Observed Y
y <- 20 + D * 100 + 1.5 * x - 0.015 * x^2 + rnorm(length(x), 0, 100)


# Simple plot of raw data
plot(x, y, pch = 19, cex = 0.2)
abline(v = 50, col = "red")

# Binned plot (using 50 bins left and right of c_0)
# with p=2 polynomial fitted below and above c_0
# (h = full span of data)
rdplot(y = y, x = x, c = 50, p = 4, nbins = c(50, 50))

# automatic bin N selection
# (selects MSE optimal choice)
rdplot(y = y, x = x, c = 50, binselect = "esmv")

# Illustration: the role of h
rdplot(y = y, x = x, c = 50, binselect = "esmv", h = 50)
rdplot(y = y, x = x, c = 50, binselect = "esmv", h = 30)
rdplot(y = y, x = x, c = 50, binselect = "esmv", h = 10)

# BW selection
summary(rdbwselect(y = y, x = x, c = 50, bwselect = "mserd"))
summary(rdbwselect(y = y, x = x, c = 50, bwselect = "msetwo"))


# Simple RDD est, p = 1 (locally linear)
# optimal bw selected (same below and above c_0)
m_p1 <- rdrobust(y = y, x = x, c = 50, p = 1, bwselect = "mserd")
# Simple RDD est, p = 2 (locally quadratic)
m_p2 <- rdrobust(y = y, x = x, c = 50, p = 2, bwselect = "mserd")
summary(m_p1)
summary(m_p2)
 

# Construct plot (using obs within selected bw)
xc <- x - 50
m <- rdrobust(y = y, x = xc, p = 1)
rdplot(y = y, x = xc, 
    subset = -m$bws[1,1]<= xc & xc <= m$bws[1,2],
    binselect = "esmv", kernel = "triangular", 
    h = c(m$bws[1,1], m$bws[1,2]), p = 1)




# Notable options:

# choose other VCE: e.g., option vce("hc2") for the usual plug-in estimator
# the default in the package is to use the sample covariance estimator of the 
# J nearest neighbors to i (you can choose J using option nnmatch())

# Clustering: option cluster()

# add covariates: option covs()
