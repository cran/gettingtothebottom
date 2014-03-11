\name{plot_nnm_coef}
\alias{plot_nnm_coef}
\title{MM Algorithm - Plotting the NNMLS regression coefficients}
\usage{
  plot_nnm_coef(nnm_sol)
}
\arguments{
  \item{nnm_sol}{Solution object from nnls_mm function}
}
\description{
  \code{plot_nnm_coef} Function for plotting the NNMLS
  regression coefficients
}
\examples{
# Setup mixture example
n <- 1e3
p <- 10
nnm <- generate_nnm(n,p)

set.seed(12345)
X <- nnm$X
b <- double(p)
nComponents <- 3
k <- sample(1:p,nComponents,replace=FALSE)
b[k] <- matrix(runif(nComponents),ncol=1)
y <- X\%*\%b + 0.25*matrix(abs(rnorm(n)),n,1)

# Obtain solution to mixture problem
nnm_sol <- nnls_mm(y,X,runif(p))

# Plot the regression coefficients
plot_nnm_coef(nnm_sol)
}

