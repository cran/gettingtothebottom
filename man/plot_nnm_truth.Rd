\name{plot_nnm_truth}
\alias{plot_nnm_truth}
\title{MM Algorithm - Plotting the True Signal}
\usage{
  plot_nnm_truth(X, b, nnm)
}
\arguments{
  \item{nnm}{NNM object from generate_nnm function}

  \item{X}{Nonnegative design matrix}

  \item{b}{Nonnegative initial regression vector}
}
\description{
  \code{plot_nnm_truth} Function for plotting the true
  mixture signal
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

# Plot the truth
plot_nnm_truth(X,b,nnm)
}

