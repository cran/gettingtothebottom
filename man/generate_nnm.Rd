\name{generate_nnm}
\alias{generate_nnm}
\title{Generate random nonnegative mixture components}
\usage{
  generate_nnm(n, p, seed = 12345)
}
\arguments{
  \item{n}{Number of samples}

  \item{p}{Number of components}

  \item{seed}{Random seed}
}
\description{
  \code{generate_nnm} Function to random nonnegative
  mixture components
}
\examples{
n <- 1e3
p <- 10
nnm <- generate_nnm(n,p)
}

