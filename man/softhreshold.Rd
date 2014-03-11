\name{softhreshold}
\alias{softhreshold}
\title{MM Algorithm - Softhreshold Function}
\usage{
  softhreshold(x, lambda)
}
\arguments{
  \item{x}{Vector of values to be softhresholded}

  \item{lambda}{Softhreshold parameter}
}
\description{
  \code{softthreshold} Function for computing the
  softhreshold
}
\examples{
x <- seq(-10,10,1)
softhreshold(x,lambda=3)
}
\author{
  Jocelyn T. Chi
}

