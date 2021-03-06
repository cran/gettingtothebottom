% Generated by roxygen2 (4.0.2): do not edit by hand
\name{makeY}
\alias{makeY}
\title{MM Algorithm - Make Y}
\usage{
makeY(X, Z, omega)
}
\arguments{
\item{X}{Matrix containing observed entries}

\item{Z}{Matrix containing last iterates}

\item{omega}{Vector containing indices of unobserved entries (by column)}
}
\description{
\code{makeY} Function for making the Y matrix
}
\examples{
n <- 5
A <- matrix(rnorm(n^2),n,n)
omega <- c(1,5,8,10,16,23)
Z <- Matrix(0,n,n,sparse=TRUE)
makeY(A,Z,omega)
}
\author{
Jocelyn T. Chi
}

