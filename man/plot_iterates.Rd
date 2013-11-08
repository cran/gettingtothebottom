\name{plot_iterates}
\alias{plot_iterates}
\title{Gradient Descent Algorithm - Plotting the Iterates}
\usage{
  plot_iterates(obj)
}
\arguments{
  \item{obj}{Object containing the results of a gradient
  descent implementation}
}
\description{
  \code{plot_iterates} Plots the iterates of an object
  containing the results of a gradient descent object
  implementation
}
\examples{
# Generate some data for a simple bivariate example
set.seed(12345)
x <- sample(seq(from = -1, to = 1, by = 0.1), size = 50, replace = TRUE)
y <- 2*x + rnorm(50)

# Components required for gradient descent
X <- as.matrix(x)
y <- as.vector(y)
f <- function(X,y,b) {
   (1/2)*norm(y-X\%*\%b,"F")^{2}
}
grad_f <- function(X,y,b) {
   t(X)\%*\%(X\%*\%b - y)
}

# Run a simple gradient descent example
simple_ex <- gdescent(f,grad_f,X,y,alpha=0.01)

# Plot the iterates
plot_iterates(simple_ex)
}
\author{
  Jocelyn T. Chi
}

