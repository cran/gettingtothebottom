\name{plot_loss}
\alias{plot_loss}
\title{Gradient Descent Algorithm - Plotting the Loss Function}
\usage{
  plot_loss(obj)
}
\arguments{
  \item{obj}{Object containing the results of a gradient
  descent implementation}
}
\description{
  \code{plot_loss} Plots the loss function of an object
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

# Plot the loss function
plot_loss(simple_ex)
}
\author{
  Jocelyn T. Chi
}

