#' Linear Programming - Least absolute deviations (LAD) regression
#' 
#' \code{ladlp} Function for performing least absolute deviations using linear programming
#' 
#' @param formula Formula for LAD regression (e.g., y~X)
#' 
#' @import lpSolve
#' 
#' @export
#' 
#' @examples
#' # Simulate data for LAD regression
#' m <- 1000 # number of observations
#' n <- 100 # number of variables
#' reg <- generate_data(m,n)
#' X <- reg$X
#' y <- reg$y
#' 
#' # Obtain LAD solution (not run)
#' # sol_lad <- ladlp(y~X)
#' # sol_lad
#' 
#' @author Jocelyn T. Chi
#' 
ladlp <- function(formula){
  
  # Process formula and call
  model_frame <- match.call(expand.dots = FALSE)
  temp <- match(c("formula"), names(model_frame), 0)
  model_frame <- model_frame[c(1, temp)]
  model_frame$drop.unused.levels <- TRUE
  model_frame[[1]] <- as.name("model.frame")
  model_frame <- eval.parent(model_frame)
  model_terms <- attr(model_frame, "terms")
  y <- model.response(model_frame)
  X <- model.matrix(model_terms, model_frame, contrasts)
  
  m <- nrow(X)
  n <- ncol(X)
  
  # Set up linear program for LAD
  f.obj <- c(rep(1,2*m),rep(0,2*n))
  A_top <- cbind(diag(m),-diag(m),X,-X)
  A_bottom <- diag((2*m+2*n))
  A <- rbind(A_top,A_bottom)
  f.con <- A
  f.rhs <- c(y,t(rep(0,2*m+2*n)))
  f.dir <- c(rep("=",length(y)),rep(">=",length(f.rhs)-length(y)))
  
  # Run the lpSolve
  lad_sol <- lp("min", f.obj, f.con, f.dir, f.rhs)
  
  # Extracting r and theta from the solution
  r_plus <- lad_sol$solution[1:m]
  r_minus <- lad_sol$solution[(m+1):(2*m)]
  r <- r_plus - r_minus
  theta_plus <- lad_sol$solution[((2*m)+1):((2*m)+n)]
  theta_minus <- lad_sol$solution[((2*m)+n+1):length(lad_sol$solution)]
  theta <- theta_plus - theta_minus
  
  return(list(coefficients=theta,residuals=r))
}

#' Linear Programming - Linear programming solver for quantile regression
#' 
#' \code{quantreglp} Function for performing quantile regression using linear programming
#' 
#' @param y An mx1 vector containing the response variables in the model
#' @param X An nxp matrix containing the predicting variables in the model
#' @param tau (optional) quantile
#' @param lambda (optional) regularization parameter
#' 
#' @import lpSolve
#' 
#' @export
#' 
#' @examples
#' set.seed(12345)
#' n <- 20
#' p <- 20
#' X <- matrix(rnorm(n*p),n,p)
#' b0 <- double(p)
#' k <- 4
#' b0[sample(1:p,k,replace=FALSE)] <- 10*rnorm(k)
#' y <- X%*%b0 + 0.1*rnorm(n)
#' 
#' lambda <- 2
#' tau <- 0.5
#' sol <- quantreglp(y,X,tau,lambda)
#'
#' @author Jocelyn T. Chi
#' 
quantreglp <- function(y,X,tau=0.5,lambda=0){
  
  # Set up LP
  n <- nrow(X)
  p <- ncol(X)
  I <- diag(1,n)
  A_top <- cbind(I,-I,X,-X)
  A_bottom <- diag(ncol(A_top))
  A <- rbind(A_top,A_bottom)
  
  # Setup lp for lpSolve
  f.obj <- c(rep(tau,n),rep(1-tau,n),rep(lambda,2*p))
  f.con <- A
  f.dir <- c(rep("=",nrow(A_top)),rep(">=",ncol(A_bottom)))
  f.rhs <- c(y,rep(0,ncol(A_bottom)))
  temp <- lp("min", f.obj, f.con, f.dir, f.rhs)
  
  # Extract solution
  r_plus <- temp$solution[1:n]
  r_minus <- temp$solution[(n+1):(2*n)]
  r <- r_plus - r_minus
  theta_plus <- temp$solution[((2*n)+1):((2*n)+p)]
  theta_minus <- temp$solution[((2*n)+p+1):(2*(n+p))]
  theta <- theta_plus - theta_minus
  
  return(list(coefficients=theta, residuals=r))
  
#   if (intercept == TRUE){
#     intercept <- theta[1]
#     theta <- theta[2:length(theta)]
#     
#     # Return solution
#     return(list(intercept=intercept,coefficients=theta, residuals=r))
#   } 
#   else {
#     return(list(coefficients=theta, residuals=r))
#   }
  
}

#' Linear Programming - Quantile regression 
#' 
#' \code{quantreg} quantreg is used to fit quantile regression models.
#' 
#' @param formula An object of class "formula" (e.g., y~X)
#' @param tau (optional) Quantile for quantile regression.  Default is 0.5 but can also be a vector (e.g., tau = c(0.25,0.5,0.75)).
#' @param lambda (optional) Regularization parameter.  Default is lambda=0 (i.e., no regularization).
#' 
#' @export
#' 
#' @examples
#' set.seed(12345)
#' n <- 20
#' p <- 20
#' X <- matrix(rnorm(n*p),n,p)
#' b0 <- double(p)
#' k <- 4
#' b0[sample(1:p,k,replace=FALSE)] <- 10*rnorm(k)
#' y <- X%*%b0 + 0.1*rnorm(n)
#' 
#' lambda <- 0
#' tau <- c(0.05,0.25,0.5,0.75,0.95)
#' sol <- quantreg(y~X,tau,lambda)
#' coef(sol)
#'
#' @author Jocelyn T. Chi
#' 
quantreg <- function(formula, tau=0.5, lambda=0) { 
  
  # Process formula and call
  model_frame <- match.call(expand.dots = FALSE)
  temp <- match(c("formula"), names(model_frame), 0)
  model_frame <- model_frame[c(1, temp)]
  model_frame$drop.unused.levels <- TRUE
  model_frame[[1]] <- as.name("model.frame")
  model_frame <- eval.parent(model_frame)
  model_terms <- attr(model_frame, "terms")
  y <- model.response(model_frame)
  X <- model.matrix(model_terms, model_frame, contrasts)
  
  # Begin quantile regression
  t <- length(tau)
  if (t > 1) {
    if (any(tau <= 0) || any(tau >= 1)) {
      stop("Values for tau must be between 0 and 1.") 
    }
    coef <- matrix(0, ncol(X), t)
    resids <- matrix(0, nrow(X), t)
    for (i in 1:t) {
      sol <- quantreglp(y, X, tau=tau[i], lambda=lambda)
      coef[, i] <- sol$coefficients
      resids[, i] <- sol$residuals
      tlabels <- paste("quantile:", format(round(tau, 3)))
      dimnames(coef) <- list(c(dimnames(X)[[2]]), tlabels)
    }
  }
  else {
    if ( tau <= 0 || tau >= 1 ) {
      stop("Values for tau must be between 0 and 1.") 
    }
    sol <- quantreglp(y, X, tau=tau, lambda=lambda)
    coef <- sol$coefficients
    resids <- sol$residuals
  }
  
  return( list(coefficients=coef, residuals=resids) )
  
}

#' Linear Programming - Generate simulated data for LAD regression
#' 
#' \code{generate_data} Function for generating simulated data for LAD regression
#' 
#' @param n Number of data points
#' @param p Number of total covariates
#' @param sigma Standard deviance
#' @param seed Random seed
#' 
#' @export
#' 
#' @examples
#' generate_data(n=1000,p=100)
#' 
#' @author Jocelyn T. Chi
#'
generate_data <- function(n=100,p=1,sigma=1,seed=12345) {
  set.seed(seed)
  X <- matrix(rnorm(n*p),n,p)
  b <- rnorm(p)
  y <- X%*%b + sigma*rnorm(n)
  return(list(y=y,X=X,b=b))
}

#' Linear Programming - Generate simulated data for sparse quantile regression
#' 
#' \code{generate_sparse_data} Function for generating simulated data for sparse quantile regression
#' 
#' @param n Number of data points
#' @param p Number of total covariates
#' @param k Number of true covariates
#' @param seed Random seed
#' 
#' @export
#' 
#' @examples
#' data <- generate_sparse_data(n=100,p=20,k=3)
#'
#' @author Jocelyn T. Chi
#'
generate_sparse_data <- function(n=100,p=20,k=3,seed=12345) {
  set.seed(seed)
  X <- matrix(rnorm(n*p),n,p)
  b <- double(p)
  b[1:k] <- 0.5
  y <- X%*%b + rnorm(n)
  return(list(y=y,X=X,b=b))
}

#' Linear Programming - Function for plotting results of bivariate quantile regression.
#' 
#' \code{plot_quantreg} Function for plotting results of bivariate quantile regression.
#' 
#' @param x Variable on the x-axis.
#' @param y Variable on the y-axis.
#' @param taus (optional) Quantiles for plotting, default is 'c(0.05,0.25,0.5,0.75,0.95)'.
#' @param xlab (optional)  Default is x.
#' @param ylab (optional)  Default is y.
#' 
#' @import ggplot2
#' 
#' @export
#' 
#' @examples
#' data(baltimoreyouth)
#' x <- baltimoreyouth$teenbir10 # Teen births
#' y <- baltimoreyouth$compl10 # Highschool completion rate
#' plot_quantreg(x,y,xlab="Teen Births per 1000 Females (aged 15-19)", 
#'      ylab="High School Completion Rates")
#' 
#' @author Jocelyn T. Chi
#' 
plot_quantreg <- function(x,y,xlab='x',ylab='y',taus=c(0.05,0.25,0.5,0.75,0.95)){
  int <- Quantile <- NULL
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  if (ncol(x) > 1 | ncol(y) > 1){
    stop("Inputs must be vectors or matrices of dimension 1.")
  }
  
  # Make data frames for plotting
  obj <- t(coef(quantreg(y~x,tau=taus)))
  obj <- cbind(obj,taus)
  colnames(obj) <- c('int', 'x', 'Quantile')
  obj <- data.frame(obj)
  
  lad <- data.frame(t(coef(quantreg(y~x,tau=0.5))))
  colnames(lad) <- c('int','x')
  
  # Plot the data
  df1 <- data.frame(xlab=x, ylab=y)
  plot <- ggplot(df1, aes(x=xlab, y=ylab))
  plot <- plot + geom_point(colour="black") + theme_bw() + theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()) +  theme(axis.line = element_line(color = 'black')) + xlab(xlab) + ylab(ylab)
  plot <- plot + geom_abline(aes(intercept=int,slope=x,colour=Quantile),data=obj) + scale_colour_gradient2(high="#000000")
  plot <- plot + stat_smooth(method="lm", se=FALSE, colour="#31a354") 
  plot <- plot + geom_abline(aes(intercept=int,slope=x),colour="#4eb3d3",data=lad)
  plot
}

#' Linear Programming - Function for plotting results of sparse quantile regression
#' 
#' \code{plot_noisyquantreg} Function for plotting results of bivariate quantile regression with noise.
#' 
#' @param y Response variable
#' @param X Noisy matrix (includes x)
#' @param taus (optional) Quantiles for plotting, default is 'c(0.05,0.25,0.5,0.75,0.95)'
#' @param lambdas (optional) Regularization parameters (Default is lambdas=c(0,0.5,1,1.5).)
#' 
#' @import ggplot2
#' @import reshape2
#' 
#' @export
#' 
#' @examples
#' data <- generate_sparse_data(n=100,p=20,k=3)
#' y <- data$y
#' X <- data$X
#' 
#' lambdas=c(0,1,5,7)
#' plot_quantreg_noisy(y,X,lambdas=lambdas) # Not run
#' 
#' @author Jocelyn T. Chi
#' 
plot_quantreg_noisy <- function(y,X,taus=c(0.05,0.25,0.5,0.75,0.95),lambdas=c(0,0.5,1,1.5)){
  Variables <- Coefficient <- NULL
  
  X <- data.frame(X)
  nameslist <- colnames(X)
  
  # Make lists
  l <- length(lambdas)
  coef_noisy <- vector("list", length=length(lambdas))
  names(coef_noisy) <- lambdas
  
  # Run quantregs
  X <- as.matrix(X)
  coef_qr <- coef(quantreg(y~X,taus))
  for (i in 1:l) {
    # Get data
    coef_noisy[[i]] <- coef(quantreg(y~X,taus,lambda=lambdas[i]))
    coef_noisy[[i]] <- cbind(rep(lambdas[i],nrow(coef_noisy[[i]])),coef_noisy[[i]])
    
    # Add column for lambda & fix row names
    colnames(coef_noisy[[i]])[1] <- "lambda"
    rownames(coef_noisy[[i]]) <- c('intercept',nameslist)
    
    # Add column for row.names (because duplicate rownames not allowed later in merging)
    coef_noisy[[i]] <- cbind(c('intercept',nameslist),coef_noisy[[i]])
    colnames(coef_noisy[[i]])[1] <- "Variables"
  }
  
  # Combine quantreg data into dataframe for ggplot2
  suppressWarnings(coef_noisy <- do.call(rbind.data.frame, coef_noisy))
  dfnoisy <- subset(coef_noisy, Variables!='intercept')
  dfnoisy$Variables <- factor(nameslist, levels = nameslist, ordered=TRUE)
  suppressWarnings(dfnoisy2 <- melt(dfnoisy, id.vars=c('Variables','lambda')))
  colnames(dfnoisy2) <- c("Variables","Lambda","Quantile","Coefficient")
  class(dfnoisy2$Coefficient) <- "numeric"
  
  # Plot the data
  p <- ggplot(dfnoisy2, aes(Variables, Coefficient)) + geom_point()
  p <- p + facet_grid(Lambda ~ Quantile) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p
}

#' Linear Programming - Generate matrix of noisy data for sparse quantile regression
#' 
#' \code{make_noise} Function for generating matrix of noisy data for sparse quantile regression
#' 
#' @param y Response vector
#' @param X True predictor
#' @param p Number of desired columns of noise
#' @param seed (optional) Random seed (Default seed=12345)
#' 
#' @export
#' 
#' @examples
#' y <- engel$foodexp
#' X <- engel$income
#' noisydata <- make_noise(y,X,10)
#' head(noisydata)
#' 
#' @author Jocelyn T. Chi
#' 
make_noise <- function(y,X,p,seed=12345) {
  X <- as.matrix(X)
  m <- nrow(X)
  truth <- X
  set.seed(seed)
  noise <- abs(matrix(sample(truth,m*p,replace=TRUE),length(truth),p))
  noise <- matrix(runif(m*p,min=min(truth),max=max(truth)),m,p)
  colnames(noise) <- c(1:p)
  data_noise <- as.matrix(cbind(truth,noise))
  if (ncol(X) == 1){
    colnames(data_noise)[1] <- 'truth'
  }
  return(data_noise)
}

#' Linear Programming - Check function
#' 
#' \code{check_func} Check function
#' 
#' @param r A vector to feed into the check function
#' @param tau Tau values to be used in the check function
#' 
#' @export
#' 
#' @examples
#' r <- seq(-1,1,length.out=1e4)
#' tau <- 0.5
#' check_func(r,tau)
#' 
check_func <- function(r,tau) {
  rho <- tau*r
  ix <- which(r < 0)
  rho[ix] <- -(1-tau)*r[ix]
  return(rho)
}

#' Linear Programming - Plot Check Function
#' 
#' \code{plot_check} Plot Check Function
#' 
#' @param r A vector to feed into the check function
#' @param taus A vector of tau values to be used in the check function
#' 
#' @export
#' 
#' @import ggplot2
#' 
#' @examples
#' r <- seq(-1,1,length.out=1e4)
#' taus <- c(0.25,0.5,0.75)
#' plot_check(r,taus)
#'
plot_check <- function(r,taus=taus){
  x <- y <- tau <- NULL
  data <- data.frame(x=c(),y=c(),tau=c())
  for (i in 1:length(taus)) {
    temp <- data.frame(x=r,y=check_func(r,taus[i]),tau=taus[i])
    data <- rbind(data,temp)
  }
  data$tau <- factor(data$tau)
  q <- ggplot(data=data,aes(x=x,y=y,group=tau)) + geom_line(aes(colour=tau,linetype=tau))
  q <- q + xlab("r") + ylab(expression(rho[tau](r))) +
    scale_colour_manual(name=expression(tau),values = c("#a8ddb5","#4eb3d3", "#08589e")) +
    scale_linetype_manual(name=expression(tau),values=c(1:3))
  q + theme_bw()
}