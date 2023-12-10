#' @title Reversible Jump Markov Chain Monte Carlo
#' @description When modeling count data, one question that is often of interest is whether the data is too dispersed for a Poisson distribution. If it is too dispersed, it may be more appropriate to use a negative binomial distribution to fit, so for a given data y, the reversible jump MCMC method is used for model selection. That's what this function rjmc does.
#' @param y Counting data
#' @param model Type of model
#' @param lambda The mean of Poisson distribution and negative binomial distribution
#' @param theta Prior distribution parameter
#' @param kkappa Negative binomial distribution parameter
#' @param alpha Prior distribution parameter
#' @param beta Prior distribution parameter
#' @param alpha.kappa Prior distribution parameter
#' @param beta.kappa Prior distribution parameter
#' @param sigma Normal distribution parameters constructed during model transfer
#' @param mu Invariants during model transfer
#' @param run Length of chain
#' @import bootstrap
#' @import boot
#' @import DAAG
#' @import microbenchmark
#' @return Model selection result
#' @export
rjmc <- function(y, model, lambda, theta, kkappa, alpha, beta, alpha.kappa, beta.kappa, sigma, mu, run){
  output <- matrix(0,nrow = run, ncol = 3)
  
  for (i in 1:run) {
    if(model==1){
      u <- rnorm(1, sd=sigma)
      kkappa <- mu*exp(u)
      lik1 <- prod(dpois(y,lambda))
      lik2 <- prod(dnbinom(y, 1/kkappa, 1-kkappa*lambda/(1+kkappa*lambda)))
      post1 <- dgamma(lambda, alpha, beta)*lik1
      post2 <- dgamma(lambda, alpha, beta)*dgamma(kkappa, alpha.kappa, beta.kappa)*lik2
      ratio <- post2/post1*mu*exp(u)/dnorm(u, sd=sigma)
      if(runif(1) < ratio){
        mod <- 2
        theta <- lambda
        kkappa <- kkappa
      }
    }
    
    if(model==2){
      u <- log(kkappa)/mu
      lik1 <- prod(dpois(y,theta))
      lik2 <- prod(dnbinom(y, 1/kkappa, 1-kkappa*theta/(1+kkappa*lambda)))
      post1 <- dgamma(theta, alpha, beta)*lik1
      post2 <- dgamma(theta, alpha, beta)*dgamma(kkappa, alpha.kappa, beta.kappa)*lik2
      ratio <- post1/post2*dnorm(log(kkappa/mu), sd=sigma)/kkappa
      if(runif(1) < ratio){
        mod <- 1
        lambda <- theta
      }
    }
    
    model <- mod
    if(model==1) output[i,] <- c(1,lambda,0)
    if(model==2) output[i,] <- c(2,theta,kkappa)
  }
  output
}
