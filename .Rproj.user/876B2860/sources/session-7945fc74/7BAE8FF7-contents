## functions for data generation
x.gen <- function(seed, n, p){
  set.seed(seed)
  #return(matrix(rep(1, n), ncol = 1))

  return(matrix(runif(n*p, -3, 3), ncol = p))
}





#' Generate  Data for the Bayesian Mixed-Effects Model for Multilevel Two-Way Functional Data (BMEF)
#'
#' This function generates synthetic data to illustrate and test the Bayesian Mixed-Effects Model
#' for Multilevel Two-Way Functional Data (BMEF). The simulated data include functional fixed effects and random effects, and noise
#' components over a time-frequency grid. Natural cubic B-splines are used as marginal bases for both the time and frequency dimensions.

#'
#' @param seed An integer seed for reproducibility.
#' @param J Number of experimental conditions.
#' @param TT Number of time points (equally spaced between 0 and 1).
#' @param FF Number of frequency points (equally spaced between 0 and 1).
#' @param K_T Number of marginal basis functions for the time dimension. (integer)
#' @param K_F Number of marginal basis functions for the frequency dimension. (integer)
#' @param U An orthonormal factor matrix of dimension \code{K_F × R}
#' @param V  An orthonormal factor matrix of dimension \code{K_T × R}
#' @param Delta A 3D array of dimension \code{J × R × p}, representing covariate coefficients across conditions and components.
#' @param X A numeric matrix of subject-level covariates of dimension \code{n × p}.
#' The intercept should be included explicitly (e.g., for intercept-only, \code{p = 1}).
#' @param Sigma2S A list of variance components including:
#' \itemize{
#'   \item{\code{Sigma2_gamma}}: {A scalar representing \eqn{\sigma_{\gamma}^2}.}
#'   \item{\code{Sigma2_omega}}: {A numeric vector of length \code{n}, representing \eqn{\sigma_{\omega_i}^2}.}
#'   \item{\code{Sigma2_epsilon}}: {A scalar representing \eqn{\sigma_{\epsilon}^2}.}
#' }
#' @return A list containing the following components:
#'  \item{tt}{A numeric vector of length \code{TT}, representing the time grid (equally spaced from 0 to 1).}
#'  \item{ff}{A numeric vector of length \code{FF}, representing the frequency grid (equally spaced from 0 to 1).}
#'  \item{O_tilde}{A matrix of dimension \code{(TT× FF) × K}, representing orthogonalized basis evaluations.}
#'  \item{A}{A 4D array of dimension \code{n × J × K_T × K_F}, storing the fixed effect coefficients \eqn{\alpha_j(x_i)} (for each [i,j] pair, arranged columnwise as a \code{K_T × K_F} matrix).}
#'  \item{B}{A 3D array of dimension \code{n × K_T × K_F}, storing subject-level random effects coefficients \eqn{\gamma_i} (for each [i], arranged ccolumnwise as a \code{K_T × K_F} matrix).}
#'  \item{C}{A 4D array of dimension \code{n × J × K_T × K_F}, storing subject-by-condition random effects coefficients \eqn{\omega_{i,j}} (for each [i,j] pair, arranged columnwise as a \code{K_T × K_F} matrix).}
#'  \item{Y}{A 3D array of dimension \code{n × J × (K_T × K_F)}, representing the observed two-way functional observations (with noise)}
#'
#' @export
#'

dat.gen <- function(seed, n, J, TT, FF,  K_T, K_F, U, V, Delta, X, Sigma2s){

  set.seed(seed)

  K <- K_T * K_F
  tt <- seq(0, 1, length.out = TT)  # time points
  ff <- seq(0, 1, length.out = FF) # frequency points

  # generate the basis evaluations
  bsMat_tt <- ns(tt,  knots = seq(0, 1, length.out = K_T)[2:(K_T-1)], intercept = TRUE)
  bsMat_ff <- ns(ff,  knots = seq(0, 1, length.out = K_F)[2:(K_F-1)], intercept = TRUE)

  O <- kronecker(bsMat_ff, bsMat_tt)


  C_t <- svd(bsMat_tt)$v
  C_f <- svd(bsMat_ff)$v

  O_tilde <- kronecker(bsMat_ff%*%C_f, bsMat_tt%*%C_t)
  #O%*%kronecker(C_f, C_t)  # same result as the previous line


  # generate the basis coefficients
  A <- array(0, dim = c(n, J, K_T, K_F))
  B <- array(NA, dim = c(n, K_T, K_F))
  C <- array(NA, dim = c(n, J,  K_T, K_F))
  R <- dim(Delta)[2]

  for(i in 1:n){
    for(j in 1:J){
      for(r in 1:R){
        A[i, j,,] <-   A[i, j,,] + as.numeric(t(Delta[j,r,])%*%t(matrix(X[i,], ncol = ncol(X))))* U[, r]%*%t(V[,r])
      }
    }
  }

  for(i in 1:n){
    B[i,,] <- matrix(rnorm(K, 0, sd = sqrt(Sigma2s$Sigma2_gamma)), K_T, K_F)
  }


  for(i in 1:n){
    for(j in 1:J){
    C[i,j,, ] <- matrix(rnorm(K, 0, sd = sqrt(Sigma2s$Sigma2_omega)[i]), K_T, K_F)
    }
  }


  # generate the responses
  Y <- Y_clean <- array(NA, dim = c(n, J, TT*FF))
  for(i in 1:n){
    for(j in 1:J){
      Y_clean[i,j,] <-  O_tilde %*% c(A[i,j,,]+ B[i,,] + C[i,j,,])
      Y[i,j,] <- Y_clean[i,j,] + rnorm(TT*FF, 0, sd =  sqrt(Sigma2s$Sigma2_epsilon))
    }
  }


  dat <- list(tt = tt, ff = ff, O_tilde = O_tilde, Y = Y, A = A, B = B, C = C)

  return(dat)
}

