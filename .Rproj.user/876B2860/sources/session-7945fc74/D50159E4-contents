#' Fit the Bayesian Mixed-Effects Model for Multilevel Two-Way Functional Data (BMEF)
#'
#' This function fits the proposed Bayesian mixed-effects models for multilevel two-way functional data (BMEF).
#' Two frameworks are implemented: \code{"BMEF-1"}, which accounts for heterogeneous variances, and \code{"BMEF-2"}, which assumes homogeneous variance.
#'
#' The inputs are multilevel two-way functional data evaluated for \code{n} subjects,
#' \code{J} conditions, \code{TT} time points, and \code{FF} frequency points. This implementation assumes uniform priors on the factor matrices \code{U} and \code{V}.
#'
#' @param Y A 3D array of functional responses of dimension \code{n × J × (TT × FF)}.
#' The third dimension represents the column-vectorized time-frequency matrix of dimension \code{TT × FF} for each subject-condition pair.
#' @param X A numeric matrix of subject-level covariates of dimension \code{n × p}.
#' The intercept should be included explicitly (e.g., for intercept-only, \code{p = 1}).
#' @param JJ A binary matrix of dimension \code{n × J}, where \code{JJ[i, j] = 1} if subject \code{i} has data observed for condition \code{j}, and \code{0} otherwise.
#' @param tt A numeric vector of time points defining the temporal grid.
#' @param ff A numeric vector of frequency points defining the frequency grid.
#' @param R The maximum CP rank.
#' @param K_T Number of marginal basis functions for the time dimension.
#' @param K_F Number of marginal basis functions for the frequency dimension.
#'
#' @param n_burn Number of burn-in MCMC iterations.
#'
#' @param n_sample Number of post-burn-in MCMC samples to retain.
#'
#' @param params A list specifying prior hyperparameters. The list includes:
#' \itemize{
#'   \item{\code{Sigma2_delta}}: {A numeric vector of length \code{p}, specifying the diagonal elements of the covariance matrix \eqn{\Sigma_{\delta}} for the prior for \eqn{\delta_{j,r}}. Defaults to \code{rep(10, p)}.}
#'   \item{\code{a_gamma}}: {\eqn{a_{\gamma}}  for the Inverse-Gamma prior on \eqn{\sigma_{\gamma}^2}. Defaults to  \code{1}.}
#'   \item{\code{b_gamma}}: {\eqn{b_{\gamma}}  for the Inverse-Gamma prior on \eqn{\sigma_{\gamma}^2}. Defaults to  \code{1}.}
#'   \item{\code{a_omega}}: {\eqn{a_{\omega}}  for the Inverse-Gamma prior on \eqn{\sigma_{\omega_i}^2}. Defaults to  \code{1}.}
#'   \item{\code{b_omega}}: {\eqn{b_{\omega}}  for the Inverse-Gamma prior on \eqn{\sigma_{\omega_i}^2}. Defaults to  \code{1}.}
#'   \item{\code{a_delta}}: {\eqn{a_{\delta}}  for the Beta prior on \eqn{\delta_{j,r}}. Defaults to  \code{1}.}
#'   \item{\code{b_delta}}: {\eqn{b_{\delta}}  for the Beta prior on \eqn{\delta_{j,r}}. Defaults to  \code{1}.}
#' }
#' @param alg_type A character string specifying the algorithm type. Need to be either \code{"BMEF-1"} or \code{"BMEF-2"}.
#' Default is \code{"BMEF-1"}.
#'
#' @param threshold A numeric threshold used to truncate \eqn{\tau_r} for rank selection. Defaults to \code{0.05}.
#'
#' @param save_all Logical. If \code{TRUE}, all MCMC samples including burn-in are saved. Defaults to \code{FALSE}.
#'
#' @return A list containing the following components (with \code{K = K_T × K_F}, \code{Rs} denoting the selected rank, and \code{ns} denoting the number of saved MCMC iterations, where \code{ns = n_sample} if \code{save_all = FALSE}, and \code{ns = n_burn + n_sample} otherwise):
#'
#' \item{alpha_s}{An 4D array of dimension \code{ns × K × J × n}, storing samples of the fixed effect coefficients.}
#' \item{gamma_s}{An 3D array of dimension \code{ns × K × n}, storing samples of the subject-level random effects (shared across conditions).}
#' \item{omega_s}{An 4D array of dimension \code{ns × K × J × n}, storing samples of the subject-by-condition-level random effects.}
#' \item{u_s}{A 3D array of dimension \code{ns × Rs × K_T}, storing samples of the factor matrix \code{U}.}
#' \item{v_s}{A 3D array of dimension \code{ns × Rs × K_F}, storing samples of the factor matrix \code{V}.}
#' \item{delta_s}{A 4D array of dimension \code{ns × J × Rs × p}, storing samples of the covariate coefficients.}
#' \item{Sigma2_s}{A list containing posterior samples of variance components:}
#'   \itemize{
#'     \item{\code{Sigma2_epsilon_s}}: {A numeric vector of length \code{ns}, storing samples of  \eqn{\sigma^2_\epsilon}.}
#'     \item{\code{Sigma2_gamma_s}}: {A numeric vector of length \code{ns}, storing samples of  \eqn{\sigma^2_\gamma}.}
#'     \item{\code{Sigma2_omega_s}}: {A matrix of dimension \code{ns × n}, storing samples of  \eqn{\sigma^2_{\omega_i}}.}
#'  }
#' \item{tt}{A numeric vector of time points defining the temporal grid (from the input).}
#' \item{ff}{A numeric vector of frequency points defining the frequency grid (from the input).}
#' \item{Y}{A 3D array of functional responses of dimension \code{n × J × (TT × FF)} (from the input).}
#' \item{X}{A numeric matrix of subject-level covariates of dimension \code{n × p}(from the input).}
#'
#' @export
#' @importFrom splines ns
#' @importFrom MASS mvrnorm ginv
#' @importFrom VGAM dlaplace
#' @importFrom mcompanion null_complement
bmef <- function(Y, X, JJ,  tt, ff, R, K_T, K_F,  n_burn, n_sample, params = NULL, alg_type = "BMEF-1", threshold = 0.05, save_all = FALSE){

  # control the randomness of sampling
  set.seed(123)

  # extract dimensions
  n <- dim(Y)[1]
  J <- dim(Y)[2] # maximum number of conditions
  p <- dim(X)[2]
  K <- K_T * K_F
  TT <- length(tt)
  FF <- length(ff)

  ## if no J_missing, then do the fast algorithm
  idx_no_missing <- 1:n
  if(sum(JJ) == n*J){
    J_missing <- FALSE
  }else{
    J_missing <- TRUE
    idx_missing <- which(apply(JJ, 1, sum) <J)
    idx_no_missing <- setdiff(1:n, idx_missing)
  }

  J_n <- sum(JJ)


  R_num <- 1 # start from the simplest model

  K_T_tilde <- K_T - R_num + 1
  K_F_tilde <- K_F - R_num + 1

  # preprocess the data (assume evenly spaced in the time and frequency domain)
  tmp_prep <- prepocess_new(tt, ff, K_T, K_F)
  O <- tmp_prep$O
  O_tilde <- tmp_prep$O_tilde
  d_ls <- diag(t(O_tilde)%*%O_tilde) # a vector of length K
  C_t <- tmp_prep$C_t
  C_f <- tmp_prep$C_f

  # initialize objects to save samples
  n_total <- n_sample+n_burn # total number of iterations

  # number of saved iterations
  if(save_all){
    n_save <- n_total
  }else{
    n_save <- n_sample
  }

  # --- save posterior samples
  omega_s <- array(NA, dim = c(n_save, K, J, n)) # for each l = 1,..., K, one value per (i, j) combination
  gamma_s <- array(NA, dim = c(n_save, K, n)) # for each l = 1,..., K, one value per i
  alpha_s <- array(NA, dim = c(n_save, K, J, n)) # for each l = 1,..., K, one value per j
  u_s <- array(NA, dim = c(n_save, R, K_T))
  v_s <- array(NA, dim = c(n_save, R, K_F))
  delta_s <- array(NA, dim = c(n_save,J, R, p))
  Sigma2_epsilon_s <- rep(NA, n_save)
  Sigma2_gamma_s <- rep(NA, n_save)
  Sigma2_omega_s <- matrix(NA, n_save, n)
  tau_s <- matrix(0, n_total, R)


  # --- initialization
  tmp_init <- init_params(inits = NULL,  K_T, K_F, J, R, n, p)
  Sigma2_epsilon_current <- tmp_init$Sigma2_epsilon
  Sigma2_gamma_current <- tmp_init$Sigma2_gamma
  Sigma2_omega_current <- tmp_init$Sigma2_omega
  u_current <- tmp_init$u
  v_current <- tmp_init$v
  delta_current <- tmp_init$delta
  pi_current <- 0.5
  tau_current <- rep(NA, R)
  tau_current[1] <- 1

  # --- extract hyper-parameters
  tmp_params <- set_params(params)
  Sigma2_delta <- tmp_params$Sigma2_delta
  a_gamma <- tmp_params$a_gamma # for the Inverse Gamma distribution for the variance of gamma
  b_gamma <- tmp_params$b_gamma  # for the Inverse Gamma distribution for variance of gamma
  a_omega <- tmp_params$a_omega # for the Inverse Gamma distribution for the variance of omega
  b_omega <- tmp_params$b_omega  # for the Inverse Gamma distribution for variance of omega
  a_delta <- tmp_params$a_delta # for the beta distribution for pi
  b_delta <- tmp_params$b_delta # for the beta distribution for pi
  h0 <- tmp_params$h0 # for the mixture Laplace  distribution
  h1 <- tmp_params$h1 # for the mixture Laplace  distribution


  # --- orthogonalize basis evaluations
  Y_tilde <- array(NA, dim = c(n, J, K))  # this is of dimension n, J, K (corresponding to \tilde{y}_{i,j}'s)
  O_tilde_tmp <- diag(1/d_ls) %*% t(O_tilde)

  Y_tilde <- aperm(apply(Y, c(1, 2), function(vec) O_tilde_tmp %*% vec), c(2, 3, 1))


  Z <- generate_Z1(n, J) # an (nJ) \times n matrix

  UV_patch_vec <- matrix(NA, R, K) # need to be updated

  for(r in 1:R){
    UV_patch_vec[r, ] <- kronecker(v_current[r, ], u_current[r, ]) # form the previous iteration
  }

  Y_ss <- sum(apply(Y, c(1, 2), function(x) t(x) %*% x), na.rm = TRUE) # used to fast compute the residuals, removed missing


  ff_theta_eta <- function(Q, g){
    Q_inv <- solve(Q)
    tmp <- mvrnorm(1, mu = Q_inv%*%g, Sigma = Q_inv)
    return(tmp/sqrt(sum(tmp^2)))
  }

  remove_attmpt <- 0  #once have attempted removal, then no more adding attempt
  R_indices <- 1
  R_num <- length(R_indices)

  for(s in 1:n_total){

    if(s%%100 == 0){
      print(paste(s, "-th iteration", sep = ""))
      #print(paste( "Active ranks: ", R_indices, sep = ""))
      #print("Tau_r")
      #print(tau_current)
    }

    if(s%%100 == 0){

      if(length(R_indices)>1){
        rank_remove <- R_indices[which(apply(tau_s[(s-100+1):s, R_indices], 2, mean) < threshold)]
      }else{
        rank_remove <- R_indices[mean(tau_s[(s-100+1):s, R_indices]) < threshold]
      }


      if(length(rank_remove) >= 1){
        remove_attmpt <- 1
        R_indices <- setdiff(R_indices, rank_remove)
        R_num <- length(R_indices)
        K_T_tilde <- K_T - R_num + 1
        K_F_tilde <- K_F - R_num + 1
      }else{
        if(remove_attmpt == 0){
          if(length(R_indices)<R){

            u_current[setdiff(1:R, R_indices),] <- t(null_complement(matrix(t(u_current[R_indices, ]), ncol = R_num), universe = NULL, na.allow = TRUE))[1:length(setdiff(1:R, R_indices)),]
            v_current[setdiff(1:R, R_indices),] <- t(null_complement(matrix(t(v_current[R_indices, ]), ncol = R_num), universe = NULL, na.allow = TRUE))[1:length(setdiff(1:R, R_indices)),]

            tau_current[setdiff(1:R, R_indices)[1]] <- 1
            R_indices <- c(R_indices, setdiff(1:R, R_indices)[1]) # add one rank

            R_num <- length(R_indices)
            K_T_tilde <- K_T - R_num + 1
            K_F_tilde <- K_F - R_num + 1

          }else{
            print("R might be too small!")
          }
        }
      }
    }


    # Block 1
    # Block 1.1: sample the fixed effects
    # Step (a):
    if(alg_type == "BMEF-2"){ # homogeneous variance

      Sigma_inv <- Sigma_inv_tilde <- matrix(0, nrow = J*K, ncol = J*K)
      Sigma2_epsilon_rep <- rep(Sigma2_epsilon_current/d_ls, J)
      Sigma2_epsilon_inv_rep <- rep(Sigma2_epsilon_current/c(t(matrix(d_ls,nrow = K_T,ncol = K_F))), J)


      i_inv_diagonal <- 1/(Sigma2_omega_current[1] + Sigma2_epsilon_rep)
      scaling_fac <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_inv_diagonal))
      Sigma_inv <- diag(i_inv_diagonal) - scaling_fac *  outer(i_inv_diagonal, i_inv_diagonal, "*")

      i_tilde_inv_diagonal <- 1/(Sigma2_omega_current[1] +  Sigma2_epsilon_inv_rep)
      scaling_tilde_fac <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_tilde_inv_diagonal))
      Sigma_inv_tilde <- diag(i_tilde_inv_diagonal)  - scaling_tilde_fac* outer(i_tilde_inv_diagonal,i_tilde_inv_diagonal, "*")


    }else{

      Sigma_inv_i_s <- Sigma_inv_i_tilde_s <- lapply(1:n, function(x) matrix(0, nrow = J*K, ncol = J*K))
      Sigma2_epsilon_rep <- rep(Sigma2_epsilon_current/d_ls, J)
      #Sigma2_epsilon_inv_rep <- rep(Sigma2_epsilon_current/d_ls, each = J)
      Sigma2_epsilon_inv_rep <- rep(Sigma2_epsilon_current/c(t(matrix(d_ls, nrow = K_T,ncol = K_F))), J)

      for(i in 1:n){

        i_inv_diagonal <- 1/(Sigma2_omega_current[i] + Sigma2_epsilon_rep)
        scaling_fac <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_inv_diagonal))
        Sigma_inv_i_s[[i]] <-  diag(i_inv_diagonal) - scaling_fac *  outer(i_inv_diagonal, i_inv_diagonal, "*")

        i_tilde_inv_diagonal <- 1/(Sigma2_omega_current[i] +  Sigma2_epsilon_inv_rep)
        scaling_tilde_fac <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_tilde_inv_diagonal))
        Sigma_inv_i_tilde_s[[i]] <- diag(i_tilde_inv_diagonal)  - scaling_tilde_fac* outer(i_tilde_inv_diagonal,i_tilde_inv_diagonal, "*")

      }
    }

    if(J_missing){ # compute the missing covariates
      Sigma_inv_i_s_missing <- Sigma_inv_i_tilde_s_missing <- lapply(1:length(idx_missing), function(x) matrix(0, nrow = J*K, ncol = J*K))

      for(qq in 1:length(idx_missing)){ # go through all missing individuals

        J_indices <- which(JJ[idx_missing[qq],]!=0) # index for non-missing conditions
        existing_cov_indices <- NULL

        for(kk in J_indices){
          existing_cov_indices  <- c(existing_cov_indices, ((kk-1)*K + 1): (kk*K))
        }

        Sigma2_epsilon_rep_tmp <- rep(Sigma2_epsilon_current/d_ls, length(J_indices))

        i_inv_diagonal_tmp <- 1/(Sigma2_omega_current[idx_missing[qq]] + Sigma2_epsilon_rep_tmp)
        scaling_fac_tmp <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_inv_diagonal_tmp))
        Sigma_inv_tmp <- diag(i_inv_diagonal_tmp) - scaling_fac_tmp *  outer(i_inv_diagonal_tmp, i_inv_diagonal_tmp, "*")
        Sigma_inv_i_s_missing[[qq]][existing_cov_indices, existing_cov_indices] <- Sigma_inv_tmp

        Sigma2_epsilon_inv_rep_tmp <- rep(Sigma2_epsilon_current/c(t(matrix(d_ls,nrow = K_T,ncol = K_F))), length(J_indices))
        i_tilde_inv_diagonal_tmp <- 1/(Sigma2_omega_current[idx_missing[qq]] +  Sigma2_epsilon_inv_rep_tmp)
        scaling_tilde_fac_tmp <- Sigma2_gamma_current/(1 + Sigma2_gamma_current* sum(i_tilde_inv_diagonal_tmp))
        Sigma_inv_i_tilde_tmp <- diag(i_tilde_inv_diagonal_tmp)  - scaling_tilde_fac_tmp* outer(i_tilde_inv_diagonal_tmp,i_tilde_inv_diagonal_tmp, "*")
        Sigma_inv_i_tilde_s_missing[[qq]][existing_cov_indices, existing_cov_indices] <- Sigma_inv_i_tilde_tmp

      }
    }


    for(r in R_indices){

      if((s == 1) & (r == 1)){
        q_r_j_s <- lapply(1:n, function(x) matrix(0, nrow = J, ncol = K)) # computed for this specific r, initialization
      }

      R_indices_neg <- setdiff(R_indices, r)

      for(i in 1:n){
        for(j in 1:J){
          if(R_num > 2){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] -  (t(X[i,]) %*% t(delta_current[j,R_indices_neg, ])) %*% UV_patch_vec[R_indices_neg, ] # n by k
          }
          if(R_num == 2){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] -  (t(X[i,]) %*% delta_current[j,R_indices_neg, ]) %*% UV_patch_vec[R_indices_neg, ] # n by k
          }
          if(R_num == 1){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,]
          }
        }
      }

      if(R_num == 1){
        B_r <- diag(K_T)
      }else{
        B_r <- null_complement(matrix(t(u_current[R_indices_neg, ]), ncol = R_num-1), universe = NULL, na.allow = TRUE) # need to double check for larger R, can be replaced for Gram-Schmidt later
      }

      lambda_tilde_r <- X%*% t(delta_current[ ,r, ]) # dim = c(n, J)
      v_tilde_r <- t(apply(lambda_tilde_r, 1, function(row) kronecker(row, v_current[r, ])))

      H_r <- array(0, c(n, K_T, K_T))
      w_r <- array(0, c(n, K_T))

      v_outer <- array(0, c(n, J*K_F, J*K_F))  # Store outer products of v_tilde_r

      for (l1 in 1:(J*K_F)){
        for (l2 in 1:(J*K_F)) {
          v_outer[, l1, l2] <- v_tilde_r[, l1] * v_tilde_r[, l2]
        }
      }

      for(l2 in 1:(J*K_F)){
        K_F_idx <- (l2 - 1) %% K_F + 1
        J_idx <- (l2 - 1) %/% K_F + 1
        q_indices <-  ((K_F_idx-1)*K_T+1): (K_F_idx*K_T)

        for(l1 in 1:(J*K_F)){
          block_row <- l1; block_col = l2; block_size = K_T
          row_indices <- ((block_row - 1) * block_size + 1):(block_row * block_size)
          col_indices <- ((block_col - 1) * block_size + 1):(block_col * block_size)

          if(alg_type == "BMEF-2"){
            inv_block_mat  <-  Sigma_inv[row_indices, col_indices]
            H_r[idx_no_missing,,] <- H_r[idx_no_missing,,] + array(t(sapply(v_outer[idx_no_missing, l1, l2], FUN = function(a){inv_block_mat*a})), dim = c(length(idx_no_missing), K_T, K_T))
            if(J_missing){
              for(qq in 1:length(idx_missing)){
                if(JJ[idx_missing[qq],J_idx]!=0){ # this condition exist
                  inv_block_mat <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
                  H_r[idx_missing[qq],,] <- H_r[idx_missing[qq],,] + v_outer[idx_missing[qq], l1, l2]*inv_block_mat
                }
              }
            }
            for(i in 1:n){
              if(i %in% idx_no_missing){
                w_r[i, ] <- w_r[i, ] +  v_tilde_r[i, l1] * (inv_block_mat %*% q_r_j_s[[i]][J_idx, q_indices])
              }else{
                if(JJ[i,J_idx]!=0){ # this condition exist
                  qq <- which(idx_missing == i)
                  inv_block_mat <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
                  w_r[i, ] <- w_r[i, ] +  v_tilde_r[i, l1] * (inv_block_mat %*% q_r_j_s[[i]][J_idx, q_indices])
                }
              }
            }
          }else{
            for(i in 1:n){
              if(i %in% idx_no_missing){
                inv_block_mat  <-  Sigma_inv_i_s[[i]][row_indices, col_indices]
                H_r[i,,] <- H_r[i,,] + v_outer[i, l1, l2] * inv_block_mat
                w_r[i, ] <- w_r[i, ] +  v_tilde_r[i, l1] * (inv_block_mat %*% q_r_j_s[[i]][J_idx, q_indices])
              }else{
                if(JJ[i,J_idx]!=0){ # this condition exist
                  qq <- which(idx_missing == i)
                  inv_block_mat <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
                  H_r[i,,] <- H_r[i,,] + v_outer[i, l1, l2] * inv_block_mat
                  w_r[i, ] <- w_r[i, ] +  v_tilde_r[i, l1] * (inv_block_mat %*% q_r_j_s[[i]][J_idx, q_indices])
                }
              }
            }
          }
        }
      }



      #Q_u  <- diag(rep(1/Sigma2_theta, K_T_tilde)) +
      Q_u  <-  t(B_r)%*%apply(H_r, c(2,3), sum)%*%B_r
      g_u  <- t(B_r) %*% apply(w_r, 2, sum)
      theta_r_current <- tryCatch({c(rFisherBingham(1, mu = g_u, Aplus = -Q_u/2, mtop = 10000))}, error = function(e){ff_theta_eta(Q_u, g_u)})
      u_current[r,] <- B_r%*%theta_r_current


      # Step a.2
      if(R_num == 1){
        S_r <- diag(K_F)
      }else{
        S_r <- null_complement(matrix(t(v_current[R_indices_neg, ]), ncol = R_num-1), universe = NULL, na.allow = TRUE) # can be replaced for Gram-Schmidt later
      }

      u_tilde_r <- t(apply(lambda_tilde_r, 1, function(row) kronecker(row, u_current[r, ]))) # dim = c(n, J*K_T)

      L_r <- array(0, c(n, K_F, K_F))
      h_r <- array(0, c(n, K_F))

      u_outer <- array(0, c(n, J*K_T, J*K_T))  # Store outer products of v_tilde_r

      for (l1 in 1:(J*K_T)){
        for (l2 in 1:(J*K_T)) {
          u_outer[, l1, l2] <- u_tilde_r[, l1] * u_tilde_r[, l2]
        }
      }

      for(l2 in 1:(J*K_T)){
        K_T_idx <- (l2 - 1) %% K_T + 1
        J_idx <- (l2 - 1) %/% K_T + 1
        q_indices <-  ((K_T_idx-1)*K_F+1): (K_T_idx*K_F)

        for(l1 in 1:(J*K_T)){
          block_row <- l1; block_col = l2; block_size = K_F
          row_indices <- ((block_row - 1) * block_size + 1):(block_row * block_size)
          col_indices <- ((block_col - 1) * block_size + 1):(block_col * block_size)

          if(alg_type == "BMEF-2"){
            inv_block_mat  <-  Sigma_inv_tilde[row_indices, col_indices]
            L_r[idx_no_missing,,] <- L_r[idx_no_missing,,] + array(t(sapply(u_outer[idx_no_missing, l1, l2], FUN = function(a){inv_block_mat*a})), dim = c(length(idx_no_missing), K_F, K_F))

            if(J_missing){
              for(qq in 1:length(idx_missing)){
                if(JJ[idx_missing[qq],J_idx]!=0){ # this condition exist
                  inv_block_mat <- Sigma_inv_i_tilde_s_missing[[qq]][row_indices, col_indices]
                  L_r[idx_missing[qq],,] <- L_r[idx_missing[qq],,] + u_outer[idx_missing[qq], l1, l2]*inv_block_mat
                }
              }
            }
            for(i in 1:n){
              if(i %in% idx_no_missing){
                h_r[i, ] <- h_r[i, ] +  u_tilde_r[i, l1] * (inv_block_mat %*% c(t(matrix(q_r_j_s[[i]][J_idx, ], nrow = K_T, ncol = K_F)))[q_indices])
              }else{
                if(JJ[i,J_idx]!=0){ # this condition exist
                  qq <- which(idx_missing == i)
                  inv_block_mat <- Sigma_inv_i_tilde_s_missing[[qq]][row_indices, col_indices]
                  h_r[i, ] <- h_r[i, ] +  u_tilde_r[i, l1] * (inv_block_mat %*% c(t(matrix(q_r_j_s[[i]][J_idx, ], nrow = K_T, ncol = K_F)))[q_indices])
                }
              }
            }
          }else{

            for(i in 1:n){
              if(i %in% idx_no_missing){
                inv_block_mat  <-  Sigma_inv_i_tilde_s[[i]][row_indices, col_indices]
                L_r[i,,] <- L_r[i,,] + u_outer[i, l1, l2] * inv_block_mat
                h_r[i, ] <- h_r[i, ] +  u_tilde_r[i, l1] * (inv_block_mat %*% c(t(matrix(q_r_j_s[[i]][J_idx, ], nrow = K_T, ncol = K_F)))[q_indices])
              }else{
                if(JJ[i,J_idx]!=0){ # this condition exist
                  qq <- which(idx_missing == i)
                  inv_block_mat <- Sigma_inv_i_tilde_s_missing[[qq]][row_indices, col_indices]
                  L_r[i,,] <- L_r[i,,] + u_outer[i, l1, l2] * inv_block_mat
                  h_r[i, ] <- h_r[i, ] +  u_tilde_r[i, l1] * (inv_block_mat %*% c(t(matrix(q_r_j_s[[i]][J_idx, ], nrow = K_T, ncol = K_F)))[q_indices])
                }
              }
            }
          }
        }
      }

      Q_v  <-  t(S_r)%*%apply(L_r, c(2,3), sum)%*%S_r
      g_v  <- t(S_r) %*% apply(h_r, 2, sum)
      eta_r_current <- tryCatch({c(rFisherBingham(1, mu = g_v, Aplus = -Q_v/2, mtop = 10000))}, error = function(e){ff_theta_eta(Q_v, g_v)})

      v_current[r,] <- S_r%*%eta_r_current
      UV_patch_vec[r, ] <- kronecker(v_current[r, ], u_current[r, ]) # does affect the partial residuals
    }


    # Step (b)
    for(r in R_indices){

      R_indices_neg <- setdiff(R_indices, r)

      for(i in 1:n){
        for(j in 1:J){
          if(R_num > 2){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] -  (t(X[i,]) %*% t(delta_current[j,R_indices_neg, ])) %*% UV_patch_vec[R_indices_neg, ] # n by k
          }
          if(R_num == 2){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,] -  (t(X[i,]) %*% delta_current[j,R_indices_neg, ]) %*% UV_patch_vec[R_indices_neg, ] # n by k
          }
          if(R_num == 1){
            q_r_j_s[[i]][j,] <- Y_tilde[i,j,]
          }
        }
      }

      Delta_r_s <- array(NA, c(n, J, p, p)) # for a specific R
      Chi_r_s <- array(NA, c(n, J, p))


      for(i in 1:n){

        for(j in 1:J){
          block_row <- j; block_col = j; block_size = K
          row_indices <- col_indices <- ((block_row - 1) * block_size + 1):(block_row * block_size)

          if(alg_type == "BMEF-2"){
            if(i %in% idx_no_missing){
              inv_block_mat  <- Sigma_inv[row_indices, col_indices]
            }else{
              qq <- which(idx_missing == i)
              inv_block_mat  <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
            }
          }else{
            if(i %in% idx_no_missing){
              inv_block_mat  <- Sigma_inv_i_s[[i]][row_indices, col_indices]
            }else{
              qq <- which(idx_missing == i)
              inv_block_mat  <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
            }
          }

          if(p == 1){
            Delta_r_s[i,j,,] <- X[i,]^2 * (t(UV_patch_vec[r,])%*%inv_block_mat%*% UV_patch_vec[r,])
          }else{
            Delta_r_s[i,j,,] <- (X[i,]%*%t(X[i,])) * c(t(UV_patch_vec[r,])%*%inv_block_mat%*% UV_patch_vec[r,])
          }

          tmp_sum <- rep(0, K)
          for(jj in setdiff(1:J, j)){
            block_row <- j; block_col = jj; block_size = K
            col_indices <- ((block_col - 1) * block_size + 1):(block_col * block_size)
            if(alg_type == "BMEF-2"){
              if(i %in% idx_no_missing){
                inv_block_mat_tmp  <- Sigma_inv[row_indices, col_indices]
              }else{
                qq <- which(idx_missing == i)
                inv_block_mat_tmp  <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
              }
            }else{
              if(i %in% idx_no_missing){
                inv_block_mat_tmp  <- Sigma_inv_i_s[[i]][row_indices, col_indices]
              }else{
                qq <- which(idx_missing == i)
                inv_block_mat_tmp  <- Sigma_inv_i_s_missing[[qq]][row_indices, col_indices]
              }
            }
            if(JJ[i,jj]==1){
              tmp_sum <- tmp_sum + inv_block_mat_tmp %*% (q_r_j_s[[i]][jj, ] - c(t(X[i,]) %*% delta_current[jj,r, ])* UV_patch_vec[r,])
            }
          }
          Chi_r_s[i,j, ] <- X[i,] * c(t(UV_patch_vec[r,]) %*% (inv_block_mat %*% q_r_j_s[[i]][j, ] + tmp_sum))
        }
      }

      for(j in 1:J){
        if(p==1){ # only has one covariate
          if(tau_current[r] == 0){
            delta_current[j,r,] <- 0
          }else{
            Q_delta <- (1/tau_current[r])*(1/Sigma2_delta) + sum(Delta_r_s[,j,,])
            g_delta <- sum(Chi_r_s[,j, ])
            delta_current[j,r,] <- rnorm(1, mean = ginv(Q_delta)%*%g_delta, sd = ginv(Q_delta))
          }
        }else{
          if(tau_current[r] == 0){
            delta_current[j,r,] <- rep(0,p)
          }else{
            Q_delta <-  (1/tau_current[r])*diag(1/Sigma2_delta) + apply(Delta_r_s[,j,,], c(2,3), sum)
            g_delta <- apply(Chi_r_s[which(JJ[,j]==1),j, ], 2, sum)
            delta_current[j,r,] <- mvrnorm(1, mu = ginv(Q_delta)%*%g_delta, Sigma = ginv(Q_delta))
          }
        }
      }
    }

    delta_sigma_inv_delta <- matrix(0, J, R)

    if(p >1){
      for(jj in 1:J){
        for(rr in R_indices){
          delta_sigma_inv_delta[jj,rr] <- t(delta_current[jj,rr,])%*%diag(1/Sigma2_delta)%*%delta_current[jj,rr,]
        }
      }
    }else{
      for(jj in 1:J){
        for(rr in R_indices){
          delta_sigma_inv_delta[jj,rr] <- delta_current[jj,rr,]^2*(1/Sigma2_delta)
        }
      }
    }

    b_current <- rep(0, R)
    for(r in R_indices){

      p_laplace <- -p*J/2+ 1
      b_laplace <-  sum(delta_sigma_inv_delta[,r])

      if(b_current[r]==1){
        a_laplace <- 2/h0
      }else{
        a_laplace <- 2/h1
      }

      tau_current[r] <- (p_laplace-1 + sqrt((p_laplace-1)^2 + a_laplace*b_laplace))/a_laplace
      p_bernouli_tmp_1 <- dlaplace(tau_current[r], 0, h0)*pi_current
      p_bernouli_tmp_2 <- dlaplace(tau_current[r], 0, h1)*(1-pi_current)

      p_bernouli <- p_bernouli_tmp_1/(p_bernouli_tmp_1+ p_bernouli_tmp_2)

      b_current[r] <- rbinom(1, 1, p_bernouli)
    }

    tau_s[s,] <- tau_current

    pi_current <- rbeta(1, a_delta + sum(b_current), b_delta + R_num  - sum(b_current))


    alpha_current <- array(0, dim = c(K, J, n))

    for(r in R_indices){
      X_delta <- X %*% t(delta_current[, r, ])  # Resulting in an n x J matrix
      for (k in 1:K) {
        alpha_current[k, , ] <- alpha_current[k, , ] + t(X_delta)*UV_patch_vec[r, k]
      }
    }

    gamma_current <- matrix(NA, K, n)
    for(l in 1:K){
      Sigma_1_vec <- rep(Sigma2_omega_current, each = J) + Sigma2_epsilon_current/d_ls[l]
      qq <- ((c(t(Y_tilde[,,l])) - c(alpha_current[l,,]))/Sigma_1_vec)
      idx_tmp <- which(!is.na(qq))
      g_gamma <- t(Z)[, idx_tmp]%*%qq[idx_tmp]
      Q_gamma_inv_vec <- 1/(1/Sigma2_gamma_current +   diag(t(Z[idx_tmp,]) %*%diag(1/Sigma_1_vec[idx_tmp]) %*%Z[idx_tmp,]))
      gamma_current[l, ] <- rnorm(n, mean = Q_gamma_inv_vec * g_gamma, sd = sqrt(Q_gamma_inv_vec))
    }


    D_omega_vec <- rep(Sigma2_omega_current, each = J)
    omega_current <- array(NA, dim = c(K, J, n))

    for(l in 1:K){
      g_omega <- (d_ls[l]/Sigma2_epsilon_current)* ((c(t(Y_tilde[,,l])) - c(alpha_current[l,,]) -  Z%*%gamma_current[l,]))
      idx_tmp <- which(is.na(g_omega))
      if(length(idx_tmp)!=0){
        g_omega[idx_tmp] <- 0
      }
      Q_omega_inv_vec <- 1/(1/D_omega_vec + d_ls[l]/Sigma2_epsilon_current)
      omagea_tmp <- rnorm(n * J, mean = Q_omega_inv_vec * g_omega, sd = sqrt(Q_omega_inv_vec))
      if(length(idx_tmp)!=0){
        omagea_tmp[idx_tmp] <- NA
      }
      omega_current[l,,] <- array(omagea_tmp, dim = c(J, n))
    }


    res_tmp <- matrix(NA, nrow = n, ncol = J) # save the sum of squared residuals for i, j combination (except for the 1st ter, y^Ty)
    for(i in 1:n){
      for(j in 1:J){
        beta_tmp <-alpha_current[,j,i] + gamma_current[,i] + omega_current[,j, i]
        beta_tmp_star <- beta_tmp*d_ls
        res_tmp[i,j] <- -2* t(beta_tmp_star)%*%Y_tilde[i,j,] + t(beta_tmp_star)%*%beta_tmp
      }
    }


    Sigma2_epsilon_current <- 1/rgamma(1, shape = (J_n*TT*FF)/2, rate = (sum(res_tmp, na.rm = TRUE) + Y_ss)/2)

    Sigma2_gamma_current <- 1/rgamma(1, a_gamma + n*K/2,  rate = b_gamma + sum(gamma_current^2)/2)

    if(alg_type == "BMEF-2"){
      Sigma2_omega_current <- rep(1 / rgamma(1, shape =  a_omega + J_n* K / 2, rate =  b_omega + sum(omega_current^2, na.rm = TRUE) / 2), n)
    }else{
      Sigma2_omega_current <- 1 / rgamma(n, shape =  a_omega + apply(JJ, 1, sum) * K / 2, rate =  b_omega + colSums(omega_current^2, dims = 2, na.rm = TRUE) / 2)
    }

    # save the posterior samples
    if(save_all){
      Sigma2_epsilon_s[s]  <- Sigma2_epsilon_current
      Sigma2_gamma_s[s] <- Sigma2_gamma_current
      Sigma2_omega_s[s,] <-Sigma2_omega_current
      v_s[s,,] <- v_current
      u_s[s,,] <- u_current
      alpha_s[s,,,] <- alpha_current
      delta_s[s,,, ] <- delta_current
      gamma_s[s,,] <- gamma_current
      omega_s[s,,,]<- omega_current
      sample_idx <- (n_burn+1):n_save

    }else{
      if(s > n_burn){
        Sigma2_epsilon_s[s - n_burn]  <- Sigma2_epsilon_current
        Sigma2_gamma_s[s - n_burn] <- Sigma2_gamma_current
        Sigma2_omega_s[s - n_burn,] <-Sigma2_omega_current
        v_s[s - n_burn,,] <- v_current
        u_s[s - n_burn,,] <- u_current
        alpha_s[s - n_burn,,,] <- alpha_current
        delta_s[s - n_burn,,, ] <- delta_current
        gamma_s[s - n_burn,,] <- gamma_current
        omega_s[s - n_burn,,,] <- omega_current
        sample_idx <- 1:n_save
      }
    }
  }

  Sigma2_s <- list(Sigma2_epsilon_s = Sigma2_epsilon_s, Sigma2_gamma_s = Sigma2_gamma_s,  Sigma2_omega_s  =  Sigma2_omega_s)
  dat2return <- list(Y = Y, X = X, tt = tt, ff = ff, alpha_s = alpha_s, omega_s = omega_s, gamma_s = gamma_s,  v_s = v_s[,R_indices,], u_s = u_s[,R_indices,], delta_s = delta_s[,,R_indices,], Sigma2_s = Sigma2_s)

  return(dat2return)
}


#' Posterior Inference for the Decomposition of Fixed Effects from a Fitted BMEF Model
#'
#' This function generates posterior summaries of two-way functional fixed effects using a fitted
#' Bayesian Mixed-Effects Model for Multilevel Two-Way Functional Data (BMEF). It returns inference results for the estimated "base" patterns,
#' covariate-dependent weights, and marginal "principal" functions over the time and frequency domains.
#'
#' @param bmef_obj A fitted BMEF object returned by the function \code{bmef()}.
#' @param cred_level A numeric value specifying the credible level. Defaults to 0.95.
#'
#' @return A list containing the following components (with \code{TT} and \code{FF} denoting the number of time and frequency grid points used to fit BMEF, and \code{Rs} the selected CP rank of the fixed effects):
#' \item{base_patterns}{A list with elements \code{m}, \code{l}, and \code{u}, representing the posterior mean, lower bound, and upper bound of the credible intervals. Each element is a 3D array of size   \code{Rs} × \code{TT} × \code{FF} summarizing the two-way "base" patterns.}
#' \item{marginal_t}{A list with elements \code{m}, \code{l}, and \code{u}, representing the posterior mean, lower bound, and upper bound of the credible intervals. Each element is a matrix of size  \code{Rs} × \code{TT} , summarizing the marginal "principal" functions over time.}
#' \item{marginal_f}{A list with elements \code{m}, \code{l}, and \code{u}, representing the posterior mean, lower bound, and upper bound of the credible intervals. Each element is a matrix of size \code{Rs} × \code{FF}, summarizing the marginal "principal" functions over frequency.}
#' \item{weights_vals}{A list with elements \code{m}, \code{l}, and \code{u},  representing the posterior mean, lower bound, and upper bound of the credible intervals. Each element is a 3D array of size   \code{Rs} ×  \code{n} × \code{J}, summarizing the covariate-dependent weights for each subject ×  condition ×  rank combination.}
#'
#' @export


inference.decompose.bmef <- function(bmef_obj, cred_level){


  Rs <- dim(bmef_obj$v_s)[2]
  K_T <- dim(bmef_obj$u_s)[3]
  K_F <- dim(bmef_obj$v_s)[3]
  sample_idx <-  1:dim(bmef_obj$u_s)[1]
  K <- K_T*K_F
  J <- dim(bmef_obj$delta_s)[2]
  X <- bmef_obj$X
  tt <- bmef_obj$tt
  ff <- bmef_obj$ff
  tmp_prep <- prepocess_new(tt, ff, K_T, K_F)
  O_tilde <- tmp_prep$O_tilde
  d_ls <- diag(t(O_tilde)%*%O_tilde) # a vector of length K
  C_t <- tmp_prep$C_t
  C_f <- tmp_prep$C_f
  bsMat_tt <-  tmp_prep$bsMat_tt
  bsMat_ff <-  tmp_prep$bsMat_ff


  base_patterns_samples <- array(NA,  dim = c(length(sample_idx), Rs, TT, FF))

  marginal_t_samples <- array(NA,  dim = c(length(sample_idx), Rs, TT))
  marginal_f_samples <- array(NA,  dim = c(length(sample_idx), Rs, FF))

  weights_samples <- array(NA, dim = c(length(sample_idx), Rs, n, J))
  UV_patch_vec_s <- array(NA, dim = c(length(sample_idx), Rs, K))

  for(ss in 1:length(sample_idx)){
    for(rr in  1:Rs){
      UV_patch_vec_s[ss ,rr,] <- kronecker(bmef_obj$v_s[ss ,rr,],bmef_obj$u_s[ss, rr,])
      base_patterns_samples[ss, rr, ,] <- O_tilde%*%UV_patch_vec_s[ss ,rr,]
      marginal_t_samples[ss, rr, ] <- bsMat_tt%*%C_t%*%bmef_obj$u_s[ss,rr,]
      marginal_f_samples[ss, rr, ] <- bsMat_ff%*%C_f%*%bmef_obj$v_s[ss,rr,]

      for(jj in 1:J){
        weights_samples[ss, rr, , jj] <- X %*%bmef_obj$delta_s[ss,jj,rr, ]
      }
    }
  }

  base_patterns <- list()
  base_patterns$m <- apply(base_patterns_samples, c(2,3,4), mean)
  base_patterns$l <- apply(base_patterns_samples, c(2,3,4), function(u){quantile(u, (1-cred_level)/2)})
  base_patterns$u <- apply(base_patterns_samples, c(2,3,4), function(u){quantile(u, 1-(1-cred_level)/2)})

  weights_vals <- list()
  weights_vals$m <- apply(weights_samples, c(2,3,4), mean)
  weights_vals$l <- apply(weights_samples, c(2,3,4), function(u){quantile(u, (1-cred_level)/2)})
  weights_vals$u <- apply(weights_samples, c(2,3,4), function(u){quantile(u, 1-(1-cred_level)/2)})


  marginal_t <- list()
  marginal_t$m <- apply(marginal_t_samples, c(2,3), mean)
  marginal_t$l <- apply(marginal_t_samples, c(2,3), function(u){quantile(u, (1-cred_level)/2)})
  marginal_t$u <- apply(marginal_t_samples, c(2,3), function(u){quantile(u, 1-(1-cred_level)/2)})

  marginal_f <- list()
  marginal_f$m <- apply(marginal_f_samples, c(2,3), mean)
  marginal_f$l <- apply(marginal_f_samples, c(2,3), function(u){quantile(u, (1-cred_level)/2)})
  marginal_f$u <- apply(marginal_f_samples, c(2,3), function(u){quantile(u, 1-(1-cred_level)/2)})

  return(list(base_patterns = base_patterns, weights_vals = weights_vals,
              marginal_t = marginal_t, marginal_f = marginal_f))

}



#' Posterior Inference for the two-way functional effects from a Fitted BMEF Model
#'
#' This function generates posterior summaries of two-way functional fixed and random effects
#' using a fitted Bayesian Mixed-Effects Model for Multilevel Two-Way Functional Data (BMEF).
#'
#' @param bmef_obj A fitted BMEF object returned by the function \code{bmef()}.
#'
#' @return A list containing posterior predictive samples for the estimated two-way functional fixed and random effects using the posterior mean
#' (with \code{J} denoting the number of conditions):
#' \itemize{
#'   \item{\code{A_s}}{A 5D array of dimension \code{n × J × TT_pred × FF_pred}, storing posterior mean of the two-way functional fixed effects.}
#'   \item{\code{B_s}}{A 4D array of dimension \code{n × TT_pred × FF_pred}, storing posterior mean of subject-level two-way functional  random effects shared across conditions.}
#'   \item{\code{C_s}}{A 5D array of dimension \code{n × J × TT_pred × FF_pred}, storing posterior mean of subject-by-condition-level two-way functional random effects.}
#' }
#' @export
#'


inference.mixed.bmef <- function(bmef_obj){

  n <- dim(bmef_obj$Y)[1]
  Rs <- dim(bmef_obj$v_s)[2]
  K_T <- dim(bmef_obj$u_s)[3]
  K_F <- dim(bmef_obj$v_s)[3]
  sample_idx <-  1:dim(bmef_obj$u_s)[1]
  K <- K_T*K_F
  J <- dim(bmef_obj$delta_s)[2]
  X <- bmef_obj$X
  tt <- bmef_obj$tt
  ff <- bmef_obj$ff
  TT <- length(tt)
  FF <- length(ff)
  tmp_prep <- prepocess_new(tt, ff, K_T, K_F)
  O_tilde <- tmp_prep$O_tilde # original O_tilde
  C_t <- tmp_prep$C_t
  C_f <- tmp_prep$C_f
  bsMat_tt <-  tmp_prep$bsMat_tt
  bsMat_ff <-  tmp_prep$bsMat_ff

  A_s <- array(NA, dim = c(n, J, TT, FF))
  B_s <- array(NA,dim = c(n, TT, FF))
  C_s <- array(NA,dim = c(n,J, TT, FF))

  for(ii in 1:n){
    B_s[ii,,] <- matrix(apply(O_tilde%*%t(bmef_obj$gamma_s[sample_idx, , ii]),1, mean), TT, FF)
  }

  for(ii in 1:n){
    for(jj in 1:J){
      A_s[ii,jj,, ] <-  matrix( apply( O_tilde%*% t(bmef_obj$alpha_s[sample_idx, ,jj, ii]),1, mean), TT, FF)
      C_s[ii,jj,, ] <- matrix(apply( O_tilde%*%  t(bmef_obj$omega_s[sample_idx, ,jj, ii]),1, mean), TT, FF)
    }
  }

  return(list(A_s = A_s, B_s = B_s, C_s = C_s))
}



