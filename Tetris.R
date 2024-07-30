library(MASS)
library(statmod)
library(matlab)
library(combinat)
library(R.utils)
library(matrixcalc)
library(clue)
library(nloptr)

# Sets parameters and prior hyperparameters for posterior sampling
new_control <- function(nrun = 10000, burn = 8000, thin = 1,
                        nu = 3, 
                        a1 = 2.1, b1 = 1,
                        a2 = 3.1, b2 = 1,
                        apsi = 1, bpsi = 0.3,
                        alpha = 10, beta = 0.5)
{
  return(list(nrun = nrun, burn = burn, thin = thin,
              nu = nu, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
              apsi = apsi, bpsi = bpsi, alpha = alpha, beta = beta))
}

# Gibbs sampler for Tetris
# Inputs: X_s = list of S matrices in samples x features format
#         alpha, beta = IBP hyperparameters
#         trace = TRUE if iterations should be printed
#         nprint = how often to print iterations if trace is TRUE
#         initial = optional list of initial parameters 
#         fixed = FALSE if the factor indicator matrix A should be sampled 
#         A_fixed = what A should be fixed to, if fixed is TRUE
# Outputs: out = list of post-burn-in chains for the loadings matrix Lambda, the errors Psi, the latent factors l_s
#                for each study s, and the factor indicator matrix A
tetris <- function(X_s,  alpha, beta, trace = TRUE, nprint = 1000, initial = 0, fixed = FALSE, A_fixed = 0,
                         control = list(...), ...)
{
  ## Read in data
  S <- length(X_s)                      # number of studies
  P <- dim(X_s[[1]])[2]                 # number of variables
  Y_s <- lapply(X_s, scale,scale=FALSE) # centered data
  n_s <- sapply(X_s, nrow)              # sample sizes
  
  ## Set hyperparameters (see new_control() documentation)
  control <- new_control()
  nrun <- control$nrun
  thin <- control$thin
  burn <- control$burn
  sp  <- (nrun - burn) / thin
  apsi <- control$apsi 
  bpsi<- control$bpsi 
  nu <- control$nu
  a1 <- control$a1
  b1 <- control$b1 
  a2 <- control$a2
  b2 <-  control$b2 
  
  ## Initialize model parameters 
  # If no initial values provided, draw Lambda, A, Psi, and l_s from prior distributions;
  # otherwise, use the initial values that were given
  if (!is.list(initial)) {
    # Inverse error variances
    psi_s <- lapply(1:S, function(x) rgamma(P, shape = apsi, scale = 1 / bpsi))
    # Inverse error matrix
    psi_inv_s <- lapply(1:S, function(x) diag(psi_s[[x]]))
    # To hold the factor indicator matrix A
    A <- array(0, dim=c(S,S*S*alpha))
    # Draw the number of factors 
    num_factors <- rpois(1, alpha)
    while (num_factors==0) {
      num_factors <- rpois(1,alpha)
    }
    # Set the first study to have each of the instantiated factors
    if (num_factors > 0) {
      A[1,1:num_factors] <- 1
    }
    # For the remaining studies...
    for (i in 2:S) {
      # Assign each of the existing factors with a scaled probability
      if (num_factors > 0) {
        probs <- colSums(A)/(beta+i-1)
        realization <- rbinom(n=num_factors,size=1,prob=probs[1:num_factors])
      } else {
        realization <- c()
      }
      # Draw the number of new factors to add
      new_factors <- rpois(1,alpha*beta/(beta+i-1))
      while (new_factors==0) {
        new_factors <- rpois(1,alpha*beta/(beta+i-1))
      }
      num_factors <- num_factors + new_factors
      # Assign each of these new factors to the study being considered
      if (num_factors > 0) {
        A[i,1:num_factors] <- c(realization,rep(1,new_factors))
      }
    }
    A <- matrix(A[1:S,1:num_factors],S,num_factors)
    colnames(A) <- NULL
    # Update the number of factors
    K <- num_factors
    # If we are fixing A, establish those parameters
    if (fixed) {
      A <- A_fixed
      K <- dim(A)[2]
    }
  } else {
    # If we are starting with initial parameters, fill them in here
    Lambda <- initial[[1]]
    A <- initial[[2]]
    K <- dim(A)[2]
    psi_s <- initial[[3]]
    psi_inv_s <- lapply(1:S, function(x) diag(c(psi_s[[x]])))
    l_s <- initial[[4]]
  }
  # Draw omega, delta, and tau from prior distributions
  omega <- matrix(rgamma(P * K, shape = nu/2, scale = 2/nu),P,K)
  delta <- rgamma(K, shape = c(a1, rep(a2, K-1)),
                  scale = c(1/b1, rep(1/b2, K-1)))
  tau <- cumprod(delta)
  D <- matvec(omega, tau)
  # If no initialization, draw Lambda and l_s from priors
  if (!is.list(initial)) {
    Lambda <- t(sapply(1:P, function(x) sapply(1:K, function(y) rnorm(1,0,1/(omega[x,y]*tau[y])))))
    Lambda <- matrix(Lambda,P,K)
    l_s <- lapply(1:S, function(x) mvrnorm(n_s[x], rep(0, K), diag(rep(1, K))))
  }
  
  ## Initialize output
  Lambdaout <- Aout <- list()
  Psiout <- l_out <- lapply(1:S, function(x) list())

  ## Start posterior sampling
  for(iter in 1:nrun)
  {
    # Only update factor indicator matrix if not fixed
    if (!fixed) {
      # Step 1: factor indicator matrix (A)
      if (colSums(A)[1]!=0) {
        # For each study...
        for(s in 1:S) {
  	        # For each factor...
            for(k in 1:K) {
              curr <- A[s,k]
              # Ak1 = A where study s has factor k
              Ak1 <- A
              Ak1[s,k] <- 1
              # Ak0 = A where study s doesn't have factor k
              Ak0 <- A
              Ak0[s,k] <- 0
              # Hold terms for the likelihood ratios
              r_l_term1 <- list()
              # For each sample in a given study...
              for(n in 1:n_s[s]) {
        	      if (K>2) {
                  # Compute the marginalized likelihood term
                  lbar <- (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)^(-1)*
                        (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%
                           (t(Y_s[[s]][n,,drop=F])-Lambda[,-k,drop=F]%*%diag(Ak1[s,-k])%*%t(l_s[[s]][n,-k,drop=F])))
        	      }
        	      if (K==1) {
        		      lbar <- (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)^(-1)*
                        (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%
                           (t(Y_s[[s]][n,,drop=F])-0))
                }
        	      if (K==2) {
        		      lbar <- (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)^(-1)*
                        (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%
                           (t(Y_s[[s]][n,,drop=F])-Lambda[,-k,drop=F]%*%diag(Ak1[s,-k,drop=F])%*%t(l_s[[s]][n,-k,drop=F])))
        	      }
                r_l_term1[[n]] <- log((t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)^(-0.5))+
                        0.5*(t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)*lbar^2
              }
              # Sum the terms to get log-likelihood of each case across samples
              r_l_term1 <- Reduce('+',r_l_term1)
              # Find the ratio and convert out of log scale
              r_l <- exp(r_l_term1)
              # Prior contribution
              r_p <- colSums(Ak0)[k]/(beta+S-1-colSums(Ak0)[k])
              # Draw whether or not study s should have factor k
              A[s,k] <- rbinom(1,1,(r_p*r_l)/(r_p*r_l+1))
              if (is.na(A[s,k])) {
                if (r_p==0) {
                  A[s,k] <- 0
                } else {
                  A[s,k] <- 1
                }
              }
              # Sample l if A[s,k] changed from 0->1
              if (A[s,k]==1&curr==0) {
                if ((K == 1)&(sum(A[s,])==0)) {
                  As <- matrix(0,1,1)
                } else {
                  As <- diag(A[s,])
                }
                # Compute inverse variance term
                l_s_term <- diag(K) + As%*%t(Lambda)%*%psi_inv_s[[s]]%*%Lambda%*%As
                # Compute standard deviation, then variance 
                l_s_sd <- solve(qr.R(qr(chol(l_s_term))))
                l_s_var <- tcrossprod(l_s_sd)
                # Compute mean term
                l_s_mean <- t(l_s_var%*%As%*%t(Lambda)%*%psi_inv_s[[s]]%*%t(Y_s[[s]]))
                # Put together: mean + sd*rnorm
                l_s[[s]][,k] <- l_s_mean[,k] + (matrix(rnorm(n_s[s] * K), nrow = n_s[s], ncol = K)%*%t(l_s_sd))[,k]
              }
           }
        }
        # Delete any factors that became all zeros
        factors_to_keep <- which(colSums(A)>0)
        if (length(factors_to_keep)==0) {
          factors_to_keep <- c(1)
        }
        K <- length(factors_to_keep)
        A <- matrix(A[1:S,factors_to_keep],S,K)
        Lambda <- matrix(Lambda[1:P,factors_to_keep],P,K)
        omega <- matrix(omega[1:P,factors_to_keep],P,K)
        l_s <- lapply(1:S, function(x) matrix(l_s[[x]][1:n_s[x],factors_to_keep],n_s[x],K))
        delta <- delta[factors_to_keep]
        tau <- cumprod(delta)
        D <- matrix(matvec(omega,tau),P,K)
      }

      # Step 2: new factors (A)
      # For each study...
      for(s in 1:S) {
        # Draw number of new factors to add
        k_new <- rpois(1,(alpha*beta)/(beta+S-1))
        # If adding a non-zero number of factors...
        if (k_new > 0) {
          # Draw new parameters from priors
          omega_new <- matrix(rgamma(P * k_new, shape = nu/2, scale = 2/nu),P,k_new)
          if (colSums(A)[1] != 0) {
            delta_new <- c(delta,rgamma(k_new, shape = c(rep(a2, k_new)),
                                        scale = c(rep(b2, k_new))))
          } else {
            delta_new <- rgamma(k_new, shape=c(a1,rep(a2, k_new-1)),
                                scale = c(b1,rep(b2, k_new-1)))
          }
          tau_new <- cumprod(delta_new)
          tau_new <- tau_new[(length(delta_new)-k_new+1):length(delta_new)]
          l_s_new <- lapply(1:S, function(x) mvrnorm(n_s[x], rep(0, k_new), diag(rep(1, k_new))))
          Lambda_new <- matrix(sapply(1:P, function(x) 
            sapply(1:k_new, function(y) rnorm(1,0,1/(omega_new[x,y]*tau_new[y])))),P,k_new)
          # The new factors being added are automatically individual ones
          A_new <- array(0, dim=c(S,k_new))
          A_new[s,] <- 1
          As_new <- diag(k_new)
          Dpnew <- matrix(matvec(omega_new,tau_new),P,k_new)
          if ((K == 1)&(A[s,1]==0)) {
            As <- matrix(0,1,1)
          } else {
            As <- diag(A[s,])
          }
          # Compute the acceptance ratio
          r_prod <- 1
          # For each variable...
          for(p in 1:P) {
            if (k_new<=1) {
              Dpnew_r <- Dpnew[p,]
              Dpnew_inv <- 1/Dpnew[p,]
            } else {
              Dpnew_r <- diag(Dpnew[p,])
              Dpnew_inv <- diag(1/Dpnew[p,])
            }
            # Construct variance term 
            var_term <- Dpnew_r + psi_s[[s]][p]*As_new%*%t(l_s_new[[s]])%*%l_s_new[[s]]%*%As_new
            # Find inverse of variance
            inv_var_term <- tcrossprod(solve(qr.R(qr(chol(var_term)))))
            # Compute Lambda bar
            Lambda_bar <- t(psi_s[[s]][p]*inv_var_term%*%As_new%*%t(l_s_new[[s]])%*%
                              (matrix(Y_s[[s]][,p],n_s[s],1)-l_s[[s]]%*%As%*%matrix(Lambda[p,],K,1)))
            # Multiply this variable's term onto the acceptance ratio
            r_prod <- r_prod * det(matrix(2*pi*Dpnew_inv,k_new,k_new))^(-0.5) * det(2*pi*inv_var_term)^(0.5) *
              exp(0.5*Lambda_bar%*%var_term%*%t(Lambda_bar))
          }
          # Find acceptance probability
          r <- min(1,r_prod)
          # If we choose to accept...
          if (runif(1)<r) {
            # Append new parameters onto existing ones
            if (colSums(A)[1] != 0) {
              l_s_new <- lapply(1:S, function(x) cbind(l_s[[x]],l_s_new[[x]]))
              Lambda_new <- cbind(Lambda,Lambda_new)
              A_new <- cbind(A,A_new)
              omega_new <- cbind(omega,omega_new)
            } 
            K <- dim(A_new)[2]
            omega <- omega_new
            delta <- delta_new
            tau <- cumprod(delta)
            D <- matvec(omega, tau)
            l_s <- l_s_new
            A <- A_new
            Lambda <- Lambda_new
          }
        }
      }
    }
    
    # Step 3: factor loadings (Lambda)
    # For each variable...
    for(p in 1:P) {
      Lambda_term1 <- Lambda_term2 <- list()
      # For each study...
      for(s in 1:S) {
        if ((K == 1)&(sum(A[s,])==0)) {
          As <- matrix(0,1,1)
        } else {
          As <- diag(A[s,])
        }
        # Term 1 is part of the inverse variance
        Lambda_term1[[s]] <- psi_s[[s]][p]*As%*%t(l_s[[s]])%*%l_s[[s]]%*%As
        # Term 2 is part of the mean
        Lambda_term2[[s]] <- psi_s[[s]][p]*As%*%t(l_s[[s]])%*%Y_s[[s]][,p]
      }
      # Sum up over all studies
      Lambda_term1 <- Reduce('+',Lambda_term1)
      Lambda_term2 <- Reduce('+',Lambda_term2)
      if (K<=1) {
        Dp <- D[p,]
      } else {
        Dp <- diag(D[p,])
      }
      Lambda_var <- t(chol(Dp + Lambda_term1))
      Lambda_sd <- solve(qr.R(qr(chol(Dp + Lambda_term1))))
      # Put together: sd*rnorm + mean 
      Lambda[p,] <- matrix(rnorm(K),1,K)%*%t(Lambda_sd)+
                        t(backsolve(t(Lambda_var),forwardsolve(Lambda_var,Lambda_term2)))
    }
    
    # Step 4: latent factors (l_s)
    # For each study...
    for(s in 1:S) {
      if ((K == 1)&(sum(A[s,])==0)) {
        As <- matrix(0,1,1)
      } else {
        As <- diag(A[s,])
      }
      # Compute inverse variance term
      l_s_term <- diag(K) + As%*%t(Lambda)%*%psi_inv_s[[s]]%*%Lambda%*%As
      # Compute standard deviation, then variance 
      l_s_sd <- solve(qr.R(qr(chol(l_s_term))))
      l_s_var <- tcrossprod(l_s_sd)
      # Compute mean term
      l_s_mean <- t(l_s_var%*%As%*%t(Lambda)%*%psi_inv_s[[s]]%*%t(Y_s[[s]]))
      # Put together: mean + sd*rnorm
      l_s[[s]] <- l_s_mean + matrix(rnorm(n_s[s] * K), nrow = n_s[s], ncol = K)%*%t(l_s_sd)
    }
    
    # Steps 5-7: factor loading priors (omega, delta, tau, D)
    if (K<=1) {
      dtau <- tau
    } else {
      dtau <- diag(tau)
    }
    # Update omega
    tau_omega_prod <- Lambda^2 %*% dtau
    omega <- matrix(rgamma(P*K, shape = (nu + 1) / 2,
                           rate = (nu + tau_omega_prod)/2),P,K)
    # Update first delta value
    delta[1] <- rgamma(1, shape = a1 + (P*K)/2, scale = 1/(b1 + 0.5 * 
                                                             sum(tau * colSums(omega * Lambda^2))/delta[1]))
    # Update tau
    tau <- cumprod(delta)
    # Iteratively update subsequent delta values and tau
    if (K > 1) {
      for (l in 2:K) {
        delta[l] <- rgamma(1, shape = a2 + 0.5*P*(K-l+1), 
                           scale = 1 / (b2 + 0.5*sum(tau[l:K]*colSums(as.matrix(omega[,l:K]*Lambda[,l:K]^2))/
                                                       delta[l])))
        tau <- cumprod(delta)
      }
    }
    D <- matvec(omega, tau)
    
    # Step 8: error terms (psi_s, psi_inv_s)
    psi_param <- list()
    # For each study...
    for (s in 1:S) {
      if ((K == 1)&(sum(A[s,])==0)) {
        As <- matrix(0,1,1)
      } else {
        As <- diag(A[s,])
      }
      # Study-dependent gamma parameter in posterior
      psi_param[[s]] <- Y_s[[s]] - l_s[[s]]%*%As%*%t(Lambda)
    }
    # Update psi_s (inverse errors) and psi_inv_s (matrix of inverse errors)
    psi_s <- lapply(1:S, function(x)
      rgamma(P, shape = apsi + (n_s[x])/2, 
             rate = bpsi + 0.5*colSums(psi_param[[x]]^2)))  
    psi_inv_s <- lapply(1:S, function(x) diag(psi_s[[x]]))
    
    # Store for output
    if(iter > burn){
      neff <- (iter - burn) / thin
      # Store loadings
      Lambdaout[[neff]] <- Lambda
      # Store factor indicator matrix
      Aout[[neff]] <- A
      # For each study...
      for(s in 1:S) {
        # Store errors 
        Psiout[[s]][[neff]] <- 1/psi_s[[s]]
        # Store latent factors
        l_out[[s]][[neff]] <- l_s[[s]]
      }
    }
    if (trace & iter %% nprint == 0) cat("iter=",iter,"\n")
  }
  
  ## Save and exit
  out <- list(Lambda = Lambdaout, Psi = Psiout, l_s = l_out, A = Aout)
  return(structure(out,  class="sp_msfa"))
}

# Gibbs sampler for Tetris with clustering
# Inputs: X = matrix in samples x features format
#         S = number of studies to cluster into
#         alpha, beta = IBP hyperparameters
#         trace = TRUE if iterations should be printed
#         nprint = how often to print iterations if trace is TRUE
#         initial = optional list of initial parameters 
#         fixed = FALSE if the factor indicator matrix A should be sampled 
#         A_fixed = what A should be fixed to, if fixed is TRUE
# Outputs: out = list of post-burn-in chains for the loadings matrix Lambda, the errors Psi, the latent factors l_s
#                for each study s, the factor indicator matrix A, and the study labels 
tetris_clustering <- function(X, S, alpha, beta, trace = TRUE, nprint = 1000, initial = 0, fixed = FALSE, A_fixed = 0,
                             control = list(...), total_iters = NULL, total_burn=NULL, ...)
{
  ## Read in data
  P <- ncol(X)               # number of variables
  Y <- scale(X,scale=F)      # centered data
  
  ## Set hyperparameters (see new_control() documentation)
  control <- new_control()
  nrun <- control$nrun
  thin <- control$thin
  burn <- control$burn
  sp  <- (nrun - burn) / thin
  apsi <- control$apsi 
  bpsi<- control$bpsi 
  nu <- control$nu
  a1 <- control$a1
  b1 <- control$b1 
  a2 <- control$a2
  b2 <-  control$b2
  
  if (!is.null(total_iters)) {
    nrun <- total_iters
    burn <- total_burn
  }
  
  ## Initialize studies
  idents <- sample(1:S,nrow(X),replace=T)
  Y_s <- lapply(1:S,function(s) X[idents==s,,drop=F])
  n_s <- sapply(Y_s,nrow)
  
  ## Initialize model parameters 
  # If no initial values provided, draw Lambda, A, Psi, and l_s from prior distributions;
  # otherwise, use the initial values that were given
  if (!is.list(initial)) {
    # Inverse error variances
    psi_s <- lapply(1:S, function(x) rgamma(P, shape = apsi, scale = 1 / bpsi))
    # Inverse error matrix
    psi_inv_s <- lapply(1:S, function(x) diag(psi_s[[x]]))
    # To hold the factor indicator matrix A
    A <- array(0, dim=c(S,S*S*alpha))
    # Draw the number of factors 
    num_factors <- rpois(1, alpha)
    while (num_factors==0) {
      num_factors <- rpois(1,alpha)
    }
    # Set the first study to have each of the instantiated factors
    if (num_factors > 0) {
      A[1,1:num_factors] <- 1
    }
    # For the remaining studies...
    for (i in 2:S) {
      # Assign each of the existing factors with a scaled probability
      if (num_factors > 0) {
        probs <- colSums(A)/(beta+i-1)
        realization <- rbinom(n=num_factors,size=1,prob=probs[1:num_factors])
      } else {
        realization <- c()
      }
      # Draw the number of new factors to add
      new_factors <- rpois(1,alpha*beta/(beta+i-1))
      while (new_factors==0) {
        new_factors <- rpois(1,alpha*beta/(beta+i-1))
      }
      num_factors <- num_factors + new_factors
      # Assign each of these new factors to the study being considered
      if (num_factors > 0) {
        A[i,1:num_factors] <- c(realization,rep(1,new_factors))
      }
    }
    A <- matrix(A[1:S,1:num_factors],S,num_factors)
    colnames(A) <- NULL
    # Update the number of factors
    K <- num_factors
    # If we are fixing A, establish those parameters
    if (fixed) {
      A <- A_fixed
      K <- dim(A)[2]
    }
  } else {
    # If we are starting with initial parameters, fill them in here
    Lambda <- initial[[1]]
    A <- initial[[2]]
    K <- dim(A)[2]
    psi_s <- initial[[3]]
    psi_inv_s <- lapply(1:S, function(x) diag(c(psi_s[[x]])))
    l_s <- initial[[4]]
  }
  # Draw omega, delta, and tau from prior distributions
  omega <- matrix(rgamma(P * K, shape = nu/2, scale = 2/nu),P,K)
  delta <- rgamma(K, shape = c(a1, rep(a2, K-1)),
                  scale = c(1/b1, rep(1/b2, K-1)))
  tau <- cumprod(delta)
  D <- matvec(omega, tau)
  # If no initialization, draw Lambda and l_s from priors
  if (!is.list(initial)) {
    Lambda <- t(sapply(1:P, function(x) sapply(1:K, function(y) rnorm(1,0,1/(omega[x,y]*tau[y])))))
    Lambda <- matrix(Lambda,P,K)
    l_s <- lapply(1:S, function(x) mvrnorm(n_s[x], rep(0, K), diag(rep(1, K))))
  }
  
  ## Initialize output
  Lambdaout <- Aout <- identsout <- list()
  Psiout <- l_out <- lapply(1:S, function(x) list())

  ## Start posterior sampling
  for(iter in 1:nrun)
  {
    # Only update factor indicator matrix if not fixed
    if (!fixed) {
      # Step 1: factor indicator matrix (A)
      if (colSums(A)[1]!=0) {
        # For each study...
        for(s in 1:S) {
          if (n_s[s]==0) {
            for (k in 1:K) {
              Ak0 <- A
              Ak0[s,k] <- 0
              r <- colSums(Ak0)[k]/(beta+S-1-colSums(Ak0)[k]) 
              A[s,k] <- rbinom(1,1,r/(r+1))
            }
          } else {
            # For each factor...
            for(k in 1:K) {
              curr <- A[s,k]
              # Ak1 = A where study s has factor k
              Ak1 <- A
              Ak1[s,k] <- 1
              # Ak0 = A where study s doesn't have factor k
              Ak0 <- A
              Ak0[s,k] <- 0
              # Hold terms for the likelihood ratios
              r_l_term1 <- list()
              # For each sample in a given study...
              for(n in 1:n_s[s]) {
                if (K>2) {
                  # Compute the marginalized likelihood term
                  lbar <- (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)^(-1)*
                    (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%
                       (t(Y_s[[s]][n,,drop=F])-Lambda[,-k,drop=F]%*%diag(Ak1[s,-k])%*%t(l_s[[s]][n,-k,drop=F])))
                }
                if (K==1) {
                  lbar <- (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)^(-1)*
                    (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%
                       (t(Y_s[[s]][n,,drop=F])-0))
                }
                if (K==2) {
                  lbar <- (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)^(-1)*
                    (t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%
                       (t(Y_s[[s]][n,,drop=F])-Lambda[,-k,drop=F]%*%diag(Ak1[s,-k,drop=F])%*%t(l_s[[s]][n,-k,drop=F])))
                }
                r_l_term1[[n]] <- log((t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)^(-0.5))+
                  0.5*(t(Lambda[,k,drop=F])%*%psi_inv_s[[s]]%*%Lambda[,k,drop=F]+1)*lbar^2
              }
              # Sum the terms to get log-likelihood of each case across samples
              r_l_term1 <- Reduce('+',r_l_term1)
              # Find the ratio and convert out of log scale
              r_l <- exp(r_l_term1)
              # Prior contribution
              r_p <- colSums(Ak0)[k]/(beta+S-1-colSums(Ak0)[k])
              # Draw whether or not study s should have factor k
              A[s,k] <- rbinom(1,1,(r_p*r_l)/(r_p*r_l+1))
              if (is.na(A[s,k])) {
                if (r_p==0) {
                  A[s,k] <- 0
                } else {
                  A[s,k] <- 1
                }
              }
              # Sample l if A[s,k] changed from 0->1
              if (A[s,k]==1&curr==0) {
                if ((K == 1)&(sum(A[s,])==0)) {
                  As <- matrix(0,1,1)
                } else {
                  As <- diag(A[s,])
                }
                # Compute inverse variance term
                l_s_term <- diag(K) + As%*%t(Lambda)%*%psi_inv_s[[s]]%*%Lambda%*%As
                # Compute standard deviation, then variance 
                l_s_sd <- solve(qr.R(qr(chol(l_s_term))))
                l_s_var <- tcrossprod(l_s_sd)
                # Compute mean term
                l_s_mean <- t(l_s_var%*%As%*%t(Lambda)%*%psi_inv_s[[s]]%*%t(Y_s[[s]]))
                # Put together: mean + sd*rnorm
                l_s[[s]][,k] <- l_s_mean[,k] + (matrix(rnorm(n_s[s] * K), nrow = n_s[s], ncol = K)%*%t(l_s_sd))[,k]
              }
            }
          }
        }
        # Delete any factors that became all zeros
        factors_to_keep <- which(colSums(A)>0)
        if (length(factors_to_keep)==0) {
          factors_to_keep <- c(1)
        }
        K <- length(factors_to_keep)
        A <- matrix(A[1:S,factors_to_keep],S,K)
        Lambda <- matrix(Lambda[1:P,factors_to_keep],P,K)
        omega <- matrix(omega[1:P,factors_to_keep],P,K)
        l_s <- lapply(1:S, function(x) matrix(l_s[[x]][1:max(1,n_s[x]),factors_to_keep],max(1,n_s[x]),K))
        delta <- delta[factors_to_keep]
        tau <- cumprod(delta)
        D <- matrix(matvec(omega,tau),P,K)
      }
      
      # Step 2: new factors (A)
      # For each study...
      for(s in 1:S) {
        # Draw number of new factors to add
        k_new <- rpois(1,(alpha*beta)/(beta+S-1))
        # If adding a non-zero number of factors...
        if (k_new > 0) {
          # Draw new parameters from priors
          omega_new <- matrix(rgamma(P * k_new, shape = nu/2, scale = 2/nu),P,k_new)
          if (colSums(A)[1] != 0) {
            delta_new <- c(delta,rgamma(k_new, shape = c(rep(a2, k_new)),
                                        scale = c(rep(b2, k_new))))
          } else {
            delta_new <- rgamma(k_new, shape=c(a1,rep(a2, k_new-1)),
                                scale = c(b1,rep(b2, k_new-1)))
          }
          tau_new <- cumprod(delta_new)
          tau_new <- tau_new[(length(delta_new)-k_new+1):length(delta_new)]
          l_s_new <- lapply(1:S, function(x) 
            array(mvrnorm(max(1,n_s[x]), rep(0, k_new), diag(rep(1, k_new))),dim=c(max(n_s[x],1),k_new)))
          Lambda_new <- matrix(sapply(1:P, function(x) 
            sapply(1:k_new, function(y) rnorm(1,0,1/(omega_new[x,y]*tau_new[y])))),P,k_new)
          # The new factors being added are automatically individual ones
          A_new <- array(0, dim=c(S,k_new))
          A_new[s,] <- 1
          As_new <- diag(k_new)
          Dpnew <- matrix(matvec(omega_new,tau_new),P,k_new)
          if ((K == 1)&(A[s,1]==0)) {
            As <- matrix(0,1,1)
          } else {
            As <- diag(A[s,])
          }
          # Compute the acceptance ratio
          r_prod <- 1
          if (n_s[s]>0) {
            # For each variable...
            for(p in 1:P) {
              if (k_new<=1) {
                Dpnew_r <- Dpnew[p,]
                Dpnew_inv <- 1/Dpnew[p,]
              } else {
                Dpnew_r <- diag(Dpnew[p,])
                Dpnew_inv <- diag(1/Dpnew[p,])
              }
              # Construct variance term 
              var_term <- Dpnew_r + psi_s[[s]][p]*As_new%*%t(l_s_new[[s]])%*%l_s_new[[s]]%*%As_new
              # Find inverse of variance
              inv_var_term <- tcrossprod(solve(qr.R(qr(chol(var_term)))))
              # Compute Lambda bar
              Lambda_bar <- t(psi_s[[s]][p]*inv_var_term%*%As_new%*%t(l_s_new[[s]])%*%
                                (matrix(Y_s[[s]][,p],n_s[s],1)-l_s[[s]]%*%As%*%matrix(Lambda[p,],K,1)))
              # Multiply this variable's term onto the acceptance ratio
              r_prod <- r_prod * det(matrix(2*pi*Dpnew_inv,k_new,k_new))^(-0.5) * det(2*pi*inv_var_term)^(0.5) *
                exp(0.5*Lambda_bar%*%var_term%*%t(Lambda_bar))
            }
          }
          # Find acceptance probability
          r <- min(1,r_prod)
          # If we choose to accept...
          if (runif(1)<r) {
            # Append new parameters onto existing ones
            if (colSums(A)[1] != 0) {
              l_s_new <- lapply(1:S, function(x) cbind(l_s[[x]],l_s_new[[x]]))
              Lambda_new <- cbind(Lambda,Lambda_new)
              A_new <- cbind(A,A_new)
              omega_new <- cbind(omega,omega_new)
            } 
            K <- dim(A_new)[2]
            omega <- omega_new
            delta <- delta_new
            tau <- cumprod(delta)
            D <- matvec(omega, tau)
            l_s <- l_s_new
            A <- A_new
            Lambda <- Lambda_new
          }
        }
      }
    }
    
    # Step 3: factor loadings (Lambda)
    # For each variable...
    for(p in 1:P) {
      Lambda_term1 <- Lambda_term2 <- list()
      # For each study...
      for(s in 1:S) {
        if (n_s[s]>0) {
          if ((K == 1)&(sum(A[s,])==0)) {
            As <- matrix(0,1,1)
          } else {
            As <- diag(A[s,])
          }
          # Term 1 is part of the inverse variance
          Lambda_term1[[s]] <- psi_s[[s]][p]*As%*%t(l_s[[s]])%*%l_s[[s]]%*%As
          # Term 2 is part of the mean
          Lambda_term2[[s]] <- psi_s[[s]][p]*As%*%t(l_s[[s]])%*%Y_s[[s]][,p]
        } else {
          Lambda_term1[[s]] <- 0
          Lambda_term2[[s]] <- 0
        }
      }
      # Sum up over all studies
      Lambda_term1 <- Reduce('+',Lambda_term1)
      Lambda_term2 <- Reduce('+',Lambda_term2)
      if (K<=1) {
        Dp <- D[p,]
      } else {
        Dp <- diag(D[p,])
      }
      Lambda_var <- t(chol(Dp + Lambda_term1))
      Lambda_sd <- solve(qr.R(qr(chol(Dp + Lambda_term1))))
      # Put together: sd*rnorm + mean 
      Lambda[p,] <- matrix(rnorm(K),1,K)%*%t(Lambda_sd)+
        t(backsolve(t(Lambda_var),forwardsolve(Lambda_var,Lambda_term2)))
    }
    
    # Step 4: latent factors (l_s)
    # For each study...
    for(s in 1:S) {
      if (n_s[s]>0) {
        if ((K == 1)&(sum(A[s,])==0)) {
          As <- matrix(0,1,1)
        } else {
          As <- diag(A[s,])
        }
        # Compute inverse variance term
        l_s_term <- diag(K) + As%*%t(Lambda)%*%psi_inv_s[[s]]%*%Lambda%*%As
        # Compute standard deviation, then variance 
        l_s_sd <- solve(qr.R(qr(chol(l_s_term))))
        l_s_var <- tcrossprod(l_s_sd)
        # Compute mean term
        l_s_mean <- t(l_s_var%*%As%*%t(Lambda)%*%psi_inv_s[[s]]%*%t(Y_s[[s]]))
        # Put together: mean + sd*rnorm
        l_s[[s]] <- l_s_mean + matrix(rnorm(n_s[s] * K), nrow = n_s[s], ncol = K)%*%t(l_s_sd)
      }
    }
    
    # Steps 5-7: factor loading priors (omega, delta, tau, D)
    if (K<=1) {
      dtau <- tau
    } else {
      dtau <- diag(tau)
    }
    # Update omega
    tau_omega_prod <- Lambda^2 %*% dtau
    omega <- matrix(rgamma(P*K, shape = (nu + 1) / 2,
                           rate = (nu + tau_omega_prod)/2),P,K)
    # Update first delta value
    delta[1] <- rgamma(1, shape = a1 + (P*K)/2, scale = 1/(b1 + 0.5 * 
                                                             sum(tau * colSums(omega * Lambda^2))/delta[1]))
    # Update tau
    tau <- cumprod(delta)
    # Iteratively update subsequent delta values and tau
    if (K > 1) {
      for (l in 2:K) {
        delta[l] <- rgamma(1, shape = a2 + 0.5*P*(K-l+1), 
                           scale = 1 / (b2 + 0.5*sum(tau[l:K]*colSums(as.matrix(omega[,l:K]*Lambda[,l:K]^2))/
                                                       delta[l])))
        tau <- cumprod(delta)
      }
    }
    D <- matvec(omega, tau)
    
    # Step 8: error terms (psi_s, psi_inv_s)
    psi_param <- list()
    # For each study...
    for (s in 1:S) {
      if ((K == 1)&(sum(A[s,])==0)) {
        As <- matrix(0,1,1)
      } else {
        As <- diag(A[s,])
      }
      # Study-dependent gamma parameter in posterior
      if (n_s[s]>0) {
        psi_param[[s]] <- Y_s[[s]] - l_s[[s]]%*%As%*%t(Lambda)
      } else {
        psi_param[[s]] <- array(0,dim=c(P,1))
      }
    }
    # Update psi_s (inverse errors) and psi_inv_s (matrix of inverse errors)
    psi_s <- lapply(1:S, function(x)
      rgamma(P, shape = apsi + (n_s[x])/2, 
             rate = bpsi + 0.5*colSums(psi_param[[x]]^2)))  
    psi_inv_s <- lapply(1:S, function(x) diag(psi_s[[x]]))
    
    # Step 9: identities 
    for (n in 1:nrow(X)) {
      marg_liks <- rep(0,S)
      for (s in 1:S) {
        marg_liks[s] <- mvtnorm::dmvnorm(X[n,],rep(0,P),Lambda%*%diag(A[s,])%*%t(Lambda)+diag(1/psi_s[[s]]),log=T)
      }
      idents[n] <- sample(1:S,1,prob=exp(marg_liks-max(marg_liks)))
    }
    Y_s <- lapply(1:S,function(s) X[idents==s,,drop=F])
    n_s <- sapply(Y_s,nrow)
    l_s <- do.call(rbind,l_s)
    l_s <- lapply(1:S,function(s) l_s[idents==s,,drop=F])
    
    for(s in 1:S) {
      if ((K == 1)&(sum(A[s,])==0)) {
        As <- matrix(0,1,1)
      } else {
        As <- diag(A[s,])
      }
      # Compute inverse variance term
      l_s_term <- diag(K) + As%*%t(Lambda)%*%psi_inv_s[[s]]%*%Lambda%*%As
      # Compute standard deviation, then variance 
      l_s_sd <- solve(qr.R(qr(chol(l_s_term))))
      l_s_var <- tcrossprod(l_s_sd)
      # Compute mean term
      l_s_mean <- t(l_s_var%*%As%*%t(Lambda)%*%psi_inv_s[[s]]%*%t(Y_s[[s]]))
      if (n_s[s]==0) {
        l_s_mean <- 0
      }
      # Put together: mean + sd*rnorm
      l_s[[s]] <- l_s_mean + matrix(rnorm(max(1,n_s[s]) * K), nrow = max(1,n_s[s]), ncol = K)%*%t(l_s_sd)
      if (max(1,n_s[s])==1) {
        l_s[[s]] <- array(l_s[[s]],dim=c(1,K))
      }
    }
    
    # Store for output
    if(iter > burn){
      neff <- (iter - burn) / thin
      # Store loadings
      Lambdaout[[neff]] <- Lambda
      # Store factor indicator matrix
      Aout[[neff]] <- A
      # Store clustering identities
      identsout[[neff]] <- idents
      # For each study...
      for(s in 1:S) {
        # Store errors 
        Psiout[[s]][[neff]] <- 1/psi_s[[s]]
        # Store latent factors
        l_out[[s]][[neff]] <- l_s[[s]]
      }
    }
    if (trace & iter %% nprint == 0) cat("iter=",iter,"\n")
  }
  
  ## Save and exit
  out <- list(Lambda = Lambdaout, Psi = Psiout, l_s = l_out, A = Aout, idents = identsout)
  return(structure(out,  class="sp_msfa"))
}

# Helper function: convert binary matrix to integer representation
bitsToInt<-function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}

# Helper function: compute density under IBP
IBP.prob <- function(A,alpha_IBP) {
  if (is.null(dim(A))) {
    if (sum(A) == 0) { return(dpois(0,alpha_IBP)) }
  } else {
    if (sum(colSums(A))==0) { return(dpois(0,alpha_IBP)) }
  }
  A <- A[,which(colSums(A)>0),drop=F]
  Kplus <- ncol(A)
  N <- nrow(A)
  Kh <- table(factor(apply(A,2,bitsToInt),levels=1:(2^N-1)))[1:(2^N-1)]
  return(Kplus*log(alpha_IBP)-sum(log(sapply(Kh,function(k) factorial(k))))+sum(log(sapply(colSums(A),function(mk) factorial(N-mk)*factorial(mk-1)/factorial(N)))))
}

# Helper function: compute distance between two A matrices
dist.A <- function(A1,A2) {
  if (ncol(A1) < ncol(A2)) {
    temp <- A2
    A2 <- A1
    A1 <- temp
  }
  if (ncol(A1)!=ncol(A2)) {
    A2 <- cbind(A2,array(0,dim=c(nrow(A2),ncol(A1)-ncol(A2))))
  }
  M <- array(0,dim=c(ncol(A1),ncol(A1)))
  for (k in 1:ncol(A1)) {
    for (j in 1:ncol(A2)) {
      M[k,j] <- sum(abs(A1[,k]-A2[,j]))
    }
  }
  s <- solve_LSAP(M)
  return(sum(M[cbind(seq_along(s), s)]))
}

# Choose A from a chain
# Input:  out = posterior output from running tetris()
#         alpha_IBP = IBP hyperparameter used when running tetris()
#         S = number of studies
# Output: point estimate of A
choose.A <- function(out,alpha_IBP,S) {
  A.chain <- out[[4]]
  Lambda.chain <- out[[1]]
  iters <- length(A.chain)
  dists <- array(0,dim=c(iters,iters))
  for (i in 1:iters) {
    for (j in i:iters) {
      dists[i,j] <- dist.A(A.chain[[i]],A.chain[[j]])
    }
  }
  for (i in 1:iters) {
    for (j in 1:(i-1)) {
      dists[i,j] <- dists[j,i]
    }
  }
  r <- rowSums(dists)
  thresh <- max(quantile(c(dists),0.05),S)
  r2 <- sapply(1:nrow(dists),function(x) sum(dists[x,]<=thresh))
  r2.max <- which(r2==max(r2))
  if (length(r2.max)==1) {
    radius <- A.chain[[r2.max]]
  } else {
    facs <- sapply(r2.max,function(x) ncol(A.chain[[x]]))
    if (length(which(facs==min(facs)))==1) {
      radius <- A.chain[[r2.max[which.min(facs)]]]
    } else {
      minfacs <- which(facs==min(facs))
      priors <- sapply(minfacs,function(x) IBP.prob(A.chain[[r2.max[x]]],alpha_IBP))
      radius <- A.chain[[r2.max[minfacs[which.max(priors)]]]]
    }
  }
  return(radius)
}

# Helper functions for Lambda recovery
eval_f <- function(Lambda,Sigmas,A) {
  Lambda <- array(Lambda,dim=c(nrow(Sigmas[[1]]),ncol(A)))
  return(sum(sapply(1:nrow(A),function(s) norm(Sigmas[[s]]-Lambda%*%diag(A[s,])%*%t(Lambda),type='f')^2)))
} 

eval_g <- function(Lambda,Sigmas,A) {
  Lambda <- array(Lambda,dim=c(nrow(Sigmas[[1]]),ncol(A)))
  return(c(Reduce('+',lapply(1:nrow(A),function(s) 
    (-4)*(Sigmas[[s]]-Lambda%*%diag(A[s,])%*%t(Lambda))%*%Lambda%*%diag(A[s,])))))
}

# Recover Lambda from a chain with fixed A
# Input: out = posterior output from running tetris() with fixed=T
#        A = fixed A when running tetris()
# Output: Lambda point estimate
getLambda <- function(out,A) {
  P <- nrow(out$Lambda[[1]])
  K <- ncol(A)
  S <- nrow(A)
  
  Sigmas <- list()
  for (s in 1:S) {
    LLTs <- array(0,dim=c(P,P,2000)) 
    for (i in 1:2000) {
      LLTs[,,i] <- out$Lambda[[i]]%*%diag(A[s,])%*%t(out$Lambda[[i]])
    }
    Sigmas[[s]] <- apply(LLTs,c(1,2),mean)
  }
  
  res <- nloptr(c(out$Lambda[[2000]]),eval_f=eval_f,eval_grad_f=eval_g,
                opts=list(algorithm='NLOPT_LD_LBFGS',maxeval=1000),
                Sigmas=Sigmas,A=A) 
  return(array(res$solution,dim=c(P,K)))
}	
