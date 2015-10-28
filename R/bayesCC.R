#' Bayesian Consensus Clustering
#' 
#' @param X             a list of data matrices, each with D_i rows & N columns.
#' @param K             integer, maximum number of clusters for K-means
#' @param a             numeric, hyperparameter for Alpha ~ Beta(a, b)
#' @param b             numeric, hyperparameter for Alpha ~ Beta(a, b)
#' @param IndivAlpha    boolean, whether to fit individual random effects
#' @param mu0           list of initial mean parameters for the Normal-Gamma
#' @param a0            list of initial shape parameters for the Normal-Gamma
#' @param b0            list of initial rate parameters for the Normal-Gamma
#' @param Concentration initial concentration parameter for Dirichlet process
#' @param maxiter       how many iterations of the MCMC sampler should be run?
#' 
#' @return a list with elements (Alpha, AlphaBounds, Cbest, Lbest, AlphaVec)
#'
#' @details
#' 
#' Reference: 
#' Lock EF and Dunson DB, 
#' "Bayesian Consensus Clustering", Bioinformatics, 29(20), 2013. 
#' 
#' The output of bayesCC(...) has several pieces:
#'
#' \itemize{ 
#'   \item Alpha. the average adherence (by data source, if IndivAlpha==T).
#'   \item AlphaBounds. the 95 percent credible interval for Alpha.
#'   \item Cbest. the "hard" overall clustering, as a binary matrix.
#'   \item Lbest. a list of the separate clusterings by data source.
#'   \item AlphaVec. a vector of alpha values over MCMC draws to assess mixing.
#' }
#'
#' Data matrices in X should have the same number of columns (one per subject),
#' but may have different numbers of rows.  If a subject is missing for a data
#' source, a nice improvement would be to marginalize over the remaining 
#' columns, perhaps after determining their overall cluster membership(s). 
#' If a row is missing for a data source, k-NN imputation should suffice. 
#' 
#' It would be nice to parallelize the runs over all candidate values for K. 
#' Similarly, PAM or NMF can be more robust than K-means in some situations. 
#' Expect the next point release of the package to support either or both. 
#'
#' Note that the first (maxiter / 2) iterations are used as burn-in for MCMC.
#' 
#' Implementation details are given in the PDF found at 
#' \url{http://www.tc.umn.edu/~elock/software/BCC.pdf}
#' This is more extensive than the Bioinformatics paper. 
#'
#' FIXME (maybe): 
#' Might be nice to use PAM and/or NMF clustering instead of K-means.
#'
#' FIXME (maybe):
#' See if it's possible to do matrix completion aided by cluster assignments,
#' for the case when entire columns are NA or mostly-NA (*cough* TARGET *cough*)
#'
#' @examples
#' 
#' \dontrun{ 
#'
#'   # try a few
#'   Ks <- 2:5
#'   names(Ks) <- paste0("K", Ks)
#'
#'   # can take a while...
#'   data(BRCAData)
#'   runK <- function(k) bayesCC(BRCAData, K=k, IndivAlpha=T, maxiter=10000)
#'   Results <- mclapply(Ks, runK)
#'  
#'   # ?alphaStar
#'   alphaStarDist <- data.frame(lapply(Results, alphaStar))
#'   boxplot(alphaStarDist, main="Mean-adjusted adherence by K (optimal: K=3)")
#' 
#' }
#' 
#' @seealso alphaStar
#'
#' @export
bayesCC <- function(X, K=2, a=1, b=1, IndivAlpha=FALSE, mu0=list(), a0=list(), 
                    b0 = list(), Concentration = 1, maxiter = 1000, ...) {

  # Initialize parameters
  Gamma <- rep(1, K)  #Posterior Dirichlet concentration
  Lbest <- list()  #List of best clustering for each data source
  L <- list() #Current clustering for each source (updates each MCMC draw)
  Llist <- list()  #saves all source clustering realizations
  Clist <- list()  #saves all overall clustering realizations
  mu <- list()  #list of mean vector for each data source (updates each draw)
  logF <- list()  #log likeihood matrix for each data source
  logL <- list()  #log clustering probabilities for each data source
  Lp <- list()  #clustering probabilities for each data source
  M <- length(X)  #Number of data sources
  N <- dim(X[[1]])[2]  #sample size
  d <- sapply(X, nrow)  #Vector of dimension for each data source
  S <- list()  #sample variance for each data source
  A <- list()  #Posterior gamma shape parameter
  B <- list()  #Posterior gamma rate parameter
  Tau <- list()  #1/Tau = posterior variance
  Sigma <- list()  #Posterior variance

  # Backwards compatibility with the original BCC.r
  if ("NumDraws" %in% names(list(...))) maxiter <- list(...)$NumDraws 

  # overall clustering probabilities
  Cprobs <- matrix(nrow = N, ncol = K)

  # for each data source...  
  for (m in 1:M) {
    Llist[[m]] <- list()
    S[[m]] <- matrix(nrow = d[m], ncol = K)
    A[[m]] <- matrix(nrow = d[m], ncol = K)
    B[[m]] <- matrix(nrow = d[m], ncol = K)
    Tau[[m]] <- matrix(nrow = d[m], ncol = K)
    Sigma[[m]] <- matrix(nrow = d[m], ncol = K)

    # Determine a0, b0 based on overall sample variance
    StDev <- apply(X[[m]], 1, sd) 

    a0[[m]] <- rep(1, d[m])
    b0[[m]] <- StDev ^ 2
    if (d[m] > 1) mu0[[m]] <- rowMeans(X[[m]]) # Find mu0 from overall mean
    if (d[m] == 1) mu0[[m]] <- mean(X[[m]])

    # Initialize via K-means clustering of each data source
    # FIXME: perhaps use PAM and/or consensus clustering instead
    InitL <- kmeans(t(X[[m]]), K) 

    if (m == 1) prevClusters <- InitL$cluster
    if (m > 1) { 
      # Quick, imperfect step that helps to align cluster indices
      InitL$cluster <- alignClusters(prevClusters, InitL$cluster) 
      prevClusters <- InitL$cluster
    }
    mu[[m]] <- matrix(nrow = d[m], ncol = K)
    for (k in 1:K) {
      if (d[m] > 1) {
        Sigma[[m]][,k] <- apply(X[[m]][, InitL$cluster == k], 1, sd)
        mu[[m]][,k] <- apply(X[[m]][, InitL$cluster == k], 1, mean)
      }
      if (d[m] == 1) {
        Sigma[[m]][,k] <- sd(X[[m]][, InitL$cluster == k])
        mu[[m]][,k] <- mean(X[[m]][, InitL$cluster == k])
      }
    }
    logF[[m]] <- matrix(nrow = N, ncol = K)
    logL[[m]] <- matrix(nrow = N, ncol = K)
    L[[m]] <- matrix(nrow = N, ncol = K)
  }

  # Vector of alpha (adherence) parameters
  alphaVec <- c()

  # overall clustering
  C <- matrix(nrow = N, ncol = K)

  # source-specific cluster probabiities (given C)
  nu <- array(rep(1 / K, N * M * K), dim = c(N, M, K))

  # cluster probabilites
  Pi <- rep(1 / K, K)

  # size of each source-specific cluster
  n <- matrix(nrow = M, ncol = K)

  alpha <- rep(0, M)
  if (IndivAlpha == TRUE) {
    for (m in 1:M) {
      while (alpha[m] < 1 / K) {
        alpha[m] <- rbeta(1, a, b)
      }
    }
  } else { 
    while (alpha[1] <= 1 / K) {
      alpha[] <- rbeta(1, a, b)
    }
  }

  # MCMC loop
  for (w in 1:maxiter) {
    # w is the current MCMC iteration
    for (m in 1:M) {
      for (k in 1:K) {
        # Determine logF from normal density
        logF[[m]][,k] <- log(nu[,m,k]) - 
                          sum(log(Sigma[[m]][,k])) + 
                          (-d[m] * log(2 * pi) - 
                           colSums((((X[[m]] - mu[[m]][,k]) / 
                                    Sigma[[m]][,k]) ^ 2))) / 2 
      }
      logL[[m]] <- logF[[m]] - apply(logF[[m]], 1, logSum)  #Normalize
    }

    # log-sum-exp FTW?
    Lp <- lapply(logL, exp)

    for (m in 1:M) {

      for (i in 1:N) L[[m]][i, ] <- rmultinom(1, 1, Lp[[m]][i, ])
      if (w > 1) L[[m]] <- alignClusters(C, L[[m]], type = "mat")
      n[m,] <- colSums(L[[m]])

      for (k in 1:K) {

        # Update cluster parameters based on normal-gamma distribution
        if (d[m] == 1 & n[m,k] > 1) {
          S[[m]][,k] <- sd(X[[m]][, L[[m]][,k] == 1]) ^ 2
          PostMean <- sum(X[[m]][, L[[m]][,k] == 1]) / (n[m,k] + 1)
          B[[m]][,k] <- b0[[m]] + 0.5 * (n[m,k] * S[[m]][,k] + n[m,k] * 
                                          (mean(X[[m]][, L[[m]][,k] == 1]) - 
                                           mu0[[m]]) ^ 2 / (1 + n[m,k]))
        }

        if (d[m] > 1 & n[m,k] > 1) {

          # posterior mean
          PostMean <- (mu0[[m]] + rowSums(X[[m]][,L[[m]][,k] == 1])) / 
                      (n[m,k] + 1)

          # variance
          S[[m]][,k] <- apply(X[[m]][,L[[m]][,k] == 1], 1, sd) ^ 2

          # posterior gamma rate parameter
          B[[m]][,k] <- b0[[m]] + 0.5 * (n[m,k] * S[[m]][,k] + n[m,k] * 
          (rowMeans(X[[m]][,L[[m]][,k] == 1]) - mu0[[m]]) ^ 2 / (1 + n[m,k]))

        }

        if (n[m,k] == 1) {

          # posterior mean
          PostMean <- (mu0[[m]] + X[[m]][, L[[m]][,k] == 1]) / 2
          
          # posterior gamma rate parameter
          B[[m]][,k] <- b0[[m]] + 
                         0.5 * (X[[m]][, L[[m]][,k] == 1] - mu0[[m]]) ^ 2 / 2
        
        }
        
        if (n[m,k] == 0) {
          PostMean <- mu0[[m]]
          B[[m]][,k] <- b0[[m]]
        }

        Lambda <- 1 + n[m,k]
        A[[m]][,k] <- a0[[m]] + n[m,k] / 2
        Tau[[m]][,k] <- rgamma(d[m], shape = A[[m]][,k], rate = B[[m]][,k])
        mu[[m]][,k] <- rnorm(d[m], PostMean, sqrt(1 / (Tau[[m]][,k] * Lambda)))
        Sigma[[m]][,k] <- sqrt(1 / Tau[[m]][,k])
      }
    }
    
    if (w > 1) {

      ### Update alpha
      alpha <- rep(1 / K, M)
      if (!IndivAlpha) {
        NumEq <- 0
        for (m in 1:M) {
          NumEq <- NumEq + sum(L[[m]] == 1 & C == 1)
        }
        for (Count in 1:10) {
          alphaTemp <- rbeta(1, a + NumEq, b + M * N - NumEq)
          
          # generate from beta dist until result >1/K 
          # (or set to 1/K after 10 trys)
          if (alphaTemp > 1 / K) {
            alpha[] <- alphaTemp
            break
          }
        }
      }

      if (IndivAlpha) {
        alpha <- rep(1 / K, M)
        for (m in 1:M) {
          NumEq <- sum(L[[m]] == 1 & C == 1)
          for (Count in 1:10) {
            alphaTemp <- rbeta(1, a + NumEq, b + N - NumEq)
            # generate from beta dist until result >1/K 
            # (or set to 1/K after 10 trys)
            if (alphaTemp > 1 / K) {
              alpha[m] <- alphaTemp
              break
            }
          }
        }
      }
    }

    ## Update C
    for (k in 1:K) {
      Cprobs[,k] <- Pi[k]
      for (m in 1:M) {
        Cprobs[,k] <- Cprobs[,k] * alpha[m] ^ (L[[m]][,k]) * 
                       ((1 - alpha[m]) / (K - 1)) ^ (1 - L[[m]][,k])
      }
    }

    Denom <- rowSums(Cprobs)
    for (i in 1:N) {
      Cprobs[i,] <- Cprobs[i,] / Denom[i]
      C[i,] <- rmultinom(1, 1, Cprobs[i,])
    }

    # update Pi
    Gamma <- Concentration + colSums(C)
    Pi <- rdirichlet(1, Gamma)
    
    if (IndivAlpha) {
      alphaVec <- cbind(alphaVec, alpha)
    } else { 
      alphaVec[w] <- alpha[1]
    }   

    # update tau
    for (m in 1:M) {
      for (k in 1:K) {
        nu[as.logical(C[,k]),m,k] <- alpha[m]
        nu[!as.logical(C[,k]),m,k] <- (1 - alpha[m]) / (K - 1)
      }
    }
    
    for (m in 1:M) Llist[[m]][[w]] <- L[[m]]
    Clist[[w]] <- C
  }

  ## Find Alpha by averaging over iteration (and 95% cred interval)
  if (IndivAlpha) {
    AlphaBounds <- matrix(nrow = length(alpha), ncol = 2)
    for (m in 1:M) {
      AlphaBounds[m,1] <- quantile(alphaVec[m,floor(maxiter / 5):maxiter],0.025)
      AlphaBounds[m,2] <- quantile(alphaVec[m,floor(maxiter / 5):maxiter],0.975)
    }
    Alpha <- rowMeans(alphaVec[, floor(maxiter / 5):maxiter])
  } else { 
    AlphaBounds <- c()
    AlphaBounds[1] <- quantile(alphaVec[floor(maxiter / 5):maxiter], 0.025)
    AlphaBounds[2] <- quantile(alphaVec[floor(maxiter / 5):maxiter], 0.975)
    Alpha <- mean(alphaVec[floor(maxiter / 5):maxiter])
  }
  
  # Choose hard clustering by least squares as in (Dahl,2006)
  for (m in 1:M) {
    Lkern <- Llist[[m]][[floor(maxiter / 5)]] %*% 
             t(Llist[[m]][[floor(maxiter / 5)]])
    for (w in floor(maxiter / 5 + 1):maxiter) {
      Lkern <- Lkern + Llist[[m]][[w]] %*% t(Llist[[m]][[w]])
    }
    Lkern <- Lkern / (maxiter - floor(maxiter / 5) + 1)
    CountLbest <- N ^ 2 + 1
    for (w in floor(maxiter / 5):maxiter) {
      CountL <- norm(Lkern - Llist[[m]][[w]] %*% t(Llist[[m]][[w]]), "F") ^ 2
      if (CountL < CountLbest) {
        Lbest[[m]] <- Llist[[m]][[w]]
        CountLbest <- CountL
      }
    }
  }
  
  # Choose overall clustering
  Ckern <- Clist[[floor(maxiter / 5)]] %*% t(Clist[[floor(maxiter / 5)]])
  for (w in floor(maxiter / 5 + 1):maxiter) {
    Ckern <- Ckern + Clist[[w]] %*% t(Clist[[w]])
  }
  Ckern <- Ckern / (maxiter - floor(maxiter / 5) + 1)
  CountCbest <- N ^ 2 + 1
  for (w in floor(maxiter / 5):maxiter) {
    CountC <- norm(Ckern - Clist[[w]] %*% t(Clist[[w]]), "F") ^ 2
    if (CountC < CountCbest) {
      Cbest <- Clist[[w]]
      CountCbest <- CountC
    }
  }
  return(list(Alpha=Alpha, AlphaBounds=AlphaBounds, Cbest=Cbest, 
              Lbest=Lbest, AlphaVec=alphaVec))
}
