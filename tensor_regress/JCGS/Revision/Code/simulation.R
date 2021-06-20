# This file collect the bricks for the simulations
library(pracma)

# massive GLM to estimate the covariate effects
massive_glm <- function(response, X, dist) {
  d1 <- dim(response)[1]
  d2 <- dim(response)[2]
  coe <- array(0, dim = c(d1, d2, dim(X)[2]))
  if (dist == "normal") {
    for (i in 1:d1) {
      for (j in 1:d2) {
        fit <- lm(response[i, j, ] ~ -1 + X)
        coe[i, j, ] <- coef(fit)
      }
    }
  }
  else if (dist == "binary") {
    for (i in 1:d1) {
      for (j in 1:d2) {
        fit <- suppressWarnings(glm(response[i, j, ] ~ -1 + X, family = binomial("logit")))
        coe[i, j, ] <- coef(fit)
      }
    }
  }
  else if (dist == "poisson") {
    for (i in 1:d1) {
      for (j in 1:d2) {
        fit <- suppressWarnings(glm(response[i, j, ] ~ -1 + X, family = poisson("log")))
        coe[i, j, ] <- coef(fit)
      }
    }
  }
  return(coe)
}

# valid rank range for core tensor ## the maximum rank should be no smaller than the product of other
# Input: rank_range = expand.grid(...)
valid_rank <- function(rank_range) {
  N <- dim(rank_range)[1]
  output <- NULL
  for (i in 1:N) {
    rank <- rank_range[i, ]
    rank <- sort(rank)
    if (rank[1] * rank[2] >= rank[3]) {
      output <- rbind(output, rank_range[i, ])
    }
  }
  return(output)
}

# compute loglikelihood under different models
loglike <- function(data, linearpre, dist) {
  if (dist == "binary") {
    p <- plogis(linearpre) # log-likelihood
    L <- dbinom(data, 1, p, log = TRUE)
    return(sum(L))
  }
  else if (dist == "normal") {
    sigma_est <- mean((data - linearpre)^2)
    L <- dnorm(data, linearpre, sqrt(sigma_est), log = TRUE)
    sum(L) # log-likelihood
  }
  else if (dist == "poisson") {
    lambda <- exp(linearpre)
    L <- dpois(data, lambda, log = TRUE) # log-likelihood for Poisson
    return(sum(L))
  }
}

# wrapped simulation function to estimate covariate effects using tensor methods or massive glm methods
conv_rate <- function(seed, signal = 10, niter = 20, cons = "non", lambda = 1, alpha = 10, solver = "GC", c_range, dist = "binary", dup = 10, d_range, p_range, naive = TRUE, initial) {
  n <- m <- nrow(d_range)
  s <- nrow(c_range)
  error <- error_naive <- array(0, dim = c(n, s, dup))
  for (i in 1:n) {
    for (k in 1:s) {
      data <- sim_data(seed, d_range[i, ], c_range[k, ], p_range[i, ], dist, dup, signal)
      X_covar1 <- data$X_covar1
      X_covar2 <- data$X_covar2
      X_covar3 <- data$X_covar3
      for (l in 1:dup) {
        result <- tensor_regress(data$tsr[[l]], X_covar1, X_covar2, X_covar3, c_range[k, ], niter, cons, lambda, alpha, solver, dist, initial = initial)
        if (naive == TRUE) {
          naive_C <- massive_glm(data$tsr[[l]], X_covar3, dist)
          error_naive[i, k, l] <- mean((naive_C - data$C_ts)^2)
        }
        error[i, k, l] <- mean((result$C_ts - data$C_ts)^2) ## use mean not the total F-norm because X has been rescaled. See sim_data(...)
        print(paste(l, "-th replicate---- when dimension is ", d_range[i, 1], d_range[i, 2], d_range[i, 3], "-- covariate is", p_range[i, 1], p_range[i, 2], p_range[i, 3], " --------core is ", c_range[k, 1], c_range[k, 2], c_range[k, 3]))
      }
    }
  }
  return(list(error, error_naive))
}


# This function is to select the lambda in the penalty constrain version
sele_lambda <- function(seed, lambda, ...) {
  re <- lapply(lambda, FUN = conv_rate, seed = seed, ...)
  re <- lapply(seq(length(re)), function(x) re[[x]]$RMSE)
  return(re)
}

# This function transports the indices in tensor
tensor_trans <- function(tsr) {
  # Y \in R^{m x l x n} --> Y \in R^{n x m x l}

  mode <- dim(tsr)
  new_tsr <- array(0, c(mode[3], mode[1:2]))
  for (i in 1:mode[3]) {
    new_tsr[i, , ] <- tsr[, , i]
  }
  return(new_tsr)
}


# simulate data from modified GLSNet model
sparse_sim <- function(d, d3, p, R, r3, s, signal, dup) {

  # mode 3 covariates
  X <- cbind(rep(1, d3), randortho(d3)[, 1:p])

  # intercept matrix
  A <- matrix(rnorm(d * R), ncol = R, nrow = d)
  B0 <- A %*% t(A)

  G <- array(runif(R^2 * r3, min = -1, max = 1), dim = c(R, R, r3))
  W1 <- randortho(d)[, 1:R]
  W2 <- randortho(d)[, 1:R]
  W3 <- randortho(p)[, 1:r3]

  B <- array(0, dim = c(d, d, (p + 1)))
  B[, , 1] <- B0

  B1 <- ttl(as.tensor(G), list(W1, W2, W3), c(1, 2, 3))@data

  U1 <- ttm(as.tensor(B1), X[, 2:(p + 1)], 3)@data

  # rescale subject to entrywise constraint
  G <- G / max(abs(U1)) * signal
  U1 <- U1 / max(abs(U1)) * signal

  B1 <- ttl(as.tensor(G), list(W1, W2, W3), ms = c(1, 2, 3))@data

  a <- 1:(d^2 * p)
  b <- round(((1 - s) * d^2 * p))
  sparse_index <- sample(a, b, replace = FALSE)
  B1[sparse_index] <- 0

  B[, , 2:(p + 1)] <- B1

  U <- ttm(as.tensor(B), X, 3)@data

  tsr <- lapply(seq(dup), function(x) array(rbinom(d * d * d3, 1, prob = as.vector(1 / (1 + exp(-U)))), dim = c(d, d, d3)))

  return(list(tsr = tsr, X = X, U = U, B = B, G = G))
}



# simulate data from modified Envelope model
env_sim <- function(D, u, p, n, signal_level, cor_level, dup) {
  C <- as.tensor(array(runif(u[1] * u[2] * p, min = -signal_level, max = signal_level), dim = c(u[1], u[2], p))) # core tensor

  gamma <- list()
  # orthogonal factor matrices
  gamma0 <- list()
  Omega <- list()
  # list of Omega
  Omega0 <- list()
  # list of covariance matrices on different modes, Sigma = Gamma Omega Gamma^T
  Sigma <- list()
  # Sigma^{1/2}
  Sigmasqrtm <- list()

  for (i in 1:2) {
    d <- D[i]
    r <- u[i]

    g_othro <- randortho(d, "orthonormal")
    gamma[[i]] <- g_othro[, 1:r]
    gamma0[[i]] <- g_othro[, (r + 1):d]

    A <- array(runif(r * r, min = -cor_level, max = cor_level), dim = c(r, r))
    Omega[[i]] <- A %*% t(A)

    A <- array(runif((d - r) * (d - r), min = -cor_level, max = cor_level), dim = c(d - r, d - r))
    Omega0[[i]] <- A %*% t(A)

    Sigma[[i]] <- gamma[[i]] %*% Omega[[i]] %*% t(gamma[[i]]) + gamma0[[i]] %*% Omega0[[i]] %*% t(gamma0[[i]])
    Sigma[[i]] <- round(Sigma[[i]], digits = 3) + diag(rep(1, d)) # sigma is non identity matrix

    Sigmasqrtm[[i]] <- sqrtm(Sigma[[i]])$B
  }

  B <- ttl(C, gamma, c(1, 2)) # coefficient tensor
  X <- array(rnorm(n * p), dim = c(n, p)) # covariates
  U <- ttm(B, X, 3)

  Y <- list()
  for (k in 1:dup) {
    # generate noise
    Epsilon <- as.tensor(array(rnorm(D[1] * D[2] * n), dim = c(D[1], D[2], n)))
    Epsilon <- ttl(Epsilon, Sigmasqrtm, c(1, 2)) # tE = Sigma \otimes Sigma

    Y[[k]] <- Epsilon + ttm(B, X, 3)
  }

  return(list(tsr = Y, U = U, X = X, B = B, Sigma = Sigma, gamma = gamma))
}
