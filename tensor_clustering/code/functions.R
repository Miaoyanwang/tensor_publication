#### main functions for hDCBM ####
# Sep 2. Jiaxin

library(rTensor)
library(tensorregress)
library(WeightedCluster)
library(EnvStats)
library(plot.matrix)
source("bricks.R")


################### Initialization ##################
# two-step SVD + weighted k-median
# tensor_unfold(,1) = unfold(,1,c(3,2))
# tensor_unfold(,2) = unfold(,2,c(1,3))
# tensor_unfold(,3) = unfold(,3,c(1,2))

wkmedian <- function(Y, r) {
  z <- rep(0, dim(Y)[1])
  sc <- 1:length(z)

  ### two-step SVD
  Y <- as.tensor(Y)
  # first SVD
  u1 <- svd(unfold(Y, 1, c(3, 2))@data)$u[, 1:r]
  u2 <- svd(unfold(Y, 2, c(1, 3))@data)$u[, 1:r]
  u3 <- svd(unfold(Y, 3, c(1, 2))@data)$u[, 1:r]
  # second SVD
  hu1 <- svd(unfold(ttl(Y, list(t(u2), t(u3)), ms = c(2, 3)), 1, c(3, 2))@data)$u[, 1:r]
  hu2 <- svd(unfold(ttl(Y, list(t(u1), t(u3)), ms = c(1, 3)), 2, c(1, 3))@data)$u[, 1:r]
  hu3 <- svd(unfold(ttl(Y, list(t(u1), t(u2)), ms = c(1, 2)), 3, c(1, 2))@data)$u[, 1:r]

  ### estimated X
  X1 <- hu1 %*% t(hu1) %*% unfold(ttl(Y, list(hu2 %*% t(hu2), hu3 %*% t(hu3)), ms = c(2, 3)), 1, c(3, 2))@data

  ### normalization
  l1 <- apply(X1, 1, function(x) sum(abs(x)))

  if (length(which(l1 == 0) > 0)) {
    sc <- which(l1 != 0)
    X1 <- X1[sc, ]
    l1 <- l1[sc]
  }

  Xs <- diag(l1^(-1)) %*% X1

  ### weighted k-median
  diss <- dist(Xs, method = "minkowski", p = 1, upper = T, diag = T)

  z[sc] <- wcKMedoids(diss, k = r, weights = l1, method = "PAMonce", cluster.only = T)
  z[-sc] <- sample(unique(z[sc]), length(z[-sc]), replace = T)

  s0 <- setdiff(1:length(z), sc)
  z <- as.factor(z)
  levels(z) <- 1:r

  return(list(z0 = as.numeric(z), s0 = s0))
}


################### Refinement ##################
# estimate S + angle distance cluster

hALloyd <- function(Y, z0, max_iter, alpha1 = 0.1) {
  z <- renumber(z0) 
  for (iter in 1:max_iter) {
    cat("iter = ", iter, "\n")

    z_new <- rep(0, length(z))
    sc <- 1:length(z)

    # estimate S
    S_unfold <- unfold(as.tensor(updateS(Y, z)), 1, c(3, 2))@data
    lS <- apply(S_unfold, 1, function(x) sum(x^2))
    index <- which(lS == 0)
    if (length(index) > 0) {
      for (i in index) {
        S_unfold[i, i] <- alpha1
      }
    }

    # calculate Y1 for clustering
    Y1_unfold <- unfold(as.tensor(Cal_Y1(Y, z)), 1, c(3, 2))@data

    l2 <- apply(Y1_unfold, 1, function(x) sum(x^2))
    if (length(which(l2 == 0) > 0)) {
      sc <- which(l2 != 0)
      Y1_unfold <- Y1_unfold[sc, ]
    }

    # update cluster via anlge distance
    z_new[sc] <- updatez(Y1_unfold, S_unfold)
    z_new[-sc] <- sample(unique(z_new[sc]), length(z_new[-sc]), replace = T)

    z_new <- renumber(z_new)

    if (identical(z_new, z)) {
      break
    }

    z <- z_new
  }

  return(z)
}


################### Generator ##################

sim_hDCBM_network <- function(seed = NA, p, r, core_conrtol = c("signal_m","signal_a2"),
                              delta = NULL, s_min = NULL, s_max = NULL, dist = c("normal", "binary"), sigma = 1,
                              theta_dist = c("abs_normal", "pareto", "non"), alpha = NULL, beta = NULL) {
  # cospi(angle) = cos(angle*pi)
  # s_min, s_max control the angle, i.e., signal
  # usually, s_max = c s_min for some c>1 and s_min > 0.


  if (is.na(seed) == FALSE) set.seed(seed)

  # generate degree heterogenity
  if (theta_dist == "abs_normal") {
    theta <- abs(rnorm(p, 0, 0.5)) + 1 - sqrt(1 / (2 * pi))
  } else if (theta_dist == "pareto") {
    theta <- rpareto(p, location = beta, shape = alpha) # choose alpha*beta = alpha - 1
  } else if (theta_dist == "non"){
    theta = rep(1, p)
  }

  # generate S with control signal
  if(core_conrtol == "signal_a2"){
    S <- sim_S_network(r, s_min, s_max)
  }else if(core_conrtol == "signal_m"){
    S <- sim_S_magnitude(r, delta)
  }
  

  # generate assignment
  z <- sample(1:r, p, replace = T)
  z <- renumber(z)

  # generate mean tensor
  X <- ttl(as.tensor(S[z, z, z]), list(diag(theta), diag(theta), diag(theta)), ms = c(1, 2, 3))@data

  # generate data
  if (dist == "normal") {
    Y <- X + array(rnorm(p^3, 0, sigma), dim = c(p, p, p))
  } else if (dist == "binary") {
    X[which(X > 1)] <- 1
    Y <- array(rbinom(p^3, 1, as.vector(X)), dim = c(p, p, p))
  }

  return(list(Y = Y, X = X, S = S, theta = theta, z = z))
}


