########## Bricks for hDCBM ############

renumber <- function(z) {
  z_re <- rep(0, length(z))
  uniq <- unique(z)

  for (i in 1:length(uniq)) {
    z_re[z == uniq[i]] <- i
  }
  return(z_re)
}

updateS <- function(Y, z) {
  r <- length(unique(z))
  S <- array(0, dim = c(r, r, r))

  for (a in 1:r) {
    for (b in 1:r) {
      for (c in 1:r) {
        S[a, b, c] <- mean(Y[z == a, z == b, z == c], na.rm = T)
      }
    }
  }
  return(S)
}

Cal_Y1 <- function(Y, z) {
  r <- length(unique(z))
  p <- dim(Y)[1]

  Y1 <- array(0, dim = c(p, r, r))
  for (b in 1:r) {
    for (c in 1:r) {
      Y1[, b, c] <- apply(Y[, z == b, z == c], 1, mean)
    }
  }
  return(Y1)
}

cosin <- function(v1, v2) {
  v1_norm <- sqrt(sum(v1^2))
  v2_norm <- sqrt(sum(v2^2))
  if (v1_norm == 0 | v2_norm == 0) {
    return(cosin = 1)
  } else {
    return(cosin = t(v1) %*% v2 / (v1_norm * v2_norm))
  }
}

updatez <- function(Y1_unfold, S_unfold) {
  z_new <- rep(0, dim(Y1_unfold)[1])

  for (j in 1:length(z_new)) {
    dist <- c()
    for (a in 1:dim(S_unfold)[1]) {
      dist <- c(dist, cosin(Y1_unfold[j, ], S_unfold[a, ]))
    }
    z_new[j] <- which.max(dist)
  }

  return(z_new)
}

# pos_randothrox <- function(n, r) {
#   mat <- array(0, dim = c(n, r))
#   sc <- 1:n
# 
#   for (i in 1:r) {
#     mat[-sc, i] <- 0
#     mat[sc, i] <- rbinom(length(sc), 1, prob = n / ((r + 0.5) * length(sc)))
# 
#     sc <- setdiff(sc, which(mat[, i] == 1))
#   }
# 
#   for (i in 1:r) {
#     mat[, i] <- mat[, i] / sqrt(sum(mat[, i]^2))
#   }
# 
#   return(mat)
# }

sim_S_network <- function(r, s_min, s_max) {
  # assortative s
  S <- array(s_min, dim = c(r, r, r))
  for (i in 1:r) {
    S[i, i, i] <- s_max
  }
  return(S)
}


sim_S_magnitude <- function(r, delta) {
  # assortative s
  S <- array(0.5, dim = c(r, r, r))
  for (i in 1:r) {
    S[i, i, i] <- 1
  }
  
  delta_min = sqrt( sum( (S[1,,] - S[2,,])^2 ))
  S = S * delta / delta_min
  
  return(S)
}
