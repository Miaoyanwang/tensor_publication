#' Weighted higher-order initialization
#'
#' Weighted higher-order initialization for multiway spherical clustering under degree-corrected tensor block model.
#' This function takes the tensor/matrix observation, the number of clusters, and a logic variable indicating the symmetry
#' of the tensor as input. The output is the estimated clustering assignment.
#'
#'
#' @param Y     order-3 array or matrix, order-3 tensor or matrix observation
#' @param r     vector, the number of clusters on each model; see "details"
#' @param asymm logic variable, if "TRUE", assume the clustering assignment differs in different modes; if "FALSE", assume all the modes share the same clustering assignment
#' @return a list containing the following:
#'
#' \code{z0} { vector ( if \code{asymm = F} ) or a list of vectors ( if \code{asymm = T} ) recording the estimated clustering assignment }
#'
#' \code{s0} { vector ( if \code{asymm = F} ) or a list of vectors ( if \code{asymm = T} ) recording the index of degenerate entities with random clustering assignment}
#'
#' @details   For the number of clusters \code{r}, the length of \code{r} is no longer than 3 with integers larger than 1;
#'            the input \code{r = r0} is equal to input \code{r = rep(r0,3)} for any integer \code{r0 > 1};
#'            matrix observation only accounts the frist two elements and requires all the elements in \code{r} are identicial.
#'
#' @export
#' @examples
#' seed = 1
#' p = 30
#' r = 3
#' delta = 0.5
#' s_min = 0.05
#' s_max = 1
#' dist = "normal"
#' theta_dist = "pareto"
#' alpha = 4
#' beta = 3/4
#' sigma = 0.2
#'
#' data = sim_hDCBM_network(seed = seed, p, r,  delta = delta,
#'  s_min = s_min, s_max = s_max, dist =  dist, sigma = sigma,
#'  theta_dist = theta_dist, alpha = alpha, beta = beta, imat = F)
#'
#' initialization = wkmeans(data$Y,r = 3, asymm = F)



wkmeans <- function(Y, r, asymm = F) {
  imat <- F

  if (length(dim(Y)) == 2) {
    cat("matrix case \n")
    dim(Y) <- c(dim(Y), 1)
    imat <- T
  }

  # r should be a vector
  if (length(r) == 1) {
    r1 <- r2 <- r3 <- r
    if (imat == T) {
      r3 <- 1
    }
  } else if (length(r) > 1) {
    r1 <- r[1]
    r2 <- r[2]

    if (imat == T) {
      r3 <- 1

      if (r1 != r2) {
        warning("matrix case requires the same number of clusters on two modes, \n here use the minimal element in r.")
        r2 <- r1 <- min(r1, r2)
      }
    } else if (imat == F) {
      r3 <- r[3]

      if (asymm == F) {
        if (r1 != r2 | r2 != r3 | r1 != r3) {
          warning("symmetric case requires the same number of clusters on every mode, \n here use the minimal element in r.")
        }
        r3 <- r2 <- r1 <- min(c(r1, r2, r3))
      }
    }
  }

  ### two-step SVD
  Y <- as.tensor(Y)

  # first SVD
  u1 <- svd(unfold(Y, 1, c(3, 2))@data)$u[, 1:r1]
  u2 <- svd(unfold(Y, 2, c(1, 3))@data)$u[, 1:r2]
  if (imat == F) {
    u3 <- svd(unfold(Y, 3, c(1, 2))@data)$u[, 1:r3]
  } else if (imat == T) {
    u3 <- as.matrix(1)
  }

  # second SVD

  hu1 <- svd(unfold(ttl(Y, list(t(u2), t(u3)), ms = c(2, 3)), 1, c(3, 2))@data)$u[, 1:r1]
  hu2 <- svd(unfold(ttl(Y, list(t(u1), t(u3)), ms = c(1, 3)), 2, c(1, 3))@data)$u[, 1:r2]
  if (imat == F) {
    hu3 <- svd(unfold(ttl(Y, list(t(u1), t(u2)), ms = c(1, 2)), 3, c(1, 2))@data)$u[, 1:r3]
  } else if (imat == T) {
    hu3 <- as.matrix(1)
  }

  ### estimated X
  X1 <- hu1 %*% t(hu1) %*% unfold(ttl(Y, list(hu2 %*% t(hu2), hu3 %*% t(hu3)), ms = c(2, 3)), 1, c(3, 2))@data
  X2 <- hu2 %*% t(hu2) %*% unfold(ttl(Y, list(hu1 %*% t(hu1), hu3 %*% t(hu3)), ms = c(1, 3)), 2, c(1, 3))@data
  X3 <- hu3 %*% t(hu3) %*% unfold(ttl(Y, list(hu1 %*% t(hu1), hu2 %*% t(hu2)), ms = c(1, 2)), 3, c(1, 2))@data

  if (asymm == T) {
    res1 <- single_wkmeans(X1, r1)
    res2 <- single_wkmeans(X2, r2)
    if (imat == F) {
      res3 <- single_wkmeans(X3, r3)
    } else if (imat == T) {
      res3 <- list(z = 1, s0 = 0)
    }

    return(list(
      z0 = list(z1 = as.numeric(res1$z), z2 = as.numeric(res2$z), z3 = as.numeric(res3$z)),
      s0 = list(s01 = res1$s0, s02 = res2$s0, s03 = res3$s0)
    ))
  } else if (asymm == F) {
    z <- rep(0, dim(Y)[1])
    sc <- 1:length(z)
    ### normalization
    l2 <- apply(X1, 1, function(x) sqrt(sum(x^2))) # l2 norm

    if (length(which(l2 == 0) > 0)) {
      sc <- which(l2 != 0)
      X1 <- X1[sc, ]
      l2 <- l2[sc]
    }

    Xs <- diag(l2^(-1)) %*% X1

    ## weighted kmeans
    diss <- dist(Xs, method = "euclidean", p = 2, upper = T, diag = T)

    z[sc] <- wcKMedoids(diss^2, k = r1, weights = l2^2, method = "PAMonce", cluster.only = T)
    z[-sc] <- sample(unique(z[sc]), length(z[-sc]), replace = T)


    s0 <- setdiff(1:length(z), sc)
    z <- as.factor(z)
    levels(z) <- 1:r1

    return(list(z0 = as.numeric(z), s0 = s0))
  }
}


#' Angle-based iteration
#'
#' Angle-based iteration for multiway spherical clustering under degree-corrected tensor block model.
#' This function takes the tensor/matrix observation, initial clustering assignment, and a logic variable indicating the symmetry
#' of the tensor as input. The output is the refined clustering assignment.
#'
#'
#'@param Y         order-3 array or matrix, order-3 tensor or matrix observation
#'@param z0        vector or list of vectors, initial clustering assignment; see "details"
#'@param max_iter  integer, max number of iterations if update does not converge
#'@param alpha1    number, substitution of degenerate core tensor; see "details"
#'@param asymm     logic variable, if "TRUE", assume the clustering assignment differs in different modes; if "FALSE", assume all the modes share the same clustering assignment
#'
#'@return a list containing the following:
#'
#'\code{z} {vector ( if \code{asymm = F} ) or a list of vectors ( if \code{asymm = T} ) recording the estimated clustering assignment}
#'
#'\code{s_deg} {logic variable, if "TRUE", degenerate estimated core tensor occurs during the iteration; if "FALSE", otherwise}
#'
#'@details For initial clustering assignment \code{z0}, if \code{z0} is a vector, assume all the modes share the same initialization \code{z0},
#'         if \code{z0} is a vector, \code{z0} specifies the initialization on each mode; matrix observation only accounts the first two vectors in the list.
#'         When the estimated core tensor has a degenerate slice, i.e., a slice with all zero elements, randomly pick an entry in the degenerate slice with value \code{alpha1}.
#'
#'@export
#'@examples
#' seed = 1
#' p = 30
#' r = 3
#' delta = 0.5
#' s_min = 0.05
#' s_max = 1
#' dist = "normal"
#' theta_dist = "pareto"
#' alpha = 4
#' beta = 3/4
#' sigma = 0.2
#'
#' data = sim_hDCBM_network(seed = seed, p, r,  delta = delta,
#'  s_min = s_min, s_max = s_max, dist =  dist, sigma = sigma,
#'  theta_dist = theta_dist, alpha = alpha, beta = beta, imat = F)
#'
#' initialization = wkmeans(test_dat$Y,r = 3, asymm = F)
#'
#' iteration = hALloyd(data$Y, initialization$z,max_iter = 20,asymm = F)


hALloyd <- function(Y, z0, max_iter, alpha1 = 0.01, asymm) {
  imat <- F
  s_deg <- F

  # z0 should be a list
  if (length(dim(Y)) == 2) {
    cat("matrix case \n")
    dim(Y) <- c(dim(Y), 1)
    imat <- T
  }

  if (is.list(z0) == T) {
    z <- lapply(z0, renumber)
    if (imat == T) {
      z$z3 <- 1
    }
  } else if (is.list(z0) == F) {
    if (imat == F) {
      z <- list(z1 = renumber(z0), z2 = renumber(z0), z3 = renumber(z0))
    } else if (imat == T) {
      z <- list(z1 = renumber(z0), z2 = renumber(z0), z3 = 1)
    }
  }


  if (asymm == T) {
    for (iter in 1:max_iter) {
      cat("iter = ", iter, "\n")

      # estimate S
      est_S <- updateS(Y, z, imat)

      # update z1
      # reducde Y
      Y1 <- Cal_Y1(Y, z)
      re1 <- single_Aiteration(unfold(as.tensor(est_S), 1, c(3, 2))@data, unfold(as.tensor(Y1), 1, c(3, 2))@data, alpha1)
      z1_new = re1$z
      z1_new <- renumber(z1_new)


      # update z2
      # reducde Y
      Y2 <- Cal_Y2(Y, z)
      re2 <- single_Aiteration(unfold(as.tensor(est_S), 2, c(1, 3))@data, unfold(as.tensor(Y2), 2, c(1, 3))@data, alpha1)
      z2_new = re2$z
      z2_new <- renumber(z2_new)


      # update z3
      # reducde Y
      if (imat == T) {
        z3_new <- 1

        if(re1$s_deg == T | re2$s_deg == T){
          s_deg = T
        }

      } else if (imat == F) {
        Y3 <- Cal_Y3(Y, z)
        re3 <- single_Aiteration(unfold(as.tensor(est_S), 3, c(1, 2))@data, unfold(as.tensor(Y3), 3, c(1, 2))@data, alpha1)
        z3_new = re3$z
        z3_new <- renumber(z3_new)

        if(re1$s_deg == T | re2$s_deg == T|re3$s_deg == T){
          s_deg = T
        }

      }

      z_new_list <- list(z1 = z1_new, z2 = z2_new, z3 = z3_new)

      if (identical(z_new_list, z)) {
        break
      }

      z <- z_new_list
    }

    return(list(z = z, s_deg = s_deg))
  } else if (asymm == F) {
    z <- z[[1]]
    for (iter in 1:max_iter) {
      cat("iter = ", iter, "\n")

      S_unfold <- unfold(as.tensor(updateS_old(Y, z, imat)), 1, c(3, 2))@data
      lS <- apply(S_unfold, 1, function(x) sum(x^2))
      index <- which(lS == 0)
      if (length(index) > 0) {
        for (i in index) {
          S_unfold[i, sample(1:dim(S_unfold)[2],1)] <- alpha1
          s_deg <- T
        }
      }

      z_new <- rep(0, length(z))
      sc <- 1:length(z)


      # calculate Y1 for clustering
      Y1_unfold <- unfold(as.tensor(Cal_Y1_old(Y, z)), 1, c(3, 2))@data

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

    return(list(z = z, s_deg = s_deg))
  }
}
