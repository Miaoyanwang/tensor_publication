# Code for Table 3. Accuracy of BIC selection for normal and Poisson models under different (d,r) combinations

# load dependencies ---------------------------------------------------------------------------------------------
# for simulation
install.packages("software/tensorregress.tar.gz", repos = NULL, type = "source")
library(tensorregress)

library(pracma)
source("simulation.R")

# report simulation statistics from prior simulation -------------------------------------------------------------
# if a new run of simulation is desired, please go to line #58 to run the code ----------------------------------

# report rank estimates under the normal model
# load pre-saved data
load("presaved/Table3_normal.RData")

# d = 20
apply(final[1, 1, , ], 2, mean) # ground truth (3,3,3)
apply(final[2, 1, , ], 2, mean) # ground truth (4,4,6)
apply(final[3, 1, , ], 2, mean) # ground truth (6,8,8)
apply(final[1, 1, , ], 2, sd) # ground truth (3,3,3)
apply(final[2, 1, , ], 2, sd) # ground truth (4,4,6)
apply(final[3, 1, , ], 2, sd) # ground truth (6,8,8)

# d = 40
apply(final[1, 2, , ], 2, mean) # sd for estimate under ground truth (3,3,3)
apply(final[2, 2, , ], 2, mean) # sd for estimate under ground truth (4,4,6)
apply(final[3, 2, , ], 2, mean) # sd for estimate under ground truth (6,8,8)
apply(final[1, 2, , ], 2, sd) # sd for estimate under ground truth (3,3,3)
apply(final[2, 2, , ], 2, sd) # sd for estimate under ground truth (4,4,6)
apply(final[3, 2, , ], 2, sd) # sd for estimate under ground truth (6,8,8)

# report rank estimates under poisson model
# load pre-saved data
load("presaved/Table3_poisson.RData")

# d = 20
apply(final[1, 1, , ], 2, mean) # ground truth (3,3,3)
apply(final[2, 1, , ], 2, mean) # ground truth (4,4,6)
apply(final[3, 1, , ], 2, mean) # ground truth (6,8,8)
apply(final[1, 1, , ], 2, sd) # ground truth (3,3,3)
apply(final[2, 1, , ], 2, sd) # ground truth (4,4,6)
apply(final[3, 1, , ], 2, sd) # ground truth (6,8,8)

# d = 40
apply(final[1, 2, , ], 2, mean) # sd for estimate under ground truth (3,3,3)
apply(final[2, 2, , ], 2, mean) # sd for estimate under ground truth (4,4,6)
apply(final[3, 2, , ], 2, mean) # sd for estimate under ground truth (6,8,8)
apply(final[1, 2, , ], 2, sd) # sd for estimate under ground truth (3,3,3)
apply(final[2, 2, , ], 2, sd) # sd for estimate under ground truth (4,4,6)
apply(final[3, 2, , ], 2, sd) # sd for estimate under ground truth (6,8,8)

# end of report -------------------------------------------------------------------------------------------------



# If a new run of simulation is desired, please run the code from here and save the results as .RData -----------
# Then run the above code to generate numbers in table ----------------------------------------------------------

set.seed(0) # set seed for this simulation
seed <- 0

# all covariates included 
d <- c(20, 40)
rank_range <- array(0, dim = c(3, 7, 3))
rank_range[1, , ] <- rbind(c(3, 3, 3), c(3, 3, 2), c(3, 2, 2), c(2, 2, 2), c(3, 3, 4), c(3, 4, 4), c(4, 4, 4))
rank_range[2, , ] <- rbind(c(4, 4, 6), c(4, 4, 5), c(3, 3, 6), c(4, 5, 7), c(3, 4, 6), c(4, 5, 6), c(5, 5, 5))
rank_range[3, , ] <- rbind(c(6, 8, 8), c(6, 7, 8), c(5, 8, 8), c(7, 8, 8), c(5, 7, 7), c(5, 7, 8), c(6, 7, 7))
core_shape <- rbind(c(3, 3, 3), c(4, 4, 6), c(6, 8, 8))
dup <- 30


# Normal model 
dist <- "normal"
final <- array(0, dim = c(3, 2, dup, 3))
for (s in 1:3) { # three possible ranks
  for (i in 1:2) { # two possible dimensions
    whole_shape <- c(d[i], d[i], d[i])
    p <- 0.4 * whole_shape
    data <- sim_data(seed, whole_shape = whole_shape, core_shape = core_shape[s, ], p = p, dist = dist, dup = dup, signal = 4)
    table <- NULL
    for (j in 1:10) { # simulation replicates
      res <- sele_rank(data$tsr[[j]], X_covar1 = data$X_covar1, X_covar2 = data$X_covar2, X_covar3 = data$X_covar3, rank_range = rank_range[s, , ], niter = 10, cons = "non", dist = dist, initial = "random")
      table <- rbind(table, res$rank)
    }
    final[s, i, , ] <- table
  }
}

save(final, file = "presaved/Table3_normal.RData")


# Poisson model 
set.seed(0) ## set seed for this simulation
seed <- 0
dist <- "poisson"
final <- array(0, dim = c(3, 2, dup, 3))
for (s in 1:3) { # three possible ranks
  for (i in 1:2) { # two possible dimensions
    whole_shape <- c(d[i], d[i], d[i])
    p <- 0.4 * whole_shape
    data <- sim_data(seed, whole_shape = whole_shape, core_shape = core_shape[s, ], p = p, dist = dist, dup = dup, signal = 4)
    table <- NULL
    for (j in 1:dup) { # simulation replicates
      res <- sele_rank(data$tsr[[j]], X_covar1 = data$X_covar1, X_covar2 = data$X_covar2, X_covar3 = data$X_covar3, rank_range = rank_range[s, , ], niter = 10, cons = "non", dist = dist, initial = "random")
      table <- rbind(table, res$rank)
    }
    final[s, i, , ] <- table
  }
}
save(final, file = "presaved/Table3_poisson.RData")
