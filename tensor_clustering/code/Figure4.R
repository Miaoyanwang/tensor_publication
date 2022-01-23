# Code for Figure 4. Verification of SNR phase transition

# load dependencies ---------------------------------------------------------------------------------------------
install.packages("software/dTBM.tar.gz", repos = NULL, type = "source")
library(dTBM)
library(RSKC)

# for plotting
library(ggplot2)
library(patchwork)


# start plotting with saved data from prior simulation ---------------------------------------------------------------
# if a new run of simulation is desired, please go to line # to run the code --------------------------------------

# tensor case
load("presaved/Figure4_tensor100.RData")

mis1 <- mis
mis_sd1 <- mis_sd

load("presaved/Figure4_tensor80.RData")

mis2 <- mis
mis_sd2 <- mis_sd

sel_color <- c("#78A153", "#E4930A")

gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4)
data <- data.frame(
  MCR = c(mis1[1, ], mis1[2, ], mis2[1, ], mis2[2, ]), sd = c(mis_sd1[1, ], mis_sd1[2, ], mis_sd2[1, ], mis_sd2[2, ]) / sqrt(30),
  Algorithm = rep(c("Ours", "Oracle"), each = 8, times = 2), gamma = rep(gamma_range, times = 4),
  p = rep(c(100, 80), each = 16)
)
data[, 3] <- factor(data[, 3], levels = c("Ours", "Oracle"))
data[, 5] <- factor(data[, 5], levels = c(100, 80))
p1 <- ggplot(data = data, aes(x = gamma, y = MCR)) +
  geom_line(aes(color = Algorithm, linetype = p), size = 1.2) +
  scale_colour_manual(values = sel_color) +
  geom_point(size = 3, aes(shape = Algorithm)) +
  scale_shape_manual(values = c(16, 5)) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  labs(title = "Tensor (K = 3)") +
  xlab("Gamma") +
  ylab("Clustering Error Rate (CER)") +
  coord_cartesian(ylim = c(0, 0.35))
p1 <- p1 + geom_errorbar(aes(ymin = MCR - sd, ymax = MCR + sd), width = 0.03, position = position_dodge(0.05)) + theme_light() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 14))


# matrix case
load("presaved/Figure4_matrix100.RData")

mis1 <- mis
mis_sd1 <- mis_sd

load("presaved/Figure4_matrix80.RData")

mis2 <- mis
mis_sd2 <- mis_sd

sel_color <- c("#78A153", "#E4930A")

gamma_range <- c(-1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4)
data <- data.frame(
  MCR = c(mis1[1, ], mis1[2, ], mis2[1, ], mis2[2, ]), sd = c(mis_sd1[1, ], mis_sd1[2, ], mis_sd2[1, ], mis_sd2[2, ]) / sqrt(30),
  Algorithm = rep(c("Ours", "Oracle"), each = 9, times = 2), gamma = rep(gamma_range, times = 4),
  p = rep(c(100, 80), each = 18)
)
data[, 3] <- factor(data[, 3], levels = c("Ours", "Oracle"))
data[, 5] <- factor(data[, 5], levels = c(100, 80))
p2 <- ggplot(data = data, aes(x = gamma, y = MCR)) +
  geom_line(aes(color = Algorithm, linetype = p), size = 1.2) +
  scale_colour_manual(values = sel_color) +
  geom_point(size = 3, aes(shape = Algorithm)) +
  scale_shape_manual(values = c(16, 5)) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  labs(title = "Matrix (K = 2)") +
  xlab("Gamma") +
  ylab("Clustering Error Rate (CER)") +
  coord_cartesian(ylim = c(0, 0.35))
p2 <- p2 + geom_errorbar(aes(ymin = MCR - sd, ymax = MCR + sd), width = 0.03, position = position_dodge(0.05)) + theme_light() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
  theme(legend.position = "none")

pdf("figure/Figure4.pdf", height = 4, width = 10)
(p2 | p1)
dev.off()



# end of plotting -----------------------------------------------------------------------------------------------



# If a new run of simulation is desired, please run the code from here and save the results as .RData -----------
# Then run the above code to generate figures -------------------------------------------------------------------

# tensor case with p = 100
p <- 100
r <- 5

gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4)
alpha_range <- 1 + 180 * p^(gamma_range / 2)

s_min <- 0.5
delta <- 0.5

dist <- "normal"
theta_dist <- "abs_normal"

mis <- mis_sd <- array(0, dim = c(2, 8))
dup <- 30

for (g in 1:length(gamma_range)) { # gamma

  s_max <- s_min * alpha_range[g]

  comp_err <- c()
  stat_err <- c()

  for (n in 1:dup) { # dup

    seed <- 1000 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = F, asymm = F, p = rep(p, 3), r = rep(r, 3),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 1, theta_dist = theta_dist
    )

    # comp
    comp_res <- dtbm(dat$Y, r = rep(r, 3), max_iter = 50, asymm = F)
    comp_err <- c(comp_err, CER(comp_res$z[[1]], dat$z[[1]]))

    # stats
    stat_res <- angle_iteration(dat$Y, dat$z, max_iter = 50, asymm = F)
    stat_err <- c(stat_err, CER(stat_res$z[[1]], dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(comp_err)
  mis[2, g] <- mean(stat_err)

  mis_sd[1, g] <- sd(comp_err)
  mis_sd[2, g] <- sd(stat_err)
} # end gamma

save(mis, mis_sd, file = "presaved/Figure4_tensor100.RData")



# tensor case with p = 80
p <- 80
r <- 5

gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4)
alpha_range <- 1 + 160 * p^(gamma_range / 2)

s_min <- 0.5
delta <- 0.5

dist <- "normal"
theta_dist <- "abs_normal"

mis <- mis_sd <- array(0, dim = c(2, 8))
dup <- 30

for (g in 1:length(gamma_range)) { # gamma

  s_max <- s_min * alpha_range[g]

  comp_err <- c()
  stat_err <- c()

  for (n in 1:dup) { # dup

    seed <- 100 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = F, asymm = F, p = rep(p, 3), r = rep(r, 3),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 1, theta_dist = theta_dist
    )

    # comp
    comp_res <- dtbm(dat$Y, r = rep(r, 3), max_iter = 50, asymm = F)
    comp_err <- c(comp_err, CER(comp_res$z[[1]], dat$z[[1]]))

    # stats
    stat_res <- angle_iteration(dat$Y, dat$z, max_iter = 50, asymm = F)
    stat_err <- c(stat_err, CER(stat_res$z[[1]], dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(comp_err)
  mis[2, g] <- mean(stat_err)

  mis_sd[1, g] <- sd(comp_err)
  mis_sd[2, g] <- sd(stat_err)
} # end gamma

save(mis, mis_sd, file = "presaved/Figure4_tensor80.RData")


# matrix case with p = 100
p <- 100
r <- 5

gamma_range <- c(-1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4)
alpha_range <- 1 + 30 * p^(gamma_range / 2)

s_min <- 0.5
delta <- 0.8

dist <- "normal"
theta_dist <- "abs_normal"

mis <- mis_sd <- array(0, dim = c(2, 9))
dup <- 30

for (g in 1:length(gamma_range)) { # gamma

  s_max <- s_min * alpha_range[g]

  comp_err <- c()
  stat_err <- c()

  for (n in 1:dup) { # dup

    seed <- 1000 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = T, asymm = F, p = rep(p, 2), r = rep(r, 2),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 1, theta_dist = theta_dist
    )

    # comp
    comp_res <- dtbm(dat$Y, r = rep(r, 2), max_iter = 50, asymm = F)
    comp_err <- c(comp_err, CER(comp_res$z[[1]], dat$z[[1]]))

    # stats
    stat_res <- angle_iteration(dat$Y, dat$z, max_iter = 50, asymm = F)
    stat_err <- c(stat_err, CER(stat_res$z[[1]], dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(comp_err)
  mis[2, g] <- mean(stat_err)

  mis_sd[1, g] <- sd(comp_err)
  mis_sd[2, g] <- sd(stat_err)
} # end gamma

save(mis, mis_sd, file = "presaved/Figure4_matrix100.RData")


# matrix case with p = 80
p <- 80
r <- 5

gamma_range <- c(-1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4)
alpha_range <- 1 + 25 * p^(gamma_range / 2)

s_min <- 0.5
delta <- 0.8

dist <- "normal"
theta_dist <- "abs_normal"

mis <- mis_sd <- array(0, dim = c(2, 9))
dup <- 30

for (g in 1:length(gamma_range)) { # gamma

  g <- 5
  n <- 5
  s_max <- s_min * alpha_range[g]

  comp_err <- c()
  stat_err <- c()

  for (n in 1:dup) { # dup

    seed <- 100 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = T, asymm = F, p = rep(p, 2), r = rep(r, 2),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 1, theta_dist = theta_dist
    )

    # comp
    comp_res <- dtbm(dat$Y, r = rep(r, 2), max_iter = 50, asymm = F)
    comp_err <- c(comp_err, CER(comp_res$z[[1]], dat$z[[1]]))

    # stats
    stat_res <- angle_iteration(dat$Y, dat$z, max_iter = 50, asymm = F)
    stat_err <- c(stat_err, CER(stat_res$z[[1]], dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(comp_err)
  mis[2, g] <- mean(stat_err)

  mis_sd[1, g] <- sd(comp_err)
  mis_sd[2, g] <- sd(stat_err)
} # end gamma

save(mis, mis_sd, file = "presaved/Figure4_matrix80.RData")
