# Code for Figure 5. Initialization vs Combined algorithm

# load dependencies ---------------------------------------------------------------------------------------------
install.packages("software/dTBM.tar.gz", repos = NULL, type = "source")
library(dTBM)
library(RSKC)

# for plotting
library(ggplot2)
library(patchwork)

# start plotting with saved data from prior simulation ---------------------------------------------------------------
# if a new run of simulation is desired, please go to line #102 to run the code --------------------------------------

# normal case


load("presaved/Figure5_normal100.RData")
mis1 <- mis
mis_sd1 <- mis_sd

load("presaved/Figure5_normal80.RData")
mis2 <- mis
mis_sd2 <- mis_sd

gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4)

sel_color <- c("#F3C937", "#2C7CE6")

data <- data.frame(
  MCR = c(mis1[1, ], mis1[2, ], mis2[1, ], mis2[2, ]), sd = c(mis_sd1[1, ], mis_sd1[2, ], mis_sd2[1, ], mis_sd2[2, ]) / sqrt(30),
  Algorithm = rep(c("Wkmeans", "Wkmeans + hALloyd"), each = 8, times = 2), gamma = rep(gamma_range, times = 4),
  p = rep(c(100, 80), each = 16)
)
data[, 3] <- factor(data[, 3], levels = c("Wkmeans", "Wkmeans + hALloyd"))
data[, 5] <- factor(data[, 5], levels = c(100, 80))
p21 <- ggplot(data = data, aes(x = gamma, y = MCR)) +
  geom_line(aes(color = Algorithm, linetype = p), size = 1.2) +
  scale_colour_manual(values = sel_color) +
  geom_point(size = 3, aes(shape = Algorithm)) +
  scale_shape_manual(values = c(16, 5)) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  # theme(axis.text=element_text(size=12),axis.title=element_text(size=10),plot.title = element_text(hjust = 0.5,size = 11)) +
  xlab("Gamma") +
  ylab("Clustering Error Rate (CER)") +
  labs(title = "Gaussian") +
  coord_cartesian(ylim = c(0, 0.35)) +
  guides(colour = guide_legend(order = 2), shape = guide_legend(order = 1))
p21 <- p21 + geom_errorbar(aes(ymin = MCR - sd, ymax = MCR + sd), width = 0.02, position = position_dodge(0.05)) + theme_light() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
  theme(legend.position = "none")


# binary case



load("presaved/Figure5_binary100.RData")
mis1 <- mis
mis_sd1 <- mis_sd

load("presaved/Figure5_binary80.RData")
mis2 <- mis
mis_sd2 <- mis_sd

gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4)

sel_color <- c("#F3C937", "#2C7CE6")

data <- data.frame(
  MCR = c(mis1[1, ], mis1[2, ], mis2[1, ], mis2[2, ]), sd = c(mis_sd1[1, ], mis_sd1[2, ], mis_sd2[1, ], mis_sd2[2, ]) / sqrt(30),
  Algorithm = rep(c("Wkmeans", "Wkmeans + hALloyd"), each = 8, times = 2), gamma = rep(gamma_range, times = 4),
  p = rep(c(100, 80), each = 16)
)
data[, 3] <- factor(data[, 3], levels = c("Wkmeans", "Wkmeans + hALloyd"))
data[, 5] <- factor(data[, 5], levels = c(100, 80))
p22 <- ggplot(data = data, aes(x = gamma, y = MCR)) +
  geom_line(aes(color = Algorithm, linetype = p), size = 1.2) +
  scale_colour_manual(values = sel_color) +
  geom_point(size = 3, aes(shape = Algorithm)) +
  scale_shape_manual(values = c(16, 5)) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  # theme(axis.text=element_text(size=12),axis.title=element_text(size=10),plot.title = element_text(hjust = 0.5,size = 11)) +
  xlab("Gamma") +
  ylab("Clustering Error Rate (CER)") +
  labs(title = "Binary") +
  coord_cartesian(ylim = c(0, 0.35)) +
  guides(linetype = guide_legend(order = 2), shape = guide_legend(order = 1), color = guide_legend(order = 1))
p22 <- p22 + geom_errorbar(aes(ymin = MCR - sd, ymax = MCR + sd), width = 0.02, position = position_dodge(0.05)) + theme_light() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 14))

pdf("figure/Figure5.pdf", width = 10, height = 4)
(p21 | p22)
dev.off()



# end of plotting -----------------------------------------------------------------------------------------------



# If a new run of simulation is desired, please run the code from here and save the results as .RData -----------
# Then run the above code to generate figures -------------------------------------------------------------------

# normal case with p = 100

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

  ini_err <- c()
  re_err <- c()

  for (n in 1:dup) { # dup

    seed <- 100 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = F, asymm = F, p = rep(p, 3), r = rep(r, 3),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 1, theta_dist = theta_dist
    )

    # initial
    ini_res <- wkmeans(dat$Y, r = rep(r, 3), asymm = F)
    ini_err <- c(ini_err, CER(ini_res$z0[[1]], dat$z[[1]]))

    # refinement
    re_res <- angle_iteration(dat$Y, ini_res$z0, max_iter = 50, asymm = F)
    re_err <- c(re_err, CER(re_res$z[[1]], dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(ini_err)
  mis[2, g] <- mean(re_err)

  mis_sd[1, g] <- sd(ini_err)
  mis_sd[2, g] <- sd(re_err)
} # end gamma


save(mis, mis_sd, file = "presaved/Figure5_normal100.RData")

# normal case with p = 80

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

  ini_err <- c()
  re_err <- c()

  for (n in 1:dup) { # dup

    seed <- 100 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = F, asymm = F, p = rep(p, 3), r = rep(r, 3),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 1, theta_dist = theta_dist
    )

    # initial
    ini_res <- wkmeans(dat$Y, r = rep(r, 3), asymm = F)
    ini_err <- c(ini_err, CER(ini_res$z0[[1]], dat$z[[1]]))

    # refinement
    re_res <- angle_iteration(dat$Y, ini_res$z0, max_iter = 50, asymm = F)
    re_err <- c(re_err, CER(re_res$z[[1]], dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(ini_err)
  mis[2, g] <- mean(re_err)

  mis_sd[1, g] <- sd(ini_err)
  mis_sd[2, g] <- sd(re_err)
} # end gamma

save(mis, mis_sd, file = "presaved/Figure5_normal80.RData")


# binary case with p = 100

p <- 100
r <- 5
gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4)
alpha_range <- 1 + 80 * p^(gamma_range / 2)

s_min <- 0.5
delta <- 0.25

dist <- "binary"
theta_dist <- "abs_normal"

mis <- mis_sd <- array(0, dim = c(2, 8))
dup <- 30

for (g in 1:length(gamma_range)) { # gamma

  s_max <- s_min * alpha_range[g]

  ini_err <- c()
  re_err <- c()

  for (n in 1:dup) { # dup

    seed <- 100 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = F, asymm = F, p = rep(p, 3), r = rep(r, 3),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 1, theta_dist = theta_dist
    )

    # initial
    ini_res <- wkmeans(dat$Y, r = rep(r, 3), asymm = F)
    ini_err <- c(ini_err, CER(ini_res$z0[[1]], dat$z[[1]]))

    # refinement
    re_res <- angle_iteration(dat$Y, ini_res$z0, max_iter = 50, asymm = F)
    re_err <- c(re_err, CER(re_res$z[[1]], dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(ini_err)
  mis[2, g] <- mean(re_err)

  mis_sd[1, g] <- sd(ini_err)
  mis_sd[2, g] <- sd(re_err)
} # end gamma


save(mis, mis_sd, file = "presaved/Figure5_binary100.RData")


# binary case with p = 80

p <- 80
r <- 5
gamma_range <- c(-2.1, -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4)
alpha_range <- 1 + 60 * p^(gamma_range / 2)

s_min <- 0.5
delta <- 0.25

dist <- "binary"
theta_dist <- "abs_normal"

mis <- mis_sd <- array(0, dim = c(2, 8))
dup <- 30

for (g in 1:length(gamma_range)) { # gamma

  s_max <- s_min * alpha_range[g]

  ini_err <- c()
  re_err <- c()

  for (n in 1:dup) { # dup

    seed <- 100 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = F, asymm = F, p = rep(p, 3), r = rep(r, 3),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 1, theta_dist = theta_dist
    )

    # initial
    ini_res <- wkmeans(dat$Y, r = rep(r, 3), asymm = F)
    ini_err <- c(ini_err, CER(ini_res$z0[[1]], dat$z[[1]]))

    # refinement
    re_res <- angle_iteration(dat$Y, ini_res$z0, max_iter = 50, asymm = F)
    re_err <- c(re_err, CER(re_res$z[[1]], dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(ini_err)
  mis[2, g] <- mean(re_err)

  mis_sd[1, g] <- sd(ini_err)
  mis_sd[2, g] <- sd(re_err)
} # end gamma


save(mis, mis_sd, file = "presaved/Figure5_binary80.RData")
