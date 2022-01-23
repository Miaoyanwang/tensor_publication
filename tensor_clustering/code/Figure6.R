# Code for Figure 6. Comparison among HLloyd, HOSVD, HOSVD+, SCORE with varying signal exponent

# load dependencies ---------------------------------------------------------------------------------------------
install.packages("software/dTBM.tar.gz", repos = NULL, type = "source")
library(dTBM)

source("software/HOLloyd.R")
source("software/tensor_score.R")
source("software/hosvd_related.R")

library(RSKC)

# for plotting
library(ggplot2)
library(patchwork)

# start plotting with saved data from prior simulation ---------------------------------------------------------------
# if a new run of simulation is desired, please go to line #87 to run the code --------------------------------------


# normal case

load("presaved/Figure6_normal.RData")

gamma_range <- c(-1.5, -1.4, -1.3, -1.2, -1.1)
sel_color <- c("#D4613E", "#F3A002", "#B3A86A", "#99B6BD", "#F3C937")

data <- data.frame(
  MCR = c(mis[1, ], mis[2, ], mis[3, ], mis[4, ], mis[5, ]),
  sd = c(mis_sd[1, ], mis_sd[2, ], mis_sd[3, ], mis_sd[4, ], mis_sd[5, ]) / sqrt(30),
  Algorithm = rep(c("Ours", "HLloyd", "HOSVD", "HOSVD+", "SCORE"), each = 5),
  gamma = rep(gamma_range, times = 5)
)
data[, 3] <- factor(data[, 3], levels = c("Ours", "HLloyd", "HOSVD", "HOSVD+", "SCORE"))
p3 <- ggplot(data = data, aes(x = gamma, y = MCR)) +
  geom_line(aes(color = Algorithm), size = 1.2) +
  scale_colour_manual(values = sel_color) +
  geom_point(size = 3, aes(shape = Algorithm)) +
  scale_shape_manual(values = c(16, 17, 5, 2, 4)) +
  xlab("Gamma") +
  ylab("Clustering Error Rate (CER)") +
  labs(title = "Gaussian") +
  coord_cartesian(ylim = c(0, 0.25)) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 1))
p3 <- p3 + geom_errorbar(aes(ymin = MCR - sd, ymax = MCR + sd), width = 0.02, position = position_dodge(0.05)) + theme_light() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 14)) +
  theme(legend.position = "none")


# binary case

load("presaved/Figure6_binary.RData")

gamma_range <- c(-1.5, -1.4, -1.3, -1.2, -1.1)
sel_color <- c("#D4613E", "#F3A002", "#B3A86A", "#99B6BD", "#F3C937")


data <- data.frame(
  MCR = c(mis[1, ], mis[2, ], mis[3, ], mis[4, ], mis[5, ]),
  sd = c(mis_sd[1, ], mis_sd[2, ], mis_sd[3, ], mis_sd[4, ], mis_sd[5, ]) / sqrt(30),
  Algorithm = rep(c("Ours", "HLloyd", "HOSVD", "HOSVD++", "SCORE"), each = 5),
  gamma = rep(gamma_range, times = 5)
)
data[, 3] <- factor(data[, 3], levels = c("Ours", "HLloyd", "HOSVD", "HOSVD++", "SCORE"))
p4 <- ggplot(data = data, aes(x = gamma, y = MCR)) +
  geom_line(aes(color = Algorithm), size = 1.2) +
  scale_colour_manual(values = sel_color) +
  geom_point(size = 3, aes(shape = Algorithm)) +
  scale_shape_manual(values = c(16, 17, 5, 2, 4)) +
  xlab("Gamma") +
  ylab("Clustering Error Rate (CER)") +
  labs(title = "Binary") +
  coord_cartesian(ylim = c(0, 0.6)) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 1))
p4 <- p4 + geom_errorbar(aes(ymin = MCR - sd, ymax = MCR + sd), width = 0.02, position = position_dodge(0.05)) + theme_light() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 14))

pdf("figure/Figure6.pdf", height = 4, width = 10)
(p3 | p4)
dev.off()


# end of plotting -----------------------------------------------------------------------------------------------



# If a new run of simulation is desired, please run the code from here and save the results as .RData -----------
# Then run the above code to generate figures -------------------------------------------------------------------

# normal case

p <- 100
r <- 5

gamma_range <- c(-1.5, -1.4, -1.3, -1.2, -1.1)
alpha_range <- 1 + 100 * p^(gamma_range / 2)

s_min <- 0.5
delta <- 0.5

dist <- "normal"
theta_dist <- "abs_normal"

mis <- mis_sd <- array(0, dim = c(5, 5))
dup <- 30

for (g in 1:length(gamma_range)) { # gamma

  s_max <- s_min * alpha_range[g]

  our_err <- c()
  han_err <- c() # HLloyd
  hosvd_err1 <- c() # HOSVD
  hosvd_err2 <- c() # HOSVD+
  ke_err <- c() # SCORE

  for (n in 1:dup) { # dup

    seed <- 100 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = F, asymm = F, p = rep(p, 3), r = rep(r, 3),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 0.75, theta_dist = theta_dist
    )

    # ours
    our_res <- dtbm(dat$Y, rep(r, 3), max_iter = 20, asymm = F)
    our_err <- c(our_err, CER(our_res$z[[1]], dat$z[[1]]))

    # HLloyd
    z.HOSC <- HO.SC(as.tensor(dat$Y), r) # high-order initialization
    z.Lloyd.HOSC <- HO.Lloyd(as.tensor(dat$Y), z.HOSC) # HLloyd iteration
    han_err <- c(han_err, min(CER(z.Lloyd.HOSC[[1]], dat$z[[1]]), CER(z.Lloyd.HOSC[[2]], dat$z[[1]]), CER(z.Lloyd.HOSC[[3]], dat$z[[1]])))

    # HOSVD
    hosvd1 <- HOSVD_kmeans(dat$Y, r)
    hosvd_err1 <- c(hosvd_err1, CER(hosvd1, dat$z[[1]]))

    # HOSVD
    hosvd2 <- HOSVD_score(dat$Y, r)
    hosvd_err2 <- c(hosvd_err2, CER(hosvd2, dat$z[[1]]))

    # SCORE
    ke_result <- tensor_score_adj(dat$Y, r, rm_diag = T, hooi = T, reg_hooi = T, score_reg = T, max_iter = 20)
    ke_err <- c(ke_err, CER(ke_result, dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(our_err)
  mis[2, g] <- mean(han_err)
  mis[3, g] <- mean(hosvd_err1)
  mis[4, g] <- mean(hosvd_err2)
  mis[5, g] <- mean(ke_err)

  mis_sd[1, g] <- sd(our_err)
  mis_sd[2, g] <- sd(han_err)
  mis_sd[3, g] <- sd(hosvd_err1)
  mis_sd[4, g] <- sd(hosvd_err2)
  mis_sd[5, g] <- sd(ke_err)
} # end gamma

save(mis, mis_sd, file = "presaved/Figure6_normal.RData")



# binary case

p <- 100
r <- 5

gamma_range <- c(-1.5, -1.4, -1.3, -1.2, -1.1)
alpha_range <- 1 + 50 * p^(gamma_range / 2)

s_min <- 0.5
delta <- 0.5

dist <- "binary"
theta_dist <- "abs_normal"

mis <- mis_sd <- array(0, dim = c(5, 5))
dup <- 30

for (g in 1:length(gamma_range)) { # gamma

  s_max <- s_min * alpha_range[g]

  our_err <- c()
  han_err <- c() # HLloyd
  hosvd_err1 <- c() # HOSVD
  hosvd_err2 <- c() # HOSVD+
  ke_err <- c() # SCORE

  for (n in 1:dup) { # dup

    seed <- 100 * g + n
    set.seed(seed)

    dat <- sim_dTBM(
      seed = seed, imat = F, asymm = F, p = rep(p, 3), r = rep(r, 3),
      core_control = "control", delta = delta, s_min = s_min, s_max = s_max,
      dist = dist, sigma = 0.75, theta_dist = theta_dist
    )

    # ours
    our_res <- dtbm(dat$Y, rep(r, 3), max_iter = 20, asymm = F)
    our_err <- c(our_err, CER(our_res$z[[1]], dat$z[[1]]))

    # HLloyd
    z.HOSC <- HO.SC(as.tensor(dat$Y), r) # high-order initialization
    z.Lloyd.HOSC <- HO.Lloyd(as.tensor(dat$Y), z.HOSC) # HLloyd iteration
    han_err <- c(han_err, min(CER(z.Lloyd.HOSC[[1]], dat$z[[1]]), CER(z.Lloyd.HOSC[[2]], dat$z[[1]]), CER(z.Lloyd.HOSC[[3]], dat$z[[1]])))

    # HOSVD
    hosvd1 <- HOSVD_kmeans(dat$Y, r)
    hosvd_err1 <- c(hosvd_err1, CER(hosvd1, dat$z[[1]]))

    # HOSVD
    hosvd2 <- HOSVD_score(dat$Y, r)
    hosvd_err2 <- c(hosvd_err2, CER(hosvd2, dat$z[[1]]))

    # SCORE
    ke_result <- tensor_score_adj(dat$Y, r, rm_diag = T, hooi = T, reg_hooi = T, score_reg = T, max_iter = 20)
    ke_err <- c(ke_err, CER(ke_result, dat$z[[1]]))
  } # end dup

  mis[1, g] <- mean(our_err)
  mis[2, g] <- mean(han_err)
  mis[3, g] <- mean(hosvd_err1)
  mis[4, g] <- mean(hosvd_err2)
  mis[5, g] <- mean(ke_err)

  mis_sd[1, g] <- sd(our_err)
  mis_sd[2, g] <- sd(han_err)
  mis_sd[3, g] <- sd(hosvd_err1)
  mis_sd[4, g] <- sd(hosvd_err2)
  mis_sd[5, g] <- sd(ke_err)
} # end gamma

save(mis, mis_sd, file = "presaved/Figure6_binary.RData")
