# Code for Figure 5: ME comparison among STD, Envelope, GLSNet, mRRR, and SupCP with binary data 

# load dependencies ---------------------------------------------------------------------------------------------
# for plotting
library(ggplot2)
library(patchwork)

# for simulation
install.packages("software/tensorregress.tar.gz", repos = NULL, type = "source")
library(tensorregress)

library(TRES) 
library(rrpack)
library(rmatio)

library(pracma)
library(rTensor)
source("software/netglm_code_adjust.R")
source("simulation.R")


# start plotting with saved data from prior simulation ---------------------------------------------------------------
# if a new run of simulation is desired, please go to line #200 to run the code --------------------------------------

new_color <- c("#069AA0", "#7B533E", "#BCA455", "#CCC591", "#D6CFC4")

# ME vs Number of modes with available features ----------------------------

# load pre-saved data
load("presaved/Figure5_mode.RData")

final <- 1 - final_mode_cor[, , , c(1, 5, 2, 4, 3)]
finalsd <- final_sd_mode_cor[, , , c(1, 5, 2, 4, 3)]

s <- 1
r <- 1
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ]) / sqrt(30), Method = rep(c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"), 3), Category = c(rep(1, 5), rep(2, 5), rep(3, 5)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"))
p1 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Number of modes with available features", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Low Effect, Low Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10))
# p1

s <- 1
r <- 2
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ]) / sqrt(30), Method = rep(c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"), 3), Category = c(rep(1, 5), rep(2, 5), rep(3, 5)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"))
p2 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Number of modes with available features", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Low Effect, High Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  theme(legend.position = "none")
# p2

s <- 2
r <- 1
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ]) / sqrt(30), Method = rep(c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"), 3), Category = c(rep(1, 5), rep(2, 5), rep(3, 5)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"))
p3 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Number of modes with available features", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "High Effect, Low Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10))
# p3

s <- 2
r <- 2
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ]) / sqrt(30), Method = rep(c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"), 3), Category = c(rep(1, 5), rep(2, 5), rep(3, 5)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"))
p4 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Number of modes with available features", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "High Effect, High Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10))
# p4


# ME vs Sample size --------------------------------------------------------

# load pre-saved data
load("presaved/Figure5_sample.RData")

d <- c(20, 40, 60, 80, 100)

final <- 1 - final_sample_cor[, , , c(1, 5, 2, 4, 3)]
finalsd <- final_sd_sample_cor[, , , c(1, 5, 2, 4, 3)]

s <- 1
r <- 1
data <- data.frame(
  d = rep(d, 5), Method = c(rep("STD (Our method)", length(d)), rep("SupCP", length(d)), rep("Envelope", length(d)), rep("mRRR", length(d)), rep("GLSNet", length(d))),
  PMSE = c(final[s, r, , 1], final[s, r, , 2], final[s, r, , 3], final[s, r, , 4], final[s, r, , 5]), sd = c(finalsd[s, r, , 1], finalsd[s, r, , 2], finalsd[s, r, , 3], finalsd[s, r, , 4], finalsd[s, r, , 5]) / sqrt(30)
)
data[, 2] <- factor(data[, 2], levels = c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"))
p5 <- ggplot(data, aes(x = d * 400, y = PMSE)) +
  geom_line(aes(color = Method), size = 1.2) +
  scale_colour_manual(values = new_color) +
  geom_point(size = 2, aes(shape = Method)) +
  scale_shape_manual(values = c(16, 5, 17, 4, 1, 2)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  xlab("Sample Size") +
  labs(title = "Low Effect, Low Rank") +
  ylab("Misalignment   error") +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  coord_cartesian(ylim = c(0, 1))
p5 <- p5 + geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = 0.5, position = position_dodge(0.05))
# p5

s <- 1
r <- 2
data <- data.frame(
  d = rep(d, 5), Method = c(rep("STD (Our method)", length(d)), rep("SupCP", length(d)), rep("Envelope", length(d)), rep("mRRR", length(d)), rep("GLSNet", length(d))),
  PMSE = c(final[s, r, , 1], final[s, r, , 2], final[s, r, , 3], final[s, r, , 4], final[s, r, , 5]), sd = c(finalsd[s, r, , 1], finalsd[s, r, , 2], finalsd[s, r, , 3], finalsd[s, r, , 4], finalsd[s, r, , 5]) / sqrt(30)
)
data[, 2] <- factor(data[, 2], levels = c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"))
p6 <- ggplot(data, aes(x = d * 400, y = PMSE)) +
  geom_line(aes(color = Method), size = 1.2) +
  scale_colour_manual(values = new_color) +
  geom_point(size = 2, aes(shape = Method)) +
  scale_shape_manual(values = c(16, 5, 17, 4, 1, 2)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  xlab("Sample Size") +
  labs(title = "Low Effect, High Rank") +
  ylab("Misalignment   error") +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  coord_cartesian(ylim = c(0, 1))
p6 <- p6 + geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = 0.5, position = position_dodge(0.05)) + theme(legend.position = "none")
# p6

s <- 2
r <- 1
data <- data.frame(
  d = rep(d, 5), Method = c(rep("STD (Our method)", length(d)), rep("SupCP", length(d)), rep("Envelope", length(d)), rep("mRRR", length(d)), rep("GLSNet", length(d))),
  PMSE = c(final[s, r, , 1], final[s, r, , 2], final[s, r, , 3], final[s, r, , 4], final[s, r, , 5]), sd = c(finalsd[s, r, , 1], finalsd[s, r, , 2], finalsd[s, r, , 3], finalsd[s, r, , 4], finalsd[s, r, , 5]) / sqrt(30)
)
data[, 2] <- factor(data[, 2], levels = c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"))
p7 <- ggplot(data, aes(x = d * 400, y = PMSE)) +
  geom_line(aes(color = Method), size = 1.2) +
  scale_colour_manual(values = new_color) +
  geom_point(size = 2, aes(shape = Method)) +
  scale_shape_manual(values = c(16, 5, 17, 4, 1, 2)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  xlab("Sample Size") +
  labs(title = "High Effect, Low Rank") +
  ylab("Misalignment   error") +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  coord_cartesian(ylim = c(0, 1))
p7 <- p7 + geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = 0.5, position = position_dodge(0.05))
# p7

s <- 2
r <- 2
data <- data.frame(
  d = rep(d, 5), Method = c(rep("STD (Our method)", length(d)), rep("SupCP", length(d)), rep("Envelope", length(d)), rep("mRRR", length(d)), rep("GLSNet", length(d))),
  PMSE = c(final[s, r, , 1], final[s, r, , 2], final[s, r, , 3], final[s, r, , 4], final[s, r, , 5]), sd = c(finalsd[s, r, , 1], finalsd[s, r, , 2], finalsd[s, r, , 3], finalsd[s, r, , 4], finalsd[s, r, , 5]) / sqrt(30)
)
data[, 2] <- factor(data[, 2], levels = c("STD (Our method)", "SupCP", "Envelope", "mRRR", "GLSNet"))
p8 <- ggplot(data, aes(x = d * 400, y = PMSE)) +
  geom_line(aes(color = Method), size = 1.2) +
  scale_colour_manual(values = new_color) +
  geom_point(size = 2, aes(shape = Method)) +
  scale_shape_manual(values = c(16, 5, 17, 4, 1, 2)) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  xlab("Sample Size") +
  labs(title = "High Effect, High Rank") +
  ylab("Misalignment   error") +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  coord_cartesian(ylim = c(0, 1))
p8 <- p8 + geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = 0.5, position = position_dodge(0.05))
# p8

pdf("Figures/Figure5.pdf", width = 9, height = 6)
(p2 | p3) /
  (p6 | p7)
dev.off()


# end of plotting -----------------------------------------------------------------------------------------------



# If a new run of simulation is desired, please run the code from here and save the results as .RData -----------
# Then run the above code to generate figures -------------------------------------------------------------------

# ME vs Number of modes with available features ----------------------------

signal_range <- c(3, 6)
core_range <- rbind(c(3, 3, 3), c(4, 5, 6))
dup <- 10
d <- 20
whole_shape <- rep(d, 3)

final_mode <- final_sd_mode <- array(0, dim = c(2, 2, 3, 5)) 
err_res <- err_env <- err_rrr <- err_mreg1 <- array(0, dim = c(2, 2, 3, dup))

final_mode_cor <- final_sd_mode_cor <- array(0, dim = c(2, 2, 3, 5))
cor_res <- cor_env <- cor_rrr <- cor_mreg1 <- array(0, dim = c(2, 2, 3, dup))

dist <- "binary"


tsr_list_s <- list()
x_list_s <- list()
u_list_s <- list()
para_list_s <- list()
for (s in 1:2) { # signal
  tsr_list_r <- list()
  x_list_r <- list()
  u_list_r <- list()
  para_list_r <- list()

  for (r in 1:2) { # rank
    tsr_list_i <- list()
    x_list_i <- list()
    u_list_i <- list()
    para_list_i <- list()

    for (i in 1:3) { # # of informative modes
      tsr_list <- list()
      x_list <- list()
      u_list <- list()
      para_list <- list()
      for (n in 1:dup) {
        seed <- 10000 * s + 1000 * r + 100 * i + n
        set.seed(seed)

        core_shape <- core_range[r, ]
        signal <- signal_range[s]
        cat("signal = ", signal, ", rank = ", core_shape, ", info number = ", i, ", dup = ", n, "\n")

        if (i == 1) { # one mode side information
          data <- sim_data(seed, whole_shape = whole_shape, core_shape = core_shape, p = c(0, 0, 0.4 * d), dist = dist, dup = dup, signal = signal, ortho = T)
          res <- tensor_regress(data$tsr[[n]], X_covar3 = data$X_covar3, core_shape = core_shape, cons = "non", dist = dist, initial = "QR_tucker") ### tensor_regress
        } else if (i == 2) { # two mode side information
          data <- sim_data(seed, whole_shape = whole_shape, core_shape = core_shape, p = c(0, 0.4 * d, 0.4 * d), dist = dist, dup = dup, signal = signal, ortho = T)
          res <- tensor_regress(data$tsr[[n]], X_covar2 = data$X_covar2, X_covar3 = data$X_covar3, core_shape = core_shape, cons = "non", dist = dist, initial = "QR_tucker") ### tensor_regress
        } else if (i == 3) { # three model side information
          data <- sim_data(seed, whole_shape = whole_shape, core_shape = core_shape, p = c(0.4 * d, 0.4 * d, 0.4 * d), dist = dist, dup = dup, signal = signal, ortho = T)
          res <- tensor_regress(data$tsr[[n]], X_covar1 = data$X_covar1, X_covar2 = data$X_covar2, X_covar3 = data$X_covar3, core_shape = core_shape, cons = "non", dist = dist, initial = "QR_tucker") ### tensor_regress
        }

        # fit STD
        err_res[s, r, i, n] <- mean((res$U - data$U)^2)
        cor_res[s, r, i, n] <- cor(as.vector(res$U), as.vector(data$U))

        # fit GLSNet
        rank_range <- c(2, 3, 4, 5, 6)
        sparsity_range <- (d^3) * 0.4 * c(0.1, 0.3, 0.5, 0.7, 0.9)

        mreg_res <- search_GLM(as.tensor(data$tsr[[n]]), data$X_covar3, data$U, rank_range, sparsity_range, dist = "binary")

        err_mreg1[s, r, i, n] <- mreg_res$err_opt
        cor_mreg1[s, r, i, n] <- mreg_res$cor_opt

        # fit Envelope
        dat <- data$tsr[[n]]
        dat[which(dat == 0)] <- -1
        
        env_fit <- TRR.fit(t(data$X_covar3), dat, core_shape[1:2], method = "FG")

        Ey <- -1 + 2 * logistic(data$U)
        err_env[s, r, i, n] <- mean((env_fit$fitted.values@data - Ey)^2)
        cor_env[s, r, i, n] <- cor(as.vector(env_fit$fitted.values@data), as.vector(Ey))

        # fit mRRR
        data0 <- unfold(as.tensor(data$tsr[[n]]), row_idx = 3, col_idx = c(1, 2))@data
        U0 <- unfold(as.tensor(data$U), row_idx = 3, col_idx = c(1, 2))@data

        rrr_fit <- mrrr(data0, X = data$X_covar3, family = list(binomial()), penstr = list(penaltySVD = "rankCon", lambdaSVD = core_shape[3]))

        err_rrr[s, r, i, n] <- mean((rrr_fit$mu - logistic(U0))^2)
        cor_rrr[s, r, i, n] <- cor(as.vector(rrr_fit$mu), as.vector(logistic(U0)))

        # transpose the data and save it in the list
        tsr_list[[n]] <- tensor_trans(dat)
        x_list[[n]] <- data$X_covar3
        u_list[[n]] <- tensor_trans(Ey)
        para_list[[n]] <- c("signal", signal_range[s], "rank", core_shape, "info", i, "dup", n)
      }

      final_mode[s, r, i, 1] <- mean(err_res[s, r, i, ])
      final_mode[s, r, i, 2] <- mean(err_env[s, r, i, ])
      final_mode[s, r, i, 3] <- mean(err_mreg1[s, r, i, ])
      final_mode[s, r, i, 4] <- mean(err_rrr[s, r, i, ])

      final_sd_mode[s, r, i, 1] <- sd(err_res[s, r, i, ])
      final_sd_mode[s, r, i, 2] <- sd(err_env[s, r, i, ])
      final_sd_mode[s, r, i, 3] <- sd(err_mreg1[s, r, i, ])
      final_sd_mode[s, r, i, 4] <- sd(err_rrr[s, r, i, ])

      final_mode_cor[s, r, i, 1] <- mean(cor_res[s, r, i, ])
      final_mode_cor[s, r, i, 2] <- mean(cor_env[s, r, i, ])
      final_mode_cor[s, r, i, 3] <- mean(cor_mreg1[s, r, i, ])
      final_mode_cor[s, r, i, 4] <- mean(cor_rrr[s, r, i, ])

      final_sd_mode_cor[s, r, i, 1] <- sd(cor_res[s, r, i, ])
      final_sd_mode_cor[s, r, i, 2] <- sd(cor_env[s, r, i, ])
      final_sd_mode_cor[s, r, i, 3] <- sd(cor_mreg1[s, r, i, ])
      final_sd_mode_cor[s, r, i, 4] <- sd(cor_rrr[s, r, i, ])

      
      tsr_list_i[[i]] <- tsr_list
      x_list_i[[i]] <- x_list
      u_list_i[[i]] <- u_list
      para_list_i[[i]] <- para_list
    } 
    tsr_list_r[[r]] <- tsr_list_i
    x_list_r[[r]] <- x_list_i
    u_list_r[[r]] <- u_list_i
    para_list_r[[r]] <- para_list_i
  } 
  tsr_list_s[[s]] <- tsr_list_r
  x_list_s[[s]] <- x_list_r
  u_list_s[[s]] <- u_list_r
  para_list_s[[s]] <- para_list_r
}

# save data in .mat file
sample_list <- list(tsr = tsr_list_s, y = x_list_s, x_true = u_list_s, para = para_list_s)
write.mat(sample_list, "Figure5_mode_data.mat")

# read output from Figure5.m
supcp_result <- read.mat("mat_output/Figure5_mode.mat")

final_mode[, , , 5] <- supcp_result$final_mode_sup
final_sd_mode[, , , 5] <- supcp_result$final_sd_mode_sup

final_mode_cor[, , , 5] <- supcp_result$final_mode_cor_sup
final_sd_mode_cor[, , , 5] <- supcp_result$final_sd_mode_cor_sup

save(final_mode, final_sd_mode, final_mode_cor, final_sd_mode_cor, file = "presaved/Figure5_mode.RData")



# ME vs Sample size --------------------------------------------------------

d <- c(20, 40, 60, 80, 100)
signal_range <- c(3, 6)
core_range <- rbind(c(3, 3, 3), c(4, 5, 6))
dup <- 30

final_sample <- final_sd_sample <- array(0, dim = c(2, 2, 5, 5)) 
err_res <- err_env <- err_rrr <- err_mreg1 <- array(0, dim = c(2, 2, 5, dup))


final_sample_cor <- final_sd_sample_cor <- array(0, dim = c(2, 2, 5, 5))
cor_res <- cor_env <- cor_rrr <- cor_mreg1 <- array(0, dim = c(2, 2, 5, dup))

dist <- "binary"


tsr_list_s <- list()
x_list_s <- list()
u_list_s <- list()
para_list_s <- list()
for (s in 1:2) { # signal
  tsr_list_r <- list()
  x_list_r <- list()
  u_list_r <- list()
  para_list_r <- list()

  for (r in 1:2) { # rank
    tsr_list_i <- list()
    x_list_i <- list()
    u_list_i <- list()
    para_list_i <- list()

    for (i in 1:5) { # dimension
      tsr_list <- list()
      x_list <- list()
      u_list <- list()
      para_list <- list()

      whole_shape <- c(20, 20, d[i])
      for (n in 1:dup) {
        seed <- 10000 * s + 1000 * r + 100 * i + n
        set.seed(seed)

        core_shape <- core_range[r, ]
        signal <- signal_range[s]
        cat("signal = ", signal, ", rank = ", c(core_shape), ", dim = ", d[i], " dist = ", dist, ", dup = ", n, "\n")

        data <- sim_data(seed, whole_shape = whole_shape, core_shape = core_shape, p = c(0, 0, 0.4 * d[i]), dist = dist, dup = dup, signal = signal, ortho = T)
        res <- tensor_regress(data$tsr[[n]], X_covar3 = data$X_covar3, core_shape = core_shape, niter = 10, cons = "non", dist = dist, initial = "QR_tucker") ### tensor_regress

        # fit STD
        err_res[s, r, i, n] <- mean((res$U - data$U)^2)
        cor_res[s, r, i, n] <- cor(as.vector(res$U), as.vector(data$U))

        # fit GLSNet
        rank_range <- c(2, 3, 4, 5, 6)
        sparsity_range <- (d^3) * 0.4 * c(0.1, 0.3, 0.5, 0.7, 0.9)

        mreg_res <- search_GLM(as.tensor(data$tsr[[n]]), data$X_covar3, data$U, rank_range, sparsity_range, dist = "binary")

        err_mreg1[s, r, i, n] <- mreg_res$err_opt
        cor_mreg1[s, r, i, n] <- mreg_res$cor_opt

        # fit Envelope
        dat <- data$tsr[[n]]
        dat[which(dat == 0)] <- -1
        
        env_fit <- TRR.fit(t(data$X_covar3), dat, core_shape[1:2], method = "FG")

        Ey <- -1 + 2 * logistic(data$U)
        err_env[s, r, i, n] <- mean((env_fit$fitted.values@data - Ey)^2)
        cor_env[s, r, i, n] <- cor(as.vector(env_fit$fitted.values@data), as.vector(Ey))

        # fit mRRR
        data0 <- unfold(as.tensor(data$tsr[[n]]), row_idx = 3, col_idx = c(1, 2))@data
        U0 <- unfold(as.tensor(data$U), row_idx = 3, col_idx = c(1, 2))@data

        rrr_fit <- mrrr(data0, X = data$X_covar3, family = list(binomial()), penstr = list(penaltySVD = "rankCon", lambdaSVD = core_shape[3]))

        err_rrr[s, r, i, n] <- mean((rrr_fit$mu - logistic(U0))^2)
        cor_rrr[s, r, i, n] <- cor(as.vector(rrr_fit$mu), as.vector(logistic(U0)))

        # transpose the data and save it in the list
        tsr_list[[n]] <- tensor_trans(dat)
        x_list[[n]] <- data$X_covar3
        u_list[[n]] <- tensor_trans(Ey)
        para_list[[n]] <- c("signal", signal, "rank", core_shape, "dim", d[i], "dup", n)
      } 

      final_sample[s, r, i, 1] <- mean(err_res[s, r, i, ])
      final_sample[s, r, i, 2] <- mean(err_env[s, r, i, ])
      final_sample[s, r, i, 3] <- mean(err_mreg1[s, r, i, ])
      final_sample[s, r, i, 4] <- mean(err_rrr[s, r, i, ])

      final_sd_sample[s, r, i, 1] <- sd(err_res[s, r, i, ])
      final_sd_sample[s, r, i, 2] <- sd(err_env[s, r, i, ])
      final_sd_sample[s, r, i, 3] <- sd(err_mreg1[s, r, i, ])
      final_sd_sample[s, r, i, 4] <- sd(err_rrr[s, r, i, ])

      final_sample_cor[s, r, i, 1] <- mean(cor_res[s, r, i, ])
      final_sample_cor[s, r, i, 2] <- mean(cor_env[s, r, i, ])
      final_sample_cor[s, r, i, 3] <- mean(cor_mreg1[s, r, i, ])
      final_sample_cor[s, r, i, 4] <- mean(cor_rrr[s, r, i, ])

      final_sd_sample_cor[s, r, i, 1] <- sd(cor_res[s, r, i, ])
      final_sd_sample_cor[s, r, i, 2] <- sd(cor_env[s, r, i, ])
      final_sd_sample_cor[s, r, i, 3] <- sd(cor_mreg1[s, r, i, ])
      final_sd_sample_cor[s, r, i, 4] <- sd(cor_rrr[s, r, i, ])

      tsr_list_i[[i]] <- tsr_list
      x_list_i[[i]] <- x_list
      u_list_i[[i]] <- u_list
      para_list_i[[i]] <- para_list
    } 
    tsr_list_r[[r]] <- tsr_list_i
    x_list_r[[r]] <- x_list_i
    u_list_r[[r]] <- u_list_i
    para_list_r[[r]] <- para_list_i
  } 
  tsr_list_s[[s]] <- tsr_list_r
  x_list_s[[s]] <- x_list_r
  u_list_s[[s]] <- u_list_r
  para_list_s[[s]] <- para_list_r
}

# save data in .mat file
sample_list <- list(tsr = tsr_list_s, y = x_list_s, x_true = u_list_s, para = para_list_s)
write.mat(sample_list, "Figure5_sample_data.mat")

# read output from Figure5.m
supcp_result <- read.mat("mat_output/Figure5_sample.mat")

final_sample[, , , 5] <- supcp_result$final_sample_sup
final_sd_sample[, , , 5] <- supcp_result$final_sd_sample_sup

final_sample_cor[, , , 5] <- supcp_result$final_sample_cor_sup
final_sd_sample_cor[, , , 5] <- supcp_result$final_sd_sample_cor_sup

save(final_sample, final_sd_sample, final_sample_cor, final_sd_sample_cor, file = "presaved/Figure5_sample.RData")
