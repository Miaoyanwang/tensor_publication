# Code for Figure 6: ME comparison of STD and other methods under model misspecification
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
# if a new run of simulation is desired, please go to line #184 to run the code --------------------------------------

# model misspecification with non-i.i.d. noise ----------------------------

# load pre-saved data
load("presaved/Figure6_noniid.RData")

# ME vs Correlation_level
new_color <- c("#069AA0", "#BCA455", "#7B533E")

final <- 1 - final_c_cor
finalsd <- final_sd_c_cor

s <- 1 # rank
r <- 1 # signal

data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ], final[s, r, 4, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ], finalsd[s, r, 4, ]) / sqrt(30), Method = rep(c("STD (Our method)", "Envelope", "SupCP"), 4), Category = c(rep(0, 3), rep(0.2, 3), rep(0.3, 3), rep(0.4, 3)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "Envelope", "SupCP"))
p1 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Correlation level", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(title = "Low Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  theme(legend.position = "none")
# p1

s <- 1
r <- 2
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ], final[s, r, 4, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ], finalsd[s, r, 4, ]) / sqrt(30), Method = rep(c("STD (Our method)", "Envelope", "SupCP"), 4), Category = c(rep(0, 3), rep(0.2, 3), rep(0.3, 3), rep(0.4, 3)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "Envelope", "SupCP"))
p2 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Correlation level", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 0.15)) +
  labs(title = "High Signal, Low Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  theme(legend.position = "none")
# p2

s <- 2
r <- 1
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ], final[s, r, 4, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ], finalsd[s, r, 4, ]) / sqrt(30), Method = rep(c("STD (Our method)", "Envelope", "SupCP"), 4), Category = c(rep(0, 3), rep(0.2, 3), rep(0.3, 3), rep(0.4, 3)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "Envelope", "SupCP"))
p3 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Correlation level", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(title = "High Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10))
# p3

s <- 2
r <- 2
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ], final[s, r, 4, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ], finalsd[s, r, 4, ]) / sqrt(30), Method = rep(c("STD (Our method)", "Envelope", "SupCP"), 4), Category = c(rep(0, 3), rep(0.2, 3), rep(0.3, 3), rep(0.4, 3)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "Envelope", "SupCP"))
p4 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Correlation level", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 0.15)) +
  labs(title = "High Signal, High Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10))
# p4


pdf("Figures/Figure6_noniid.pdf", width = 8, height = 3)
(p1 | p3)
dev.off()

# model misspecification with sprasity ----------------------------

# load pre-saved data
load("presaved/Figure6_sparse.RData")

# ME vs Sparsity_level
new_color <- c("#069AA0", "#D6CFC4")

final <- 1 - final_s_cor
finalsd <- final_sd_s_cor

s <- 1 # signal
r <- 1 # rank
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ]) / sqrt(30), Method = rep(c("STD (Our method)", "GLSNet"), 3), Category = c(rep(0.1, 2), rep(0.5, 2), rep(0.9, 2)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "GLSNet"))
p1 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Sparsity level", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 0.1)) +
  labs(title = "Low Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  theme(legend.position = "none")
# p1

s <- 1
r <- 2
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ]) / sqrt(30), Method = rep(c("STD (Our method)", "GLSNet"), 3), Category = c(rep(0.1, 2), rep(0.5, 2), rep(0.9, 2)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "GLSNet"))
p2 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Sparsity level", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 0.1)) +
  labs(title = "High Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10))
# p2

s <- 2
r <- 1
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ]) / sqrt(30), Method = rep(c("STD (Our method)", "GLSNet"), 3), Category = c(rep(0.1, 2), rep(0.5, 2), rep(0.9, 2)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "GLSNet"))
p3 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Sparsity level", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 0.1)) +
  labs(title = "High Signal, Low Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10)) +
  theme(legend.position = "none")
# p3

s <- 2
r <- 2
data <- data.frame(PMSE = c(final[s, r, 1, ], final[s, r, 2, ], final[s, r, 3, ]), sd = c(finalsd[s, r, 1, ], finalsd[s, r, 2, ], finalsd[s, r, 3, ]) / sqrt(30), Method = rep(c("STD (Our method)", "GLSNet"), 3), Category = c(rep(0.1, 2), rep(0.5, 2), rep(0.9, 2)))
data[, 3] <- factor(data[, 3], levels = c("STD (Our method)", "GLSNet"))
p4 <- ggplot(data = data, aes(x = as.factor(Category), y = PMSE, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = PMSE - sd, ymax = PMSE + sd), width = .2, position = position_dodge(.9)) +
  labs(x = "Sparsity level", y = "Misalignment error") +
  coord_cartesian(ylim = c(0, 0.1)) +
  labs(title = "High Signal, High Rank", size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  scale_fill_manual(values = new_color) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 10))
# p4

pdf("Figures/Figure6_sparse.pdf", width = 8, height = 3)
(p1 | p2)
dev.off()

# end of plotting -----------------------------------------------------------------------------------------------



# If a new run of simulation is desired, please run the code from here and save the results as .RData -----------
# Then run the above code to generate figures -------------------------------------------------------------------

# model misspecification with non-i.i.d. noise ----------------------------
D <- c(20, 20)
n <- 20
p <- 5
dup <- 30
u_range <- rbind(c(3, 3), c(4, 5))
signal_level_range <- c(3, 6)
cor_level_range <- c(0, 0.2, 0.3, 0.4)

final_c <- final_sd_c <- array(0, dim = c(2, 2, 4, 3))
err_res <- err_env <- array(0, dim = c(2, 2, 4, dup))

final_c_cor <- final_sd_c_cor <- array(0, dim = c(2, 2, 4, 3))
cor_res <- cor_env <- array(0, dim = c(2, 2, 4, dup))

tsr_list_i <- list()
x_list_i <- list()
u_list_i <- list()
para_list_i <- list()
for (i in 1:2) { # u
  tsr_list_j <- list()
  x_list_j <- list()
  u_list_j <- list()
  para_list_j <- list()
  for (j in 1:2) { # signal
    tsr_list_k <- list()
    x_list_k <- list()
    u_list_k <- list()
    para_list_k <- list()
    for (k in 1:4) { # cor_level
      tsr_list <- list()
      x_list <- list()
      u_list <- list()
      para_list <- list()
      for (m in 1:dup) { # dup

        seed <- 10000 * i + 1000 * j + 100 * k + m
        set.seed(seed)

        u <- u_range[i, ]
        signal_level <- signal_level_range[j]
        cor_level <- cor_level_range[k]

        # generate data
        dat <- env_sim(D, u, p, n, signal_level, cor_level, dup)

        # fit STD
        std_fit <- tensor_regress(dat$tsr[[m]]@data, X_covar3 = dat$X, core_shape = c(u, p), niter = 10, cons = "non", dist = "normal", initial = "QR_tucker")
        err_res[i, j, k, m] <- mean((std_fit$U - dat$U@data)^2)
        cor_res[i, j, k, m] <- cor(as.vector(std_fit$U), as.vector(dat$U@data))

        # fit Envelope
        env_fit <- TRR.fit(t(dat$X), dat$tsr[[m]]@data, u, method = "FG")
        err_env[i, j, k, m] <- mean((env_fit$fitted.values@data - dat$U@data)^2)
        cor_env[i, j, k, m] <- cor(as.vector(env_fit$fitted.values@data), as.vector(dat$U@data))

        # transpose the data and save it in the list
        tsr_list[[m]] <- tensor_trans(dat$tsr[[m]]@data)
        x_list[[m]] <- dat$X
        u_list[[m]] <- tensor_trans(dat$U@data)
        para_list[[m]] <- c("u", u, "signal", signal_level, "cor_level", cor_level, "dup", m)
      }

      final_c[i, j, k, 1] <- mean(err_res[i, j, k, ])
      final_c[i, j, k, 2] <- mean(err_env[i, j, k, ])

      final_sd_c[i, j, k, 1] <- sd(err_res[i, j, k, ])
      final_sd_c[i, j, k, 2] <- sd(err_env[i, j, k, ])

      final_c_cor[i, j, k, 1] <- mean(cor_res[i, j, k, ])
      final_c_cor[i, j, k, 2] <- mean(cor_env[i, j, k, ])

      final_sd_c_cor[i, j, k, 1] <- sd(cor_res[i, j, k, ])
      final_sd_c_cor[i, j, k, 2] <- sd(cor_env[i, j, k, ])


      tsr_list_k[[k]] <- tsr_list
      x_list_k[[k]] <- x_list
      u_list_k[[k]] <- u_list
      para_list_k[[k]] <- para_list
    }

    tsr_list_j[[j]] <- tsr_list_k
    x_list_j[[j]] <- x_list_k
    u_list_j[[j]] <- u_list_k
    para_list_j[[j]] <- para_list_k
  }

  tsr_list_i[[i]] <- tsr_list_j
  x_list_i[[i]] <- x_list_j
  u_list_i[[i]] <- u_list_j
  para_list_i[[i]] <- para_list_j
}

# save data in .mat file
sample_list <- list(tsr = tsr_list_i, y = x_list_i, x_true = u_list_i, para = para_list_i)
write.mat(sample_list, "Figure6_noniid_data.mat")

# read output from Figure6.m
supcp_result <- read.mat("mat_output/Figure6_noniid.mat")

final_c[, , , 3] <- supcp_result$final_c_sup
final_sd_c[, , , 3] <- supcp_result$final_sd_c_sup

final_c_cor[, , , 3] <- supcp_result$final_c_cor_sup
final_sd_c_cor[, , , 3] <- supcp_result$final_sd_c_cor_sup

save(final_c, final_sd_c, final_c_cor, final_sd_c_cor, file = "presaved/Figure6_noniid.RData")


# model misspecification with sprasity ----------------------------
d <- 20
d3 <- 20
p <- 8
dup <- 30
signal_range <- c(2, 4)
r_range <- c(2, 4)
s_range <- c(0.1, 0.5, 0.9)

final_s <- final_sd_s <- array(0, dim = c(2, 2, 3, 2))
err_res <- err_mreg <- array(0, dim = c(2, 2, 3, dup))

final_s_cor <- final_sd_s_cor <- array(0, dim = c(2, 2, 3, 2))
cor_res <- cor_mreg <- array(0, dim = c(2, 2, 3, dup))

for (i in 1:2) { # signal
  for (j in 1:2) { # rank
    for (k in 1:3) { # sparsity

      signal <- signal_range[i]
      R <- r3 <- r_range[j]
      s <- s_range[k]

      for (m in 1:dup) {
        seed <- 10000 * i + 1000 * j + 100 * k + m
        set.seed(seed)

        dat <- sparse_sim(d, d3, p, R, r3, s, signal, dup)

        # fit STD
        std_fit <- tensor_regress(dat$tsr[[m]], X_covar3 = dat$X, core_shape = c(R, R, r3), niter = 20, cons = "non", dist = "binary", initial = "QR_tucker")
        err_res[i, j, k, m] <- mean((std_fit$U - dat$U)^2)
        cor_res[i, j, k, m] <- cor(as.vector(std_fit$U), as.vector(dat$U))

        # fit GLSNet
        mreg_res <- mySymGLM(as.tensor(dat$tsr[[m]]), dat$X, R, sparsity = d * d * p * s, niter = 50)

        Theta <- mreg_res$A %*% diag(mreg_res$w) %*% t(mreg_res$A)
        mreg_B <- ttm(mreg_res$BB, dat$X, m = 3)@data
        mreg_U <- mreg_B
        for (l in 1:dim(mreg_B)[3]) {
          mreg_U[, , l] <- mreg_B[, , l] + Theta
        }
        err_mreg[i, j, k, m] <- mean((mreg_U - dat$U)^2)
        cor_mreg[i, j, k, m] <- cor(as.vector(mreg_U), as.vector(dat$U))
      }

      final_s[i, j, k, 1] <- mean(err_res[i, j, k, ])
      final_s[i, j, k, 2] <- mean(err_mreg[i, j, k, ])

      final_sd_s[i, j, k, 1] <- sd(err_res[i, j, k, ])
      final_sd_s[i, j, k, 2] <- sd(err_mreg[i, j, k, ])

      final_s_cor[i, j, k, 1] <- mean(cor_res[i, j, k, ])
      final_s_cor[i, j, k, 2] <- mean(cor_mreg[i, j, k, ])

      final_sd_s_cor[i, j, k, 1] <- sd(cor_res[i, j, k, ])
      final_sd_s_cor[i, j, k, 2] <- sd(cor_mreg[i, j, k, ])
    }
  }
}

save(final_s, final_sd_s, final_s_cor, final_sd_s_cor, file = "presaved/Figure6_sparse.RData")
