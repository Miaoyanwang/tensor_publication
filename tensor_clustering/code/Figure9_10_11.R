# Code for Figures 9,10,11. HCP data analysis

# load dependencies & data ---------------------------------------------------------------------------------------------
install.packages("software/dTBM.tar.gz", repos = NULL, type = "source")
library(dTBM)

# for plotting
library(RColorBrewer)
library(ggplot2)
library(patchwork)

# start analysis ------------------------------------------------------------------------------------------------

r <- c(4, 4, 3)
result <- dtbm(HCP$tensor, r = r, max_iter = 50, asymm = T)

# start plotting Figure 10 ------------------------------------------------------------------------------------------------

# obtain estimated S
est_signal <- array(0, dim = r)
for (a in 1:r[1]) {
  for (b in 1:r[2]) {
    for (c in 1:r[3]) {
      est_signal[a, b, c] <- mean(HCP$tensor[result$z[[1]] == a, result$z[[2]] == b, result$z[[3]] == c])
    }
  }
}

ave_s <- (est_signal[, , 1] * sum(result$z[[3]] == 1) +
  est_signal[, , 2] * sum(result$z[[3]] == 2) +
  est_signal[, , 3] * sum(result$z[[3]] == 3)) / 136

# Average

mat_s <- data.frame(
  x = rep(c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"), times = 4),
  y = rep(c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"), each = 4),
  value = as.vector(ave_s)
)
mat_s$x <- factor(mat_s$x, levels = c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"))
mat_s$y <- factor(mat_s$y, levels = c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"))

gg_ave_s <- ggplot(data = mat_s, aes(x, y)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2)), color = "black") +
  scale_fill_distiller(palette = "Purples", direction = 1, limits = c(0, 1)) +
  labs(title = "Weighted Average Brain Connection") +
  xlab("Brain node clusters") +
  ylab("Brain node clusters") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, angle = 45, hjust = 1),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

# Group 1

mat_s1 <- data.frame(
  x = rep(c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"), times = 4),
  y = rep(c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"), each = 4),
  value = as.vector(est_signal[, , 1] - ave_s)
)
mat_s1$x <- factor(mat_s1$x, levels = c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"))
mat_s1$y <- factor(mat_s1$y, levels = c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"))


gg_s1 <- ggplot(data = mat_s1, aes(x, y)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-0.15, 0.15)) +
  theme_light() +
  labs(title = "Group 1 Enrichment") +
  xlab("Brain node clusters") +
  ylab("Brain node clusters") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, angle = 45, hjust = 1),
    plot.title = element_text(size = 16),
    legend.position = "none"
  )

# Group 2

mat_s2 <- data.frame(
  x = rep(c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"), times = 4),
  y = rep(c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"), each = 4),
  value = as.vector(est_signal[, , 2] - ave_s)
)
mat_s2$x <- factor(mat_s2$x, levels = c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"))
mat_s2$y <- factor(mat_s2$y, levels = c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"))


gg_s2 <- ggplot(data = mat_s2, aes(x, y)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-0.15, 0.15)) +
  theme_light() +
  labs(title = "Group 2 Enrichment") +
  xlab("Brain node clusters") +
  ylab("Brain node clusters") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 16),
    # axis.title.y=element_text(size=16),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, angle = 45, hjust = 1),
    plot.title = element_text(size = 16),
    legend.position = "none"
  )

# Group 3
mat_s3 <- data.frame(
  x = rep(c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"), times = 4),
  y = rep(c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"), each = 4),
  value = as.vector(est_signal[, , 3] - ave_s)
)
mat_s3$x <- factor(mat_s3$x, levels = c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"))
mat_s3$y <- factor(mat_s3$y, levels = c("L.Hemis", "R.Temporal", "R.Occiptial", "R.Supra"))


gg_s3 <- ggplot(data = mat_s3, aes(x, y)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-0.15, 0.15)) +
  theme_light() +
  labs(title = "Group 3 Enrichment") +
  xlab("Brain node clusters") +
  ylab("Brain node clusters") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 16),
    # axis.title.y=element_text(size=16),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, angle = 45, hjust = 1),
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

pdf("figure/Figure10.pdf", width = 24, height = 6)
(gg_ave_s | gg_s1 | gg_s2 | gg_s3)
dev.off()

# preparing csv files for Figure 9 ------------------------------------------------------------------------------------------------

node_file <- read.table("presaved/for_Figure9_10/Desikan_dti.node.txt", header = F)
node_file[, 6] <- colnames(HCP$tensor[, , 1])
node_file[, 6] <- gsub("\\_.*", "", node_file[, 6])
node_file[which(node_file[, 6] == "RMF"), 6] <- "SupF"
node_file[, 4] <- result$z[[1]]
write.table(node_file, "presaved/for_Figure9_11/hpc_node.txt", row.names = F, col.names = F, quote = F)

edge_file <- diag(rep(1, 68))
write.table(edge_file, "presaved/for_Figure9_11/hpc_edge_null.txt", row.names = F, col.names = F)


# preparing csv files for Figure 11 ------------------------------------------------------------------------------------------------

# Average

p <- 1:136
ave_X <- matrix(0, nrow = 68, ncol = 68)

for (c in p) {
  ave_X <- ave_X + HCP$tensor[, , c]
}
ave_X <- ave_X / length(p)

index <- which(ave_X < 1 & 0.01 < ave_X, arr.ind = T)
edge_file <- ave_X
edge_file[index] <- 0

write.table(edge_file, "presaved/for_Figure9_11/hpc_edge.txt", row.names = F, col.names = F)


# Group 1

p <- 1:136
p1 <- p[result$z[[3]] == 1]
ave_X1 <- matrix(0, nrow = 68, ncol = 68)
for (c in p1) {
  ave_X1 <- ave_X1 + HCP$tensor[, , c]
}
ave_X1 <- ave_X1 / length(p1)

edge_file <- ave_X1 - ave_X

write.table(edge_file, "presaved/for_Figure9_11/hpc_edge1.txt", row.names = F, col.names = F)


# Group 2

p <- 1:136
p2 <- p[result$z[[3]] == 2]
ave_X2 <- matrix(0, nrow = 68, ncol = 68)
for (c in p2) {
  ave_X2 <- ave_X2 + HCP$tensor[, , c]
}
ave_X2 <- ave_X2 / length(p2)

edge_file <- ave_X2 - ave_X

write.table(edge_file, "presaved/for_Figure9_11/hpc_edge2.txt", row.names = F, col.names = F)


# Group 3

p <- 1:136
p3 <- p[result$z[[3]] == 3]
ave_X3 <- matrix(0, nrow = 68, ncol = 68)
for (c in p3) {
  ave_X3 <- ave_X3 + HCP$tensor[, , c]
}
ave_X3 <- ave_X3 / length(p3)

edge_file <- ave_X3 - ave_X

write.table(edge_file, "presaved/for_Figure9_11/hpc_edge3.txt", row.names = F, col.names = F)
