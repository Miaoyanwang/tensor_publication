# Code for Figure 8. Nations data analysis

# load dependencies ---------------------------------------------------------------------------------------------
install.packages("software/tensorregress.tar.gz", repos = NULL, type = "source")
library(tensorregress)

library(lattice)
library(RColorBrewer)

library(pracma)
source("simulation.R")

# load data -----------------------------------------------------------------------------------------------------
data <- load("rawdata/nations.RData")

# start analysis ------------------------------------------------------------------------------------------------
set.seed(0)
seed <- 0

tsr <- R
tsr[is.na(tsr)] <- 0
X <- cov

# If a new run of rank selection is desired, run the following two lines
# rank_range=valid_rank(expand.grid(c(3:5),c(3:5),c(3:5)))
# res=sele_rank(tsr,X_covar1=cov,X_covar2=cov,X_covar3=NULL,rank_range=rank_range,niter=10,cons = 'non',dist="binary",initial = "random")

core_shape <- c(4, 4, 4)
result <- tensor_regress(tsr, X_covar1 = cov, X_covar2 = cov, X_covar3 = NULL, core_shape, niter = 10, cons = "penalty", lambda = 0.1, alpha = 10, solver = "CG", dist = "binary", initial = "random")

save(result, file = "presaved/output_nations.RData")
clustering <- kmeans(result$W[[3]], 4, nstart = 5)
relnames <- names(tsr[1, 1, ])
relnames[clustering$cluster == 1]
relnames[clustering$cluster == 2]
relnames[clustering$cluster == 3]
relnames[clustering$cluster == 4]


# start plotting ------------------------------------------------------------------------------------------------
my.palette <- brewer.pal(n = 11, name = "RdBu")
myPanel <- function(x, y, z, ...) {
  panel.levelplot(x, y, z, ...)
  # use handy matrix indexing
  panel.text(x, y, round(M[cbind(x, y)], 2)) 
}

pdf("Figure8_1.pdf", width = 6, height = 6)
M <- result$C_ts[, , which(relnames == "warning")]
levelplot(M, panel = myPanel, col.regions = my.palette, at = seq(-1.6, 1.6, length.out = 12), main = "Relation: warning")
dev.off()

pdf("Figure8_2.pdf", width = 6, height = 6)
M <- result$C_ts[, , which(relnames == "violentactions")]
levelplot(M, panel = myPanel, col.regions = my.palette, at = seq(-1.6, 1.6, length.out = 12), main = "Relation: violentactions")
dev.off()

pdf("Figure8_3.pdf", width = 6, height = 6)
M <- result$C_ts[, , which(relnames == "treaties")]
levelplot(M, panel = myPanel, col.regions = my.palette, at = seq(-1.6, 1.6, length.out = 12), main = "Relation: treaties")
dev.off()

pdf("Figure8_4.pdf", width = 6, height = 6)
M <- result$C_ts[, , which(relnames == "aidenemy")]
levelplot(M, panel = myPanel, col.regions = my.palette, at = seq(-1.6, 1.6, length.out = 12), main = "Relation: aidenemy")
dev.off()
