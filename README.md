# Flexible Analysis of Spatial Transcriptomics (FAST): A Deconvolution Approach
## Authors: Meng Zhang; Joel Parker; Xiaoxiao Sun

### Example codes
load("Simulation.RData")
load("Adjacent.RData")
X <- output[[1]]
true_props <-  output[[2]]


config <- dmain_config
config$r <- 2
config$lambda_1 <- 0.05
config$lambda_2 <- 0.1
config$n_iter <- 1000
X <- log(X+1)
X <- X/max(X)
set.seed(1)
system.time(res<-dmain(X,adj_final,config))
cor(res$H[,1], true_props[,1])