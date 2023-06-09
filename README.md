# FAST

load("Simulation.RData")
load("Adjacent.RData")
X <- simu_3ct[[1]]
true_props <-  simu_3ct[[2]]

set.seed(2023)
config <- dmain_config
config$r <- 3
config$lambda_1 <- 0.5
config$lambda_2 <- 1.0
config$n_iter <- 2000
system.time(res<-dmain(log(X+1),adj_final,config))

