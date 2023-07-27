# Flexible Analysis of Spatial Transcriptomics (FAST): A Deconvolution Approach

## Authors: Meng Zhang; Joel Parker; Xiaoxiao Sun

## Example codes

### Install
    library("devtools")
    install_github("shawnstat/FAST/FAST")
    library(FAST)
    data(simulate)
### Run FAST
    X <- simulated_data$X
    true_props <-  simulated_data$true_proportions
    config <- dmain_config
    config$r <- 2
    config$lambda_1 <- 0.05
    config$lambda_2 <- 0.1
    config$n_iter <- 1000
    X <- log(X+1)
    X <- X/max(X)
    set.seed(1)
    system.time(res<-dmain(X,simulated_data$adjacent,config))
    cor(res$H[,1], true_props[,1])
