library(parallel)
library(doParallel)
library(foreach)

tryCatch({parallel::stopCluster(cl)}, error = function(e) {})
closeAllConnections()

total_seeds <- 100 
lambda_reg <- 0.001 

K_gauss <- function(d, bw) dnorm(d, mean = 0, sd = bw)

num_cores <- min(64, parallel::detectCores() - 1)
cl <- parallel::makeCluster(num_cores)
registerDoParallel(cl)

  sim_results <- foreach(s = 1:total_seeds, .packages=c("stats")) %dopar% {
  tryCatch({
    
    lambda_reg <- 0.001 
    max_rho_soft <- 10.0  
    K_gauss <- function(d, bw) dnorm(d, mean = 0, sd = bw)
    
    set.seed(s)
    
    N <- 1000
    r <- rbinom(N, 1, 0.5)
    
    n0 <- sum(r == 1) 
    n1 <- N - n0      
    pi_val <- n0 / N
    
    y_s_orig <- rnorm(n0, 0, sqrt(2))
    y_t_orig <- rnorm(n1, 1, 1)
    Y_orig <- c(y_s_orig, y_t_orig)
    
    x1_orig <- -0.5 * Y_orig + rnorm(N, 0, 1)
    x2_orig <-  0.5 * Y_orig + rnorm(N, 0, 1)
    x3_orig <-  Y_orig + rnorm(N, 0, 1)
    X_orig <- cbind(x1_orig, x2_orig, x3_orig)
    
    x_s_orig <- X_orig[1:n0, ]
    x_t_orig <- X_orig[(n0+1):N, ]
    