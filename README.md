# unsupervied domain adaption under label shift

## Background 
*   In real-world, distribution shifts frequently occur between source and target domains. If the target labels ($Y$) are observed, the label shift problem can be addressed simply by calculating the density ratio $\rho$. However, this repository focuses on Unsupervised Domain Adaptation (UDA), a much more practical and challenging setting where the target labels remain completely unobserved.

## Simulation Setting
```R
set.seed(1)
N <- 1000

# 1. Indicator (Source: r=1, Target: r=0)
r <- rbinom(N, 1, 0.5)
n0 <- sum(r == 1)
n1 <- N - n0
pi_val <- n0 / N

# 2. DATA Generation
y_s_orig <- rnorm(n0, mean = 0, sd = sqrt(2))  # Source Y
y_t_orig <- rnorm(n1, mean = 1, sd = 1)        # Target Y
Y_orig <- c(y_s_orig, y_t_orig)

# 3. Feature Generation
x1_orig <- -0.5 * Y_orig + rnorm(N, 0, 1)
x2_orig <-  0.5 * Y_orig + rnorm(N, 0, 1)
x3_orig <-        Y_orig + rnorm(N, 0, 1)
X_orig <- cbind(x1_orig, x2_orig, x3_orig)
```

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
