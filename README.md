# unsupervied domain adaption under label shift

## Background 
*   In real-world, distribution shifts frequently occur between source and target domains. If the target labels ($Y$) are observed, the label shift problem can be addressed simply by calculating the density ratio $\rho$. However, this repository focuses on Unsupervised Domain Adaptation (UDA), a much more practical and challenging setting where the target labels remain completely unobserved.

## Simulation Setting
```R
# 0. Set Seed & Sample Size
set.seed(42)
N <- 1000

# 1. Domain Indicator (Source: r=1, Target: r=0)
r <- rbinom(N, 1, 0.5)
n0 <- sum(r == 1)
n1 <- N - n0
pi_val <- n0 / N

# 2. Label Shift Generation (Target Outcome has different Mean/Var)
y_s_orig <- rnorm(n0, mean = 0, sd = sqrt(2))  # Source Y
y_t_orig <- rnorm(n1, mean = 1, sd = 1)        # Target Y
Y_orig <- c(y_s_orig, y_t_orig)

# 3. Feature Generation (P(X|Y) remains invariant)
x1_orig <- -0.5 * Y_orig + rnorm(N, 0, 1)
x2_orig <-  0.5 * Y_orig + rnorm(N, 0, 1)
x3_orig <-        Y_orig + rnorm(N, 0, 1)
X_orig <- cbind(x1_orig, x2_orig, x3_orig)

# 4. Target Task (3D Multiple Linear Regression)
# y = b0 + b1*x1 + b2*x2 + b3*x3
x_s_orig <- X_orig[1:n0, ]
x_t_orig <- X_orig[(n0+1):N, ]
```

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
