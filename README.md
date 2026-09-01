# DEAL: Debiased learning for unsupervised domain adaptation under label shift

## Overview
*   In real-world, distribution shifts frequently occur between source and target domains. If the target labels ($Y$) are observed, the label shift problem can be addressed simply by calculating the density ratio $\rho$. However, this repository focuses on Unsupervised Domain Adaptation (UDA), a much more practical and challenging setting where the target labels remain completely unobserved.
*   The objective of DEAL is to discover the optimal model parameters for the target distribution even when target labels ($Y$) are completely unobserved by leveraging the Efficient Influence Function based on a doubly flexible estimator to robustly combine an imperfect prediction model and the density ratio ($\rho$)
  
## Method
This simulation study compares the following methodologies
*   Naive
*   IPW(importance wieghted estimator
*   Deal
*   Oracle
  
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
## Method

- Transitions the Doubly flexible EIF from simple scalar estimation into an estimating equation for model parameter ($\theta$) optimization.
- Injects the loss gradient of the prediction model into the shift-corrected EIF to construct a pseudo-gradient ($\psi$) for parameter updates.
- Solves $\sum \psi(\theta) = 0$ via optimization, ultimately converging to the optimal parameter $\theta$ of the unobserved target distribution.

## Simulation Results

- To ensure stable convergence, the initial parameter $\theta$ for the DEAL optimization was initialized using the estimates derived from the baseline IPW (Inverse Probability Weighting) regression model.
- The figure below compares the parameter estimation accuracy and target prediction errors (RMSE, MAE) across the three denominator stabilization strategies (PureH, Clip0.01, GLM).

![Simulation Results](./image_a9bda6.png)


## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
