# unsupervided domain adaption under label shift

## Background 
*   **Label Shift and the Instability of IPW:** Under $P(Y)$ shift (Label Shift) where the label distribution differs between the source and target domains, Inverse Probability Weighting (IPW) is traditionally utilized. However, it suffers from a critical vulnerability: even a minuscule error in estimating the density ratio can cause the entire prediction performance to collapse.
*   **Limitations of Traditional Doubly Robust (DR) Frameworks:** Existing DR frameworks are limited to estimating a **single scalar statistic** (e.g., population mean or median) using the Efficient Influence Function (EIF). They possess a clear limitation: they cannot be used to directly train the parameters (weights) of a prediction model that takes new data $X$ as input to predict $Y$.
*   **The Paradigm Shift of the DEAL Algorithm:** DEAL maintains the structural skeleton of the DR equation but replaces the target variable (traditionally used to find a simple mean) with the gradient of the prediction model's loss function. In other words, it is a universal framework that re-engineers the EIF equation into a **pseudo-gradient ($\psi$) engine for error correction** to update machine learning parameters.
*   **Safely Controlling the Fredholm Integral Equation:** The Fredholm integral equation of the first kind, which is used to profile the nuisance function, is a notoriously ill-posed problem where even minor noise can cause the solution to explode to infinity. DEAL fundamentally prevents parameter collapse by not tying this equation directly to the deep learning optimizer. Instead, it securely stabilizes the solution in advance through matrix profiling and Ridge regularization.

## 🚀 Core Contribution of this Repository
This repository conducts parallel simulations to comparatively verify three strategies for perfectly controlling the numerical explosion caused when the denominator approaches zero during the internal matrix ($H$) inversion process of the DEAL algorithm.
*   **PureH:** Completely excludes safety mechanisms (clipping) to observe the pure algorithm's collapse threshold.
*   **Clip0.01:** Forcibly applies hard clipping to analyze the limits of bias generation.
*   **GLM (Quasipoisson):** Proposes an optimal numerical stabilization technique that maximizes prediction performance by introducing a Log-link function to inherently prevent negative values and smooth the weight curve.

## ⚙️ Simulation Setting & Data Generating Process (DGP)
To simulate a Label Shift scenario, data is generated as follows. The invariant assumption $P(X \mid Y)$ is strictly maintained across both domains, while the distribution (mean, variance) of $P(Y)$ is explicitly altered to induce a domain mismatch.

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

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
