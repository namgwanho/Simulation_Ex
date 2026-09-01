# DEAL: Debiased learning for unsupervised domain adaptation under label shift

## Overview
 Distribution shifts are occured between source and target domains. If the target labels ($Y$) are observed, the label shift problem can be addressed simply by calculating the density ratio $\rho$.
 
  However, in Unsupervised Domain Adaptation scenarios where target labels remain completely unobserved, handling label shift becomes a highly difficult problem.
  
 The objective of DEAL is to discover the optimal model parameters for the target distribution even when target labels ($Y$) are completely unobserved by leveraging the Efficient Influence Function based on a doubly flexible estimator to robustly combine an imperfect prediction model and the density ratio ($\rho$)
  
## Method
This simulation study compares the following methodologies
*   Naive
*   IPW(importance wieghted estimator
*   Deal
*   Oracle
  
\-   Transitions the Doubly flexible EIF from simple scalar estimation into an estimating equation for model parameter ($\theta$) optimization. <br>
\-   Injects the loss gradient of the prediction model into the shift-corrected EIF to construct a pseudo-gradient ($\psi$) for parameter updates. <br>
\-   Solves $\sum \psi(\theta) = 0$ via optimization, ultimately converging to the optimal parameter $\theta$ of the unobserved target distribution. 

## Simulation Settings

| Parameter | Value / Distribution | Description |
| :--- | :--- | :--- |
| **Total Seeds** | `100` | Number of independent simulation iterations |
| **$N$** | `1000` | Total sample size |
| **$n_0$** | $r \sim \text{Binomial}(N, 0.5)$ | Number of Source domain data ($\approx 500$) |
| **$n_1$** | $N - n_0$ | Number of Target domain data ($\approx 500$) |
| **$Y_s$** | $\mathcal{N}(0, 2)$ | Source outcome distribution (sd = $\sqrt{2}$) |
| **$Y_t$** | $\mathcal{N}(1, 1)$ | Target outcome distribution (sd = $1$) |
| **$X_1$** | $-0.5Y + \epsilon$ | Feature 1 ($\epsilon \sim \mathcal{N}(0, 1)$) |
| **$X_2$** | $0.5Y + \epsilon$ | Feature 2 ($\epsilon \sim \mathcal{N}(0, 1)$) |
| **$X_3$** | $Y + \epsilon$ | Feature 3 ($\epsilon \sim \mathcal{N}(0, 1)$) |
| **$\theta$** | $(b_0, b_1, b_2, b_3)^\top$ | Target parameters for 3D Multiple Linear Regression |

## Simulation Results

- To ensure stable convergence, the initial parameter $\theta$ for the DEAL optimization was initialized using the estimates derived from the baseline IPW (Inverse Probability Weighting) regression model.
- The figure below compares the parameter estimation accuracy and target prediction errors (RMSE, MAE) across the three denominator stabilization strategies (PureH, Clip0.01, GLM).

![Simulation Results](./image_a9bda6.png)


## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
