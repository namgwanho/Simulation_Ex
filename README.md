# Unsupervised domain adaptation with Label Shift

## Overview
 Distribution shifts are occured between source and target domains. If the target labels ($Y$) are observed, the label shift problem can be addressed simply by calculating the density ratio $\rho$.
 
  However, in Unsupervised Domain Adaptation scenarios where target labels remain completely unobserved, handling label shift becomes a highly difficult problem.
  
 The objective of my algoruthm is to discover the optimal model parameters for the target distribution even when target labels $Y$ are completely unobserved by leveraging the Efficient Influence Function based on a doubly flexible estimator to robustly combine an imperfect prediction model and the density ratio $\rho$
  
## Method
This simulation study compares the following methodologies
*   Naive
*   IPW(importance wieghted estimator)
*   Mine
*   Oracle
  
\-   Transitions the Doubly flexible EIF from simple scalar estimation into an estimating equation for model parameter $\theta$ optimization. <br>
\-   Injects the loss gradient of the prediction model into the shift-corrected EIF to construct a $\psi$ for parameter updates. <br>
\-   Solves $\sum \psi(\theta) = 0$ via optimization, ultimately converging to the optimal parameter $\theta$ of the unobserved target distribution. 

## Simulation Settings

| Parameter | Value / Distribution | Description |
| :--- | :--- | :--- |
| **Total Seeds** | `100` | Number of independent simulation iterations |
| **$N$** | `1000` | Total sample size |
| **$n_0$** | $r \sim \text{Binomial}(N, 0.5)$ | Number of Source domain data ($\approx 500$) |
| **$n_1$** | $N - n_0$ | Number of Target domain data ($\approx 500$) |
| **$Y_s$** | $\mathcal{N}(0, 2)$ | Source outcome distribution |
| **$Y_t$** | $\mathcal{N}(1, 1)$ | Target outcome distribution |
| **$X_1$** | $-0.5Y + \epsilon$ | Feature 1 ($\epsilon$ ~ N(0, 1)) |
| **$X_2$** | $0.5Y + \epsilon$ | Feature 2 ($\epsilon$ ~ N(0, 1)) |
| **$X_3$** | $Y + \epsilon$ | Feature 3 ($\epsilon$ ~ N(0, 1)) |
| **$\theta$** | $(b_0, b_1, b_2, b_3)^\top$ | Target parameters by Linear Regression |

## Simulation Results

| Metric | Naive | IPW | Mine | ORACLE |
| :--- | :--- | :--- | :--- | :--- |
| **beta 0** | 0.0001 | 0.3371 | 0.4014 | 0.4024 |
| **beta 1** | -0.2474 | -0.2159 | -0.2089 | -0.2063 |
| **beta 2** | 0.2477 | 0.2124 | 0.1925 | 0.1992 |
| **beta 3** | 0.5031 | 0.4294 | 0.4033 | 0.4024 |
| **RMSE** | 0.7045 | 0.6357 | 0.6650 | 0.6240 |
| **MAE** | 0.5637 | 0.5068 | 0.5312 | 0.4950 |

## Repository Structure
Mine sim.R - simulation code for this algorithm

## Requiremets
R

## Research Status

This repository contains ongoing and unpublished research.

The code is provided to document the implementation and experimental development of the project. Full theoretical derivations, detailed algorithms, simulation settings, and complete numerical results are intentionally omitted at this stage.

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/
