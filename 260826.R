# =====================================================================
# 🔥 [최종 완벽본] 1D Label Shift: 3가지 분모 방식 완전 분리 실험 및 종합
#    - 상황 1: PureH (어떠한 클리핑도 없는 완벽한 순수 H행렬)
#    - 상황 2: Clip001 (클리핑 0.01)
#    - 상황 3: GLM (Quasipoisson Log-link)
#    - 평가: Naive, IPW(S/T/H), DEAL(S/T/H), Oracle
# =====================================================================
library(parallel)
library(doParallel)
library(foreach)

tryCatch({parallel::stopCluster(cl)}, error = function(e) {})
closeAllConnections()

total_seeds <- 100 
num_cores <- min(64, parallel::detectCores() - 1)

cat("=============================================================\n")
cat("🔥 3가지 상황(PureH, Clip0.01, GLM) 분리 실험을 시작합니다...\n")
cat("=============================================================\n")


# =====================================================================
# 🟩 [1단계] 순수 H행렬 (PureH) 단독 실험 - 💡 모든 클리핑 완전 제거
# =====================================================================
cat("\n▶ [1/3] PureH (순수 H행렬) 실험 진행 중...\n")
cl <- parallel::makeCluster(num_cores)
registerDoParallel(cl)

sim_pureh <- foreach(s = 1:total_seeds, .packages=c("stats")) %dopar% {
  tryCatch({
    lambda_reg <- 0.001; max_rho_soft <- 10.0
    K_gauss <- function(d, bw) dnorm(d, mean = 0, sd = bw)
    set.seed(s)
    
    # 0. 데이터 생성 (원본 스케일)
    N <- 1000; r <- rbinom(N, 1, 0.5); n0 <- sum(r == 1); n1 <- N - n0; pi_val <- n0 / N
    y_s_orig <- rnorm(n0, 0, sqrt(2)); y_t_orig <- rnorm(n1, 1, 1); Y_orig <- c(y_s_orig, y_t_orig)
    x1_orig <- -0.5 * Y_orig + rnorm(N, 0, 1); x2_orig <- 0.5 * Y_orig + rnorm(N, 0, 1); x3_orig <- Y_orig + rnorm(N, 0, 1)
    X_orig <- cbind(x1_orig, x2_orig, x3_orig)
    x_s_orig <- X_orig[1:n0, ]; x_t_orig <- X_orig[(n0+1):N, ]; y_grid <- seq(-5, 5, length.out = 100) 
    
    # 1. 3-Stage Rho Evolution
    phi_temp <- cbind(1, x_s_orig)
    beta_temp <- solve(t(phi_temp) %*% phi_temp + 1e-4 * diag(4)) %*% t(phi_temp) %*% y_s_orig
    y_t_hat_orig <- cbind(1, x_t_orig) %*% beta_temp
    
    f_s <- approxfun(density(y_s_orig, kernel="gaussian")$x, density(y_s_orig, kernel="gaussian")$y, rule=2)
    f_t <- approxfun(density(y_t_hat_orig, kernel="gaussian")$x, density(y_t_hat_orig, kernel="gaussian")$y, rule=2)
    
    rho_star_grid <- pmax(f_t(y_grid) / (f_s(y_grid) + 1e-6), 0.01)
    rho_star_grid <- max_rho_soft * tanh(rho_star_grid / max_rho_soft); rho_star_grid <- rho_star_grid / mean(rho_star_grid)
    rho_star_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_star_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_star_final <- rho_star_final / mean(rho_star_final)
    
    h_ml <- 3 * n0^(-1/7); l_bw <- 1.5 * n0^(-1/3); h_bw <- 0.5 * n0^(-1/16)
    W_ML_unnorm <- matrix(0, N, n0)
    for (i in 1:N) for (j in 1:n0) W_ML_unnorm[i, j] <- K_gauss(X_orig[i,1]-x_s_orig[j,1], h_ml) * K_gauss(X_orig[i,2]-x_s_orig[j,2], h_ml) * K_gauss(X_orig[i,3]-x_s_orig[j,3], h_ml)
    W_ML <- W_ML_unnorm / (rowSums(W_ML_unnorm) + 1e-10)
    W_Y_unnorm <- outer(y_s_orig, y_s_orig, function(y1, y2) K_gauss(y1 - y2, l_bw)); W_Y <- W_Y_unnorm / (rowSums(W_Y_unnorm) + 1e-10)
    V_mat_grid <- outer(y_s_orig, y_grid, function(y_i, y_0) K_gauss(y_i - y_0, h_bw))
    
    solve_rho_grid <- function(rho_guess_grid) {
      rgs <- approx(y_grid, rho_guess_grid, xout=y_s_orig, rule=2)$y
      w <- 1 / as.vector(W_ML %*% (rgs^2 + (pi_val/(1-pi_val))*rgs))
      M <- W_Y %*% diag(w[1:n0]) %*% W_ML[1:n0, ] %*% diag(rgs)
      A <- solve(t(M)%*%M + 1e-3*diag(n0)) %*% t(M) %*% V_mat_grid
      E <- W_ML %*% diag(rgs) %*% A
      t1 <- colMeans(w[(n0+1):N] * E[(n0+1):N, , drop=FALSE]); t2 <- colMeans(rgs * (V_mat_grid - w[1:n0]*E[1:n0, , drop=FALSE]))
      (t1 + t2) / (colMeans(V_mat_grid) + 1e-8)
    }
    
    rho_tilde_grid <- pmax(0.01, pmin(solve_rho_grid(rho_star_grid), 20.0))
    rho_tilde_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_tilde_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_tilde_final <- rho_tilde_final / mean(rho_tilde_final)
    
    rho_hat_grid <- solve_rho_grid(rho_tilde_grid); rho_hat_grid[is.na(rho_hat_grid) | is.infinite(rho_hat_grid)] <- 1.0
    rho_hat_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_hat_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_hat_final <- rho_hat_final / mean(rho_hat_final)
    
    # 2. Base Models
    theta_naive  <- as.numeric(coef(lm(y_s_orig ~ x_s_orig)))
    theta_ipw_s  <- as.numeric(coef(lm(y_s_orig ~ x_s_orig, weights = rho_star_final)))
    theta_ipw_t  <- as.numeric(coef(lm(y_s_orig ~ x_s_orig, weights = rho_tilde_final)))
    theta_ipw_h  <- as.numeric(coef(lm(y_s_orig ~ x_s_orig, weights = rho_hat_final)))
    theta_oracle <- as.numeric(coef(lm(y_t_orig ~ x_t_orig)))
    
    calc_metrics <- function(theta_est) {
      if(any(is.na(theta_est))) return(c(NA, NA, NA))
      err <- as.vector(y_t_orig - cbind(1, x_t_orig) %*% theta_est)
      c(mean(err^2), sqrt(mean(err^2)), mean(abs(err)))
    }
    
    metrics_naive <- calc_metrics(theta_naive); metrics_ipw_s <- calc_metrics(theta_ipw_s)
    metrics_ipw_t <- calc_metrics(theta_ipw_t); metrics_ipw_h <- calc_metrics(theta_ipw_h); metrics_oracle <- calc_metrics(theta_oracle)
    
    # 3. DEAL Engine (PureH: 💡 안전장치 1e-10조차 완전 제거)
    Phi_orig <- cbind(1, x_s_orig); Phi_all_orig <- cbind(1, X_orig); XtX_orig <- t(Phi_orig) %*% Phi_orig
    penalty_diag <- diag(ncol(Phi_orig)); penalty_diag[1,1] <- 0 
    
    H_orig <- Phi_orig %*% solve(XtX_orig + 1e-4 * penalty_diag) %*% t(Phi_orig)
    H_all_orig <- Phi_all_orig %*% solve(XtX_orig + 1e-4 * penalty_diag) %*% t(Phi_orig)
    K_h_deal <- exp(-(as.matrix(dist(y_s_orig))^2) / (2*(1.06*sd(y_s_orig)*(n0^(-1/5)))^2)); K_h_deal <- K_h_deal / rowSums(K_h_deal)
    
    run_deal_pureh <- function(rho_val) {
      P <- diag(rho_val)
      d <- H_orig %*% (P %*% rho_val + (pi_val / (1 - pi_val)) * rho_val)
      d_all <- H_all_orig %*% (P %*% rho_val + (pi_val / (1 - pi_val)) * rho_val)
      
      # 💡 인위적 클리핑, 최소한의 안전장치(1e-10) 모두 제거. 순수한 역수 취함.
      W <- diag(1 / as.vector(d))
      W_all <- diag(1 / as.vector(d_all))
      
      L_h_deal <- K_h_deal %*% W %*% H_orig %*% P
      M_h_lambda <- solve(t(L_h_deal) %*% L_h_deal + lambda_reg * diag(n0)) %*% t(L_h_deal)
      
      compute_psi <- function(theta) {
        U <- -1 * (as.vector(y_s_orig - Phi_orig %*% theta) * Phi_orig)  
        R_h <- K_h_deal %*% (diag(n0) - W %*% H_orig %*% (P %*% P)) %*% U
        B_h_lambda <- W_all %*% H_all_orig %*% ((P %*% P) %*% U + P %*% (M_h_lambda %*% R_h))
        psi_hat <- matrix(0, nrow=N, ncol=ncol(Phi_orig))
        psi_hat[1:n0, ] <- (1/pi_val)*rho_val*(U - B_h_lambda[1:n0, ]); psi_hat[(n0+1):N, ] <- (1/(1-pi_val))*B_h_lambda[(n0+1):N, ]
        colMeans(psi_hat)
      }
      psi_0 <- compute_psi(rep(0, ncol(Phi_orig))); J <- matrix(0, ncol(Phi_orig), ncol(Phi_orig))
      for (j in 1:ncol(Phi_orig)) { e_vec <- rep(0, ncol(Phi_orig)); e_vec[j] <- 1; J[, j] <- compute_psi(e_vec) - psi_0 }
      if (any(is.na(J)) || abs(det(J)) < 1e-6) return(rep(NA, 4))
      as.numeric(solve(J + 1e-3 * penalty_diag, -psi_0))
    }
    
    theta_deal_s <- run_deal_pureh(rho_star_final)
    theta_deal_T <- run_deal_pureh(rho_tilde_final)
    theta_deal_H <- run_deal_pureh(rho_hat_final)
    
    metrics_deal_s <- calc_metrics(theta_deal_s)
    metrics_deal_T <- calc_metrics(theta_deal_T)
    metrics_deal_H <- calc_metrics(theta_deal_H)
    
    res_df <- data.frame(
      Seed = s,
      b0_Naive = theta_naive[1], b1_Naive = theta_naive[2], b2_Naive = theta_naive[3], b3_Naive = theta_naive[4],
      b0_IPW_s = theta_ipw_s[1], b1_IPW_s = theta_ipw_s[2], b2_IPW_s = theta_ipw_s[3], b3_IPW_s = theta_ipw_s[4],
      b0_IPW_T = theta_ipw_t[1], b1_IPW_T = theta_ipw_t[2], b2_IPW_T = theta_ipw_t[3], b3_IPW_T = theta_ipw_t[4],
      b0_IPW_H = theta_ipw_h[1], b1_IPW_H = theta_ipw_h[2], b2_IPW_H = theta_ipw_h[3], b3_IPW_H = theta_ipw_h[4],
      b0_Oracle = theta_oracle[1], b1_Oracle = theta_oracle[2], b2_Oracle = theta_oracle[3], b3_Oracle = theta_oracle[4],
      b0_DEAL_s_PureH = theta_deal_s[1], b1_DEAL_s_PureH = theta_deal_s[2], b2_DEAL_s_PureH = theta_deal_s[3], b3_DEAL_s_PureH = theta_deal_s[4],
      b0_DEAL_T_PureH = theta_deal_T[1], b1_DEAL_T_PureH = theta_deal_T[2], b2_DEAL_T_PureH = theta_deal_T[3], b3_DEAL_T_PureH = theta_deal_T[4],
      b0_DEAL_H_PureH = theta_deal_H[1], b1_DEAL_H_PureH = theta_deal_H[2], b2_DEAL_H_PureH = theta_deal_H[3], b3_DEAL_H_PureH = theta_deal_H[4],
      MSE_Naive = metrics_naive[1], RMSE_Naive = metrics_naive[2], MAE_Naive = metrics_naive[3],
      MSE_IPW_s = metrics_ipw_s[1], RMSE_IPW_s = metrics_ipw_s[2], MAE_IPW_s = metrics_ipw_s[3],
      MSE_IPW_T = metrics_ipw_t[1], RMSE_IPW_T = metrics_ipw_t[2], MAE_IPW_T = metrics_ipw_t[3],
      MSE_IPW_H = metrics_ipw_h[1], RMSE_IPW_H = metrics_ipw_h[2], MAE_IPW_H = metrics_ipw_h[3],
      MSE_Oracle = metrics_oracle[1], RMSE_Oracle = metrics_oracle[2], MAE_Oracle = metrics_oracle[3],
      MSE_DEAL_s_PureH = metrics_deal_s[1], RMSE_DEAL_s_PureH = metrics_deal_s[2], MAE_DEAL_s_PureH = metrics_deal_s[3],
      MSE_DEAL_T_PureH = metrics_deal_T[1], RMSE_DEAL_T_PureH = metrics_deal_T[2], MAE_DEAL_T_PureH = metrics_deal_T[3],
      MSE_DEAL_H_PureH = metrics_deal_H[1], RMSE_DEAL_H_PureH = metrics_deal_H[2], MAE_DEAL_H_PureH = metrics_deal_H[3]
    )
    list(res = res_df)
  }, error = function(e) { NULL })
}
parallel::stopCluster(cl)

results_pureh <- do.call(rbind, lapply(sim_pureh, function(x) x$res))
cat("✅ PureH 완료!\n")


# =====================================================================
# 🟩 [2단계] H행렬 + Clip 0.01 단독 실험
# =====================================================================
cat("\n▶ [2/3] Clip0.01 실험 진행 중...\n")
cl <- parallel::makeCluster(num_cores)
registerDoParallel(cl)

sim_clip001 <- foreach(s = 1:total_seeds, .packages=c("stats")) %dopar% {
  tryCatch({
    lambda_reg <- 0.001; max_rho_soft <- 10.0
    K_gauss <- function(d, bw) dnorm(d, mean = 0, sd = bw)
    set.seed(s)
    
    N <- 1000; r <- rbinom(N, 1, 0.5); n0 <- sum(r == 1); n1 <- N - n0; pi_val <- n0 / N
    y_s_orig <- rnorm(n0, 0, sqrt(2)); y_t_orig <- rnorm(n1, 1, 1); Y_orig <- c(y_s_orig, y_t_orig)
    x1_orig <- -0.5 * Y_orig + rnorm(N, 0, 1); x2_orig <- 0.5 * Y_orig + rnorm(N, 0, 1); x3_orig <- Y_orig + rnorm(N, 0, 1)
    X_orig <- cbind(x1_orig, x2_orig, x3_orig)
    x_s_orig <- X_orig[1:n0, ]; x_t_orig <- X_orig[(n0+1):N, ]; y_grid <- seq(-5, 5, length.out = 100) 
    
    phi_temp <- cbind(1, x_s_orig); beta_temp <- solve(t(phi_temp) %*% phi_temp + 1e-4 * diag(4)) %*% t(phi_temp) %*% y_s_orig; y_t_hat_orig <- cbind(1, x_t_orig) %*% beta_temp
    f_s <- approxfun(density(y_s_orig, kernel="gaussian")$x, density(y_s_orig, kernel="gaussian")$y, rule=2)
    f_t <- approxfun(density(y_t_hat_orig, kernel="gaussian")$x, density(y_t_hat_orig, kernel="gaussian")$y, rule=2)
    rho_star_grid <- pmax(f_t(y_grid) / (f_s(y_grid) + 1e-6), 0.01)
    rho_star_grid <- max_rho_soft * tanh(rho_star_grid / max_rho_soft); rho_star_grid <- rho_star_grid / mean(rho_star_grid)
    rho_star_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_star_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_star_final <- rho_star_final / mean(rho_star_final)
    
    h_ml <- 3 * n0^(-1/7); l_bw <- 1.5 * n0^(-1/3); h_bw <- 0.5 * n0^(-1/16)
    W_ML_unnorm <- matrix(0, N, n0)
    for (i in 1:N) for (j in 1:n0) W_ML_unnorm[i, j] <- K_gauss(X_orig[i,1]-x_s_orig[j,1], h_ml) * K_gauss(X_orig[i,2]-x_s_orig[j,2], h_ml) * K_gauss(X_orig[i,3]-x_s_orig[j,3], h_ml)
    W_ML <- W_ML_unnorm / (rowSums(W_ML_unnorm) + 1e-10)
    W_Y_unnorm <- outer(y_s_orig, y_s_orig, function(y1, y2) K_gauss(y1 - y2, l_bw)); W_Y <- W_Y_unnorm / (rowSums(W_Y_unnorm) + 1e-10)
    V_mat_grid <- outer(y_s_orig, y_grid, function(y_i, y_0) K_gauss(y_i - y_0, h_bw))
    
    solve_rho_grid <- function(rho_guess_grid) {
      rgs <- approx(y_grid, rho_guess_grid, xout=y_s_orig, rule=2)$y
      w <- 1 / as.vector(W_ML %*% (rgs^2 + (pi_val/(1-pi_val))*rgs))
      M <- W_Y %*% diag(w[1:n0]) %*% W_ML[1:n0, ] %*% diag(rgs)
      A <- solve(t(M)%*%M + 1e-3*diag(n0)) %*% t(M) %*% V_mat_grid
      E <- W_ML %*% diag(rgs) %*% A
      t1 <- colMeans(w[(n0+1):N] * E[(n0+1):N, , drop=FALSE]); t2 <- colMeans(rgs * (V_mat_grid - w[1:n0]*E[1:n0, , drop=FALSE]))
      (t1 + t2) / (colMeans(V_mat_grid) + 1e-8)
    }
    
    rho_tilde_grid <- pmax(0.01, pmin(solve_rho_grid(rho_star_grid), 20.0))
    rho_tilde_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_tilde_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_tilde_final <- rho_tilde_final / mean(rho_tilde_final)
    rho_hat_grid <- solve_rho_grid(rho_tilde_grid); rho_hat_grid[is.na(rho_hat_grid) | is.infinite(rho_hat_grid)] <- 1.0
    rho_hat_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_hat_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_hat_final <- rho_hat_final / mean(rho_hat_final)
    
    # DEAL Engine (Clip001)
    Phi_orig <- cbind(1, x_s_orig); Phi_all_orig <- cbind(1, X_orig); XtX_orig <- t(Phi_orig) %*% Phi_orig
    penalty_diag <- diag(ncol(Phi_orig)); penalty_diag[1,1] <- 0 
    H_orig <- Phi_orig %*% solve(XtX_orig + 1e-4 * penalty_diag) %*% t(Phi_orig)
    H_all_orig <- Phi_all_orig %*% solve(XtX_orig + 1e-4 * penalty_diag) %*% t(Phi_orig)
    K_h_deal <- exp(-(as.matrix(dist(y_s_orig))^2) / (2*(1.06*sd(y_s_orig)*(n0^(-1/5)))^2)); K_h_deal <- K_h_deal / rowSums(K_h_deal)
    
    run_deal_clip <- function(rho_val) {
      P <- diag(rho_val)
      d <- H_orig %*% (P %*% rho_val + (pi_val / (1 - pi_val)) * rho_val)
      d_all <- H_all_orig %*% (P %*% rho_val + (pi_val / (1 - pi_val)) * rho_val)
      W <- diag(1 / pmax(as.vector(d), 0.01)); W_all <- diag(1 / pmax(as.vector(d_all), 0.01)) # 💡 Clip 0.01 적용
      
      L_h_deal <- K_h_deal %*% W %*% H_orig %*% P
      M_h_lambda <- solve(t(L_h_deal) %*% L_h_deal + lambda_reg * diag(n0)) %*% t(L_h_deal)
      
      compute_psi <- function(theta) {
        U <- -1 * (as.vector(y_s_orig - Phi_orig %*% theta) * Phi_orig)  
        R_h <- K_h_deal %*% (diag(n0) - W %*% H_orig %*% (P %*% P)) %*% U
        B_h_lambda <- W_all %*% H_all_orig %*% ((P %*% P) %*% U + P %*% (M_h_lambda %*% R_h))
        psi_hat <- matrix(0, nrow=N, ncol=ncol(Phi_orig))
        psi_hat[1:n0, ] <- (1/pi_val)*rho_val*(U - B_h_lambda[1:n0, ]); psi_hat[(n0+1):N, ] <- (1/(1-pi_val))*B_h_lambda[(n0+1):N, ]
        colMeans(psi_hat)
      }
      psi_0 <- compute_psi(rep(0, ncol(Phi_orig))); J <- matrix(0, ncol(Phi_orig), ncol(Phi_orig))
      for (j in 1:ncol(Phi_orig)) { e_vec <- rep(0, ncol(Phi_orig)); e_vec[j] <- 1; J[, j] <- compute_psi(e_vec) - psi_0 }
      if (any(is.na(J)) || abs(det(J)) < 1e-6) return(rep(NA, 4))
      as.numeric(solve(J + 1e-3 * penalty_diag, -psi_0))
    }
    
    theta_deal_s <- run_deal_clip(rho_star_final); theta_deal_T <- run_deal_clip(rho_tilde_final); theta_deal_H <- run_deal_clip(rho_hat_final)
    
    calc_metrics <- function(theta_est) {
      if(any(is.na(theta_est))) return(c(NA, NA, NA))
      err <- as.vector(y_t_orig - cbind(1, x_t_orig) %*% theta_est)
      c(mean(err^2), sqrt(mean(err^2)), mean(abs(err)))
    }
    metrics_deal_s <- calc_metrics(theta_deal_s); metrics_deal_T <- calc_metrics(theta_deal_T); metrics_deal_H <- calc_metrics(theta_deal_H)
    
    res_df <- data.frame(
      Seed = s,
      b0_DEAL_s_Clip001 = theta_deal_s[1], b1_DEAL_s_Clip001 = theta_deal_s[2], b2_DEAL_s_Clip001 = theta_deal_s[3], b3_DEAL_s_Clip001 = theta_deal_s[4],
      b0_DEAL_T_Clip001 = theta_deal_T[1], b1_DEAL_T_Clip001 = theta_deal_T[2], b2_DEAL_T_Clip001 = theta_deal_T[3], b3_DEAL_T_Clip001 = theta_deal_T[4],
      b0_DEAL_H_Clip001 = theta_deal_H[1], b1_DEAL_H_Clip001 = theta_deal_H[2], b2_DEAL_H_Clip001 = theta_deal_H[3], b3_DEAL_H_Clip001 = theta_deal_H[4],
      MSE_DEAL_s_Clip001 = metrics_deal_s[1], RMSE_DEAL_s_Clip001 = metrics_deal_s[2], MAE_DEAL_s_Clip001 = metrics_deal_s[3],
      MSE_DEAL_T_Clip001 = metrics_deal_T[1], RMSE_DEAL_T_Clip001 = metrics_deal_T[2], MAE_DEAL_T_Clip001 = metrics_deal_T[3],
      MSE_DEAL_H_Clip001 = metrics_deal_H[1], RMSE_DEAL_H_Clip001 = metrics_deal_H[2], MAE_DEAL_H_Clip001 = metrics_deal_H[3]
    )
    list(res = res_df)
  }, error = function(e) { NULL })
}
parallel::stopCluster(cl)

results_clip001 <- do.call(rbind, lapply(sim_clip001, function(x) x$res))
cat("✅ Clip0.01 완료!\n")


# =====================================================================
# 🟩 [3단계] GLM (Quasipoisson) 단독 실험
# =====================================================================
cat("\n▶ [3/3] GLM (Quasipoisson) 실험 진행 중...\n")
cl <- parallel::makeCluster(num_cores)
registerDoParallel(cl)

sim_glm <- foreach(s = 1:total_seeds, .packages=c("stats")) %dopar% {
  tryCatch({
    lambda_reg <- 0.001; max_rho_soft <- 10.0
    K_gauss <- function(d, bw) dnorm(d, mean = 0, sd = bw)
    set.seed(s)
    
    N <- 1000; r <- rbinom(N, 1, 0.5); n0 <- sum(r == 1); n1 <- N - n0; pi_val <- n0 / N
    y_s_orig <- rnorm(n0, 0, sqrt(2)); y_t_orig <- rnorm(n1, 1, 1); Y_orig <- c(y_s_orig, y_t_orig)
    x1_orig <- -0.5 * Y_orig + rnorm(N, 0, 1); x2_orig <- 0.5 * Y_orig + rnorm(N, 0, 1); x3_orig <- Y_orig + rnorm(N, 0, 1)
    X_orig <- cbind(x1_orig, x2_orig, x3_orig)
    x_s_orig <- X_orig[1:n0, ]; x_t_orig <- X_orig[(n0+1):N, ]; y_grid <- seq(-5, 5, length.out = 100) 
    
    phi_temp <- cbind(1, x_s_orig); beta_temp <- solve(t(phi_temp) %*% phi_temp + 1e-4 * diag(4)) %*% t(phi_temp) %*% y_s_orig; y_t_hat_orig <- cbind(1, x_t_orig) %*% beta_temp
    f_s <- approxfun(density(y_s_orig, kernel="gaussian")$x, density(y_s_orig, kernel="gaussian")$y, rule=2)
    f_t <- approxfun(density(y_t_hat_orig, kernel="gaussian")$x, density(y_t_hat_orig, kernel="gaussian")$y, rule=2)
    rho_star_grid <- pmax(f_t(y_grid) / (f_s(y_grid) + 1e-6), 0.01)
    rho_star_grid <- max_rho_soft * tanh(rho_star_grid / max_rho_soft); rho_star_grid <- rho_star_grid / mean(rho_star_grid)
    rho_star_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_star_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_star_final <- rho_star_final / mean(rho_star_final)
    
    h_ml <- 3 * n0^(-1/7); l_bw <- 1.5 * n0^(-1/3); h_bw <- 0.5 * n0^(-1/16)
    W_ML_unnorm <- matrix(0, N, n0)
    for (i in 1:N) for (j in 1:n0) W_ML_unnorm[i, j] <- K_gauss(X_orig[i,1]-x_s_orig[j,1], h_ml) * K_gauss(X_orig[i,2]-x_s_orig[j,2], h_ml) * K_gauss(X_orig[i,3]-x_s_orig[j,3], h_ml)
    W_ML <- W_ML_unnorm / (rowSums(W_ML_unnorm) + 1e-10)
    W_Y_unnorm <- outer(y_s_orig, y_s_orig, function(y1, y2) K_gauss(y1 - y2, l_bw)); W_Y <- W_Y_unnorm / (rowSums(W_Y_unnorm) + 1e-10)
    V_mat_grid <- outer(y_s_orig, y_grid, function(y_i, y_0) K_gauss(y_i - y_0, h_bw))
    
    solve_rho_grid <- function(rho_guess_grid) {
      rgs <- approx(y_grid, rho_guess_grid, xout=y_s_orig, rule=2)$y
      w <- 1 / as.vector(W_ML %*% (rgs^2 + (pi_val/(1-pi_val))*rgs))
      M <- W_Y %*% diag(w[1:n0]) %*% W_ML[1:n0, ] %*% diag(rgs)
      A <- solve(t(M)%*%M + 1e-3*diag(n0)) %*% t(M) %*% V_mat_grid
      E <- W_ML %*% diag(rgs) %*% A
      t1 <- colMeans(w[(n0+1):N] * E[(n0+1):N, , drop=FALSE]); t2 <- colMeans(rgs * (V_mat_grid - w[1:n0]*E[1:n0, , drop=FALSE]))
      (t1 + t2) / (colMeans(V_mat_grid) + 1e-8)
    }
    
    rho_tilde_grid <- pmax(0.01, pmin(solve_rho_grid(rho_star_grid), 20.0))
    rho_tilde_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_tilde_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_tilde_final <- rho_tilde_final / mean(rho_tilde_final)
    rho_hat_grid <- solve_rho_grid(rho_tilde_grid); rho_hat_grid[is.na(rho_hat_grid) | is.infinite(rho_hat_grid)] <- 1.0
    rho_hat_final <- pmax(0.01, max_rho_soft * tanh(approx(y_grid, rho_hat_grid, xout=y_s_orig, rule=2)$y / max_rho_soft)); rho_hat_final <- rho_hat_final / mean(rho_hat_final)
    
    # DEAL Engine (GLM)
    Phi_orig <- cbind(1, x_s_orig); Phi_all_orig <- cbind(1, X_orig); XtX_orig <- t(Phi_orig) %*% Phi_orig
    penalty_diag <- diag(ncol(Phi_orig)); penalty_diag[1,1] <- 0 
    H_orig <- Phi_orig %*% solve(XtX_orig + 1e-4 * penalty_diag) %*% t(Phi_orig)
    H_all_orig <- Phi_all_orig %*% solve(XtX_orig + 1e-4 * penalty_diag) %*% t(Phi_orig)
    K_h_deal <- exp(-(as.matrix(dist(y_s_orig))^2) / (2*(1.06*sd(y_s_orig)*(n0^(-1/5)))^2)); K_h_deal <- K_h_deal / rowSums(K_h_deal)
    
    run_deal_glm <- function(rho_val) {
      P <- diag(rho_val)
      Z_val <- as.vector((rho_val^2) + (pi_val / (1 - pi_val)) * rho_val)
      df_s <- data.frame(Z = Z_val, x_s_orig); df_all <- data.frame(X_orig); colnames(df_all) <- colnames(df_s)[-1]
      fit_d <- glm(Z ~ ., data = df_s, family = quasipoisson(link = "log"))
      d <- as.numeric(predict(fit_d, type = "response")); d_all <- as.numeric(predict(fit_d, newdata = df_all, type = "response"))
      
      # GLM 특성상 1e-10 최소 안전장치 유지 (로그 링크 때문)
      W <- diag(1 / pmax(d, 1e-10)); W_all <- diag(1 / pmax(d_all, 1e-10)) 
      
      L_h_deal <- K_h_deal %*% W %*% H_orig %*% P
      M_h_lambda <- solve(t(L_h_deal) %*% L_h_deal + lambda_reg * diag(n0)) %*% t(L_h_deal)
      
      compute_psi <- function(theta) {
        U <- -1 * (as.vector(y_s_orig - Phi_orig %*% theta) * Phi_orig)  
        R_h <- K_h_deal %*% (diag(n0) - W %*% H_orig %*% (P %*% P)) %*% U
        B_h_lambda <- W_all %*% H_all_orig %*% ((P %*% P) %*% U + P %*% (M_h_lambda %*% R_h))
        psi_hat <- matrix(0, nrow=N, ncol=ncol(Phi_orig))
        psi_hat[1:n0, ] <- (1/pi_val)*rho_val*(U - B_h_lambda[1:n0, ]); psi_hat[(n0+1):N, ] <- (1/(1-pi_val))*B_h_lambda[(n0+1):N, ]
        colMeans(psi_hat)
      }
      psi_0 <- compute_psi(rep(0, ncol(Phi_orig))); J <- matrix(0, ncol(Phi_orig), ncol(Phi_orig))
      for (j in 1:ncol(Phi_orig)) { e_vec <- rep(0, ncol(Phi_orig)); e_vec[j] <- 1; J[, j] <- compute_psi(e_vec) - psi_0 }
      if (any(is.na(J)) || abs(det(J)) < 1e-6) return(rep(NA, 4))
      as.numeric(solve(J + 1e-3 * penalty_diag, -psi_0))
    }
    
    theta_deal_s <- run_deal_glm(rho_star_final); theta_deal_T <- run_deal_glm(rho_tilde_final); theta_deal_H <- run_deal_glm(rho_hat_final)
    
    calc_metrics <- function(theta_est) {
      if(any(is.na(theta_est))) return(c(NA, NA, NA))
      err <- as.vector(y_t_orig - cbind(1, x_t_orig) %*% theta_est)
      c(mean(err^2), sqrt(mean(err^2)), mean(abs(err)))
    }
    metrics_deal_s <- calc_metrics(theta_deal_s); metrics_deal_T <- calc_metrics(theta_deal_T); metrics_deal_H <- calc_metrics(theta_deal_H)
    
    res_df <- data.frame(
      Seed = s,
      b0_DEAL_s_GLM = theta_deal_s[1], b1_DEAL_s_GLM = theta_deal_s[2], b2_DEAL_s_GLM = theta_deal_s[3], b3_DEAL_s_GLM = theta_deal_s[4],
      b0_DEAL_T_GLM = theta_deal_T[1], b1_DEAL_T_GLM = theta_deal_T[2], b2_DEAL_T_GLM = theta_deal_T[3], b3_DEAL_T_GLM = theta_deal_T[4],
      b0_DEAL_H_GLM = theta_deal_H[1], b1_DEAL_H_GLM = theta_deal_H[2], b2_DEAL_H_GLM = theta_deal_H[3], b3_DEAL_H_GLM = theta_deal_H[4],
      MSE_DEAL_s_GLM = metrics_deal_s[1], RMSE_DEAL_s_GLM = metrics_deal_s[2], MAE_DEAL_s_GLM = metrics_deal_s[3],
      MSE_DEAL_T_GLM = metrics_deal_T[1], RMSE_DEAL_T_GLM = metrics_deal_T[2], MAE_DEAL_T_GLM = metrics_deal_T[3],
      MSE_DEAL_H_GLM = metrics_deal_H[1], RMSE_DEAL_H_GLM = metrics_deal_H[2], MAE_DEAL_H_GLM = metrics_deal_H[3]
    )
    list(res = res_df)
  }, error = function(e) { NULL })
}
parallel::stopCluster(cl)

results_glm <- do.call(rbind, lapply(sim_glm, function(x) x$res))
cat("✅ GLM 완료!\n")



# =====================================================================
# 🎨 [시각화 전용] 일반 H행렬 (PureH) 박스플롯 다시 그리기
# =====================================================================
# 이미 실행 완료된 results_pureh 데이터를 사용합니다.
df_pureh <- as.data.frame(results_pureh)

# 모델 순서 및 색상 세팅 (가로형 플롯을 위해 역순 정렬)
models_orig <- c("Naive", "IPW_s", "IPW_T", "IPW_H", "DEAL_s_PureH", "DEAL_T_PureH", "DEAL_H_PureH", "Oracle")
box_names <- rev(c("Naive", "IPW(S)", "IPW(T)", "IPW(H)", "DEAL(S)", "DEAL(T)", "DEAL(H)", "Oracle"))
box_colors <- rev(c("salmon", "khaki1", "khaki3", "khaki4", "lightblue", "skyblue", "dodgerblue", "gray"))
models_plot <- rev(models_orig)

# 1. 파라미터 (b0 ~ b3) 가로형 박스플롯 (2x2 배열)
cat("\n🎨 [일반 H행렬] 파라미터 박스플롯을 그립니다...\n")
par(mfrow = c(2, 2), mar = c(4, 7, 3, 1))
for(i in 0:3) {
  cols <- paste0("b", i, "_", models_plot)
  boxplot(as.list(df_pureh[, cols, drop=FALSE]), names=box_names, col=box_colors, 
          horizontal=TRUE, main=paste("[PureH] Parameter b", i), 
          outline=FALSE, las=1, cex.axis=1.4, cex.main=1.4)
  
  # 평균값 하얀색 다이아몬드로 표시
  points(sapply(as.list(df_pureh[, cols, drop=FALSE]), mean, na.rm=TRUE), 
         1:8, pch=23, bg="white", cex=1.3)
}
# 2. 예측 지표 (RMSE, MAE) 가로형 박스플롯 (1x2 배열)
cat("\n🎨 [일반 H행렬] 예측 지표 박스플롯을 그립니다...\n")
par(mfrow = c(1, 2), mar = c(4, 7, 3, 1))
for(m in c("RMSE", "MAE")) {
  cols <- paste0(m, "_", models_plot)
  boxplot(as.list(df_pureh[, cols, drop=FALSE]), names=box_names, col=box_colors, 
          horizontal=TRUE, main=paste("[PureH] Prediction:", m), 
          outline=FALSE, las=1, cex.axis=1.1, cex.main=1.4)
  
  # 평균값 하얀색 다이아몬드로 표시
  points(sapply(as.list(df_pureh[, cols, drop=FALSE]), mean, na.rm=TRUE), 
         1:8, pch=23, bg="white", cex=1.3)
}
par(mfrow = c(1, 1)) # 레이아웃 초기화 완료

# =====================================================================
# 🎨 [시각화 전용] Clip 0.01 박스플롯 다시 그리기 (에러 완벽 수정본)
# =====================================================================
# 💡 DEAL이 안 붙은 베이스 모델(Naive, IPW, Oracle)의 파라미터와 '예측 지표(RMSE 등)' 전체를 안전하게 가져옵니다.
base_cols <- grep("DEAL", names(results_pureh), invert = TRUE, value = TRUE)
df_clip <- merge(results_pureh[, base_cols], results_clip001, by="Seed")

models_orig <- c("Naive", "IPW_s", "IPW_T", "IPW_H", "DEAL_s_Clip001", "DEAL_T_Clip001", "DEAL_H_Clip001", "Oracle")
box_names <- rev(c("Naive", "IPW(S)", "IPW(T)", "IPW(H)", "DEAL(S)", "DEAL(T)", "DEAL(H)", "Oracle"))
box_colors <- rev(c("salmon", "khaki1", "khaki3", "khaki4", "lightgreen", "seagreen", "forestgreen", "gray"))
models_plot <- rev(models_orig)

# 1. 파라미터 (b0 ~ b3) 가로형 박스플롯 (2x2 배열)
cat("\n🎨 [Clip 0.01] 파라미터 박스플롯을 그립니다...\n")
par(mfrow = c(2, 2), mar = c(4, 7, 3, 1))
for(i in 0:3) {
  cols <- paste0("b", i, "_", models_plot)
  boxplot(as.list(df_clip[, cols, drop=FALSE]), names=box_names, col=box_colors, 
          horizontal=TRUE, main=paste("[Clip 0.01] Parameter b", i), 
          outline=TRUE, las=1, cex.axis=1.4, cex.main=1.4)
  
  # 평균값 하얀색 다이아몬드로 표시
  points(sapply(as.list(df_clip[, cols, drop=FALSE]), mean, na.rm=TRUE), 
         1:8, pch=23, bg="white", cex=1.3)
}

# 2. 예측 지표 (RMSE, MAE) 가로형 박스플롯 (1x2 배열)
cat("\n🎨 [Clip 0.01] 예측 지표 박스플롯을 그립니다...\n")
par(mfrow = c(1, 2), mar = c(4, 7, 3, 1))
for(m in c("RMSE", "MAE")) {
  cols <- paste0(m, "_", models_plot)
  boxplot(as.list(df_clip[, cols, drop=FALSE]), names=box_names, col=box_colors, 
          horizontal=TRUE, main=paste("[Clip 0.01] Prediction:", m), 
          outline=TRUE, las=1, cex.axis=1.1, cex.main=1.4)
  
  # 평균값 하얀색 다이아몬드로 표시
  points(sapply(as.list(df_clip[, cols, drop=FALSE]), mean, na.rm=TRUE), 
         1:8, pch=23, bg="white", cex=1.3)
}
par(mfrow = c(1, 1)) # 레이아웃 초기화 완료

# =====================================================================
# 🎨 [시각화 전용] GLM (Quasipoisson) 박스플롯 다시 그리기 (에러 완벽 수정본)
# =====================================================================
# 💡 DEAL이 안 붙은 베이스 모델(Naive, IPW, Oracle)의 파라미터와 '예측 지표(RMSE 등)' 전체를 안전하게 가져옵니다.
base_cols <- grep("DEAL", names(results_pureh), invert = TRUE, value = TRUE)
df_glm <- merge(results_pureh[, base_cols], results_glm, by="Seed")

models_orig <- c("Naive", "IPW_s", "IPW_T", "IPW_H", "DEAL_s_GLM", "DEAL_T_GLM", "DEAL_H_GLM", "Oracle")
box_names <- rev(c("Naive", "IPW(S)", "IPW(T)", "IPW(H)", "DEAL(S)", "DEAL(T)", "DEAL(H)", "Oracle"))
# GLM은 보라색 톤으로 통일하여 구분
box_colors <- rev(c("salmon", "khaki1", "khaki3", "khaki4", "plum", "orchid", "purple", "gray"))
models_plot <- rev(models_orig)

# 1. 파라미터 (b0 ~ b3) 가로형 박스플롯 (2x2 배열)
cat("\n🎨 [GLM] 파라미터 박스플롯을 그립니다...\n")
par(mfrow = c(2, 2), mar = c(4, 7, 3, 1))
for(i in 0:3) {
  cols <- paste0("b", i, "_", models_plot)
  boxplot(as.list(df_glm[, cols, drop=FALSE]), names=box_names, col=box_colors, 
          horizontal=TRUE, main=paste("[GLM] Parameter b", i), 
          outline=TRUE, las=1, cex.axis=1.4, cex.main=1.4)
  
  # 평균값 하얀색 다이아몬드로 표시
  points(sapply(as.list(df_glm[, cols, drop=FALSE]), mean, na.rm=TRUE), 
         1:8, pch=23, bg="white", cex=1.3)
}

# 2. 예측 지표 (RMSE, MAE) 가로형 박스플롯 (1x2 배열)
cat("\n🎨 [GLM] 예측 지표 박스플롯을 그립니다...\n")
par(mfrow = c(1, 2), mar = c(4, 7, 3, 1))
for(m in c("RMSE", "MAE")) {
  cols <- paste0(m, "_", models_plot)
  boxplot(as.list(df_glm[, cols, drop=FALSE]), names=box_names, col=box_colors, 
          horizontal=TRUE, main=paste("[GLM] Prediction:", m), 
          outline=TRUE, las=1, cex.axis=1.1, cex.main=1.4)
  
  # 평균값 하얀색 다이아몬드로 표시
  points(sapply(as.list(df_glm[, cols, drop=FALSE]), mean, na.rm=TRUE), 
         1:8, pch=23, bg="white", cex=1.3)
}
par(mfrow = c(1, 1)) # 레이아웃 초기화 완료

# =====================================================================
# 🎨 [최종 종합 박스플롯] Naive, IPW(H), DEAL(H) 3종, Oracle 집중 비교
# =====================================================================

# 1. 에러 원천 차단: 필요한 데이터만 직접 쏙쏙 뽑아서 새 데이터프레임 조립
df_final <- data.frame(Seed = 1:nrow(results_pureh))

# 파라미터 (b0 ~ b3) 직접 가져오기
for(i in 0:3) {
  df_final[[paste0("b", i, "_Naive")]] <- results_pureh[[paste0("b", i, "_Naive")]]
  df_final[[paste0("b", i, "_IPW_H")]] <- results_pureh[[paste0("b", i, "_IPW_H")]]
  df_final[[paste0("b", i, "_DEAL_H_PureH")]] <- results_pureh[[paste0("b", i, "_DEAL_H_PureH")]]
  df_final[[paste0("b", i, "_DEAL_H_Clip001")]] <- results_clip001[[paste0("b", i, "_DEAL_H_Clip001")]]
  df_final[[paste0("b", i, "_DEAL_H_GLM")]] <- results_glm[[paste0("b", i, "_DEAL_H_GLM")]]
  df_final[[paste0("b", i, "_Oracle")]] <- results_pureh[[paste0("b", i, "_Oracle")]]
}

# 예측 지표 (RMSE, MAE) 직접 가져오기
for(m in c("RMSE", "MAE")) {
  df_final[[paste0(m, "_Naive")]] <- results_pureh[[paste0(m, "_Naive")]]
  df_final[[paste0(m, "_IPW_H")]] <- results_pureh[[paste0(m, "_IPW_H")]]
  df_final[[paste0(m, "_DEAL_H_PureH")]] <- results_pureh[[paste0(m, "_DEAL_H_PureH")]]
  df_final[[paste0(m, "_DEAL_H_Clip001")]] <- results_clip001[[paste0(m, "_DEAL_H_Clip001")]]
  df_final[[paste0(m, "_DEAL_H_GLM")]] <- results_glm[[paste0(m, "_DEAL_H_GLM")]]
  df_final[[paste0(m, "_Oracle")]] <- results_pureh[[paste0(m, "_Oracle")]]
}

# 2. 플롯 대상 모델 및 예쁜 이름/색상 정의
models_orig <- c("Naive", "IPW_H", "DEAL_H_PureH", "DEAL_H_Clip001", "DEAL_H_GLM", "Oracle")
box_names_orig <- c("Naive", "IPW(H)", "DEAL(H) Pure", "DEAL(H) Clip", "DEAL(H) GLM", "Oracle")

# 색상: Naive(빨강) -> IPW(노랑) -> DEAL(파랑/초록/보라) -> Oracle(회색)
box_colors_orig <- c("salmon", "khaki", "lightblue", "lightgreen", "plum", "gray")

# 가로형(horizontal) 출력을 위해 순서 뒤집기
models_plot <- rev(models_orig)
box_names <- rev(box_names_orig)
box_colors <- rev(box_colors_orig)

# 3. 파라미터 (b0 ~ b3) 가로형 박스플롯 (2x2 배열)
cat("\n🎨 [최종 비교] 파라미터 b0 ~ b3 박스플롯 출력 중...\n")
par(mfrow = c(2, 2), mar = c(4, 8, 3, 1)) # y축 레이블이 잘리지 않도록 좌측 여백(8) 확보
for(i in 0:3) {
  cols <- paste0("b", i, "_", models_plot)
  boxplot(as.list(df_final[, cols, drop=FALSE]), names=box_names, col=box_colors, 
          horizontal=TRUE, main=paste("Parameter b", i), 
          outline=FALSE, las=1, cex.axis=1.3, cex.main=1.5)
  
  # 평균값 하얀색 다이아몬드로 표시
  points(sapply(as.list(df_final[, cols, drop=FALSE]), mean, na.rm=TRUE), 
         1:6, pch=23, bg="white", cex=1.5)
}
par(mfrow = c(1, 1))

# 4. 예측 지표 (RMSE, MAE) 가로형 박스플롯 (1x2 배열)
cat("\n🎨 [최종 비교] 타겟 예측 성능(RMSE, MAE) 박스플롯 출력 중...\n")
par(mfrow = c(1, 2), mar = c(4, 8, 3, 1))
for(m in c("RMSE", "MAE")) {
  cols <- paste0(m, "_", models_plot)
  boxplot(as.list(df_final[, cols, drop=FALSE]), names=box_names, col=box_colors, 
          horizontal=TRUE, main=paste("Prediction:", m),ylim =c(0,2.2),
          outline=TRUE, las=1, cex.axis=1.2, cex.main=1.5)
  
  # 평균값 하얀색 다이아몬드로 표시
  points(sapply(as.list(df_final[, cols, drop=FALSE]), mean, na.rm=TRUE), 
         1:6, pch=23, bg="white", cex=1.5)
}
par(mfrow = c(1, 1)) # 레이아웃 초기화 완료

# =====================================================================
# 📊 [에러 완전 수정본] 원본 데이터에서 직접 추출하는 중앙값 요약 표
# =====================================================================

# 1. 행(지표)과 열(모델) 기준 정의
metrics <- c("b0", "b1", "b2", "b3", "MSE", "RMSE", "MAE")
model_names <- c("Naive", "IPW(H)", "DEAL(H)_Pure", "DEAL(H)_Clip", "DEAL(H)_GLM", "Oracle")

# 2. 빈 데이터프레임 생성
median_df <- data.frame(matrix(ncol = length(model_names), nrow = length(metrics)))
colnames(median_df) <- model_names
rownames(median_df) <- metrics

# 3. 원본 데이터에서 안전하게 값을 빼오는 헬퍼 함수
get_med <- function(df, col_name) {
  if (col_name %in% names(df)) {
    return(round(median(df[[col_name]], na.rm = TRUE), 4))
  } else {
    return(NA) # 만약 열이 없으면 에러 대신 NA 반환
  }
}

# 4. 각 지표별로 정확한 원본에서 다이렉트로 값 추출
for(m in metrics) {
  median_df[m, "Naive"]        <- get_med(results_pureh, paste0(m, "_Naive"))
  median_df[m, "IPW(H)"]       <- get_med(results_pureh, paste0(m, "_IPW_H"))
  median_df[m, "DEAL(H)_Pure"] <- get_med(results_pureh, paste0(m, "_DEAL_H_PureH"))
  median_df[m, "DEAL(H)_Clip"] <- get_med(results_clip001, paste0(m, "_DEAL_H_Clip001"))
  median_df[m, "DEAL(H)_GLM"]  <- get_med(results_glm, paste0(m, "_DEAL_H_GLM"))
  median_df[m, "Oracle"]       <- get_med(results_pureh, paste0(m, "_Oracle"))
}

# 5. 결과 출력
cat("\n=============================================================\n")
cat("📊 모델별 중앙값(Median) 종합 데이터프레임\n")
cat("=============================================================\n")
print(median_df)


