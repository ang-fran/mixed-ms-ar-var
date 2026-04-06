rm(list=ls())
# Install then load my necessary packages
library(reshape2)
library(ggplot2)

set.seed(123)

# ---- Data Generation Process ----
simulate_msar_mixed_lag = function(T, K, p_vec, mu, sigma, phi, P) {
  p_max = max(p_vec)
  s = integer(T)
  y = numeric(T)
  
  # Initialize first observation
  s[1] = sample(1:K, 1)
  y[1] = rnorm(1, mu[s[1]], sigma[s[1]])
  
  if (T > 1) {
    for (t in 2:min(T, p_max)) {
      s[t] = sample(1:K, 1, prob = P[s[t-1],])
      y[t] = rnorm(1, mu[s[t]], sigma[s[t]])
    }
  }
  
  if (T > p_max) {
    for (t in (p_max+1):T) {
      s[t] = sample(1:K, 1, prob = P[s[t-1],])
      mean_t = mu[s[t]]
      for (i in 1:p_max) {
        if (i > p_vec[s[t]]) next
        mean_t = mean_t + phi[i, s[t]] * (y[t-i] - mu[s[t-i]])
      }
      y[t] = mean_t + rnorm(1, 0, sigma[s[t]])
    }
  }
  
  list(y=y, s=s)
}

# ---- Decode / Encode augmented state ----
decode_z = function(idx, K, p_max) {
  x = idx - 1L
  st = (x %% K) + 1L
  x = x %/% K
  base = K+1
  lags = integer(p_max)
  for (i in 1:p_max) {
    lags[i] = x %% base
    x = x %/% base
  }
  c(st, lags)
}

encode_z = function(z, K, p_max) {
  st = z[1]; lags = z[-1]
  x = st - 1L
  base = K+1
  mult = K; pow = 1L
  for (i in 1:p_max) {
    x = x + mult*(lags[i]*pow)
    pow = pow*base
  }
  x+1L
}

# ---- Compute emission probabilities ----
compute_g_z_allt = function(y, mu, sigma, phi, p_vec, eps=1e-300) {
  Tn = length(y); K = length(mu); p_max = nrow(phi)
  K_aug = K*(K+1)^p_max
  g = matrix(0, nrow=K_aug, ncol=Tn)
  
  for (t in 1:Tn) {
    for (z in 1:K_aug) {
      zt = decode_z(z, K, p_max)
      st = zt[1]; mean_t = mu[st]
      for (i in 1:p_max) {
        if (i > p_vec[st] || i>t-1) next
        s_lag = zt[i+1]; if (s_lag==0) next
        mean_t = mean_t + phi[i,st]*(y[t-i]-mu[s_lag])
      }
      g[z,t] = dnorm(y[t], mean=mean_t, sd=sigma[st]) + eps
    }
  }
  list(g=g, K=K_aug)
}

# ---- E-step ----
msar_e_step = function(y, mu, sigma, phi, p_vec, P_hat, pi=NULL, eps=1e-300) {
  Tn = length(y); K = length(mu); p_max = nrow(phi)
  K_aug = K*(K+1)^p_max
  if (is.null(pi)) pi = rep(1/K, K)
  
  Eg = compute_g_z_allt(y, mu, sigma, phi, p_vec, eps)
  g = Eg$g
  alpha = matrix(0,Tn,K_aug); beta = matrix(0,Tn,K_aug); csc = numeric(Tn)
  
  pi_z1 = rep(0,K_aug)
  for (j in 1:K) pi_z1[encode_z(c(j,rep(0L,p_max)),K,p_max)] = pi[j]
  
  # Forward
  a1 = pi_z1*g[,1]; c1 = sum(a1)+eps; alpha[1,] = a1/c1; csc[1] = c1
  for (t in 2:Tn) {
    anew = rep(0,K_aug)
    for (zprev in 1:K_aug) {
      ap = alpha[t-1,zprev]; if (ap==0) next
      zprev_vec = decode_z(zprev,K,p_max); i = zprev_vec[1]
      for (j in 1:K) {
        znew_vec = c(j,zprev_vec[1:p_max]); znew = encode_z(znew_vec,K,p_max)
        anew[znew] = anew[znew] + ap*P_hat[i,j]
      }
    }
    anew = anew*g[,t]; ct = sum(anew)+eps; alpha[t,] = anew/ct; csc[t] = ct
  }
  
  # Backward
  beta[Tn,] = 1
  for (t in (Tn-1):1) {
    bnew = rep(0,K_aug)
    for (zprev in 1:K_aug) {
      zprev_vec = decode_z(zprev,K,p_max); i = zprev_vec[1]; acc = 0
      for (j in 1:K) {
        znew_vec = c(j,zprev_vec[1:p_max]); znew = encode_z(znew_vec,K,p_max)
        acc = acc + P_hat[i,j]*g[znew,t+1]*beta[t+1,znew]
      }
      bnew[zprev] = acc
    }
    beta[t,] = bnew/(csc[t+1]+eps)
  }
  
  gammaZ = alpha*beta; gammaZ = gammaZ/rowSums(gammaZ)
  
  smoothed_marg = matrix(0,Tn,K)
  for (t in 1:Tn) for (z in 1:K_aug) smoothed_marg[t,decode_z(z,K,p_max)[1]] = smoothed_marg[t,decode_z(z,K,p_max)[1]]+gammaZ[t,z]
  
  smoothed_joint = vector("list",Tn); smoothed_joint[[1]] = NULL
  for (t in 2:Tn) {
    J = matrix(0,K,K)
    for (zprev in 1:K_aug) { alpha_val = alpha[t-1,zprev]; if(alpha_val==0) next
    zprev_vec = decode_z(zprev,K,p_max); i = zprev_vec[1]
    for (j in 1:K) { znew_vec = c(j,zprev_vec[1:p_max]); znew = encode_z(znew_vec,K,p_max)
    J[i,j] = J[i,j] + alpha_val*P_hat[i,j]*g[znew,t]*beta[t,znew]
    }
    }
    smoothed_joint[[t]] = J/sum(J)
  }
  
  list(alpha=alpha, beta=beta, gammaZ=gammaZ, smoothed_marg=smoothed_marg,
       smoothed_joint=smoothed_joint, loglik=sum(log(csc)))
}

# ---- M-step ----
msar_m_step = function(y, gammaZ, smoothed_joint, mu_prev, phi_prev, sigma_prev, p_vec, eps=1e-12) {
  Tn = length(y); K = length(mu_prev); p_max = nrow(phi_prev); K_aug = ncol(gammaZ)
  mu_new = numeric(K); phi_new = matrix(0,p_max,K); sigma_new = numeric(K); P_new = matrix(0,K,K)
  
  # phi update
  for (j in 1:K) {
    XtWX = matrix(0,p_max,p_max); XtWy = numeric(p_max)
    for (t in 1:Tn) for (z in 1:K_aug) {
      w = gammaZ[t,z]; if(is.na(w) || w <= 0) next
      zt = decode_z(z,K,p_max); if(zt[1]!=j) next
      x = numeric(p_max); ytil = y[t]-mu_prev[j]
      for (i in 1:p_max) { if(i>p_vec[j]||i>t-1) next
        s_lag = zt[i+1]; if(s_lag==0) next
        x[i] = y[t-i]-mu_prev[s_lag]
      }
      XtWX = XtWX + w*(x%*%t(x)); XtWy = XtWy + w*(x*ytil)
    }
    phi_new[,j] = MASS::ginv(XtWX) %*% XtWy  }
  
  # mu update
  for (j in 1:K) {
    num = 0; den = 0
    for (t in 1:Tn) for (z in 1:K_aug) {
      w = gammaZ[t,z]; if(w<=0) next
      zt = decode_z(z,K,p_max); if(zt[1]!=j) next
      pred_no_mu = 0
      for (i in 1:p_max) { if(i>p_vec[j]||i>t-1) next; s_lag=zt[i+1]; if(s_lag==0) next
      pred_no_mu = pred_no_mu + phi_new[i,j]*(y[t-i]-mu_prev[s_lag])
      }
      num = num + w*(y[t]-pred_no_mu); den = den + w
    }
    mu_new[j] = num/(den+eps)
  }
  
  # sigma update
  for (j in 1:K) {
    num = 0; den = 0
    for (t in 1:Tn) for (z in 1:K_aug) {
      w = gammaZ[t,z]; if(w<=0) next
      zt = decode_z(z,K,p_max); if(zt[1]!=j) next
      pred = mu_new[j]
      for (i in 1:p_max) { if(i>p_vec[j]||i>t-1) next; s_lag=zt[i+1]; if(s_lag==0) next
      pred = pred + phi_new[i,j]*(y[t-i]-mu_new[s_lag])
      }
      resid = y[t]-pred; num = num + w*resid^2; den = den + w
    }
    sigma_new[j] = sqrt(num/(den+eps))
  }
  
  # P update
  for (i in 1:K) for (j in 1:K) P_new[i,j] = sum(sapply(2:Tn,function(t) smoothed_joint[[t]][i,j]))
  for (i in 1:K) P_new[i,] = P_new[i,]/sum(P_new[i,])
  
  list(mu=mu_new, phi=phi_new, sigma=sigma_new, P=P_new)
}

# ---- Iteration of EM algorithm till convergence ----
msar_em_mixed_lag_runner = function(y, mu_init, phi_init, sigma_init, P_init, p_vec,
                                    max_iter=200, tol=1e-6, verbose=TRUE) {
  mu_curr = mu_init; phi_curr = phi_init; sigma_curr = sigma_init; P_curr = P_init
  loglik_history = numeric(max_iter)
  
  for (iter in 1:max_iter) {
    E = msar_e_step(y, mu_curr, sigma_curr, phi_curr, p_vec, P_curr)
    Mres = msar_m_step(y, E$gammaZ, E$smoothed_joint, mu_curr, phi_curr, sigma_curr, p_vec)
    
    mu_curr = Mres$mu; phi_curr = Mres$phi; sigma_curr = Mres$sigma; P_curr = Mres$P
    loglik_history[iter] = E$loglik
    if (iter>1 && abs(loglik_history[iter]-loglik_history[iter-1])<tol) { loglik_history = loglik_history[1:iter]; break }
  }
  
  list(mu=mu_curr, phi=phi_curr, sigma=sigma_curr, P=P_curr, loglik_history=loglik_history)
}

align_regimes = function(mu_true, sigma_true, phi_true, P_true,
                         mu_hat, sigma_hat, phi_hat, P_hat, p_max=NULL) {
  K = length(mu_true)
  
  if (is.null(p_max)) p_max = nrow(phi_hat)
  
  # Pad phi_true and phi_hat to p_max rows if needed
  phi_true_pad = matrix(0, nrow=p_max, ncol=K)
  phi_true_pad[1:nrow(phi_true),] = phi_true
  phi_hat_pad = matrix(0, nrow=p_max, ncol=K)
  phi_hat_pad[1:nrow(phi_hat),] = phi_hat
  
  if(K==2){
    cost1 = sum((mu_hat-mu_true)^2) + sum((sigma_hat-sigma_true)^2) + sum((phi_hat_pad-phi_true_pad)^2)
    cost2 = sum((mu_hat[c(2,1)]-mu_true)^2) + sum((sigma_hat[c(2,1)]-sigma_true)^2) + sum((phi_hat_pad[,c(2,1)]-phi_true_pad)^2)
    if(cost2<cost1){
      mu_hat = mu_hat[c(2,1)]; sigma_hat = sigma_hat[c(2,1)]
      phi_hat_pad = phi_hat_pad[,c(2,1),drop=FALSE]; P_hat = P_hat[c(2,1),c(2,1)]
    }
  }
  list(mu_hat=mu_hat, sigma_hat=sigma_hat, phi_hat=phi_hat_pad, P_hat=P_hat)
}

# ---- Run simulation and EM ----
T = 700; K = 2; p_vec = c(2,4)
mu_true = c(5,-5); sigma_true = c(1,2)
phi_true = matrix(0,nrow=max(p_vec),ncol=K); phi_true[1:2,1]=c(0.5,0.3); phi_true[,2]=c(0.4,0.2,0.1, 0.2)
P_true = matrix(c(0.8,0.2,0.3,0.7),2,2,byrow=TRUE)

sim = simulate_msar_mixed_lag(T,K,p_vec,mu_true,sigma_true,phi_true,P_true)
y = sim$y; s = sim$s

# Initial guesses
mu_init = c(7, -2)
phi_init = matrix(0,nrow=max(p_vec),ncol=K)
phi_init[1:2,1] = c(0.8,0.2)
phi_init[,2] = c(0.5,0.4,0.3, 0.1)
sigma_init = c(2,5)
P_init = matrix(c(0.75,0.25,0.4,0.6),2,2,byrow=TRUE)

fit = msar_em_mixed_lag_runner(y, mu_init, phi_init, sigma_init, P_init, p_vec)
aligned = align_regimes(mu_true, sigma_true, phi_true, P_true,
                        fit$mu, fit$sigma, fit$phi, fit$P, p_max = max(p_vec))

# ---- Comparison chart ----
compare_table = function(mu_true, sigma_true, phi_true, P_true,
                         mu_init, sigma_init, phi_init,
                         mu_hat, sigma_hat, phi_hat, P_hat) {
  cat("Mu comparison (true | initial | estimated):\n")
  print(round(cbind(mu_true, mu_init, mu_hat),3))
  cat("\nSigma comparison (true | initial | estimated):\n")
  print(round(cbind(sigma_true, sigma_init, sigma_hat),3))
  cat("\nPhi comparison (rows=regimes, cols=lags)\n")
  cat("True phi:\n"); print(round(phi_true,3))
  cat("Init phi:\n"); print(round(phi_init,3))
  cat("Estimated phi:\n"); print(round(phi_hat,3))
  cat("\nP comparison (true | estimated):\n")
  print(round(cbind(P_true, P_hat),3))
}

compare_table(mu_true, sigma_true, phi_true, P_true,
              mu_init, sigma_init, phi_init,
              aligned$mu_hat, aligned$sigma_hat, aligned$phi_hat, aligned$P_hat)




# Repeat simulation but for multiple values of T and different initial values

# ---- Settings ----
T_vec = c(200, 500, 700)
init_list = list(
  list(mu=c(7,-2), phi=matrix(c(0.8,0.2,0,0,0.5,0.4,0.3,0.1), nrow=4), sigma=c(2,5), P=matrix(c(0.75,0.25,0.4,0.6),2,2)),
  list(mu=c(1,-1), phi=matrix(c(0.6,0.1,0,0,0.3,0.2,0.2,0.1), nrow=4), sigma=c(1.5,3), P=matrix(c(0.6,0.4,0.2,0.8),2,2)),
  list(mu=c(15,-10), phi=matrix(c(0.7,0.4,0,0,0.4,0.3,0.1,0.2), nrow=4), sigma=c(1,4), P=matrix(c(0.4,0.6,0.1,0.9),2,2))
)
K = 2; p_vec = c(2,4)
mu_true = c(5,-5); sigma_true = c(1,2)
phi_true = matrix(0,nrow=max(p_vec),ncol=K)
phi_true[1:2,1]=c(0.5,0.3); phi_true[,2]=c(0.4,0.2,0.1,0.2)
P_true = matrix(c(0.8,0.2,0.3,0.7),2,2,byrow=TRUE)

# Loop over T for convergence to true values
results_T = list()
for(T in T_vec){
  sim = simulate_msar_mixed_lag(T,K,p_vec,mu_true,sigma_true,phi_true,P_true)
  y = sim$y
  
  fit = msar_em_mixed_lag_runner(y, init_list[[1]]$mu, init_list[[1]]$phi,
                                 init_list[[1]]$sigma, init_list[[1]]$P, p_vec)
  aligned = align_regimes(mu_true, sigma_true, phi_true, P_true,
                          fit$mu, fit$sigma, fit$phi, fit$P, p_max = max(p_vec))
  results_T[[paste0("T_",T)]] = list(fit=aligned, loglik=fit$loglik_history)
}

# Loop over initializations for local maxima
sim_full = simulate_msar_mixed_lag(500,K,p_vec,mu_true,sigma_true,phi_true,P_true)
y = sim_full$y

results_init = list()
for(i in 1:length(init_list)){
  fit = msar_em_mixed_lag_runner(y, init_list[[i]]$mu, init_list[[i]]$phi,
                                 init_list[[i]]$sigma, init_list[[i]]$P, p_vec)
  aligned = align_regimes(mu_true, sigma_true, phi_true, P_true,
                          fit$mu, fit$sigma, fit$phi, fit$P, p_max = max(p_vec))
  results_init[[paste0("init_",i)]] = list(fit=aligned, loglik=fit$loglik_history)
}

# ---- Plotting ----

plot_params = function(results_list, param="mu"){
  df_list = list()
  for(name in names(results_list)){
    res = results_list[[name]]$fit
    val = switch(param,
                 mu=res$mu_hat,
                 sigma=res$sigma_hat,
                 phi=res$phi_hat)
    if(param=="phi"){
      val = melt(val)
      colnames(val) = c("Lag","Regime","Value")
      val$Setting = name
      df_list[[name]] = val
    } else {
      df_list[[name]] = data.frame(Regime=1:K, Value=val, Setting=name)
    }
  }
  df = do.call(rbind, df_list)
  df$Setting = factor(df$Setting, levels=names(results_list))
  df
}

plot_line = function(df,param="mu"){
  if(param=="phi"){
    ggplot(df, aes(x=factor(Lag), y=Value, color=Setting, group=Setting)) +
      geom_line() + geom_point() +
      facet_wrap(~Regime, labeller=labeller(Regime=function(x) paste("Regime",x))) +
      labs(title="Phi convergence", x="Lag", y="Value")
  } else {
    ggplot(df, aes(x=Setting, y=Value, color=factor(Regime), group=Regime)) +
      geom_line() + geom_point() +
      labs(title=paste(param,"convergence"), x="Setting", y="Value", color="Regime")
  }
}


# Convergence to true values as T increases
df_mu_T = plot_params(results_T,"mu")
df_sigma_T = plot_params(results_T,"sigma")
df_phi_T = plot_params(results_T,"phi")

plot_line(df_mu_T,"mu")
plot_line(df_sigma_T,"sigma")
plot_line(df_phi_T,"phi")

# Convergence to local maxima with different initial values
df_mu_init = plot_params(results_init,"mu")
df_sigma_init = plot_params(results_init,"sigma")
df_phi_init = plot_params(results_init,"phi")

plot_line(df_mu_init,"mu")
plot_line(df_sigma_init,"sigma")
plot_line(df_phi_init,"phi")

# Plot log-likelihood history for one example
plot_loglik = function(results_list){
  df_list = list()
  for(name in names(results_list)){
    df_list[[name]] = data.frame(Iter=1:length(results_list[[name]]$loglik),
                                 LogLik=results_list[[name]]$loglik,
                                 Setting=name)
  }
  df = do.call(rbind, df_list)
  ggplot(df, aes(x=Iter, y=LogLik, color=Setting, group=Setting)) +
    geom_line() +
    labs(title="EM log-likelihood evolution", x="Iteration", y="Log-Likelihood")
}

plot_loglik(results_T)     # For T study
plot_loglik(results_init)  # For different initial values

# Mu plot
plot_params = function(param_list, true_values, param_name, filename) {
  df = do.call(rbind, lapply(names(param_list), function(nm) {
    est = param_list[[nm]]$fit[[paste0(param_name, "_hat")]]
    data.frame(Regime = 1:length(est), Value = est, Label = nm)
  }))
  df_true = data.frame(Regime = 1:length(true_values), Value = true_values, Label = "True")
  df_all = rbind(df, df_true)
  
  p = ggplot(df_all, aes(x=factor(Regime), y=Value, fill=Label)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_brewer(palette="Set2") +
    xlab("Regime") + ylab(param_name) +
    ggtitle(paste("Estimated", param_name, "vs True")) +
    theme_minimal()
  ggsave(filename, p, width=7, height=4)
}

# Phi  plot
plot_phi = function(param_list, true_phi, p_max, K, filename) {
  plot_df = do.call(rbind, lapply(names(param_list), function(nm) {
    phi_hat = param_list[[nm]]$fit$phi_hat
    df = melt(phi_hat)
    df$Iteration = nm
    colnames(df) = c("Lag", "Regime", "Value", "Iteration")
    df
  }))
  true_df = melt(true_phi)
  true_df$Iteration = "True"
  colnames(true_df) = c("Lag", "Regime", "Value", "Iteration")
  plot_df = rbind(plot_df, true_df)
  
  p = ggplot(plot_df, aes(x=factor(Lag), y=Value, fill=Iteration)) +
    geom_bar(stat="identity", position=position_dodge()) +
    facet_wrap(~Regime) +
    scale_fill_brewer(palette="Set2") +
    xlab("Lag") + ylab("phi") +
    ggtitle("Estimated phi coefficients vs True") +
    theme_minimal()
  ggsave(filename, p, width=7, height=5)
}

# Log-likelihood
plot_loglik = function(param_list, filename) {
  log_df = do.call(rbind, lapply(names(param_list), function(nm) {
    loglik = param_list[[nm]]$loglik
    data.frame(Iteration = 1:length(loglik), LogLik = loglik, Label = nm)
  }))
  
  p = ggplot(log_df, aes(x=Iteration, y=LogLik, color=Label)) +
    geom_line(linewidth =1) +
    xlab("EM Iteration") + ylab("Log-Likelihood") +
    ggtitle("Log-Likelihood Trajectories") +
    theme_minimal()
  
  ggsave(filename, p, width=7, height=4)
}

# ---- T-Varying plots ----
plot_params(results_T, mu_true, "mu", "mu_T_plot.png")
plot_params(results_T, sigma_true, "sigma", "sigma_T_plot.png")
plot_phi(results_T, phi_true, max(p_vec), K, "phi_T_plot.png")
plot_loglik(results_T, "loglik_T_plot.png")

# ---- Initialization plots ----
plot_params(results_init, mu_true, "mu", "mu_init_plot.png")
plot_params(results_init, sigma_true, "sigma", "sigma_init_plot.png")
plot_phi(results_init, phi_true, max(p_vec), K, "phi_init_plot.png")
plot_loglik(results_init, "loglik_init_plot.png")
















# ============================================================
# 9) Simple plot of series colored by true regime
# ============================================================
plot(y, type="l", col="grey50", lwd=1.5, ylab="y_t", xlab="t")
cols = ifelse(s==1,"royalblue4","violetred4")
points(y, col=cols, pch=16, cex=0.7)
legend("topleft", legend=c("Regime 1","Regime 2"), col=c("royalblue4","violetred4"), pch=16)
fit$loglik_history

par(mfrow=c(1,1), mar=c(4,4,2,1))

# 1) Mu
barplot(rbind(mu_true, aligned$mu_hat),
        beside=TRUE, col=c("skyblue","salmon"),
        main="Regime means (mu)", names.arg=paste0("Reg ",1:K))
legend("topright", legend=c("True","Estimated"), fill=c("skyblue","salmon"))

# 2) Sigma
barplot(rbind(sigma_true, aligned$sigma_hat),
        beside=TRUE, col=c("skyblue","salmon"),
        main="Regime std dev (sigma)", names.arg=paste0("Reg ",1:K))
legend("topleft", legend=c("True","Estimated"), fill=c("skyblue","salmon"))

# 3) Phi
library(reshape2); library(ggplot2)
phi_df = melt(cbind(phi_true, aligned$phi_hat))
colnames(phi_df) = c("Lag","Regime","Value")
phi_df$Type = rep(c("True","Estimated"), each=max(p_vec)*K)
ggplot(phi_df, aes(x=factor(Lag), y=Value, fill=Type)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~Regime, labeller=labeller(Regime=function(x) paste("Regime", x))) +
  labs(title="Autoregressive coefficients (phi)", x="Lag", y="Value") +
  scale_fill_manual(values=c("skyblue","salmon"))
