rm(list=ls())
library(MCMCpack)
#library(Rsolnp)
library(gsynth)
library(splines)
set.seed(12345) #for reproducibility of simulation
sim.n=100


GSCM.Mse<-rep(0,sim.n)
BayesBspline.Mse<-rep(0,sim.n)

GSCM.r<-rep(0,sim.n)
Bayes.r<-rep(0,sim.n)


BayesBspline.averagelength<-rep(0,sim.n)
GSCM.averagelength<-rep(0,sim.n)

for (M in 1:sim.n){
  set.seed(M)
  T<-100
  T0<-60
  T11<-T0+1
  T1<-T-T0 # Number of post-treatment periods
  K<-4 #number of total covariates
  p<-4 
  N<-40 #number of units
  Ntrt<-1
  Nco<-N-Ntrt
  
  Y<-matrix(0,N,T)
  
  alpha.i<-rnorm(N,0,sqrt(0.5))
  
  gamma0<-0
  gamma.sd<-sqrt(0.1)
  gamma.t<-numeric(T)
  gamma.t[1]<-gamma0
  for (t in 2:T){
    gamma.t[t]<-rnorm(1, gamma.t[t-1], gamma.sd)
  }
  c_vec<-rep(1,4)
  
  Z<-array(rnorm(N*p*T,1,sqrt(2)),c(N,p,T))
  
  r_true <- 2  # Hardcode this for the data generation
  f_t <- matrix(rnorm(r_true * T, 0, 1), nrow = r_true, ncol = T)
  lambda <- matrix(runif(N * r_true, -sqrt(3), sqrt(3)), N, r_true)
  
  
  for (i in 1:N) {
    for (t in 1:T) {
      # trt_val <- if (i==1 & t > T0) log(t/2 - 29) else 0
      trt_val <- if (i == 1 & t > T0) log(t - T0) else 0
      Y[i,t] <- alpha.i[i] + gamma.t[t] + sum(lambda[i, ] * f_t[, t])+sum(c_vec* Z[i,,t]) + trt_val + rnorm(1, 0, sqrt(0.1))
    }
  } 
  
  X<-array(0,c(T,K,N))
  for (i in 1:N){
    X[,,i]<-t(Z[i,,])
  }
  
  ########################################
  ########################################
  #Bayesian R selection ##################
  ########################################
  
  xi0<-1
  v0<-1000
  C0<-1
  a0<-1
  b0<-1
  R0<-2:9
  log.f.Y<-rep(0,length(R0))
  max.BF<--Inf
  
  ##Iteration scheme
  for(i1 in 1:length(R0)){
    
    r1=R0[i1]
    if(i1==1){
      Betahat0<-matrix(rnorm(K),K,1)
      Fhat0<-matrix(0,T0,r1)
      Fstar0<-matrix(0,T0-r1,r1)
      for (i in 1:(T0-r1)){
        for (j in 1:r1){
          Fstar0[i,j]<-rnorm(1,0,1)
        }
      }
      Fhat0<-rbind(diag(r1),Fstar0)
      
      Lambdahat0<-matrix(rnorm(N*r1),N,r1)
      
      Sigmahat0<-diag(0.1,T0) 
      C<-Lambdahat0%*%t(Fhat0)
    }else{
      Betahat0<-Betahat
      Lambdahat0<-C[,1:r1]
      Fstar<-t(solve(t(Lambdahat0)%*%Lambdahat0+0.1*diag(r1))%*%t(Lambdahat0)%*%C[,-(1:r1)])
      Fhat0<-rbind(diag(r1),Fstar) #
      Sigmahat0<-Sigmahat
    }
    
    
    esp<-0.00001
    U<-5000
    Fhat<-Fhat0
    Lambdahat<-Lambdahat0
    Sigmahat<-Sigmahat0
    
    for(j in 1:U){
      
      #inv.Sigmahat<-diag(1/diag(Sigmahat))
      Omega<-xi0*Fhat%*%t(Fhat)+Sigmahat
      inv.Omega<-solve(Omega)
      txx<-t(X[1:T0,,1])%*%inv.Omega%*%X[1:T0,,1]
      txy<-t(X[1:T0,,1])%*%inv.Omega%*%Y[1,1:T0]
      
      for (i in 2:N){
        txx<-txx+t(X[1:T0,,i])%*%inv.Omega%*%X[1:T0,,i]
        txy<-txy+t(X[1:T0,,i])%*%inv.Omega%*%Y[i,1:T0]
      }
      sd.beta<-solve(txx+diag(K)/v0)
      Betahat<-sd.beta%*%txy #MAP of Beta
      
      for (i in 1:N){
        Lambdahat[i,]<-solve(crossprod(Fhat)+Sigmahat[1,1]*diag(r1)/xi0)%*%t(Fhat)%*%(Y[i,1:T0]-X[1:T0,,i]%*%Betahat)
      }
      
      for (t in (r1+1):T0){
        sd.ft<-solve((Sigmahat[t,t])^(-1)*t(Lambdahat)%*%Lambdahat+diag(1/C0,r1))
        Fhat[t,]<-sd.ft%*%((Sigmahat[t,t])^(-1)*t(Lambdahat)%*%(Y[,t]-t(X[t,,])%*%Betahat))
      }
      
      J0<-0
      for (t in 1:T0){
        J<-Y[,t]-t(X[t,,])%*%Betahat-Lambdahat[,]%*%Fhat[t,]
        J0<-J0+sum(J^2)
      }
      
      diag(Sigmahat)<-((J0+b0)/2)/((a0+N*T0)/2+1)
      
      if(abs(Sigmahat[1,1]-Sigmahat0[1,1])<esp&max(abs(Betahat-Betahat0))<esp){break}
      if(j==U){print("No convergence")}
      Betahat0<-Betahat
      Fhat0<-Fhat
      Lambdahat0<-Lambdahat
      Sigmahat0<-Sigmahat
    }
    
    C<-Lambdahat%*%t(Fhat) #C updated
    
    #####log(y|r)##########
    Xb0<-matrix(0,N,T0)
    for(i in 1:N){
      Xb0[i,]<-X[1:T0,,i]%*%Betahat
    }
    Y.star<-Y[1:N,1:T0]-Xb0[1:N,]-C
    #matrix.Sigma<-matrix(diag(Sigmahat),N-1,T,byrow=TRUE) #corrected
    log.f.Y[i1]<-sum(dnorm(Y.star,0,sqrt(Sigmahat[1,1]),log=TRUE))-0.5*((T0-r1)*r1+N*r1+K+1)*log(N)
    #please consider the below for the revision #
    if(max.BF<log.f.Y[i1]){
      max.BF<-log.f.Y[i1]
      ini.beta.MCMC<-Betahat
      ini.Lambda.MCMC<-Lambdahat
      ini.F.MCMC<-Fhat
      ini.Sigma.MCMC<-Sigmahat}
    #  print(paste("Current r is ",r1))
    
  }
  
  r1= which.max(log.f.Y)+1 #corrected
  Bayes.r[M]<-r1
  
  
  X<-array(0,c(T,K,N))
  for (i in 1:N){
    X[,,i]<-t(Z[i,,])
  }
  ##############################
  txx <- matrix(0, K, K)
  txy <- numeric(K)
  for (i in 2:N) {
    txx <- txx + t(X[,,i]) %*% X[,,i]
    txy <- txy + t(X[,,i]) %*% Y[i,]
  }
  hat.beta <- solve(txx + diag(1e-6, K)) %*% txy # Added ridge for stability
  
  
  Fstar <- matrix(0, nrow = T, ncol = r1)
  for (t in 1:T) {
    kt <- min(t, r1)
    Fstar[t, 1:kt] <- rnorm(kt, 0, 0.1) # Initialize smaller to prevent early explosion
  }
  
  hat.Lambda <- matrix(rnorm(N * r1, 0, 0.1), N, r1)
  
  hat.Sigma <- diag(0.1, T) 
  Psi       <- diag(1, r1)
  
  # Hyper-parameters for IG priors 
  apsi <- 1; bpsi <- 1  
  
  asig <- 1; bsig <- 1 
  
  atheta <- 100; btheta <-0.1 
  
  
  # --- Set up for multiple models
  t_rel <- 1:T1
  
  # Design Matrice B
  X.spl  <- as.matrix(bs(t_rel, df = 4, intercept = TRUE, Boundary.knots = c(1, T1)))
  
  # Pre-calculate Inverses
  inv.tXX.spl  <- solve(crossprod(X.spl)  + diag(1e-6, ncol(X.spl)))
  
  # Current State Initializations
  tau_spl <- rep(0, T1); alpha_spl <- rep(0, ncol(X.spl)); thetasq_spl <- 0.5
  mu.spl  <- as.numeric(X.spl %*% alpha_spl)
  
  # Storage Matrices 
  MC<-10000
  MC.Tau.Spl  <- matrix(0, MC, T)
  
  #########################
  MC.beta  <- matrix(0, MC, K)
  MC.ATT   <- matrix(0, MC, T) 
 
  v0=1000
  
  # 1. Get initial residuals for unit 1
  initial_resid <- Y[1, (T0+1):T] - (X[(T0+1):T,,1] %*% hat.beta) - 
    (Fstar[(T0+1):T, ] %*% hat.Lambda[1, ])
  
  # 2. Project them onto the spline basis to get a 'warm' starting alpha
  alpha_spl <- solve(t(X.spl) %*% X.spl + diag(1e-4, ncol(X.spl))) %*% t(X.spl) %*% initial_resid
  mu.spl <- X.spl %*% alpha_spl
  tau_spl <- as.numeric(mu.spl)
  Tau_vec <- tau_spl
  
  for (k in 1:MC) {
    # 1. Update hat.beta 
    # Pre-calculate Omega inverse using the working factors Fstar
    Omega_star <- Fstar %*% t(Fstar) + hat.Sigma 
    inv_Omega <- solve(Omega_star + diag(1e-7, T)) # solve prevents "singular matrix" errors
    
    V_beta_inv <- diag(1/v0, K)
    M_beta_sum <- rep(0, K)
    
    for(i in 1:N) {
      y_adj <- Y[i,]
      if(i == 1) y_adj[(T0+1):T] <- y_adj[(T0+1):T] - Tau_vec
      
      Xi <- X[,,i]
      V_beta_inv <- V_beta_inv + t(Xi) %*% inv_Omega %*% Xi
      M_beta_sum <- M_beta_sum + t(Xi) %*% inv_Omega %*% y_adj
    }
    sd_beta <- solve(V_beta_inv)
    hat_beta <- as.numeric(mvrnorm(1, sd_beta %*% M_beta_sum, sd_beta))
    MC.beta[k, ] <- hat_beta
    
    # 2. Update hat.Lambda (Working Loadings)
    inv_hat_Sigma <- diag(1 / diag(hat.Sigma))
    inv_Psi <- solve(Psi)
    # Pre-calculate shared variance for all units
    V_lam <- solve(t(Fstar) %*% inv_hat_Sigma %*% Fstar + inv_Psi)
    
    for(i in 1:N) {
      y_resid <- Y[i,] - X[,,i] %*% hat_beta
      if(i == 1) y_resid[(T0+1):T] <- y_resid[(T0+1):T] - Tau_vec
      hat.Lambda[i,] <- mvrnorm(1, V_lam %*% t(Fstar) %*% inv_hat_Sigma %*% y_resid, V_lam)
    }
    
    # 3. Update Fstar (Working Factors - Lower Triangular)
    for(t in 1:T) {
      kt <- min(t, r1)
      Lam_kt <- hat.Lambda[, 1:kt, drop=FALSE] #prevent R from automatically converting 
      #a matrix into a vector when only one column/row is selected
      V_f <- solve( (1/hat.Sigma[t,t]) * t(Lam_kt) %*% Lam_kt + diag(kt) )
      
      y_resid_t <- Y[, t] - (t(X[t,,]) %*% hat_beta)
      if(t > T0) y_resid_t[1] <- y_resid_t[1] - Tau_vec[t-T0]
      
      m_f <- V_f %*% ((1/hat.Sigma[t,t]) * t(Lam_kt) %*% y_resid_t)
      Fstar[t, 1:kt] <- mvrnorm(1, m_f, V_f)
    }
    
    # 4. PX TRANSFORMATION: Update Psi and Map to Identifiable hat.F
    # Update expansion parameter Psi (Latent variance)
    for(l in 1:r1) {
      Psi[l, l] <- rinvgamma(1, (apsi + N)/2, (bpsi + sum(hat.Lambda[,l]^2))/2)
    }
    
    # Map F* to hat.F (The identifiable rotation)
    hat.F <- Fstar
    for(l in 1:r1) {
      s_l <- if(Fstar[l,l] >= 0) 1 else -1
      hat.F[,l] <- s_l * (1/sqrt(Psi[l,l])) * Fstar[,l]
    }
    
    # 5. Update hat.Sigma (Error Variance)
    for(t in 1:T) {
      y_err <- Y[,t] - (t(X[t,,]) %*% hat_beta) - (hat.Lambda %*% Fstar[t,])
      if(t > T0) y_err[1] <- y_err[1] - Tau_vec[t-T0]
      hat.Sigma[t,t] <- rinvgamma(1, (asig + N)/2, (bsig + sum(y_err^2))/2)
    }
    
####################################
    
    #B-SPLINE
    for(j in 1:T1) {
      t_idx <- T0 + j
      resid_1t <- Y[1,t_idx] - (X[t_idx,,1] %*% hat_beta) - (hat.Lambda[1,] %*% Fstar[t_idx,])
      #v_t <- 1 / (1/hat.Sigma[t_idx,t_idx] + 1/thetasq_spl)
      sig_t2 <- hat.Sigma[t_idx, t_idx]
      v_t <- sig_t2 * (1 / (1 + sig_t2 / thetasq_spl))
      tau_spl[j] <- rnorm(1, v_t*(resid_1t/hat.Sigma[t_idx,t_idx] + mu.spl[j]/thetasq_spl), sqrt(v_t))
      MC.Tau.Spl[k, t_idx] <- tau_spl[j]
    }
    alpha_spl   <- as.numeric(mvrnorm(1, inv.tXX.spl %*% t(X.spl) %*% tau_spl, thetasq_spl * inv.tXX.spl))
    mu.spl      <- X.spl %*% alpha_spl
    thetasq_spl <- rinvgamma(1, (atheta + T1)/2, (btheta + sum((tau_spl - mu.spl)^2))/2)
    
    
    # IMPORTANT: Sync back to Tau_vec for the next iteration's Factor update
    Tau_vec <- tau_spl
  }
  
  # Processing Results
  B <- floor(0.2 * MC)
  est_spl  <- colMeans(MC.Tau.Spl[-(1:B), T11:T])
  
  
  
  ########################
  ########  GSCM  ###########
  r_true1 <- 4  # Hardcode this for the data generation
  f_t <- matrix(rnorm(r_true1 * T, 0, 1), nrow = r_true, ncol = T)
  lambda <- matrix(runif(N * r_true, -sqrt(3), sqrt(3)), N, r_true)
  
  for (i in 1:N) {
    for (t in 1:T) {
      trt_val <- if (i == 1 & t > T0) log(t - T0) else 0
      Y[i,t] <- alpha.i[i] + gamma.t[t] + sum(lambda[i, ] * f_t[, t])+sum(c_vec* Z[i,,t]) + trt_val + rnorm(1, 0, sqrt(0.1))
    }
  } 
  ###### Define the treatment assignment variable####
  D<-c(rep(0,T0),rep(1,T1),rep(0,(N*T-T)))
  ###### Define the time variable######
  time <- rep(1:T, N)
  #### Define the region ######
  unit<-rep(1:N, each = T, times = 1)
  
  # Combine all variables into a data.frame
  
  data <- data.frame(Y = as.vector(t(Y)), unit=unit, time = time, D=D,
                     Z1 = as.vector(t(Z[,1,])), Z2 = as.vector(t(Z[,2,])),
                     Z3 = as.vector(t(Z[,3,])), Z4 = as.vector(t(Z[,4,])))
  
  
  out <- gsynth(Y ~ D + Z1 + Z2 + Z3 + Z4, data = data,
                index = c("unit","time"), force = "two-way",
                CV = TRUE, r = c(2,9), se = TRUE,
                inference = "parametric", nboots = 1000,
                parallel = FALSE)
  
  GSCM.r[M]<-out$r.cv 
  
  ###############################################
  # Calculate MSEs for this simulation run M
  t_post <- (T0 + 1):T
  TrueTrtEffect <- log(t_post - T0)
  
  BayesBspline.Mse[M]   <- mean((TrueTrtEffect - est_spl)^2)
  GSCM<-as.numeric((TrueTrtEffect-out$est.att[T11:T])^2)
  GSCM.Mse[M]<-mean(GSCM)
  
  # Posterior mean and pointwise 95% CI
  Post.Tau.Spl <- MC.Tau.Spl[-(1:B), T11:T]
  Bayes.95CI.spl <- apply(Post.Tau.Spl, 2, quantile, probs = c(0.025, 0.975))
  
  # ----------------------------
  # Pointwise 95% credible interval 
  # ----------------------------
  CI_pt_spl <- apply(Post.Tau.Spl, 2, quantile, probs = c(0.025, 0.975))
  
  # ----------------------------
  # Simultaneous 95% band (sup-t)
  # ----------------------------
  tau_hat <- colMeans(Post.Tau.Spl)        # posterior mean, length T1
  sd_tau  <- apply(Post.Tau.Spl, 2, sd)    # posterior sd, length T1
  sd_tau[sd_tau < 1e-12] <- 1e-12          # safety floor
  
  # standardized deviations for each draw across time
  Z <- sweep(Post.Tau.Spl, 2, tau_hat, "-")
  Z <- sweep(Z, 2, sd_tau, "/")
  
  # sup norm for each draw
  zmax <- apply(abs(Z), 1, max)
  
  # critical value for 95% simultaneous band
  c95 <- as.numeric(quantile(zmax, 0.95))
  
  # simultaneous band
  CI_sim_spl_lower <- tau_hat - c95 * sd_tau
  CI_sim_spl_upper <- tau_hat + c95 * sd_tau
  
  #######################
  Post.Tau.Spl  <- MC.Tau.Spl[-(1:B), T11:T]
  
  # Pointwise 95% CI (2 x T1)
  CI_pt_spl  <- apply(Post.Tau.Spl,  2, quantile, probs=c(0.025, 0.975))
  
  BayesBspline.averagelength[M]   <- mean(CI_pt_spl[2, ]  - CI_pt_spl[1, ])
  
  GSCM_band <- as.numeric(out$est.att[T11:T, 4] - out$est.att[T11:T, 3])
  GSCM.averagelength[M] <- mean(GSCM_band)
  
  print(paste("current iteration is",M))
}


mean(BayesBspline.Mse)
mean(GSCM.Mse)
Bayes.r
GSCM.r

mean(BayesBspline.averagelength)
mean(GSCM.averagelength)


sd(BayesBspline.Mse)/sqrt(sim.n) 
sd(GSCM.Mse)/sqrt(sim.n) 

plot(t_post, TrueTrtEffect, col=1, type="l", lwd=2, ylim=c(0,4),
     xlab = "Post-treatment time", ylab = "Treatment effect")

lines(t_post, est_spl,  col=2, lwd=2) # Bayes B-spline
lines(t_post,out$est.att[T11:T],col=5,lwd=2)
lines(t_post, CI_sim_spl_lower, col=2, lty=2, lwd=1)
lines(t_post, CI_sim_spl_upper, col=2, lty=2, lwd=1)

legend("topleft", 
       legend = c("True Effect", "Bayes B-spline", "GSCM","Spline 95% CI"),
       col = c(1, 2,5, 2),lty = c(1, 1,1, 2), lwd = c(2, 2,2, 1),cex = 0.5,                    
       bty = "n")                    


barplot(table(Bayes.r),ylim = c(0,100),ylab = "The counts r of Bayes GSCM",xlab="r")
barplot(table(GSCM.r),ylim = c(0,100),ylab = "The counts r of GSCM",xlab="r")

#####B-spline plot######
plot(t_post, TrueTrtEffect, col=1, type="l", lwd=2, ylim=c(0,4),
     xlab = "Post-treatment time", ylab = "Treatment effect")
lines(t_post, est_spl,  col=4, lwd=2) # Bayes B-spline (Red)
lines(t_post, CI_sim_spl_lower, col=4, lty=2, lwd=1)
lines(t_post, CI_sim_spl_upper, col=4, lty=2, lwd=1)
legend("topleft", 
       legend = c("True Effect",  "Bayes B-spline","Spline 95% CI"),
       col = c(1, 4, 4),lty = c(1, 1,  2), lwd = c(2, 2, 1),cex = 0.5,                    
       bty = "n") 

