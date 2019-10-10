myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}


# True purity 
KL_true = function(alpha){

  alpha = matrix(alpha,p,1)
  alpha = alpha/c(sqrt(t(alpha) %*% alpha))
  alpha = round(alpha, 2)
  
  beta1 = beta_drg
  beta2 = beta_pbo
  
  gamma1 = gamma_drg
  gamma2 = gamma_pbo
  
  D1 = bi_sigma
  D2 = bi_sigma2

  A_0 = -2 + 0.5 *  sum(diag(ginv(D1[2:3, 2:3]) %*% D2[2:3, 2:3])) + 
    0.5 * sum(diag(ginv(D2[2:3, 2:3]) %*% D1[2:3, 2:3]))
  A_1 = t(beta1[2:3] - beta2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                          solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
  A_2 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                            solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
  A_3 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                            solve(D2[2:3, 2:3])) %*% (gamma1[2:3] - gamma2[2:3]) 
  
  est_alpha = alpha
  sigma_x = diag(1,p)
  
  return_res = A_0 + A_1/2 +
    c(A_3/2) * t(est_alpha) %*% sigma_x %*% est_alpha 

  res = return_res
  
  return(res)
}

KL_true(rep(1,10)) # 3.843557

# KL with two estimation of D, with Lagrange multiplier
KL_lm = function(alpha){
  
  dat_try = dat
  
  for(wcalculation in 1){
    
    # scale alpha
    alpha = matrix(alpha,p,1)
    alpha = alpha/c(sqrt(t(alpha) %*% alpha))
    alpha = round(alpha, 2)
    
    # calculate W = alpha^T x
    covar_list = dat_try[,paste('X', 1:p, sep='')]
    covar_list = scale(covar_list)
    
    temp_W = unlist(c(covar_list))
    temp_W = matrix(temp_W, dim(dat_try)[1], p)
    dat_try$W = c(temp_W %*% alpha)
    
  }
  
  dat_pbo_est = dat_try[dat_try$trt == 1, ]
  dat_drg_est = dat_try[dat_try$trt == 2, ]
  
  fit_pbo_est = myTryCatch(lmer(y ~ tt + I(tt^2) + W + W * tt +
                                  W * I(tt^2) + (tt+I(tt^2)|subj),
                                data = dat_pbo_est, REML = FALSE))
  fit_drg_est = myTryCatch(lmer(y ~ tt + I(tt^2) + W + W * tt +
                                  W * I(tt^2) + (tt+I(tt^2)|subj),
                                data = dat_drg_est, REML = FALSE))
  
  if(is.null(fit_pbo_est$warning) + is.null(fit_drg_est$warning) == 2){
    # beta, gamma, and D matrix
    beta1 = as.matrix(fixef(fit_drg_est$value))[1:3]
    gamma1 = as.matrix(fixef(fit_drg_est$value))[4:6] 
    D1 = as.matrix(VarCorr(fit_drg_est$value)$subj)[1:3, 1:3] 
    beta2 = as.matrix(fixef(fit_pbo_est$value))[1:3] 
    gamma2 = as.matrix(fixef(fit_pbo_est$value))[4:6] 
    D2 = as.matrix(VarCorr(fit_pbo_est$value)$subj)[1:3, 1:3] 
    
    temp = covar_list 
    xx = as.numeric(unlist(c(temp)))
    xx = matrix(xx, p, dim(dat_try)[1], byrow = TRUE)
    
    # this is the function derivated from the KL-divergence
    A_0 = -2 + 0.5 *  sum(diag(ginv(D1[2:3, 2:3]) %*% D2[2:3, 2:3])) + 
      0.5 * sum(diag(ginv(D2[2:3, 2:3]) %*% D1[2:3, 2:3]))
    A_1 = t(beta1[2:3] - beta2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                            solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_2 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_3 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (gamma1[2:3] - gamma2[2:3]) 
    
    mu_x = matrix(apply(t(xx), 2, mean), p, 1)
    sigma_x = cov(t(xx))
    
    lambda0 = max(eigen(sigma_x)$values)
    n_max = which(eigen(sigma_x)$values == lambda0)
    est_alpha = eigen(sigma_x)$vectors[,n_max]
    est_alpha = matrix(est_alpha, p, 1)
    
    return_res = A_0 + A_1/2 + 
      c(A_3/2) * t(est_alpha) %*% (sigma_x) %*% est_alpha
    
    return_res =  c(-return_res)

  }else{
    print('not convergence')
    return_res = 0
  }
  return(return_res)
}
KL_lm_D = function(alpha){
  
  dat_try = dat
  
  for(wcalculation in 1){
    
    # scale alpha
    alpha = matrix(alpha,p,1)
    alpha = alpha/c(sqrt(t(alpha) %*% alpha))
    alpha = round(alpha, 2)
    
    # calculate W = alpha^T x
    covar_list = dat_try[,paste('X', 1:p, sep='')]
    covar_list = scale(covar_list)
    
    temp_W = unlist(c(covar_list))
    temp_W = matrix(temp_W, dim(dat_try)[1], p)
    dat_try$W = c(temp_W %*% alpha)
  }
  
  fit = myTryCatch(lmer(y ~ tt + I(tt^2) +
                          W + W * tt + W * I(tt^2) + 
                          trt * tt + trt *  I(tt^2) +  
                          trt * W * tt + trt * W * I(tt^2) + 
                          (tt + I(tt^2)|subj),
                        data = dat_try, REML = FALSE))
  
  if(is.null(fit$warning) == 1){
    fit_cov = as.matrix(fixef(fit$value))
    
    beta1 = fit_cov[c("(Intercept)", "tt", "I(tt^2)"),] +
      fit_cov[c("trt", "tt:trt", "I(tt^2):trt"),] * 2
    beta2 = fit_cov[c("(Intercept)", "tt", "I(tt^2)"),] + 
      fit_cov[c("trt", "tt:trt", "I(tt^2):trt"),]
    
    gamma1 = fit_cov[c("W","tt:W", "I(tt^2):W"),] + 
      fit_cov[c("W:trt", "tt:W:trt", "I(tt^2):W:trt"),] * 2
    gamma2 = fit_cov[c("W","tt:W", "I(tt^2):W"),] + 
      fit_cov[c("W:trt", "tt:W:trt", "I(tt^2):W:trt"),]
    
    D1 = as.matrix(VarCorr(fit$value)$subj)[1:3, 1:3]
    D2 = D1
    
    temp = covar_list 
    xx = as.numeric(unlist(c(temp)))
    xx = matrix(xx, p, dim(dat_try)[1], byrow = TRUE)
    
    # this is the function derivated from the KL-divergence
    A_0 = -2 + 0.5 *  sum(diag(ginv(D1[2:3, 2:3]) %*% D2[2:3, 2:3])) + 
      0.5 * sum(diag(ginv(D2[2:3, 2:3]) %*% D1[2:3, 2:3]))
    A_1 = t(beta1[2:3] - beta2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                            solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_2 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_3 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (gamma1[2:3] - gamma2[2:3]) 
    
    mu_x = matrix(apply(t(xx), 2, mean), p, 1)
    sigma_x = cov(t(xx))
    
    lambda0 = max(eigen(sigma_x)$values)
    n_max = which(eigen(sigma_x)$values == lambda0)
    
    est_alpha = eigen(sigma_x)$vectors[,n_max]
    est_alpha = matrix(est_alpha, p, 1)
    
    return_res = A_0 + A_1/2 + 
      c(A_3/2) * t(est_alpha) %*% (sigma_x) %*% est_alpha
    return_res =  c(-return_res)
 
  }else{
    print('not convergence')
    return_res = 0
  }
  return(return_res)
}


n = 200
p = 4
alpha = rep(1,p)

truekl = c()
for(theta_angle in seq(0, 180, 30)){
  angles = theta_angle
  theta_angle = theta_angle / 180 * pi
  for(parameters in 1){
    
    # # sample size: 
    # n = 100
    # 
    # # dimension:  
    # p = 4
    # 
    # # angle: 
    # theta_angle = 60
    # theta_angle = theta_angle / 180 * pi
    
    
    # others 
    ni = 7 # number of time points
    tt = as.matrix(0:6) # pt = the order of time points
    X = cbind(matrix(1, length(tt), 1), tt, tt^2)
    
    # generate beta randomly
    beta_drg = as.matrix(c(0,3.1,1),3,1)
    beta_pbo = as.matrix(c(0,3,0.9),3,1)
    
    # generate gamma randomly
    #gamma_drg=matrix(c(0,-0.2,0.1),3,1)
    #gamma_pbo=matrix(c(0,0.1,0.2),3,1)
    
    gamma1 = c(0, 1, 0)
    gamma2 = c(0, cos(theta_angle), sin(theta_angle))
    angel = sum(gamma1 * gamma2)/((t(gamma1) %*% gamma1) * (t(gamma2) %*% gamma2))
    #angel
    gamma_drg = gamma1
    gamma_pbo = gamma2
    angel = round(acos(angel) * 180/3.1415926535)
    print(angel)
    round(gamma1,3)
    round(gamma2,3)
    
    # bi
    set.seed(21)
    eign1 = diag(c(1,0.4,0.5))
    eign2 = matrix(runif(9,0,1),3,3)
    a = eign2 %*% eign1 %*% solve(eign2)
    bi_sigma = t(a) %*% (a)
    eigen(bi_sigma)
    bi_sigma2 = bi_sigma
    
    # sigma
    sigma_drg = 1
    sigma_pbo = 1
    
  }
  truekl= c(truekl, KL_true(alpha))
}






KL = function(alpha){
  
  p = length(selected.covar)
  dat_try = dat
  
  for(wcalculation in 1){
    
    # scale alpha
    alpha = matrix(alpha,p,1)
    alpha = alpha/c(sqrt(t(alpha) %*% alpha))
    alpha = round(alpha, 2)
    
    # calculate W = alpha^T x
    covar_list = dat_try[,selected.covar]
    covar_list = scale(covar_list)
    
    temp_W = unlist(c(covar_list))
    temp_W = matrix(temp_W, dim(dat_try)[1], p)
    dat_try$W = c(temp_W %*% alpha)
    
  }
  
  dat_pbo_est = dat_try[dat_try$trt == 1, ]
  dat_drg_est = dat_try[dat_try$trt == 2, ]
  
  fit_pbo_est = myTryCatch(lmer(y ~ tt + I(tt^2) + W + W * tt +
                                  W * I(tt^2) + (tt+I(tt^2)|subj),
                                data = dat_pbo_est, REML = FALSE))
  fit_drg_est = myTryCatch(lmer(y ~ tt + I(tt^2) + W + W * tt +
                                  W * I(tt^2) + (tt+I(tt^2)|subj),
                                data = dat_drg_est, REML = FALSE))
  
  if(is.null(fit_pbo_est$warning) + is.null(fit_drg_est$warning) == 2){
    # beta, gamma, and D matrix
    beta1 = as.matrix(fixef(fit_drg_est$value))[1:3]
    gamma1 = as.matrix(fixef(fit_drg_est$value))[4:6] 
    D1 = as.matrix(VarCorr(fit_drg_est$value)$subj)[1:3, 1:3] 
    beta2 = as.matrix(fixef(fit_pbo_est$value))[1:3] 
    gamma2 = as.matrix(fixef(fit_pbo_est$value))[4:6] 
    D2 = as.matrix(VarCorr(fit_pbo_est$value)$subj)[1:3, 1:3] 
    
    temp = covar_list 
    xx = as.numeric(unlist(c(temp)))
    xx = matrix(xx, p, dim(dat_try)[1], byrow = TRUE)
    
    # this is the function derivated from the KL-divergence
    A_0 = -2 + 0.5 *  sum(diag(ginv(D1[2:3, 2:3]) %*% D2[2:3, 2:3])) + 
      0.5 * sum(diag(ginv(D2[2:3, 2:3]) %*% D1[2:3, 2:3]))
    A_1 = t(beta1[2:3] - beta2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                            solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_2 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_3 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (gamma1[2:3] - gamma2[2:3]) 
    
    mu_x = matrix(apply(t(xx), 2, mean), p, 1)
    sigma_x = cov(t(xx))
    
    return_res = A_0 + A_1/2 + c(A_2) * t(mu_x) %*% alpha + 
      c(A_3/2) * t(alpha) %*% (sigma_x + mu_x %*% t(mu_x)) %*% alpha  
    
    return_res =  c(-return_res)
  }else{
    print('not convergence')
    return_res = 0
  }
  return(return_res)
}

# KL with two estimation of D, with Lagrange multiplier
KL_lm0= function(alpha){
  
  p = length(selected.covar)
  dat_try = dat
  
  for(wcalculation in 1){
    
    # scale alpha
    alpha = matrix(alpha,p,1)
    alpha = alpha/c(sqrt(t(alpha) %*% alpha))
    alpha = round(alpha, 2)
    
    # calculate W = alpha^T x
    covar_list = dat_try[,selected.covar]
    covar_list = scale(covar_list)
    
    temp_W = unlist(c(covar_list))
    temp_W = matrix(temp_W, dim(dat_try)[1], p)
    dat_try$W = c(temp_W %*% alpha)
    
  }
  
  dat_pbo_est = dat_try[dat_try$trt == 1, ]
  dat_drg_est = dat_try[dat_try$trt == 2, ]
  
  fit_pbo_est = myTryCatch(lmer(y ~ tt + I(tt^2) + W + W * tt +
                                  W * I(tt^2) + (tt+I(tt^2)|subj),
                                data = dat_pbo_est, REML = FALSE))
  fit_drg_est = myTryCatch(lmer(y ~ tt + I(tt^2) + W + W * tt +
                                  W * I(tt^2) + (tt+I(tt^2)|subj),
                                data = dat_drg_est, REML = FALSE))
  
  if(is.null(fit_pbo_est$warning) + is.null(fit_drg_est$warning) == 2){
    # beta, gamma, and D matrix
    beta1 = as.matrix(fixef(fit_drg_est$value))[1:3]
    gamma1 = as.matrix(fixef(fit_drg_est$value))[4:6] 
    D1 = as.matrix(VarCorr(fit_drg_est$value)$subj)[1:3, 1:3] 
    beta2 = as.matrix(fixef(fit_pbo_est$value))[1:3] 
    gamma2 = as.matrix(fixef(fit_pbo_est$value))[4:6] 
    D2 = as.matrix(VarCorr(fit_pbo_est$value)$subj)[1:3, 1:3] 
    
    temp = covar_list 
    xx = as.numeric(unlist(c(temp)))
    xx = matrix(xx, p, dim(dat_try)[1], byrow = TRUE)
    
    # this is the function derivated from the KL-divergence
    A_0 = -2 + 0.5 *  sum(diag(ginv(D1[2:3, 2:3]) %*% D2[2:3, 2:3])) + 
      0.5 * sum(diag(ginv(D2[2:3, 2:3]) %*% D1[2:3, 2:3]))
    A_1 = t(beta1[2:3] - beta2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                            solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_2 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_3 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (gamma1[2:3] - gamma2[2:3]) 
    
    mu_x = matrix(apply(t(xx), 2, mean), p, 1)
    sigma_x = cov(t(xx))
    
    lambda0 = max(eigen(sigma_x)$values)
    n_max = which(eigen(sigma_x)$values == lambda0)
    est_alpha = eigen(sigma_x)$vectors[,n_max]
    est_alpha = matrix(est_alpha, p, 1)
    
    return_res = A_0 + A_1/2 +
      c(A_3/2) * t(est_alpha) %*% (sigma_x) %*% est_alpha 
    return_res =  c(-return_res)
    
  }else{
    print('not convergence')
    return_res = 0
  }
  return(return_res)
}

# KL with one estimation of D, with Lagrange multiplier
KL_lm0_D = function(alpha){
  p = length(selected.covar)
  dat_try = dat
  
  for(wcalculation in 1){
    
    # scale alpha
    alpha = matrix(alpha,p,1)
    alpha = alpha/c(sqrt(t(alpha) %*% alpha))
    alpha = round(alpha, 2)
    
    # calculate W = alpha^T x
    covar_list = dat_try[,selected.covar]
    covar_list = scale(covar_list)
    
    temp_W = unlist(c(covar_list))
    temp_W = matrix(temp_W, dim(dat_try)[1], p)
    dat_try$W = c(temp_W %*% alpha)
  }
  
  fit = myTryCatch(lmer(y ~ tt + I(tt^2) +
                          W + W * tt + W * I(tt^2) + 
                          trt * tt + trt *  I(tt^2) +  
                          trt * W * tt + trt * W * I(tt^2) + 
                          (tt + I(tt^2)|subj),
                        data = dat_try, REML = FALSE))
  
  if(is.null(fit$warning) == 1){
    fit_cov = as.matrix(fixef(fit$value))
    
    beta1 = fit_cov[c("(Intercept)", "tt", "I(tt^2)"),] +
      fit_cov[c("trt", "tt:trt", "I(tt^2):trt"),] * 2
    beta2 = fit_cov[c("(Intercept)", "tt", "I(tt^2)"),] + 
      fit_cov[c("trt", "tt:trt", "I(tt^2):trt"),]
    
    gamma1 = fit_cov[c("W","tt:W", "I(tt^2):W"),] + 
      fit_cov[c("W:trt", "tt:W:trt", "I(tt^2):W:trt"),] * 2
    gamma2 = fit_cov[c("W","tt:W", "I(tt^2):W"),] + 
      fit_cov[c("W:trt", "tt:W:trt", "I(tt^2):W:trt"),]
    
    D1 = as.matrix(VarCorr(fit$value)$subj)[1:3, 1:3]
    D2 = D1
    
    temp = covar_list 
    xx = as.numeric(unlist(c(temp)))
    xx = matrix(xx, p, dim(dat_try)[1], byrow = TRUE)
    
    # this is the function derivated from the KL-divergence
    A_0 = -2 + 0.5 *  sum(diag(ginv(D1[2:3, 2:3]) %*% D2[2:3, 2:3])) + 
      0.5 * sum(diag(ginv(D2[2:3, 2:3]) %*% D1[2:3, 2:3]))
    A_1 = t(beta1[2:3] - beta2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                            solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_2 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
    A_3 = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                              solve(D2[2:3, 2:3])) %*% (gamma1[2:3] - gamma2[2:3]) 
    
    mu_x = matrix(apply(t(xx), 2, mean), p, 1)
    sigma_x = cov(t(xx))
    
    lambda0 = max(eigen(sigma_x)$values)
    n_max = which(eigen(sigma_x)$values == lambda0)
    est_alpha = eigen(sigma_x)$vectors[,n_max]
    est_alpha = matrix(est_alpha, p, 1)
    
    return_res = A_0 + A_1/2 +
      c(A_3/2) * t(est_alpha) %*% (sigma_x) %*% est_alpha 
    return_res =  c(-return_res)
 
  }else{
    print('not convergence')
    return_res = 0
  }
  return(return_res)
}
