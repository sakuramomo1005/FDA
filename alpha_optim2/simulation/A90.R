## simulation, 
## kl function made by 1005

# added lagrange multiplier
# tried one D

library(lme4)
library(mixtools)
library(MASS)
##############################################################
#*#*#*#*#*#*#*#*#*#*#*#*#*#*##*#*#*#*#*#*#**#*#*#*##*#*#*#*#*#
# change parameter here 
##############################################################

# sample size: 
n = 200
p = 4

# angle: 
theta_angle = 90
angles = theta_angle
theta_angle = theta_angle / 180 * pi


##############################################################
#*#*#*#*#*#*#*#*#*#*#*#*#*#*##*#*#*#*#*#*#**#*#*#*##*#*#*#*#*#
# Data generation 
##############################################################
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

true_generation = function(alpha, p, n, ni, tt, X, 
                           beta_drg, gamma_drg, bi_sigma, sigma_drg,
                           beta_pbo, gamma_pbo, bi_sigma2,sigma_pbo){
  # alpha
  # set.seed(123)
  alpha = as.matrix(alpha,p,1)
  dat_drg = c()
  for(i in 1:n){
    drg_temp = NULL
    drg_temp$subj = rep(paste('drg',i,sep=''),ni)
    drg_temp$subj2 = rep(i, ni)
    drg_temp$trt = rep('drg',ni)
    drg_temp = as.data.frame(drg_temp)
    baseline = as.matrix(rnorm(p,0,1),p,1)
    drg_temp = cbind(matrix(rep(baseline, each = ni),ni,p),drg_temp)
    colnames(drg_temp)[1:p] = paste('X',1:p, sep = '')
    w = rep(t(alpha) %*% baseline,ni)
    drg_temp$w = w
    drg_temp$tt = tt
    bi = mvrnorm(1, c(0,0,0), bi_sigma)
    yi = X%*%(beta_drg+bi+gamma_drg*w[1]) + sigma_drg*rnorm(ni,0,1)
    drg_temp$y = yi
    dat_drg = rbind(dat_drg, as.data.frame(drg_temp))
  }
  
  dat_pbo = c()
  for(i in 1:n){
    pbo_temp = NULL
    pbo_temp$subj = rep(paste('pbo',i,sep=''),ni)
    pbo_temp$subj2 = rep((n+i), ni)
    pbo_temp$trt = rep('pbo',ni)
    pbo_temp = as.data.frame(pbo_temp)
    baseline = as.matrix(rnorm(p),p,1)
    pbo_temp = cbind(matrix(rep(baseline, each = ni),ni,p),pbo_temp)
    colnames(pbo_temp)[1:p] = paste('X',1:p, sep = '')
    w = rep(t(alpha) %*% baseline,ni)
    pbo_temp$w = w
    pbo_temp$tt = tt
    bi = mvrnorm(1, c(0,0,0), (bi_sigma2))
    yi = X%*%(beta_pbo+bi+gamma_pbo*w[1]) + sigma_pbo*rnorm(ni,0,1)
    pbo_temp$y = yi
    dat_pbo = rbind(dat_pbo, as.data.frame(pbo_temp))
  }
  print('True data generated')
  return(list(dat_drg = dat_drg, dat_pbo = dat_pbo))
}


##############################################################
#*#*#*#*#*#*#*#*#*#*#*#*#*#*##*#*#*#*#*#*#**#*#*#*##*#*#*#*#*#
# Functions 
##############################################################

# True purity 
KL_true = function(alpha){
  
  alpha = matrix(alpha,p,1)
  alpha = alpha/c(sqrt(t(alpha) %*% alpha))
  alpha = round(alpha, 2)
  
  dat_try = dat
  
  beta1 = beta_drg
  beta2 = beta_pbo
  
  gamma1 = gamma_drg
  gamma2 = gamma_pbo
  
  D1 = bi_sigma
  D2 = bi_sigma2
  
  temp = dat[,paste('X',1:p, sep = '')]
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
  
  return_res = A_0 + A_1/2 + c(A_2) * t(mu_x) %*% est_alpha + 
    c(A_3/2) * t(est_alpha) %*% (sigma_x + mu_x %*% t(mu_x)) %*% est_alpha + 
    lambda0 * (t(est_alpha) %*%est_alpha - 1)
  
  return_res2 = A_0 + A_1/2 + c(A_2) * t(mu_x) %*% alpha + 
    c(A_3/2) * t(alpha) %*% (sigma_x + mu_x %*% t(mu_x)) %*% alpha 
  
  res = list(return_res = return_res,
             return_res2 = return_res2)
  
  return(res)
}

# Estimate purity 

KL = function(alpha){
  
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
KL_lm0 = function(alpha){
  
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
    
    return_res = A_0 + A_1/2 + c(A_2) * t(mu_x) %*% est_alpha + 
      c(A_3/2) * t(est_alpha) %*% (sigma_x + mu_x %*% t(mu_x)) %*% est_alpha + 
      lambda0 * (t(est_alpha) %*%est_alpha - 1)
    return_res =  c(-return_res)
    
  }else{
    print('not convergence')
    return_res = 0
  }
  return(return_res)
}
# KL with one estimation of D, with Lagrange multiplier
KL_lm0_D = function(alpha){
  
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
    
    return_res = A_0 + A_1/2 + c(A_2) * t(mu_x) %*% est_alpha + 
      c(A_3/2) * t(est_alpha) %*% (sigma_x + mu_x %*% t(mu_x)) %*% est_alpha + 
      lambda0 * (t(est_alpha) %*%est_alpha - 1)
    return_res =  c(-return_res)
    # 
    # return_res = A_0 + A_1/2 + c(A_2) * t(mu_x) %*% alpha + 
    #   c(A_3/2) * t(alpha) %*% (sigma_x + mu_x %*% t(mu_x)) %*% alpha - lambda0 * (t(alpha) %*% alpha - 1)
    # return_res =  c(-return_res)
    
  }else{
    print('not convergence')
    return_res = 0
  }
  return(return_res)
}


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

true_alpha = c()
true_KL1 = c(); true_KL2 = c()
intis = list()
optim_res1 = list(); optim_res2 = list(); optim_res3 = list()
optim_vaiue1 = c(); optim_vaiue2 = c(); optim_vaiue3 = c()

count = 0 
begins = Sys.time()
for(iters in 1:500){
  
  count = count + 1
  
  # true alpha 
  alpha = rep(1,p)
  alpha = alpha/sqrt(sum(alpha^2))
  sum(alpha^2)
  
  true_alpha[[count]] = alpha
  
  for(datageneration in 1){
    # set.seed(123)
    data = true_generation(alpha, p, n, ni, tt, X, 
                           beta_drg, gamma_drg, bi_sigma,sigma_drg,
                           beta_pbo, gamma_pbo, bi_sigma2,sigma_pbo)
    
    dat_drg = data$dat_drg
    dat_pbo = data$dat_pbo
    dat = rbind(dat_drg, dat_pbo)
    dat$t1 = dat$tt
    print(head(dat))
    print(dim(dat))
    
    alpha = matrix(alpha,p,1)
    alpha = alpha/c(sqrt( t(alpha) %*% alpha))
    temp = dat[,1:p]; temp = unlist(temp); temp = matrix(temp, dim(dat)[1], p)
    dat$W = temp %*% alpha
    
    # data of unique subjects 
    dat2 = unique(dat[,c('subj','subj2','trt','W')])
    print(dim(dat2)) # 200, 3
    rownames(dat2) = NULL
    print(head(dat2))
  }
  
  dat$trt = ifelse(dat$trt == 'drg', 2, 1)
  
  true_KL1 = c(true_KL1, KL_true(alpha)$return_res)
  true_KL2 = c(true_KL2, KL_true(alpha)$return_res2)
  
  inti = alpha 
  
  #intis[[count]] = inti
  
  # let's set the scenarios 
  res1 = optim(inti, KL)
  res2 = optim(inti, KL_lm0)
  res3 = optim(inti, KL_lm0_D)
  
  optim_res1[[count]] = res1$par
  optim_res2[[count]] = res2$par
  optim_res3[[count]] = res3$par
  
  optim_vaiue1 = c(optim_vaiue1, res1$value)
  optim_vaiue2 = c(optim_vaiue2, res2$value)
  optim_vaiue3 = c(optim_vaiue3, res3$value)
}

ends = Sys.time()

aaa = ends - begins

result = list(true_KL1 = true_KL1, true_KL2 = true_KL2, 
              optim_res1 = optim_res1, 
              optim_res2 = optim_res2, 
              optim_res3 = optim_res3, 
              optim_vaiue1 = optim_vaiue1, 
              optim_vaiue2 = optim_vaiue2, 
              optim_vaiue3 = optim_vaiue3, 
              time = aaa)

names = paste('D','_a_', angles, '.RData', sep='')
save(result, file = names)

