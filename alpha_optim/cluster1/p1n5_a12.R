# sample size: 
n = 5000

# angle: 
theta_angle = 120
aaa = theta_angle
theta_angle = theta_angle / 180 * pi

### the things needs to vary.
# sample size 
# dimension 
# angle

library(lme4)
library(mixtools)
library(MASS)

# dataset generation 

true_generation = function(alpha, p, n, ni, tt, X, 
                           beta_drg, gamma_drg, bi_sigma, sigma_drg,
                           beta_pbo, gamma_pbo, bi_sigma2,sigma_pbo){
  # alpha
  set.seed(123)
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
KL = function(alpha){ # the function 
  
  dat_try = dat
  
  # scale alpha
  
  alpha = matrix(alpha,p,1)
  alpha = alpha/c(sqrt(t(alpha) %*% alpha))
  alpha = round(alpha, 2)
  
  # calculate W = alpha^T x
  
  temp_W = unlist(c(dat_try[,1:p]))
  temp_W = matrix(temp_W, dim(dat_try)[1], p)
  dat_try$W = temp_W %*% alpha
  
  # estimate beta, gamma, and D
  
  dat_pbo_est = dat_try[dat_try$trt == 'pbo', ]
  dat_drg_est = dat_try[dat_try$trt == 'drg', ]
  fit_drg_est = lmer(y ~ tt + I(tt^2) + W + W * tt +
                       W * I(tt^2) + (tt+I(tt^2)|subj),
                     data = dat_drg_est, REML = FALSE)
  fit_pbo_est = lmer(y ~ tt + I(tt^2) + W + W * tt +
                       W * I(tt^2) + (tt+I(tt^2)|subj),
                     data = dat_pbo_est, REML = FALSE)
  
  # beta, gamma, and D matrix
  
  beta1 = as.matrix(fixef(fit_drg_est))[1:3]
  gamma1 = as.matrix(fixef(fit_drg_est))[4:6] 
  D1 = as.matrix(VarCorr(fit_drg_est)$subj)[1:3, 1:3] 
  beta2 = as.matrix(fixef(fit_pbo_est))[1:3] 
  gamma2 = as.matrix(fixef(fit_pbo_est))[4:6] 
  D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[1:3, 1:3] 
  
  # this is the function derivated from the KL-divergence
  
  A = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                          solve(D2[2:3, 2:3])) %*% (beta1[2:3] - beta2[2:3])
  B = t(gamma1[2:3] - gamma2[2:3]) %*% (solve(D1[2:3, 2:3]) + 
                                          solve(D2[2:3, 2:3])) %*% (gamma1[2:3] - gamma2[2:3]) 
  
  fun_res = sum(diag(ginv(D1[2:3, 2:3]) %*% D2[2:3, 2:3])) + 
    sum(diag(ginv(D2[2:3, 2:3]) %*% D1[2:3, 2:3]))
  
  temp = dat_try[,paste('X',1:p, sep='')]
  xx = as.numeric(unlist(c(temp)))
  xx = matrix(xx,p,dim(dat_try)[1], byrow = TRUE)
  fun_res1 = as.vector(c(A) * t(xx) %*% alpha + 
                         (t(xx) %*% alpha)^2* c(B))
  
  # for(ss in 1:dim(dat_try)[1]){
  #   temp = dat_try[ss,paste('X',1:p, sep='')]
  #   xx = as.numeric(temp)
  #   xx = matrix(xx,p,1)
  #   fun_res = fun_res + as.vector(c(A) * t(xx) %*% alpha + 
  #                                   t(alpha) %*% xx %*% t(xx) %*% alpha * c(B))
  # }
  
  fun_res = fun_res + sum(fun_res1)
  
  return(-fun_res/dim(dat_try)[1])
  
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
true_KL = c()
intis = list()
optim_res = list()
optim_vaiue = c()

count = 0 

for(p in c(2,4,8,16)){
  
  count = count + 1
  
  # true alpha 
  alpha = rep(1,p)
  alpha = alpha/sqrt(sum(alpha^2))
  sum(alpha^2)
  
  true_alpha[[count]] = alpha
  
  for(datageneration in 1){
    set.seed(123)
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
  
  true_KL = c(true_KL, KL(alpha))
  inti = 1:p
  
  intis[[count]] = inti
  
  # let's set the scenarios 
  res = optim(inti, KL)
  optim_res[[count]] = res$par
  optim_vaiue = c(optim_vaiue, res$value)
}

result = list(p = c(2,4,8,16), true_KL = true_KL, 
              true_alpha = true_alpha, 
              intis = intis,
              optim_res = optim_res, 
              optim_vaiue = optim_vaiue)

names = paste('p1n_',n,'_a', aaa, '.RData', sep='')
save(result, file = names)
