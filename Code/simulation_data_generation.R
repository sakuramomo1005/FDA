# Data generation for simulation
# Can work for high dimensions

# Notations:
# 1. alpha: input alpha vector
# 2. p: the length of alpha vector, which is also the dimension of baseline variables 
# 3. ni: number of subjects in each intervention arm
# 4. tt: time for longitudinal study
# 5. X: the input time matrix (1, t, t^2)
# 6. beta_drg: treated group beta value, vector, dimension = 3
# 7. gamma_drg: treated group gamma value, vector, dimension = 3
# 8. bi_sigma: treated group b matrix, the covariance matrix
# 9. sigma_drg: treated group random error sigma 
# 10. beta_pbo: placebo group beta value, vector, dimension = 3
# 11. gamma_pbo: placebo group gamma value, vector, dimension = 3
# 12. bi_sigma2: placebo group b matrix, the covariance matrix
# 13. sigma_pbo: placebo group random error sigma

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