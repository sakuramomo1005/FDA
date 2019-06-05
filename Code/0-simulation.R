### simulation 

# to check that the value make purity max is the true alpha used to generate the dataset

# this is a two dimensions example. can work for higher dimensions

setwd('')
source('simulation_data_generation.R')

#############################################
# set parameters: the parameters can change
# generate beta randomly
beta_drg = as.matrix(c(0,-1,3),3,1)
beta_pbo = as.matrix(c(0,3,1),3,1)

# generate gamma randomly
gamma_drg=matrix(c(0,-2,1),3,1)
gamma_pbo=matrix(c(0,1,2),3,1)

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

tt = as.matrix(0:6) # pt = the order of time points
ni = length(tt) # 7
X = cbind(matrix(1, length(tt), 1), tt, tt^2)

n = 100 # number of subjects in each group  
p = 2 # number of baseline covariates
ni = 7 # number of time points

##### set true alpha
theta = pi/3
alpha = c(sin(theta), cos(theta))
set.seed(123)
data = true_generation(alpha, p, n, ni, tt, X, 
                       beta_drg, gamma_drg, bi_sigma,sigma_drg,
                       beta_pbo, gamma_pbo, bi_sigma2,sigma_pbo)
dat_drg = data$dat_drg
dat_pbo = data$dat_pbo
dat = rbind(dat_drg, dat_pbo)
dat$t1 = dat$tt
head(dat)

### find the max purity
purity = c()
for(theta in seq(0, pi,0.01)){
  print(theta)
  temp= c(sin(theta), cos(theta))
  res = purity_function(temp, varname = c('X1','X2'), times = '', 
                        trt = '',
                        trtlevel = '',
                        subj = '',
                        outcome = 'y',
                        start = 0, data = dat)
  purity = c(purity, mean(res$purity))
}

seq(0, pi,0.01)[which(purity == max(purity))]
pi/3

### optim function 
f = function(x){
  temp= c(sin(x), cos(x))
  res = purity_function(temp, varname = c('X1','X2'), times = '', 
                        trt = '',
                        trtlevel = '',
                        subj = '',
                        outcome = 'y',
                        start = 0, data = dat)
  return(-mean(res$purity))
}
optim(1, f, method = 'Brent', lower = 0, upper = 3.15)

# the max value is at the true alpha value

#########################################################################
# clustering 
theta = pi/3
alpha = c(sin(theta), cos(theta))
dat$W = cbind(dat$X1, dat$X2) %*% matrix(alpha,2,1)

# 1. can do the orthogonal transformation
# but the simulated data does not need the transformation, the ellipses look good enough.
A = diag(1,3)

# 2. generate the new dataset with only unique record. 
dat2 = unique(dat[,c('subj','subj2','trt','W')])
dim(dat2) # 200, 3
rownames(dat2) = NULL
head(dat2)

# 3. fit model and calculate beta and gamma values
dat_est = dat
dat_pbo_est = dat_est[dat_est$trt == 'pbo', ]
dat_drg_est = dat_est[dat_est$trt == 'drg', ]
fit_drg_est = lmer(y ~ t1 + I(t1^2) + W + W * t1 +
                     W * I(t1^2) + (t1+I(t1^2)|subj),
                   data = dat_drg_est, REML = FALSE)
fit_pbo_est = lmer(y ~ t1 + I(t1^2) + W + W * t1 +
                     W * I(t1^2) + (t1+I(t1^2)|subj),
                   data = dat_pbo_est, REML = FALSE)

# beta, gamma, and D matrix
beta1 = as.matrix(fixef(fit_drg_est))[2:3] # -0.9127195  2.9606399
gamma1 = as.matrix(fixef(fit_drg_est))[5:6] # 
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[1:3, 1:3] 
beta2 = as.matrix(fixef(fit_pbo_est))[2:3] #  2.958115 1.081628
gamma2 = as.matrix(fixef(fit_pbo_est))[5:6] # 
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[1:3, 1:3] 

# coefficients for each subject
bis_drg = as.matrix(coef(fit_drg_est)$subj) 
dim(bis_drg) # 196
n_drg = dim(bis_drg)[1]
bis_pbo = as.matrix(coef(fit_pbo_est)$subj) 
dim(bis_pbo) # 162
n_pbo = dim(bis_pbo)[1]
bisall = as.data.frame(rbind(bis_pbo, bis_drg))
bisall$group = c(rep(1, n_pbo),rep(2, n_drg))
bisall$subj = rownames(bisall)
bisall = merge(bisall, dat2, by = 'subj')
bisall$W = bisall$W.y
head(bisall)
dim(bisall)

setwd(store_doc)
save(bisall, file = 'sim_bisall.RData')
setwd(load_doc)

# 4. calculate the X matrix (1, t, t^2) after counting the linear combination and orthogonal transformation
bisall$X0 = bisall$`(Intercept)` + bisall$W.x * bisall$W.y
bisall$X1 = bisall$t1 + bisall$`t1:W` * bisall$W.y
bisall$X2 = bisall$`I(t1^2)` + bisall$`I(t1^2):W` * bisall$W.y
bis_trans = cbind(bisall$X0, bisall$X1, bisall$X2) %*% t(solve(A))

# the dataset with coefficients of "slope" and "concavity"
data_trans = data.frame(bis_trans[,2:3])
data_trans$group = bisall$group
data_trans$W = bisall$W.y
data_trans$responder = bisall$responder
colnames(data_trans)[1:2] = c("slope", "concavity")
data_trans$subj = bisall$subj
head(data_trans)

# 5. calculate the expected MVN distributions for drug group and placebo group
# get the centor beta
beta1 = as.matrix(fixef(fit_drg_est))[1:3]
gamma1 = as.matrix(fixef(fit_drg_est))[4:6] 
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[1:3, 1:3] 
beta2 = as.matrix(fixef(fit_pbo_est))[1:3]
gamma2 = as.matrix(fixef(fit_pbo_est))[4:6] 
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[1:3, 1:3]
mu1s = matrix(rep(beta1, each = dim(bisall)[1]),dim(bisall)[1],3) + 
  matrix(bisall$W, dim(bisall)[1],1) %*% matrix(gamma1,1,3)
mu2s = matrix(rep(beta2, each = dim(bisall)[1]),dim(bisall)[1],3) + 
  matrix(bisall$W, dim(bisall)[1],1) %*% matrix(gamma2,1,3)

mu1s_trans = t(solve(A) %*% t(mu1s))
mu2s_trans = t(solve(A) %*% t(mu2s))
D1 = solve(A) %*% D1 %*% t(solve(A))
D2 = solve(A) %*% D2 %*% t(solve(A))
dim(mu1s_trans)
dim(data_trans)

# 6. run the Convexity-Based Clustering Approach algorithm
ns1 = ns2 = 100 # generate a large dataset for MC simulation
pi1 = pi2 = 0.5
# Monte Carlo simulation to calculate the boundary

XSIM0 = c()
for(i in 1:dim(mu1s_trans)[1]){
  
  miu1 = mu1s_trans[i,2:3]
  miu2 = mu2s_trans[i,2:3]
  cov1 = D1[2:3, 2:3]; cov2 = D2[2:3, 2:3]
  d = 2
  
  f1 = list(miu = as.numeric(miu1), cov = as.matrix(cov1))
  f2 = list(miu = as.numeric(miu2), cov = as.matrix(cov2))
  
  xsim_patient = data.frame(mvrnorm(ns1, f1$miu, f1$cov))
  xsim_patient$group = rep(1, ns1)
  xsim_control = data.frame(mvrnorm(ns2, f2$miu, f2$cov))
  xsim_control$group = rep(2, ns2)
  
  xsim0 = rbind(xsim_patient, xsim_control)
  xsim0$lambda = NA
  
  for (i in 1:(ns1+ns2)){
    xsim0$lambda[i] = lambda(xsim0[i,1:d], d, f1,f2,pi1,pi2)
    # calculate the lambda of each point
  }
  
  XSIM0 = rbind(XSIM0, xsim0)
}

# check the lambda values
head(XSIM0)
quantile(XSIM0$lambda)


# 7. calculate the boundary of the clusters 
# (the boundary calculated with the linear transformation W)
XSIM0 = XSIM0[order(XSIM0$group),]
rownames(XSIM0) = NULL
ns1 = sum(XSIM0$group == 1)
ns2 = sum(XSIM0$group == 2)

# calculate the boundary of the clusters 
# (the boundary calculated with the linear transformation W)
nby = which(names(XSIM0) == "group")
k = 4
niter = 100
p1 = clustering(XSIM0, k, d, niter, ns1, ns2)
x = p1$xsim #simulated observations
u = p1$bound #threshold of the clusters
u # 0.07327093 0.43855535 0.86533455
boundary = u

## 11. generate the lambda for the simulated data 
data_trans$lambda = NA
m1 = c(); m2 = c()
for(i in 1:dim(data_trans)[1]){
  temp = data_trans[i,]
  
  miu1 = matrix(fixef(fit_drg_est)[1:3],1,3) + 
    matrix(temp$W, 1,1) %*% matrix(fixef(fit_drg_est)[4:6],1,3)
  miu2 = matrix(fixef(fit_pbo_est)[1:3],1,3) + 
    matrix(temp$W, 1,1) %*% matrix( fixef(fit_pbo_est)[4:6],1,3)
  miu1 = t(solve(A) %*% t(miu1)); miu1 = miu1[2:3]
  miu2 = t(solve(A) %*% t(miu2)); miu2 = miu2[2:3]
  m1 = rbind(m1, miu1); m2 = rbind(m2, miu2)
  input = temp[,c('slope', 'concavity')]
  
  data_trans$lambda[i] = lambda(input, d = 2, f1 = list(miu=miu1, cov=D1[2:3,2:3]), 
                                f2 = list(miu=miu2, cov=D2[2:3,2:3]), 
                                pi1 = 0.5, pi2 = 0.5)
  
  if (!is.na(data_trans$lambda[i])) {
    data_trans$cluster[i] = 1
    for (j in 1:(k-1)){
      if (data_trans$lambda[i] >= u[j]) data_trans$cluster[i] = j + 1
      # u is the boundary of the clusters
    }      
  }
}

table(data_trans$cluster)
# 
# 1   4 
# 180  20 
# the simulated data is too pure. 

#########################################
# draw table1:
data = data_trans
placebo = data[data$group == 1, ]
drug = data[data$group == 2, ]

sum(drug$cluster == 1) # 100
sum(drug$cluster == 4) # 0
sum(placebo$cluster == 1) # 80
sum(placebo$cluster == 4) # 20

table1 = data.frame(drug = c(sum(drug$cluster == 1),sum(drug$cluster == 4)),
           placebo = c(sum(placebo$cluster == 1), sum(placebo$cluster == 4)))
rownames(table1) = c('cluster1','cluster4')
table1

