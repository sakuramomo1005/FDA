# code for draw the figure 5
# please run the main code first

# calculate the boudnary first and they draw the plots 

set.seed(123)
# replace W with the mean W to draw the boundary. 
bisall$W = mean(bisall$W.y)

# made new X0, X1, X2, with combination of baseline variables, which are set as the mean values
bisall$X0 = bisall$`(Intercept)` + bisall$W.x * bisall$W
bisall$X1 = bisall$t1 + bisall$`t1:W` * bisall$W
bisall$X2 = bisall$`I(t1^2)` + bisall$`I(t1^2):W` * bisall$W
bis_trans = cbind(bisall$X0, bisall$X1, bisall$X2) %*% t(solve(A))

data_trans_b = data.frame(bis_trans[,2:3])
data_trans_b$group = bisall$group
data_trans_b$W = bisall$W
data_trans_b$responder = bisall$responder
colnames(data_trans_b)[1:2] = c("slope", "concavity")
head(data_trans_b)

# get the centor beta
bisall$W = bisall$W
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

# run the cvxclustr code
dim(mu1s_trans)
dim(mu2s_trans)
dim(data_trans_b)
ns1 = ns2 = 100
pi1 = pi2 = 0.5
XSIM_scaled0 = c()
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
    xsim0$lambda[i] = lambda(xsim0[i,1:d], d, f1,f2,pi1,pi2)# calculate the lambda of each point
  }
  XSIM_scaled0 = rbind(XSIM_scaled0, xsim0)
}

# get a centored boundary
head(XSIM_scaled0)
quantile(XSIM_scaled0$lambda)

XSIM_scaled0 = XSIM_scaled0[order(XSIM_scaled0$group),]
rownames(XSIM_scaled0) = NULL
ns1 = sum(XSIM_scaled0$group == 1)
ns2 = sum(XSIM_scaled0$group == 2)

# calculate the boundary of the clusters:
nby = which(names(XSIM_scaled0) == "group")
k = 4
niter = 100
p1 = clustering(XSIM_scaled0, k, d, niter, ns1, ns2)
p1$bound
# [1] 0.3611642 0.5341726 0.7361042

x_scaledw = p1$xsim
head(x_scaledw)
table(x_scaledw$cluster)
dim(x_scaledw)
# 1     2     3     4 
# 22397 20953 16546 11704

setwd(load_doc)