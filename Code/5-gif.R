# draw gif for real data 
# please run the main code first


####################################################
# fit LMM
dat_est = dat
dat_pbo_est = dat_est[dat_est$trt == 0, ]
dat_drg_est = dat_est[dat_est$trt == 1, ]
fit_drg_est = lmer(y ~ t1 + I(t1^2) + W + W * t1 +
                     W * I(t1^2) + (t1+I(t1^2)|subj),
                   data = dat_drg_est, REML = FALSE)
fit_pbo_est = lmer(y ~ t1 + I(t1^2) + W + W * t1 +
                     W * I(t1^2) + (t1+I(t1^2)|subj),
                   data = dat_pbo_est, REML = FALSE)


####################################################
# fit the model and calculated the beta, gamma and D
beta1 = as.matrix(fixef(fit_drg_est))[1:3]
gamma1 = as.matrix(fixef(fit_drg_est))[4:6]
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[1:3, 1:3]
D1 = solve(A) %*% D1 %*% t(solve(A))
D1 = D1[2:3, 2:3]

beta2 = as.matrix(fixef(fit_pbo_est))[1:3]
gamma2 = as.matrix(fixef(fit_pbo_est))[4:6]
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[1:3, 1:3]
D2 = solve(A) %*% D2 %*% t(solve(A))
D2 = D2[2:3, 2:3]


####################################################
# calculate the line of ellipse's center // gamma lines
m1 = c(); m2 = c()
for( i in seq(-100,100,0.1)){
  mu1 = t(solve(A) %*% (beta1 + gamma1 * i))
  mu2 = t(solve(A) %*% (beta2 + gamma2 * i))
  m1 = rbind(m1, mu1)
  m2 = rbind(m2, mu2)
}
m1 = m1[,2:3]; m2 = m2[,2:3]


####################################################
# draw the series of pngs

setwd(store_doc)
nellipse = 100  # number of points for drawing ellipses
c = 4  

for(i in iseqs){
  png(paste(names,i,'.png',sep=''))
  
  wvalue = i
  
  mu1 = solve(A) %*% (beta1 + gamma1 * wvalue); mu1 = mu1[2:3,]
  mu2 = solve(A) %*% (beta2 + gamma2 * wvalue); mu2 = mu2[2:3,]
  
  points1 = (as.matrix(coef(fit_drg_est)$subj)[, 1:3] + 
               i *  as.matrix(coef(fit_drg_est)$subj)[, 4:6]) %*% t(solve(A))
  points2 = (as.matrix(coef(fit_pbo_est)$subj)[, 1:3] + 
               i *  as.matrix(coef(fit_pbo_est)$subj)[, 4:6]) %*% t(solve(A))
  points1 = points1[, 2:3]; points2 = points2[, 2:3]
  
  plot(points1[,1],points1[,2], cex = 0.6, pch = 20, asp=1, 
       col = 'orangered', xlab = "Slope", ylab = "Concavity",
       main = bquote(" " ~ alpha^T ~ x == .(i)), 
       xlim = c(-25, 10))
  points(points2[,1],points2[,2], cex = .6, pch = 18, col = 'lightblue')
  
  epbo = eigen(D1)
  eproz = eigen(D2)
  theta = seq(0,2*pi, length.out = nellipse)
  ellip1 = cbind(cos(theta), sin(theta))
  ellip2 = ellip1
  ellip1 = ellip1 %*% sqrt(diag(c*epbo$values)) %*% t(epbo$vectors)
  ellip2 = ellip2 %*% sqrt(diag(c*eproz$values)) %*% t(eproz$vectors)
  ellip1 = ellip1 + t(matrix(mu1, 2, nellipse))
  ellip2 = ellip2 + t(matrix(mu2, 2, nellipse))
  
  lines(ellip1, col = 2, lty = 2, lwd = 3)
  lines(ellip2, lwd = 3, col = 'blue', lty = 3)
  
  x1 = colMeans(points1)
  x2 = (0.2*eigen(D1)$values[1]*eigen(D1)$vectors[,1]+colMeans(points1)) 
  x3 = (0.2*eigen(D1)$values[2]*eigen(D1)$vectors[,2]+colMeans(points1))
  arrows(x1[1],x1[2],x2[1],x2[2], col = 'black', lwd = 2, length = 0.1)
  arrows(x1[1],x1[2],x3[1],x3[2], col = 'black', lwd = 2, length = 0.1)
  
  x1 = colMeans(points2)
  x2 = (0.3*eigen(D2)$values[1]*eigen(D2)$vectors[,1]+colMeans(points2)) 
  x3 = (0.3*eigen(D2)$values[2]*eigen(D2)$vectors[,2]+colMeans(points2))
  arrows(x1[1],x1[2],x2[1],x2[2], col = 'black', lwd = 2, length = 0.1)
  arrows(x1[1],x1[2],x3[1],x3[2], col = 'black', lwd = 2, length = 0.1)
  
  lines(m1[,1],m1[,2], lty = 2, col = 'darkgreen')
  lines(m2[,1],m2[,2], lty = 2, col = 'darkgreen')
  
  text(-22,5,expression(gamma[1]), cex = 0.8)
  text(-22,-10, expression(gamma[2]), adj = c(0,0), cex = 0.8)
  
  legend('topleft',legend = c('Drug','Placebo'), 
         col = c('orangered', 'lightblue'), pch = c(20, 18), cex = 0.8)
  dev.off()
}


####################################################
# plot the one without gamma
wvalue = 0

mu1 = solve(A) %*% (beta1 + gamma1 * wvalue); mu1 = mu1[2:3,]
mu2 = solve(A) %*% (beta2 + gamma2 * wvalue); mu2 = mu2[2:3,]

points1 = (as.matrix(coef(fit_drg_est)$subj)[, 1:3] + 
             i *  as.matrix(coef(fit_drg_est)$subj)[, 4:6]) %*% t(solve(A))
points2 = (as.matrix(coef(fit_pbo_est)$subj)[, 1:3] + 
             i *  as.matrix(coef(fit_pbo_est)$subj)[, 4:6]) %*% t(solve(A))
points1 = points1[, 2:3]; points2 = points2[, 2:3]

plot(points1[,1],points1[,2], cex = 0.6, pch = 20, asp=1, 
     col = 'orangered', xlab = "Slope", ylab = "Concavity",
     main = bquote(" " ~ alpha^T ~ x == .(i)), 
     xlim = c(-25, 10))
points(points2[,1],points2[,2], cex = .6, pch = 18, col = 'lightblue')

epbo = eigen(D1)
eproz = eigen(D2)
theta = seq(0,2*pi, length.out = nellipse)
ellip1 = cbind(cos(theta), sin(theta))
ellip2 = ellip1
ellip1 = ellip1 %*% sqrt(diag(c*epbo$values)) %*% t(epbo$vectors)
ellip2 = ellip2 %*% sqrt(diag(c*eproz$values)) %*% t(eproz$vectors)
ellip1 = ellip1 + t(matrix(mu1, 2, nellipse))
ellip2 = ellip2 + t(matrix(mu2, 2, nellipse))

lines(ellip1, col = 2, lty = 2, lwd = 3)
lines(ellip2, lwd = 3, col = 'blue', lty = 3)

x1 = colMeans(points1)
x2 = (0.2*eigen(D1)$values[1]*eigen(D1)$vectors[,1]+colMeans(points1)) 
x3 = (0.2*eigen(D1)$values[2]*eigen(D1)$vectors[,2]+colMeans(points1))
arrows(x1[1],x1[2],x2[1],x2[2], col = 'black', lwd = 2, length = 0.1)
arrows(x1[1],x1[2],x3[1],x3[2], col = 'black', lwd = 2, length = 0.1)

x1 = colMeans(points2)
x2 = (0.3*eigen(D2)$values[1]*eigen(D2)$vectors[,1]+colMeans(points2)) 
x3 = (0.3*eigen(D2)$values[2]*eigen(D2)$vectors[,2]+colMeans(points2))
arrows(x1[1],x1[2],x2[1],x2[2], col = 'black', lwd = 2, length = 0.1)
arrows(x1[1],x1[2],x3[1],x3[2], col = 'black', lwd = 2, length = 0.1)

lines(m1[,1],m1[,2], lty = 2, col = 'darkgreen')
lines(m2[,1],m2[,2], lty = 2, col = 'darkgreen')

text(-22,5,expression(gamma[1]), cex = 0.8)
text(-22,-10, expression(gamma[2]), adj = c(0,0), cex = 0.8)

legend('topleft',legend = c('Drug','Placebo'), 
       col = c('orangered', 'lightblue'), pch = c(20, 18), cex = 0.8)

setwd(load_doc)
