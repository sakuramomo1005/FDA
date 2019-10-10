## combine results 
## 1007

setwd('/Users/yaolanqiu/Desktop/NYU/Research/FDA/alpha_optimization_2019-10-05/cluster/result')

library('lsa')


mean_true_1 = c();
mean_true_2 = c()

kl1 = kl2 = kl3 = c()
kl_list1 = kl_list2 = kl_list3 = c()

Cos1 = Cos2 = Cos3 = c()
cos_list1 = cos_list2 = cos_list3 = list()

vec1 = c(0.5, 0.5, 0.5, 0.5)

counts = 0
for(angles in c(0, 30, 60, 90, 120,150,180)){
  counts = counts + 1
  
  names = paste('D','_a_', angles, '.RData', sep='')
  load(names)
  true_kl1 = mean(result$true_KL1)
  true_kl2 = mean(result$true_KL2)
  
  mean_true_1 = c(mean_true_1, true_kl1)
  mean_true_2 = c(mean_true_2, true_kl2)
  
  r1 = -result$optim_vaiue1
  r2 = -result$optim_vaiue2
  r3 = -result$optim_vaiue3
  
  rr1 = round(r1,2)
  rr1 = ifelse(rr1 > true_kl1 * 2, 0, rr1)
  rr1 = rr1[rr1!=0]
  
  rr2 = round(r2,2)
  rr2 = ifelse(rr2 > true_kl1 * 2, 0, rr2)
  rr2 = rr2[rr2!=0]

  rr3 = round(r3,2)
  rr3 = ifelse(rr3 > true_kl1 * 2, 0, rr3)
  rr3 = rr3[rr3!=0]
  
  kl1 = c(kl1, mean(rr1))
  kl2 = c(kl2, mean(rr2))
  kl3 = c(kl3, mean(rr3))
  
  kl_list1[[counts]] = rr1
  kl_list2[[counts]] = rr2
  kl_list3[[counts]] = rr3
  
  cos1 = c(); cos2 = c(); cos3 = c()
  for(i in 1:500){
    vec2 = result$optim_res1[[i]]
    vec2 = vec2/sqrt(sum(vec2^2))
    cos1 = c(cos1, cosine(vec1,vec2))
    
    vec2 = result$optim_res2[[i]]
    vec2 = vec2/sqrt(sum(vec2^2))
    cos2 = c(cos2, cosine(vec1,vec2))
    
    vec2 = result$optim_res3[[i]]
    vec2 = vec2/sqrt(sum(vec2^2))
    cos3 = c(cos3, cosine(vec1,vec2))
  }
  cos_list1[[counts]] = cos1
  cos_list2[[counts]] = cos2
  cos_list3[[counts]] = cos3
  
  Cos1 = c(Cos1, mean(cos1))
  Cos2 = c(Cos2, mean(cos2))
  Cos3 = c(Cos3, mean(cos3))
  
}

vec2 = vec2/sqrt(sum(vec2^2))

cosine(vec1,vec2)

result_table = data.frame(angels = c(0, 30, 60, 90, 120,150,180),
                          cos = Cos3, 
                          true_KL = mean_true_2,
                          KL = kl2,
                          KL2 = kl3)

save(result_table, file ='result_table.RData')

result_combine = list(kl1 = kl1, 
kl2 = kl2, 
kl3 = kl3, 
kl_list1 = kl_list1, 
kl_list2 = kl_list2, 
kl_list3 = kl_list3, 
Cos1 = Cos1, Cos2 = Cos2, Cos3 = Cos3, 
cos_list1 = cos_list1, 
cos_list2 = cos_list2, 
cos_list3 = cos_list3)

save(result_combine, file = 'result_combine.RData')


pdf('hist_2D.pdf', height = 6, width = 8)
par(mfrow = c(2,4))
for(i in 1:7){
  hist(kl_list2[[i]], breaks = 11,
       main = c(paste('Angle =', seq(0,180,30)[i]), 
                    paste('True purity =', round(mean_true_2[i],2))), 
       xlab = 'Estimated Purity')
  abline(v = mean_true_2[i], col = 2, lty = 2, lwd = 2)
}
dev.off()

pdf('hist_1D.pdf', height = 6, width = 8)
par(mfrow = c(2,4))
for(i in 1:7)
  hist(kl_list3[[i]], breaks = 11, 
       main = c(paste('Angle =', seq(0,180,30)[i]), 
                paste('True purity =', round(mean_true_2[i],2))), 
       xlab = 'Estimated Purity')
  abline(v = mean_true_2[i], col = 2, lty = 2, lwd = 2)
}
dev.off()


hist(kl_list2[[1]], breaks = 11)
abline(v = mean_true_2[1], col = 2, lty = 2)

i = 4
hist(kl_list1[[i]], breaks = 11)
abline(v = mean_true_1[i], col = 2, lty = 2, cex = 5)


load('one_table.RData')
one_table

load('two_table.RData')
two_table

load('four_table.RData')
four_table

load('eight_table.RData')
eight_table

load('two_table.RData')
two_table



