

theta_angle = 180
setwd("/gpfs/home/ly1192/FDA_alpha_largesample2")
n = 500; p = 4; niters = 500
angles = theta_angle; theta_angle = theta_angle / 180 * pi
source('functions0915.R')
names = paste('D4','_a_', angles, '.RData', sep='')
save(result, file = names)