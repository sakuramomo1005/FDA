# setwd('/Users/yaolanqiu/Desktop/NYU/Research/FDA/alpha_optimization_2019-09-07/cluster3-largersample')

# use estimated D

setwd("/gpfs/home/ly1192/FDA_alpha_largesample")

n = 2500
p = 4
niters = 500
# angle: 
theta_angle = 0
angles = theta_angle
theta_angle = theta_angle / 180 * pi
source('functions0915.R')

names = paste('D3','_a_', angles, '.RData', sep='')
save(result, file = names)