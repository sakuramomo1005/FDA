# setwd('/Users/yaolanqiu/Desktop/NYU/Research/FDA/alpha_optimization_2019-09-07/cluster3-largersample')

# use true D

n = 20
p = 4
niters = 2
# angle: 
theta_angle = 60
angles = theta_angle
theta_angle = theta_angle / 180 * pi
source('functions_trueD0915.R')

names = paste('trueD3','_a_', angles, '.RData', sep='')
save(result, file = names)