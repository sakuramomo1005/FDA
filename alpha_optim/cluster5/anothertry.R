setwd("/gpfs/home/ly1192/FDA_alpha_largesample")
n = 500
niters = 500

for(theta_angle in c(0, 30, 60)){
   system(paste("Rscript functions0915_2.R ", theta_angle, niters , n), intern = TRUE)
}
