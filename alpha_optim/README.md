# Alpha optimization


Files explanation

## Attemp 1 2019-09-09

* Only simulated 1 dataset, n = 1000 or 5000 

* True alpha KL calculated used the true alpha value and estimated beta. gamma, and D

* there are some singular estimation of D matrix

R file called "p1n"
  
R data file called names = paste('p1n_',n,'_a', aaa, '.RData', sep='')

Folder: cluster

## Attemp 2 2019-09-10

* Only simulated 1 dataset 

* True alpha KL calculated used the true alpha value and estimated beta. gamma, and D

* if the matrix is singular, the estimation of KL set as 0. 

R file called "n1"
  
R data file called names = paste('n_',n,'_a', theta_angle, '.RData', sep='')

Folder: cluster, cluster-changekl

## Attemp 3 2019-09-11

* Simulated 500 dataset 

* True alpha KL calculated used the true alpha value and estimated beta. gamma, and D

* if the matrix is singular, the estimation of KL set as 0. 

   + R file called "d_a": used the estimated beta, gamma, D
   + R file called "td_a": used the true beta, gamma, D
  
R data file called "D_a", "TrueD_a"

Folder: cluster2

Problem: forget to change the seed, all 500 results are the same 


## Attemp 4 2019-09-11

* Simulated 500 dataset 

* True alpha KL calculated used the true alpha value and estimated beta. gamma, and D

* if the matrix is singular, the estimation of KL set as 0. 

   + R file called "d_a": used the estimated beta, gamma, D
   + R file called "td_a": used the true beta, gamma, D
  
R data file called "D_a", "TrueD_a"

Folder: cluster2 - FDA_alpha2

Sovled Problem: forget to change the seed, all 500 results are the same 


## Attemp 5 2019-09-14

* Simulated 500 dataset 

* 5000 subjects

* True alpha KL calculated used the true alpha value and estimated beta. gamma, and D

* if the matrix is singular, the estimation of KL set as 0. 

   + R file called "d3_a": used the estimated beta, gamma, D
   + R file called "td3_a": used the true beta, gamma, D
  
R data file called "D3_a", "TrueD3_a"

Folder: cluster3_largersample

Problem: died, cannot handle this large


## Attemp 6 2019-09-15

* Simulated 500 dataset 

* 1000 subjects

* True alpha KL calculated used the true alpha value and estimated beta. gamma, and D

* if the matrix is singular, the estimation of KL set as 0. 

   + R file called "d4_a": used the estimated beta, gamma, D
   + R file called "td4_a": used the true beta, gamma, D
  
R data file called "D4_a", "TrueD4_a"

Folder: cluster4_largesample

Problem solved






