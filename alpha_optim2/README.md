
## MEETING SEP 23 . 

TO DO LIST

1. derivation
   + previous method we could get a close form of purity. however, the derivation doesn't makes a lot of sense.        maybe we could treat beta gamma and D as the function of alpha and use chain rule
   + we get another form of purity, calculate the integral
   
2. plot
   + more bars in the histogram
   
3. If the optim function runs too slow, maybe we can consider the taylor expension and get the first or second order? 

4. Embarc dataset 
   + use one covariate each time to calculate the purity
   + use 3 or 4 to calculate the purity
   + the purity by each one of the covariates should be smaller than the combination of the covariates
   
5. Iteration? talked last time
   + given the best alpha, we can get the hat beta, gamma and d
   + given the hat beta, gamma, and d, we can get the best alpha



## After meeting Oct 1. 

Results: 

##### Simulation 

+ multi.R: simulation with p = 10, try purity calculation with different number of covariates combination. 
+ multi.RData: RData file, call "multi", saved the multi table for 2,4,8,10 covariates combination, and time with each of the one covariate. 
+ one_table.RData: results are wrapped in multi.RData, not very useful. 
+ two_table.RData
+ four_table.RData
+ eight_table.RData
+ ten_table.RData

+ result_combine.RData: simulation combined results, with each KL purity and cosine similiarity. 
+ hist_2D: histogram of purity with different angles. Two different LME models were fitted. 
+ hist_2D: histogram of purity with different angles. One model was fitted to make D1 = D2


##### Embarc

+ alt_old_one_table22.RData: purity with only one covaritate included. 
+ two_table2.RData: purity with two covariates included. 
+ three_table2.RData
+ all_table2.RData: with all covaraites included

All results were wraped in the pdf file. 
