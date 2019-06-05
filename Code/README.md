# ROADMAP

## Functions

* simulation_data_generation.R: generate simulation data, can work for high dimensions (>2 covariates)
* functions.R: purity calculation functions
* cvxcluster-0513.R: convex based clustering method

(more details inside the code files)

## Code for simulation 

* 0-simulation.R: simulate data to check the alpha that max purity is the true alpha. just need to change new parameters.
* 0-sim_draw_gif.R: use simulated data to draw the ellipses gif. 

## Code for real data analysis

* 1-main code.R: main code for real data analysis, other codes nested in this one. only need to set correct file path
* 2-boundary calculation.R: code for boundary calculation based on MC simulated data
* 3-table1.R: codes to draw the table 
* 4-figures.R: codes to draw figures, fig 5 and 6 in the paper
* 5-gif.R: real data ellipses gif 
* 6-trajectory.R: fited quadratic trajectories

