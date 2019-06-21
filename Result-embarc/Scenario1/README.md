# Scenario1

Choose two covariates: age + severity (continuous, binary), scaled. 

The purity plot

When alpha = [0.9092974, -0.4161468], get the max purity. 

## Results: 

1. traj-scenario1.pdf: trajectory plots for 4 clusters and two treatment group, clustered with w = alpha * [age, severity].

2. traj-scenario1_sameW.pdf: trajectory plots for 4 clusters and two treatment group, while w = max(alpha * [age, severity]). (gamma was estimated with the really w and then used max(w) in the following steps)

3. ellipse-scenario1.pdf: ellipse plots of the slope and concavity, as different *w* values.

4. boundary plots: scenario1-boundary-diff-W-group1.pdf, scenario1-boundary-diff-W-group2.pdf: boundary were draw with different W, and points were with different W

5. boundary plots: scenario1-boundary-maxW-group1.pdf, scenario1-boundary-maxW-group2.pdf: boundary were draw with the same W (maxW), and points were with the same W (maxW)




