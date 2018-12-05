# README #

Replication files for simulations and empirical applications for "Monte Carlo Confidence Sets for Identified Sets" by Chen, Christensen and Tamer.

The replication files are written in R. The additional R packages required to run the files are: foreign, ggplot2, parallel and pbivnorm.

The replication files are set up for a multi-core environment. The makeCluster commands may need to be changed to suit the number of cores in the user's environment. Please refer to the documentation for the parallel package for further details. 

Please note: if using a different number of cores, some variation in the simulations and the output from the SMC algorithm is to be expected due to the way random seeds are handled in a multi-core environment. Please contact us only if running the programs in their original form generates results that are different from what appears in the paper.

## Simulations ##

### Example 1: Missing data ###

There are three folders: flat, curved, and cugmm. These contain replication files for the simulations using a likelihood criterion with flat prior, likelihood criterion with curved prior, and continuously-updated GMM criterion, respectively.

To replicate the results, run files ex1-100-0.R through ex1-1000-2.R. Running file ex1-nN-C.R generates the RData results file ./results/nN/C.RData for N = 100, 250, 500, 1000 and C = 0, 1, 2.

Once all of the .RData files have been generated, run tab_format.R to tabulate the simulation results.

### Example 2: Entry game ###

To replicate the results, run files ex2-100.R through ex2-1000.R. Running file ex1-N.R generates the RData results file ./results/nN/1.RData for N = 100, 250, 500, 1000.

Once all of the .RData files have been generated, run tab_format.R to tabulate the simulation results.

## Empirical Applications ##

### Application 1: Airline entry game ###

To replicate the results, run files game_0.R and game_1.R, which generate the RData results files game0.RData and game1.RData, respectively.

Once the .RData files have been generated, run results_format.R to tabulate the simulation results and generate the plots of the SMC draws.

### Application 2: Trade flows ###

Main files to be run are (in order) Trade_1.R, Trade_p2_2.R, ..., Trade_p2_10.R, Trade_2.R, Trade_3.R then Trade_4.R.

Trade_1.R computes the MLEs, generates the parameter draws from the posterior via SMC, and computes CSs based on inverting t-statistics and CSs using the percentile method.

Trade_p2_2.R, ..., Trade_p2_10.R generates procedure 2 critical values from the SMC draws. As there are a large number of numerical optimizations to be run, each has been set up to be executed in parallel in 20 jobs on a multi-core environment. Running these files generates p2_N_M.RData for N = 2, ..., 10 (each coefficient of interest) and M = 1, ..., 20 (number of jobs).

Trade_2.R computes the quantiles of the profile QLR values stored in the p2_N_M.RData files.

Trade_3.R computes procedure 2 and procedure 3 CSs for the subvectors of interest. Output is stored in subvec_2.RData, ..., subvec_10.RData. This file has been designed to run as an array job, with array index values between 2 and 10.

Trade_4.R tabulates the results.
