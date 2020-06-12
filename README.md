# misclassICS

Code to reproduce methodology and simulations from "On the Interplay Between Exposure Misclassification and Informative Cluster Size". 

## Simulations
#### sim1_2019.R,...sim12_2019.R
Run main simulations in each of 12 scenarios. In simulations 1 to 4, gamma=0; in simulations 5 to 8, gamma=-0.25; in simulations 9 to 12, gamma=0.25.In simulations 1,5, and 9, sensitivity and specificity are fixed with respect to cluster size; in simulations 2,6, and 10, sensitivity and specificity decrease in cluster size; in simulations 3,7, and 11, sensitivity and specificity increase in cluster size; in simulations 4,8, and 12, only sensitivity decreases in cluster size.

#### sim1_small.R,...sim12_small.R
Run same simulations in each of 12 scenarios, with reduced sample size for validation set: 10%.

#### sim1INDUCED_2019.R,...,sim3INDUCED_2019.R
Run simulations to show induced informativeness. In each, gamma=0, alpha_1=-1 and \beta_1=1. Sensitivity/specificity are fixed, decrease, and increase with respect to cluster size in simulations 1,2 and 3, respectively.

#### sim6COND_2019.R
Run simulation with conditional model specification (gamma=-0.25; both sensitivity/specificity decrease with cluster size).

#### inducedICS_2019.R
Run simulation to create plots showing bias due to induced informativeness.

#### simulate_sizedist_2019.R
Simulate populations of one million clusters and summarize characteristics for each of twelve scenarios.

#### make_results_tables.R
Summarize results of simulations.


## Functions
Generic functions to generate data, fit models, and correct exposure misclassification.
#### JMMICSpack
Rcpp package containing functions used in fitting joint marginalized models (and joint conditional models).
#### JMMICS.R
Wrapper functions for fitting joint marginalized models (and joint conditional models).
#### genICS.R
Generate data with informative cluster size via Poisson or ZIP model, via shared random intercepts or random slopes. Specifies either a conditionally-specified or marginally-specified model.
#### MISclassICSpack
Rcpp package containing functions used to correct for exposure misclassification under ICS.
#### MISclassICS.R
R functions to fit EEE method.
#### MISclassICS_C.R
Wrapper functions for fitting observed likelihood method.
#### genMIS.R
Function to add exposure misclassification to a dataset (generated from genICS.R).


## Data Analysis
Code for analysis of NHSII data (study of ADHD/Diethylstilbestrol). Note: data not publicly available, so code will not run.
#### 0Mis_analysis_2.R
Run data analysis. Calls "1Mis_clean_data.R" and "2Mis_EDA.R". Fits naive analyses, validation-only analyses, as well as observed likelihood and EEE correction methods.
#### 1Mis_clean_data.R
Clean NHSII data and format for model fitting.
#### 2Mis_EDA.R
Generate descriptive tables and plots.
#### 3Mis_make_results_tables_2.R
Take results from "0Mis_analysis_2.R" and report tables of results with estimates and confidence intervals.
