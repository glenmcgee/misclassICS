# misclassICS

Code to reproduce methodology and simulations from "On the Interplay Between Exposure Misclassification and Informative Cluster Size". 

## Simulations
#### sim1_2019.R,...sim12_2019.R

#### simulate_sizedist_2019.R

#### inducedICS_2019.R

#### make_results_tables.R

#### sim1INDUCED_2019.R,...,sim3INDUCED_2019.R

#### sim6COND_2019.R


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
