# misclassICS

Code to reproduce methodology and simulations from "On the Interplay Between Exposure Misclassification and Informative Cluster Size". 


## Functions
Generic functions to generate data and fit models via proposed joint marginalized model.
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

