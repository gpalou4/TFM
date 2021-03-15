rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(MOFA2)
cat("\n\n MOFA LOADED \n\n", file = stdout())
library(reticulate)
cat("\n\n reticulate LOADED \n\n", file = stdout())
#system("python3 -c 'import mofapy2; print(\"mofapy2 Version:\", mofapy2.__version__); print(\"mofapy2 Path:\", mofapy2.__path__)'")
use_python("/usr/bin/python3", required=TRUE)

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 2.MOFA TRAINING MODEL SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

cat("\n\n... 1) LOADING MOFA OBJECT ...\n\n", file = stdout())

MOFAobject.name <- load("MOFA/1.create_mofa_object/regicor/MOFAobject_regicor_ewas_nogroups_391.RData")
MOFAobject <- get(MOFAobject.name)
cat("\n\n // MOFAobject //\n", file = stdout())
MOFAobject

cat("\n\n... 2) DEFINE OPTIONS ...\n\n", file = stdout())

cat("\n\n... 2.1) DEFINE DATA OPTIONS ...\n\n", file = stdout())

### DATA OPTIONS ###

# likelihoods: likelihood per view (options are “gaussian”, “poisson”, “bernoulli”)
# scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is FALSE
# scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is FALSE

### DATA OPTIONS ###

data.opts <- get_default_data_options(MOFAobject)
data.opts

cat("\n\n... 2.2) DEFINE MODEL OPTIONS ...\n\n", file = stdout())

### MODEL OPTIONS ###

# num_factors: number of factors
# likelihods: same as in data_opts
# spikeslab_factors: use spike-slab sparsity prior in the factors? default is FALSE.
# spikeslab_weights: use spike-slab sparsity prior in the weights? default is TRUE.
# ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
# ard_weights: use ARD prior in the weights? Default is TRUE.

# Only change the default model options if you are familiar with the underlying mathematical model!

### MODEL OPTIONS ###

model.opts <- get_default_model_options(MOFAobject)
model.opts$num_factors <- 30
model.opts

cat("\n\n... 3) DEFINE TRAIN OPTIONS ...\n\n", file = stdout())

### TRAIN OPTIONS ###

# maxiter: number of iterations. Default is 1000.
# convergence_mode: “fast”, “medium”, “slow”. For exploration, the fast mode is good enough.
# startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
# gpu_mode: use GPU mode? (needs cupy installed and a functional GPU, see https://cupy.chainer.org/)
# stochastic: use stochastic inference?
#   verbose: verbose mode?
#   seed: random seed

### TRAIN OPTIONS ### 

train.opts <- get_default_training_options(MOFAobject)
train.opts$convergence_mode <- "slow"
train.opts

cat("\n\n... 4) STOCHASTIC INFERENCE OPTIONS ...\n\n", file = stdout())

# If the number of samples is very large (at the order of >1e4), you may want to try the stochastic inference scheme. 
# If combined with GPUs, it makes inference significantly faster. However, it requires some additional hyperparameters 
# that in some data sets may need to be optimised (vignette in preparation):
  
### STOCHASTIC INFERENCE OPTIONS ### 

# batch_size: numeric value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# learning_rate: learning rate (we recommend values from 0.5 to 0.75)
# forgetting_rate: forgetting rate (we recommend values from 0.25 to 0.75

### STOCHASTIC INFERENCE OPTIONS ### 

stochastic.opts <- get_default_stochastic_options(MOFAobject)
stochastic.opts

cat("\n\n... 5) BUILD THE MOFA OBJECT ...\n\n", file = stdout())

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data.opts,
  model_options = model.opts,
  training_options = train.opts
  #stochastic_options = stochastic.opts # optional
)

cat("\n// MOFAobject built //\n", file = stdout())
MOFAobject

# Find where "mofapy2" package is located
#py_config()
reticulate::py_discover_config("mofapy2")

cat("\n\n... 6) TRAIN THE MOFA OBJECT ...\n\n", file = stdout())

outfile = "MOFA/2.training_model/regicor/MOFA_trained_regicor_ewas_nogroups_391.hdf5"
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time


cat("\n\n ######################## 2.MOFA TRAINING MODEL SCRIPT ENDS ######################## \n\n", file = stdout())




