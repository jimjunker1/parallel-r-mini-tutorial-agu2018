# This is an example script to do parallel runs in R

# Why parallelize?
# E.g.
# Bootstrapping
# Cross-validation
# Fitting multiple regression models
# or
# searching large amounts of data

# There are two ways to do this: with doParallel and with lapply

#set your working directory

rm(list=ls(all=TRUE))

################# Using lapply #######################

# Regular lapply
lapply(1:3, function(x) c(x, x^2, x^3, x^4, x^5))
# or
lapply(1:6, function(x) c(2^x))

# Each element is independently calculated: "embarrasingly parallel"
# (see http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/)

library(parallel)

# detect number of cores (you can also take detectCores()-1 so you don't end up using all your memory)
num_cores <- detectCores()

# Initiate cluster
parallel_cluster <- makeCluster(num_cores)

# parallel version of lapply:
parLapply(parallel_cluster, 1:5,
          function(exponent)
            2^exponent)

# stop the cluster so it is returned to memory
stopCluster(parallel_cluster)

# Note: there is some overhead involved in sending jobs to other nodes so it does have to
# be worth it (obviously not the case in this example)

#### Mac and Windows differences: for Macs you can use "FORK", meaning that the environment also gets copied
# to another node, so it'll recognize objects defined *before* the parallelization. This is not the case for
# Windows. E.g. this works on a Mac but not on Windows:
parallel_cluster<-makeCluster(num_cores, type="FORK")
base <- 2
parLapply(parallel_cluster, 
          2:4, 
          function(exponent) 
            base^exponent)

stopCluster(parallel_cluster)

# For Windows (or on a Mac without defining "type="), you need to use clusterExport:
parallel_cluster<-makeCluster(num_cores) # without the fork
base <- 2
clusterExport(parallel_cluster, "base") # <-- THIS IS ESSENTIAL. You will need to do the same thing for packages that you need
parLapply(parallel_cluster, 
          2:4, 
          function(exponent) 
            base^exponent)

stopCluster(parallel_cluster)

##################### Using loops / foreach ##################### 
# Use the package doParallel. Has MPI in the background. Sends each loop to a different core. See:
# https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf

library(unmarked) # Package to develop occupancy models
library(tidyverse) # Data organization and manipulation
library(doParallel) # To parallelize runs (locations)

# An example without parallelization: bootstrapping!

# load data
wq_data <- read.csv("new_jx_wq.csv")
# We'll use a GLM to model dissolved oxygen from temperature
# Temperature is in column 2 and Dissolved O2 in column 6

# set number of bootstraps
trials <- 20000
# Start timer
init <- Sys.time()
# Do a regular loop
#all_coefficients <- data.frame()
for(i in 1:trials){
  ind <- sample(100, 100, replace=TRUE)
  result1 <- lm(wq_data[ind,6] ~ wq_data[ind,2])
  #all_coefficients <- data.frame(all_coefficients, coefficients(result1))
}
# get final time
Sys.time() - init
# 17.95 secs on my computer

############### Now do it parallel
num_cores <- detectCores()
# Below is an HPC approach (to get the number of cores you requested, not all of them)
#num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
# initiate cluster
myCluster <- makeCluster(num_cores, type = "FORK") 
registerDoParallel(myCluster) # But you can also say registerDoParallel(cores=2)
# The cores argument specifies the number of worker processes that doParallel will use to execute tasks,
# which will by default be equal to one-half the total number of cores on the machine.

init <- Sys.time()

# Use 'foreach', loop over the trials and use '.combine' to bind the results together
result <- foreach(icount(trials), .combine=cbind) %dopar% {
  ind <- sample(100, 100, replace=TRUE)
  result1 <- lm(wq_data[ind,6] ~ wq_data[ind,2])
  coefficients(result1)
}
# in this case cbind: result[1:2,1:15] shows what it looks like
Sys.time() - init
# 13.62 secs
stopCluster(myCluster)
