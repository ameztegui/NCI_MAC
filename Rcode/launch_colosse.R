###################################################
#### Models: supercomputer file 1
###################################################

# Import required functions
source('./Rcode/NCI_Models_DIA.R')
# source('./Rcode/NCI_Models_HEI.R')
# source('./Rcode/NCI_Models_VOL.R')

# Load libraries
library(Rmpi)
library(snow)
library(doSNOW)
library(foreach)

nCore = as.integer(Sys.getenv("MOAB_PROCCOUNT"))
if (is.na(nCore)){ nCore = 8 }

# # Multiple nodes
cl <- makeMPIcluster(nCore)
registerDoSNOW(cl)

# Single node, do not use. But could be useful.
#library(doParallel)
#cl<-makeCluster(8)
#registerDoParallel(cl)

out <- clusterEvalQ(cl, library(likelihood))
out <- clusterEvalQ(cl, library(magrittr))
out <- clusterEvalQ(cl, library(reshape2))
out <- clusterEvalQ(cl, library(forcats))
out <- clusterEvalQ(cl, library(dplyr))

out <- clusterExport(cl, 'targets')
out <- clusterExport(cl, 'neighbours')
out <- clusterExport(cl, 'par_function')


####################
# multiple model analysis
####################
test<-foreach(j=c(1:28), .packages = "likelihood") %:% foreach(i=c(1:19)) %dopar% par_function(i,j)

#####################
# delete cluster
#####################
stopCluster(cl)


