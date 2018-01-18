#
#
#### Extract data (AIC, R2, Slope and parameter values) from a set of models after performing
# maximum likelihood parameter estimation and model fittig
#
# Aitor Ameztegui CREAF - CTFC - UQAM
# October 2016


rm(list=ls()) 
library(tidyverse)


### Create a list with all the .txt files in the  directory

path <-"./Results/Model_Fits/"

all_files <-  list.files(path,pattern = "*.txt", full.names = T)

# Divide according to model
mod01a_files <- list.files(path, pattern = "Model.01a.*", full.names = T)
mod02a_files <- list.files(path, pattern = "Model.02a.*", full.names = T)
mod03a_files <- list.files(path, pattern = "Model.03a.*", full.names = T)
mod04a_files <- list.files(path, pattern = "Model.04a.*", full.names = T)
mod05a_files <- list.files(path, pattern = "Model.05a.*", full.names = T)
mod06a_files <- list.files(path, pattern = "Model.06a.*", full.names = T)
mod07a_files <- list.files(path, pattern = "Model.07a.*", full.names = T)
mod08a_files <- list.files(path, pattern = "Model.08a.*", full.names = T)
mod08b_files <- list.files(path, pattern = "Model.08b.*", full.names = T)
mod08c_files <- list.files(path, pattern = "Model.08c.*", full.names = T)
mod09a_files <- list.files(path, pattern = "Model.09a.*", full.names = T)
mod09b_files <- list.files(path, pattern = "Model.09b.*", full.names = T)
mod09c_files <- list.files(path, pattern = "Model.09c.*", full.names = T)
mod10a_files <- list.files(path, pattern = "Model.10a.*", full.names = T)
mod10b_files <- list.files(path, pattern = "Model.10b.*", full.names = T)
mod10c_files <- list.files(path, pattern = "Model.10c.*", full.names = T)
mod11a_files <- list.files(path, pattern = "Model.11a.*", full.names = T)
mod11b_files <- list.files(path, pattern = "Model.11b.*", full.names = T)
mod11c_files <- list.files(path, pattern = "Model.11c.*", full.names = T)
mod12a_files <- list.files(path, pattern = "Model.12a.*", full.names = T)
mod12b_files <- list.files(path, pattern = "Model.12b.*", full.names = T)
mod12c_files <- list.files(path, pattern = "Model.12c.*", full.names = T)
mod13a_files <- list.files(path, pattern = "Model.13a.*", full.names = T)
mod13b_files <- list.files(path, pattern = "Model.13b.*", full.names = T)
mod13c_files <- list.files(path, pattern = "Model.13c.*", full.names = T)
mod14a_files <- list.files(path, pattern = "Model.14a.*", full.names = T)
mod14b_files <- list.files(path, pattern = "Model.14b.*", full.names = T)
mod14c_files <- list.files(path, pattern = "Model.14c.*", full.names = T)



# Extract AICc, R2 and LL from models -------------------------------------

modelcomp <- map_df(all_files,function(filename)    {
        dum <- read.table(filename, skip=4, nrows=1, header=F,  sep="\t")
        dum$filename = filename
        return(dum)})

modelcomp<-modelcomp[,1:8]
        names(modelcomp)<-c("LL","Param", "AICc", "AIC", "Slope", "R2","Null","File")
        modelcomp$Variable_Size<- substr(modelcomp$File, 22, 24)
        modelcomp$Species<- substr(modelcomp$File, 26, 29)
        modelcomp$Model<- substr(modelcomp$File, 30, 38)

write.table(modelcomp, file = "./Results/ModelComparison.txt", sep="\t", row.names = F)



# Extract parameters from the models --------------------------------------

extract_parameters <- function (filename) {
        n_params <- read.table(filename, skip=4, nrows=1, header=F,  sep="\t")[[2]]
        params <- read.table(filename, skip=11, nrows=n_params, header=F,  sep="\t") %>%
                select(1:2)
        param_names <- as.character(params[,1])

        params <- params%>%
                spread(V1, V2)
        
        params <- params[,param_names]  %>%
                mutate(Variable_Size =  substr(filename, 22, 24),
                       Species = substr(filename, 26, 29),
                       Model =  substr(filename, 30, 38)) %>%
                select(-sd)
     
}

pars_model01a <- map_df(mod01a_files, extract_parameters)
pars_model02a <- map_df(mod02a_files, extract_parameters)
pars_model03a <- map_df(mod03a_files, extract_parameters)
pars_model04a <- map_df(mod01a_files, extract_parameters)
pars_model05a <- map_df(mod01a_files, extract_parameters)
pars_model06a <- map_df(mod01a_files, extract_parameters)
pars_model07a <- map_df(mod01a_files, extract_parameters)
pars_model08a <- map_df(mod01a_files, extract_parameters)
pars_model08b <- map_df(mod01a_files, extract_parameters)
pars_model01a <- map_df(mod01a_files, extract_parameters)
pars_model01a <- map_df(mod01a_files, extract_parameters)
pars_model01a <- map_df(mod01a_files, extract_parameters)
pars_model01a <- map_df(mod01a_files, extract_parameters)

        
pars_models <- map_df(all_files, extract_parameters)

write.table(pars_models, file = "./Results/ModelParameters.txt", sep="\t", col.names=NA)
