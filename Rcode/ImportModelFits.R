#
#
#### Extract data (AIC, R2, Slope and parameter values) from a set of models after performing
# maximum likelihood parameter estimation and model fittig
#
# Aitor Ameztegui CREAF - CTFC - UQAM
# October 2016


rm(list=ls()) 
library(plyr)


### Create a list with all the .txt files in the  directory

path <-"./Results/Model_Fits/"

files <-  list.files(path,pattern = "*.txt")

# Divide according to model
mod01afiles <- list.files(path, pattern = "Model.01a.*")
mod02afiles <- list.files(path, pattern = "Model.02a.*")
mod03afiles <- list.files(path, pattern = "Model.03a.*")
mod04afiles <- list.files(path, pattern = "Model.04a.*")
mod05afiles <- list.files(path, pattern = "Model.05a.*")
mod06afiles <- list.files(path, pattern = "Model.06a.*")
mod07afiles <- list.files(path, pattern = "Model.07a.*")
mod08afiles <- list.files(path, pattern = "Model.08a.*")
mod08bfiles <- list.files(path, pattern = "Model.08b.*")
mod08cfiles <- list.files(path, pattern = "Model.08c.*")
mod09afiles <- list.files(path, pattern = "Model.09a.*")
mod09bfiles <- list.files(path, pattern = "Model.09b.*")
mod09cfiles <- list.files(path, pattern = "Model.09c.*")
mod10afiles <- list.files(path, pattern = "Model.10a.*")
mod10bfiles <- list.files(path, pattern = "Model.10b.*")
mod10cfiles <- list.files(path, pattern = "Model.10c.*")
mod11afiles <- list.files(path, pattern = "Model.11a.*")
mod11bfiles <- list.files(path, pattern = "Model.11b.*")
mod11cfiles <- list.files(path, pattern = "Model.11c.*")
mod12afiles <- list.files(path, pattern = "Model.12a.*")
mod12bfiles <- list.files(path, pattern = "Model.12b.*")
mod12cfiles <- list.files(path, pattern = "Model.12c.*")
mod13afiles <- list.files(path, pattern = "Model.13a.*")
mod13bfiles <- list.files(path, pattern = "Model.13b.*")
mod13cfiles <- list.files(path, pattern = "Model.13c.*")
mod14afiles <- list.files(path, pattern = "Model.14a.*")
mod14bfiles <- list.files(path, pattern = "Model.14b.*")
mod14cfiles <- list.files(path, pattern = "Model.14c.*")

all_names <- paste0(path, files)
mod01anames <- paste0(path, mod01afiles)
mod02anames <- paste0(path, mod02afiles)
mod03anames <- paste0(path, mod03afiles)
mod04anames <- paste0(path, mod04afiles)
mod05anames <- paste0(path, mod05afiles)
mod06anames <- paste0(path, mod06afiles)
mod07anames <- paste0(path, mod07afiles)
mod08anames <- paste0(path, mod08afiles)
mod08bnames <- paste0(path, mod08bfiles)
mod08cnames <- paste0(path, mod08cfiles)
mod09anames <- paste0(path, mod09afiles)
mod09bnames <- paste0(path, mod09bfiles)
mod09cnames <- paste0(path, mod09cfiles)
mod10anames <- paste0(path, mod10afiles)
mod10bnames <- paste0(path, mod10bfiles)
mod10cnames <- paste0(path, mod10cfiles)
mod11anames <- paste0(path, mod11afiles)
mod11bnames <- paste0(path, mod11bfiles)
mod11cnames <- paste0(path, mod11cfiles)
mod12anames <- paste0(path, mod12afiles)
mod12bnames <- paste0(path, mod12bfiles)
mod12cnames <- paste0(path, mod12cfiles)
mod13anames <- paste0(path, mod13afiles)
mod13bnames <- paste0(path, mod13bfiles)
mod13cnames <- paste0(path, mod13cfiles)
mod14anames <- paste0(path, mod14afiles)
mod14bnames <- paste0(path, mod14bfiles)
mod14cnames <- paste0(path, mod14cfiles)

# Extract AICc, R2 and LL from models -------------------------------------

modelcomp <- ldply(all_names,function(filename)    {
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


par_model01 <- ldply(mod01names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],File[1])
        return(parameters) })

par_model02 <- ldply(mod02names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        SRX0 <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        SRXb <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],SRX0[,2], SRXb[,2],File[1])
        return(parameters) })

par_model03 <- ldply(mod03names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        FDX0 <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        FDXb <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],FDX0[,2], FDXb[,2],File[1])
        return(parameters) })

par_model04 <- ldply(mod04names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        C <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            C[,2],D[,2],File[1])
        return(parameters) })

par_model05 <- ldply(mod05names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        Cprim <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        f <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=18, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            Cprim[,2],f[,2],D[,2],File[1])
        return(parameters) })

par_model06 <- ldply(mod06names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        I <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        C <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=18, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            I[,2],C[,2],D[,2],File[1])
        return(parameters) })

par_model07 <- ldply(mod07names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        I <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        Cprim <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        f <-  read.table(extract, skip=18, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=19, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            I[,2],Cprim[,2],f[,2],D[,2],File[1])
        return(parameters) })

par_model08 <- ldply(mod08names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        C <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        FDX0 <-  read.table(extract, skip=18, nrows=1, header=F,  sep="\t")
        FDXb <-  read.table(extract, skip=19, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            C[,2],D[,2],FDX0[,2],FDXb[,2],File[1])
        return(parameters) })

par_model09 <- ldply(mod09names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        C <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        lambdavec1 <-  read.table(extract, skip=18, nrows=1, header=F,  sep="\t")
        lambdavec2 <-  read.table(extract, skip=19, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            C[,2],D[,2],lambdavec1[,2],lambdavec2[,2],File[1])
        return(parameters) })

par_model10 <- ldply(mod10names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        C <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        lambdavec1 <-  read.table(extract, skip=18, nrows=1, header=F,  sep="\t")
        lambdavec2 <-  read.table(extract, skip=19, nrows=1, header=F,  sep="\t")
        lambdavec3 <-  read.table(extract, skip=20, nrows=1, header=F,  sep="\t")
        lambdavec4 <-  read.table(extract, skip=21, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            C[,2],D[,2],lambdavec1[,2],lambdavec2[,2],
                            lambdavec3[,2],lambdavec4[,2],File[1])
        return(parameters) })

par_model11 <- ldply(mod11names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        C <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        lambdavec1 <-  read.table(extract, skip=18, nrows=1, header=F,  sep="\t")
        lambdavec2 <-  read.table(extract, skip=19, nrows=1, header=F,  sep="\t")
        lambdavec3 <-  read.table(extract, skip=20, nrows=1, header=F,  sep="\t")
        lambdavec4 <-  read.table(extract, skip=21, nrows=1, header=F,  sep="\t")
        lambdavec5 <-  read.table(extract, skip=22, nrows=1, header=F,  sep="\t")
        lambdavec6 <-  read.table(extract, skip=23, nrows=1, header=F,  sep="\t")
        lambdavec7 <-  read.table(extract, skip=24, nrows=1, header=F,  sep="\t")
        lambdavec8 <-  read.table(extract, skip=25, nrows=1, header=F,  sep="\t")
        lambdavec9 <-  read.table(extract, skip=26, nrows=1, header=F,  sep="\t")
        lambdavec10 <-  read.table(extract, skip=27, nrows=1, header=F,  sep="\t")
        lambdavec11 <-  read.table(extract, skip=28, nrows=1, header=F,  sep="\t")
        lambdavec12 <-  read.table(extract, skip=29, nrows=1, header=F,  sep="\t")
        lambdavec13 <-  read.table(extract, skip=30, nrows=1, header=F,  sep="\t")
        lambdavec14 <-  read.table(extract, skip=31, nrows=1, header=F,  sep="\t")
        lambdavec15 <-  read.table(extract, skip=32, nrows=1, header=F,  sep="\t")
        lambdavec16 <-  read.table(extract, skip=33, nrows=1, header=F,  sep="\t")
        lambdavec17 <-  read.table(extract, skip=34, nrows=1, header=F,  sep="\t")
        lambdavec18 <-  read.table(extract, skip=35, nrows=1, header=F,  sep="\t")
        lambdavec19 <-  read.table(extract, skip=36, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            C[,2],D[,2],lambdavec1[,2],lambdavec2[,2],
                            lambdavec3[,2],lambdavec4[,2],lambdavec5[,2],lambdavec6[,2],
                            lambdavec7[,2],lambdavec8[,2],lambdavec9[,2],lambdavec10[,2],
                            lambdavec11[,2],lambdavec12[,2],lambdavec13[,2],lambdavec14[,2],
                            lambdavec15[,2],lambdavec16[,2],lambdavec17[,2],lambdavec18[,2],
                            lambdavec19[,2],File[1])
        return(parameters) })

par_model12 <- ldply(mod12names,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        sizeX0 <-  read.table(extract, skip=12, nrows=1, header=F,  sep="\t")
        sizeXb <- read.table(extract, skip=13, nrows=1, header=F,  sep="\t")
        alpha <-  read.table(extract, skip=14, nrows=1, header=F,  sep="\t")
        beta <- read.table(extract, skip=15, nrows=1, header=F,  sep="\t")
        C <-  read.table(extract, skip=16, nrows=1, header=F,  sep="\t")
        D <-  read.table(extract, skip=17, nrows=1, header=F,  sep="\t")
        lambdavec1 <-  read.table(extract, skip=18, nrows=1, header=F,  sep="\t")
        lambdavec2 <-  read.table(extract, skip=19, nrows=1, header=F,  sep="\t")
        lambdavec3 <-  read.table(extract, skip=20, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],sizeX0[,2],sizeXb[,2],alpha[,2], beta[,2],
                            C[,2],D[,2],lambdavec1[,2],lambdavec2[,2],
                            lambdavec3[,2],File[1])
        return(parameters) })

par_pool_null <- ldply(nullnames,function(extract)    {
        PotGrowth <- read.table(extract, skip=11, nrows=1, header=F,  sep="\t")
        File <- extract
        parameters <- cbind(PotGrowth[,2],File[1])
        return(parameters) })
# 
# parameters <- rbind(par_pool_full, par_pool_noho, par_pool_nocv,
#                     par_pool_noba, par_pool_nosh,par_pool_null)
#         names(parameters)<-c("PotRichness","HOa","HOb","CVa", "CVb",
#                              "BAa", "BAb","SHa","SHb","Lambda","File")
#         parameters$Variable<- substr(parameters$File, 27, 30)
#         parameters$PDF<- substr(parameters$File, 32, 36)
#         parameters$Model<- substr(parameters$File, 43, 52)


# write.table(parameters, file = "./Results/ModelParameters_Abun.txt", sep="\t", col.names=NA)
