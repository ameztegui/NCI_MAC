###################################################
#### test models
###################################################

####################
# load data
###################

rm(list=ls())

library(likelihood)
library(forcats)
library(dplyr)
library(reshape2)

load ("./Data/NCI_MAC.Rdata")


table (neighbours$CodeSp_IN, neighbours$CodeSp_NEAR)

####################
# function
####################

# this function runs anneal to estimate parameters for species i and model j

par_function <- function(i,j) {
        SPP = as.factor(unique(targets$CodeSp))
        
        # create source files for a species i
        using.targets <- which(targets$CodeSp==levels(SPP)[i])
        using.neighbours <- which(neighbours$CodeSp_IN==levels(SPP)[i])
        wk.targets <- targets[using.targets,]
        wk.neighbours <- neighbours[using.neighbours,]
        
        ## lets refactor with forcats!!!
        n_neighbours <- nrow(wk.neighbours %>% group_by(CodeSp_NEAR) %>%
                summarise (n=n()) %>%
                filter(n > 300))
        
        wk.neighbours$CodeSp_2 <- fct_lump(wk.neighbours$CodeSp_NEAR, n = n_neighbours, other_level = "zOther")


        # Make files in format
        wk.heis <- dcast(data=wk.neighbours, Tree_ID_IN ~ voisins, value.var= "H2009")[,-1]
        wk.distances <- dcast(data=wk.neighbours, Tree_ID_IN ~voisins, value.var= "NEAR_DIST")[,-1]
        wk.species <- dcast(data=wk.neighbours, Tree_ID_IN ~voisins, value.var= "CodeSp_NEAR")[,-1]
        wk.species2 <- dcast(data=wk.neighbours, Tree_ID_IN ~voisins, value.var= "CodeSp_2")[,-1]
        wk.intra <- dcast(data=wk.neighbours, Tree_ID_IN ~voisins, value.var= "intra")[,-1]
        wk.clusters <- dcast(data=wk.neighbours, Tree_ID_IN ~voisins, value.var= "cluster")[,-1]
        wk.shades <- dcast(data=wk.neighbours, Tree_ID_IN ~voisins, value.var= "shade")[,-1]
        wk.phyl <- dcast(data=neighbours, Tree_ID_IN ~voisins, value.var= "Phyl")[,-1]
        wk.targets$sr <- apply(wk.species, 1, function(x)length(unique(x)))


        #--------------------------------------------------------------#
        # LEM
        # Create a matrix of lambda indexes. Each matrix value is the index
        # of that neighbor species' lambda in the lambda vector. These
        # indexes won't change so can be created prior to the run.
        #--------------------------------------------------------------#
        all.spp <- levels(wk.neighbours$CodeSp_2)
        neigh.lambda.indexes <- apply(wk.species2, c(1,2),
                                      function(x) {
                                        if (is.na(x)) {
                                          NA
                                        } else {
                                          which(all.spp == x)
                                        }})

        intra.spp <-levels(as.factor(as.character(unique(unlist(wk.intra)))))
        intra.spp <- intra.spp[!is.na(intra.spp)]
        intra.lambda.indexes <- apply(wk.intra, c(1,2),
                                      function(x) {
                                              if (is.na(x)) {
                                                      NA
                                              } else {
                                                      which(intra.spp == x)
                                              }})

        shades.spp <-levels(as.factor(as.character(unique(unlist(wk.shades)))))
        shades.spp <- shades.spp[!is.na(shades.spp)]
        shades.lambda.indexes <- apply(wk.shades, c(1,2),
                                      function(x) {
                                              if (is.na(x)) {
                                                      NA
                                              } else {
                                                      which(shades.spp == x)
                                              }})


        clusters.spp <-levels(as.factor(as.character(unique(unlist(wk.clusters)))))
        clusters.spp <- clusters.spp[!is.na(clusters.spp)]
        clusters.lambda.indexes <- apply(wk.clusters, c(1,2),
                                       function(x) {
                                               if (is.na(x)) {
                                                       NA
                                               } else {
                                                       which(clusters.spp == x)
                                               }})

           # Model functions ---------------------------------------------------------
        model_functions= list (
                ## Null model
                Model.01a <- function (PotGrowth) {
                        PotGrowth},

                ## Size effect
                Model.02a <- function (PotGrowth, sizeX0, sizeXb) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        PotGrowth * size.effect
                },

                ## Size effect + competition effect
                Model.03a <- function(PotGrowth, sizeX0, sizeXb, alpha, beta, C,D) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        nci <- rowSums((wk.heis^alpha)/(wk.distances^beta), na.rm=T)
                        competition.effect <- exp(-(C/100) * nci^D)
                        PotGrowth * size.effect * competition.effect
                },

                ## Size effect + competition effect (with tree size)
                Model.04a <- function(PotGrowth, sizeX0, sizeXb, alpha, beta,
                                     Cprim, D,gamma) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        nci <- rowSums((wk.heis^alpha)/(wk.distances^beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },


                ## Size effect + competition effect (with species richness)
                Model.05a <- function(PotGrowth, sizeX0, sizeXb, alpha, beta,
                                     Cprim, D,gamma, theta) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        nci <- rowSums((wk.heis^alpha)/(wk.distances^beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$shannon^theta
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },

                ## Size effect + competition effect (with FDif)
                Model.06a <- function(PotGrowth, sizeX0, sizeXb, alpha, beta,
                                     Cprim, D,gamma, theta) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        nci <- rowSums((wk.heis^alpha)/(wk.distances^beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$FDif^theta
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },

                ## Size effect + competition effect (with FDis)
                Model.07a <- function(PotGrowth, sizeX0, sizeXb, alpha, beta,
                                     Cprim, D,gamma, theta) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        nci <- rowSums((wk.heis^alpha)/(wk.distances^beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$FDis^theta
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },


                ## Size effect + competition effect (species-specific lambda)
               Model.08a <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, Cprim, D, gamma) {
                  size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                  lambda.vals <- lambdaVec[neigh.lambda.indexes]
                  dim(lambda.vals) <- dim(neigh.lambda.indexes)
                  nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                  C <- Cprim/1000 * wk.targets$DB2009^gamma
                  competition.effect <- exp(-(C) * (nci)^D)
                  PotGrowth * size.effect * competition.effect
                },

               ## Size effect + competition effect (species-specific lambda, negative lambda)
               Model.08b <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                    alpha, beta, Cprim, D, gamma) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       C <- Cprim/1000 * wk.targets$DB2009^gamma
                       competition.effect <- exp(-(C) * (nci)^D)
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (species-specific lambda, gaussian competition)
               Model.08c <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, compX0, compXb) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       competition.effect <- exp(-0.5*((nci - compX0)/compXb))
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

                ## Size effect + competition effect (intra vs. interspecific lambda)
                Model.09a <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, Cprim, D, gamma) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        lambda.vals <- lambdaVec[intra.lambda.indexes]
                        dim(lambda.vals) <- dim(intra.lambda.indexes)
                        nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },

               ## Size effect + competition effect (intra vs. interspecific lambda, negative lambda)
               Model.09b <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                    alpha, beta, Cprim, D, gamma) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[intra.lambda.indexes]
                       dim(lambda.vals) <- dim(intra.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       C <- Cprim/1000 * wk.targets$DB2009^gamma
                       competition.effect <- exp(-(C) * (nci)^D)
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (intra vs. interspecific lambda, gaussian competition)
               Model.09c <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, compX0, compXb) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[intra.lambda.indexes]
                       dim(lambda.vals) <- dim(intra.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       competition.effect <- exp(-0.5*((nci - compX0)/compXb))
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },


                ## Size effect + competition effect (group-specific lambda)
                 Model.10a <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, Cprim, D, gamma) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        lambda.vals <- lambdaVec[clusters.lambda.indexes]
                        dim(lambda.vals) <- dim(clusters.lambda.indexes)
                        nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },

               ## Size effect + competition effect (group-specific lambda, negative lambda)
               Model.10b <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                    alpha, beta, Cprim, D, gamma) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[clusters.lambda.indexes]
                       dim(lambda.vals) <- dim(clusters.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       C <- Cprim/1000 * wk.targets$DB2009^gamma
                       competition.effect <- exp(-(C) * (nci)^D)
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (group-specific lambda, Gaussian competition)
               Model.10c <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, compX0, compXb) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[clusters.lambda.indexes]
                       dim(lambda.vals) <- dim(clusters.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       competition.effect <- exp(-0.5*((nci - compX0)/compXb))
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

                ## Size effect + competition effect (shade-tolerance lambda, gaussian competition)
                Model.11a <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, Cprim, D, gamma) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        lambda.vals <- lambdaVec[shades.lambda.indexes]
                        dim(lambda.vals) <- dim(shades.lambda.indexes)
                        nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },

               ## Size effect + competition effect (shade-tolerance lambda, negative lambda)
               Model.11b <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                    alpha, beta, Cprim, D, gamma) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[shades.lambda.indexes]
                       dim(lambda.vals) <- dim(shades.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       C <- Cprim/1000 * wk.targets$DB2009^gamma
                       competition.effect <- exp(-(C) * (nci)^D)
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (shade-tolerance lambda, Gaussian competition)
               Model.11c <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, compX0, compXb) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[shades.lambda.indexes]
                       dim(lambda.vals) <- dim(shades.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       competition.effect <- exp(-0.5*((nci - compX0)/compXb))
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

                ## Size effect + competition effect (sp-spec. lambda + species richness)
                Model.12a<- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, Cprim, D, gamma, theta) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        lambda.vals <- lambdaVec[neigh.lambda.indexes]
                        dim(lambda.vals) <- dim(neigh.lambda.indexes)
                        nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$shannon^theta
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },

               ## Size effect + competition effect (sp-spec. lambda (negative lambda) + species richness)
               Model.12b <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                    alpha, beta, Cprim, D, gamma, theta) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$shannon^theta
                       competition.effect <- exp(-(C) * (nci)^D)
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (sp-spec. lambda (negative lambda) + species richness (Gaussian competition)
               Model.12c <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, compX0, compXb) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       competition.effect <- exp(-0.5*((nci - compX0)/compXb))
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (sp-spec. lambda + FDif)
               Model.13a <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, Cprim, D, gamma, theta) {
                        size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                        lambda.vals <- lambdaVec[neigh.lambda.indexes]
                        dim(lambda.vals) <- dim(neigh.lambda.indexes)
                        nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                        C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$FDif^theta
                        competition.effect <- exp(-(C) * (nci)^D)
                        PotGrowth * size.effect * competition.effect
                },

               ## Size effect + competition effect (sp-spec. lambda (negative lambda) + FDif)
               Model.13b <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                    alpha, beta, Cprim, D, gamma, theta) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$FDif^theta
                       competition.effect <- exp(-(C) * (nci)^D)
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (sp-spec. lambda (negative lambda) + FDif, Gaussian competition)
               Model.13c <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, compX0, compXb) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       competition.effect <- exp(-0.5*((nci - compX0)/compXb))
                       competition.effect <- ifelse(is.nan(competition.effect), 0, competition.effect)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (sp-spec. lambda + FDis)
               Model.14a <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                    alpha, beta, Cprim, D, gamma, theta) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$FDis^theta
                       competition.effect <- exp(-(C) * (nci)^D)
                       PotGrowth * size.effect * competition.effect
               },

               ## Size effect + competition effect (sp-spec. lambda (negative lambda) + FDis)
               Model.14b <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                    alpha, beta, Cprim, D, gamma, theta) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       C <- Cprim/1000 * wk.targets$DB2009^gamma * wk.targets$FDis^theta
                       competition.effect <- exp(-(C) * (nci)^D)
                       PotGrowth * size.effect * competition.effect
               },
               ## Size effect + competition effect (sp-spec. lambda (negative lambda) + FDis), Gaussian competition
               Model.14c <- function(PotGrowth, sizeX0, sizeXb, lambdaVec,
                                     alpha, beta, compX0, compXb) {
                       size.effect <- exp(-0.5*((log(((wk.targets$DB2009)/sizeX0)/sizeXb))^2))
                       lambda.vals <- lambdaVec[neigh.lambda.indexes]
                       dim(lambda.vals) <- dim(neigh.lambda.indexes)
                       nci <- rowSums(lambda.vals *(wk.heis ^ alpha)/(wk.distances ^ beta), na.rm=T)
                       nci <- ifelse(nci < 0, 0,nci)
                       competition.effect <- exp(-0.5*((nci - compX0)/compXb))
                       PotGrowth * size.effect * competition.effect
               }

               )


        # Get parameters ---------------------------------------------------------
                # Initial parameters ------------------------------------------------------
        par_model <- list(par_Model.01a = list(PotGrowth = 10,
                                              sd = 2),

                          par_Model.02a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              sd = 2),

                          par_Model.03a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              C = 0.00001,
                                              D =1,
                                              sd = 2),

                          par_Model.04a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              sd = 2),

                          par_Model.05a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              theta = 0.5,
                                              sd = 2),

                          par_Model.06a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              theta = 0.5,
                                              sd = 2),

                          par_Model.07a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              theta = 0.5,
                                              sd = 2),

                          par_Model.08a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              sd = 2),

                          par_Model.08b = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              sd = 2),

                          par_Model.08c = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                               compX0 = 500,
                                               compXb = 100,
                                               sd = 2),

                          par_Model.09a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              lambdaVec=rep(0.5,2),
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              sd = 2),

                          par_Model.09b = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,2),
                                               Cprim = 0.001,
                                               D =1,
                                               gamma = 0.5,
                                               sd = 2),

                          par_Model.09c = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,2),
                                               compX0 = 500,
                                               compXb = 100,
                                               sd = 2),

                          par_Model.10a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              lambdaVec=rep(0.5,nlevels(wk.targets$cluster)),
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              sd = 2),

                          par_Model.10b = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.targets$cluster)),
                                               Cprim = 0.001,
                                               D =1,
                                               gamma = 0.5,
                                               sd = 2),

                          par_Model.10c = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.targets$cluster)),
                                               compX0 = 500,
                                               compXb = 100,
                                               sd = 2),

                          par_Model.11a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              lambdaVec=rep(0.5,nlevels(wk.targets$shade)),
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              sd = 2),

                          par_Model.11b = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.targets$shade)),
                                               Cprim = 0.001,
                                               D =1,
                                               gamma = 0.5,
                                               sd = 2),

                          par_Model.11c = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.targets$shade)),
                                               compX0 = 500,
                                               compXb = 100,
                                               sd = 2),

                          par_Model.12a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              theta = 0.5,
                                              sd = 2),

                          par_Model.12b = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                               Cprim = 0.001,
                                               D =1,
                                               gamma = 0.5,
                                               theta = 0.5,
                                               sd = 2),

                          par_Model.12c = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                               compX0 = 500,
                                               compXb = 100,
                                               sd = 2),

                          par_Model.13a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                              Cprim = 0.001,
                                              D =1,
                                              gamma = 0.5,
                                              theta = 0.5,
                                              sd = 2),

                          par_Model.13b = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                               Cprim = 0.001,
                                               D =1,
                                               gamma = 0.5,
                                               theta = 0.5,
                                               sd = 2),

                          par_Model.13c = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                               compX0 = 500,
                                               compXb = 100,
                                               sd = 2),

                          par_Model.14a = list(PotGrowth = 10,
                                              sizeX0= 10,
                                              sizeXb= 1,
                                              alpha = 2,
                                              beta = 1,
                                              lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                              Cprim = 0.001,
                                              D = 1,
                                              gamma = 0.5,
                                              theta = 0.5,
                                              sd = 2),

                          par_Model.14b = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                               Cprim = 0.001,
                                               D = 1,
                                               gamma = 0.5,
                                               theta = 0.5,
                                               sd = 2),

                          par_Model.14c = list(PotGrowth = 10,
                                               sizeX0= 10,
                                               sizeXb= 1,
                                               alpha = 2,
                                               beta = 1,
                                               lambdaVec=rep(0.5,nlevels(wk.neighbours$CodeSp_2)),
                                               compX0 = 500,
                                               compXb = 100,
                                               sd = 2))

                # Low parameters ------------------------------------------------------
        par_low <- list(par_lo_Model.01a = list(PotGrowth = 0.1,
                                               sd = 0.0001),

                        par_lo_Model.02a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               sd = 0.0001),

                        par_lo_Model.03a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               C = 1e-10,
                                               D =0,
                                               sd = 0.0001),

                        par_lo_Model.04a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               Cprim =1e-10,
                                               D =0,
                                               gamma = -2,
                                               sd = 0.0001),

                        par_lo_Model.05a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               Cprim = 1e-10,
                                               D =0,
                                               gamma = -2,
                                               theta = -2,
                                               sd = 0.0001),

                        par_lo_Model.06a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               Cprim = 1e-12,
                                               D =0,
                                               gamma = -2,
                                               theta = -2,
                                               sd = 0.0001),

                        par_lo_Model.07a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               Cprim = 1e-12,
                                               D =0,
                                               gamma = -2,
                                               theta = -2,
                                               sd = 0.0001),

                        par_lo_Model.08a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               lambdaVec=rep(0,nlevels(wk.neighbours$CodeSp_2)),
                                               Cprim = 1e-12,
                                               D =0,
                                               gamma = -2,
                                               sd = 0.0001),

                        par_lo_Model.08b = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.neighbours$CodeSp_2)),
                                                Cprim = 1e-12,
                                                D =0,
                                                gamma = -2,
                                                sd = 0.0001),

                        par_lo_Model.08c = list(PotGrowth = 0.1,
                                                sizeX0 = 0.5,
                                                sizeXb = 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec = rep(-1,nlevels(wk.neighbours$CodeSp_2)),
                                                compX0  = 0.001,
                                                compXb = 0.0001,
                                                sd = 0.0001),

                        par_lo_Model.09a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               lambdaVec=rep(0,2),
                                               Cprim = 1e-12,
                                               D =0,
                                               gamma = -2,
                                               sd = 0.0001),

                        par_lo_Model.09b = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,2),
                                                Cprim = 1e-12,
                                                D =0,
                                                gamma = -2,
                                                sd = 0.0001),

                        par_lo_Model.09c = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,2),
                                                compX0  = 0.001,
                                                compXb = 0.0001,
                                                sd = 0.0001),

                        par_lo_Model.10a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               lambdaVec=rep(0,nlevels(wk.targets$cluster)),
                                               Cprim = 1e-12,
                                               D =1,
                                               gamma = -2,
                                               sd = 0.0001),

                        par_lo_Model.10b = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.targets$cluster)),
                                                Cprim = 1e-12,
                                                D =1,
                                                gamma = -2,
                                                sd = 0.0001),

                        par_lo_Model.10c = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.targets$cluster)),
                                                compX0  = 0.001,
                                                compXb = 0.0001,
                                                sd = 0.0001),

                        par_lo_Model.11a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               lambdaVec=rep(0,nlevels(wk.targets$shade)),
                                               Cprim = 1e-12,
                                               D =0,
                                               gamma = -2,
                                               sd = 0.0001),

                        par_lo_Model.11b = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.targets$shade)),
                                                Cprim = 1e-12,
                                                D =0,
                                                gamma = -2,
                                                sd = 0.0001),

                        par_lo_Model.11c = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.targets$shade)),
                                                compX0  = 0.001,
                                                compXb = 0.0001,
                                                sd = 0.0001),

                        par_lo_Model.12a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               lambdaVec=rep(0,nlevels(wk.neighbours$CodeSp_2)),
                                               Cprim = 1e-12,
                                               D =0,
                                               gamma = -2,
                                               theta = -2,
                                               sd = 0.0001),

                        par_lo_Model.12b = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.neighbours$CodeSp_2)),
                                                Cprim = 1e-12,
                                                D =0,
                                                gamma = -2,
                                                theta = -2,
                                                sd = 0.0001),

                        par_lo_Model.12c = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.neighbours$CodeSp_2)),
                                                compX0  = 0.001,
                                                compXb = 0.0001,
                                                sd = 0.0001),

                        par_lo_Model.13a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               lambdaVec=rep(0,nlevels(wk.neighbours$CodeSp_2)),
                                               Cprim = 1e-12,
                                               D =0,
                                               gamma = -2,
                                               theta = -2,
                                               sd = 0.0001),

                        par_lo_Model.13b = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.neighbours$CodeSp_2)),
                                                Cprim = 1e-12,
                                                D =0,
                                                gamma = -2,
                                                theta = -2,
                                                sd = 0.0001),

                        par_lo_Model.13c = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(-1,nlevels(wk.neighbours$CodeSp_2)),
                                                compX0  = 0.001,
                                                compXb = 0.0001,
                                                sd = 0.0001),

                        par_lo_Model.14a = list(PotGrowth = 0.1,
                                               sizeX0= 0.5,
                                               sizeXb= 0.5,
                                               alpha = 0.5,
                                               beta = 0.5,
                                               lambdaVec=rep(0,nlevels(wk.neighbours$CodeSp_2)),
                                               Cprim = 1e-12,
                                               D =0,
                                               gamma = -2,
                                               theta = -2,
                                               sd = 0.0001),

                        par_lo_Model.14b = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(0,nlevels(wk.neighbours$CodeSp_2)),
                                                Cprim = 1e-12,
                                                D =0,
                                                gamma = -2,
                                                theta = -2,
                                                sd = 0.0001),
                        par_lo_Model.14c = list(PotGrowth = 0.1,
                                                sizeX0= 0.5,
                                                sizeXb= 0.5,
                                                alpha = 0.5,
                                                beta = 0.5,
                                                lambdaVec=rep(0,nlevels(wk.neighbours$CodeSp_2)),
                                                compX0  = 0.001,
                                                compXb = 0.0001,
                                                sd = 0.0001))

                # High parameters ------------------------------------------------------
        par_high <- list(par_hi_Model.01a = list(PotGrowth = 80,
                                                sd = 10),

                         par_hi_Model.02a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                sd = 10),

                         par_hi_Model.03a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                C = 1000,
                                                D = 3,
                                                sd = 10),

                         par_hi_Model.04a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                sd = 10),

                         par_hi_Model.05a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                theta = 3,
                                                sd = 10),

                         par_hi_Model.06a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                theta = 3,
                                                sd = 10),

                         par_hi_Model.07a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                theta = 3,
                                                sd = 10),

                         par_hi_Model.08a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                sd = 10),

                         par_hi_Model.08b = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                 Cprim = 1000,
                                                 D = 3,
                                                 gamma = 3,
                                                 sd = 10),

                         par_hi_Model.08c = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                 compX0 = 5000,
                                                 compXb = 500,
                                                 sd = 10),

                         par_hi_Model.09a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                lambdaVec=rep(1,2),
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                sd = 10),

                         par_hi_Model.09b = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,2),
                                                 Cprim = 1000,
                                                 D = 3,
                                                 gamma = 3,
                                                 sd = 10),

                         par_hi_Model.09c = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,2),
                                                 compX0 = 5000,
                                                 compXb = 500,
                                                 sd = 10),

                         par_hi_Model.10a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                lambdaVec=rep(1,nlevels(wk.targets$cluster)),
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                sd = 10),

                         par_hi_Model.10b = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.targets$cluster)),
                                                 Cprim = 1000,
                                                 D = 3,
                                                 gamma = 3,
                                                 sd = 10),

                         par_hi_Model.10c = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.targets$cluster)),
                                                 compX0 = 5000,
                                                 compXb = 500,
                                                 sd = 10),

                         par_hi_Model.11a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                lambdaVec=rep(1,nlevels(wk.targets$shade)),
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                sd = 10),

                         par_hi_Model.11b = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.targets$shade)),
                                                 Cprim = 1000,
                                                 D = 3,
                                                 gamma = 3,
                                                 sd = 10),

                         par_hi_Model.11c = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.targets$shade)),
                                                 compX0 = 5000,
                                                 compXb = 500,
                                                 sd = 10),

                         par_hi_Model.12a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                theta = 3,
                                                sd = 10),

                         par_hi_Model.12b = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                 Cprim = 1000,
                                                 D = 3,
                                                 gamma = 3,
                                                 theta = 3,
                                                 sd = 10),

                         par_hi_Model.12c = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                 compX0 = 5000,
                                                 compXb = 500,
                                                 sd = 10),

                         par_hi_Model.13a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                theta = 3,
                                                sd = 10),

                         par_hi_Model.13b = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                 Cprim = 1000,
                                                 D = 3,
                                                 gamma = 3,
                                                 theta = 3,
                                                 sd = 10),

                         par_hi_Model.13c = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                 compX0 = 5000,
                                                 compXb = 500,
                                                 sd = 10),

                         par_hi_Model.14a = list(PotGrowth = 80,
                                                sizeX0= 100,
                                                sizeXb= 100,
                                                alpha = 4,
                                                beta = 4,
                                                lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                Cprim = 1000,
                                                D = 3,
                                                gamma = 3,
                                                theta = 3,
                                                sd = 10),

                         par_hi_Model.14b = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                 Cprim = 1000,
                                                 D = 3,
                                                 gamma = 3,
                                                 theta = 3,
                                                 sd = 10),

                         par_hi_Model.14c = list(PotGrowth = 80,
                                                 sizeX0= 100,
                                                 sizeXb= 100,
                                                 alpha = 4,
                                                 beta = 4,
                                                 lambdaVec=rep(1,nlevels(wk.neighbours$CodeSp_2)),
                                                 compX0 = 5000,
                                                 compXb = 500,
                                                 sd = 10)
                         )




# Annealing algorithm -----------------------------------------------------
        # Define the model to use
        # this part used to get functions from Modeles.R, which are now a in list of functions in the script
        modelname<- c('Model.01a','Model.02a','Model.03a', 'Model.04a', 'Model.05a',
                      'Model.06a', 'Model.07a','Model.08a', 'Model.08b','Model.08c',
                      'Model.09a', 'Model.09b','Model.09c','Model.10a','Model.10b', 'Model.10c',
                      'Model.11a','Model.11b', 'Model.11c','Model.12a','Model.12b','Model.12c',
                      'Model.13a','Model.13b', 'Model.13c', 'Model.14a','Model.14b', 'Model.14c' )

        model.flo = model_functions[[j]]

        # Parameter list (for effect of stage)
        par.flo = par_model[[j]]
        par_hi.flo = par_high[[j]]
        par_lo.flo  = par_low[[j]]
        
        # NOVELTY: substitute PotGrowth by actual max value                
        par_lo.flo$PotGrowth = 0.75*max(wk.targets$Dgrowth)
        par.flo$PotGrowth = max(wk.targets$Dgrowth)
        par_hi.flo$PotGrowth = 1.5*max(wk.targets$Dgrowth)
        
        #  set up the list of variables to use in annealing
        pdf.flo=dnorm
        dep.var="Dgrowth"
        var.flo <- list(mean = "predicted", x = "Dgrowth", log=T)

        #  Call the annealing algorithm, specifying the proper model
        results=anneal(model = model.flo, var = var.flo,
                        source_data =wk.targets,dep_var="Dgrowth",
                        par = par.flo, par_lo=par_lo.flo, par_hi=par_hi.flo,
                        pdf=pdf.flo, hessian = FALSE, max_iter=20000)

       write_results(results,paste("./Results/Model_Fits/","HEI_", levels(SPP)[i], modelname[j],".txt", sep=""))

        }


