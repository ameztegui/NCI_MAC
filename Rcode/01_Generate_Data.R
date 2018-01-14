rm(list=ls())


# Load libraries
library(tidyverse)
library(reshape2)
library(ggvis)
library(readr)

# Load and wrangle the datasets -------------------------------------------

## Load the database with all the trees and the neighbours
all_trees <- read_tsv('./Data/old_all_trees2009_2014.txt') 
raw_neighbours <- read_tsv('./Data/old_raw_neighbours.txt')

## Load the dataset of functional traits,  and neighbour datasets
traits <- read.table('./Data/functional_traits.dta')
std_traits<- data.frame(scale(traits))

## Determine the target trees 
  # First, select only those that have neighbours (from ArcGIS)
  targets <- filter(all_trees, Tree_ID %in% raw_neighbours$Tree_ID_IN)
  targets <- filter(targets,!is.na(Dgrowth))
  targets <- filter(targets, !is.na(H2009))
  
  # But let's keep only the neighbours of the target trees
  neighbours <- raw_neighbours[raw_neighbours$Tree_ID_IN %in% targets$Tree_ID,]

  # And let's delete the target tree as a neighbour of itself
    neighbours<- neighbours %>%
            filter(NEAR_DIST != 0)

  # Let's create a variable with the sequence of neighbours, 
  # and another one with the total number of neighbours
    neighbours$voisins <- sequence(tabulate(neighbours$IN_FID+1))
    
    n_neighbours <- neighbours %>% group_by(Tree_ID_IN) %>%
            summarise(count=n())
    targets$n_neighbours <- n_neighbours$count

    
    # Organize the files to get accumulated diameter, height and volume per each neighbour species 
    melt_sum <- neighbours %>% 
      group_by(Tree_ID_IN, CodeSp_NEAR) %>% 
      summarise (DB2009 = sum(DB2009, na.rm=T),
                 H2009 = sum(H2009, na.rm=T),
                 V2009 = sum (V2009, na.rm = T))
    
    # neig_sps_diam <- as.matrix(dcast(data=melt_sum, Tree_ID_IN~CodeSp_NEAR, value.var= "DB2009")[,-1])
    # neig_sps_heig <- as.matrix(dcast(data=melt_sum, Tree_ID_IN~CodeSp_NEAR, value.var= "H2009")[,-1])
    neig_sps_volu <-dcast(data=melt_sum, Tree_ID_IN~CodeSp_NEAR, value.var= "V2009")[,-1]

    # Organize target files to get the diam, height and volume value for each target, and in the same format as neighbours
    # targ_sps_diam <- as.matrix(dcast(targets, Tree_ID~CodeSp, value.var="DB2009"))[,-1]
    # targ_sps_heig <- as.matrix(dcast(targets, Tree_ID~CodeSp, value.var="H2009"))[,-1]
    targ_sps_volu <- dcast(targets, Tree_ID~CodeSp, value.var="V2009")[,-1]
    
    
    
    
    
    # Calculate Functional Diversity --------------------------------------------------------------
    
    library(FD)
    library(vegan)
    
    # Neighbour Traits (CWM of neighbours)
    FD_neighbours <- dbFD(std_traits,as.matrix(neig_sps_volu),
                          calc.FRic = F,calc.FGR = F,calc.FDiv = F,
                          clust.type = "ward")
    
    # Target Traits (CWM of target)
    FD_targets <- dbFD(std_traits,as.matrix(targ_sps_volu),
                       calc.FRic = F,calc.FGR = F,calc.FDiv = F,
                       clust.type = "ward")
    
    
    
    ## Functional difference ------------------------------------------- 
    
    ### (1) Based on mean of the differences between traits ###
    
    # Difference between target and neighbour CWM values
    traits_dif <- (FD_targets$CWM - FD_neighbours$CWM)
    
    # Quadratic mean
    FDif <- data.frame("FDif"= sqrt(rowSums(traits_dif^2)/ncol(traits_dif))+1)
    targets <- cbind (targets, FDif)
    
    
    ### (2) Based on the difference in functional identity ###
    # Perform PCA to determine Functional Identity
    pca_traits <- princomp(traits, cor=T)
    #summary(pca_traits)
    #loadings(pca_traits)
    #scores(pca_traits)
    PC1=data.frame(Axis1=pca_traits$scores[,1])
    
    # CWM of Axis1 for Neighbours 
    FId_Voisins=functcomp(PC1, as.matrix(neig_sps_volu))
    
    # CWM of Axis1 for Target 
    FId_Cible=functcomp(PC1, as.matrix(targ_sps_volu))
    
    # Determine the difference (in absolute value)
    FDif2<- abs(FId_Voisins - FId_Cible) +1
    colnames(FDif2) <- "FDif2"
    targets <- cbind (targets, FDif2)
    
    
    
    # Calculate Functional Dispersion ---------------------------------------------------------------
    
    FDis <- FD_neighbours$FDis +1
    targets <- cbind (targets, FDis)
    
    
    ## Calculate species richness (Shannon exponent) ---------------------
    library(vegan)
    
    sp_neighbours <- neighbours %>% group_by(Tree_ID_IN)%>%
      summarise(n_sp=length(unique(CodeSp_NEAR)))
    
    # Use the matrix of abundances based on volume, but substitute NA for 0
    neig_sps_volu_NoNA <- neig_sps_volu
    neig_sps_volu_NoNA[is.na(neig_sps_volu_NoNA)] <- 0
    shannon<- exp(diversity(neig_sps_volu_NoNA, index = "shannon"))
    targets <- cbind (targets, shannon)
  
    
    
    
## Species grouping --------------------------------------------------------------
    
    # Cluster grouping (Ward Method) --------------------------------------------------------------
    
    d<-vegdist(std_traits,method="euclidean",na.r=T)
    fit <- hclust(d, method = "ward.D2")
    
    # Plot dendrogram
    plot(fit, cex=0.6, 
         main = " Functional Groups ", 
         hang = -1)
    # Add red boxes to separate groups
    rect.hclust(fit, 4)
    
    
    # Get the cluster for each species
    cluster_groups <- data.frame(CodeSp= levels(as.factor(targets$CodeSp)),cluster=cutree(fit, 4))
    cluster_groups$cluster <- paste0("cluster",as.character(cluster_groups$cluster))
    targets <- left_join(targets, cluster_groups, by=c("CodeSp"))
    neighbours <- left_join(neighbours, cluster_groups, by=c("CodeSp_NEAR"="CodeSp"))
    
    targets$cluster <- as.factor(targets$cluster)
    neighbours$cluster <- as.factor(neighbours$cluster)
    
    # Cluster grouping (shade tolerance) ----------------------------------------------------------
    
    TolS <- as.matrix(traits$ShadeT) #cluster pour la tol?rance ? l'ombrage 
    rownames(TolS)=row.names(traits)
    
    # Ward Hierarchical Clustering
    d_TolS <- dist(TolS, method = "euclidean") # distance matrix
    fit_TolS <- hclust(d_TolS, method="ward.D2") 
    
    # Plot dendrogram
    plot(fit_TolS, cex=0.7, 
         main = " Shade Tolerance Groups ")
    # Add red boxes to separate groups
    rect.hclust(fit_TolS, 3)
    
    # Get the shade clusters for each species
    shade_groups <- data.frame(CodeSp= levels(as.factor(targets$CodeSp)),shade=cutree(fit_TolS, k=3))
    shade_groups$shade <- paste0("shade",as.character(shade_groups$shade))
    targets<- left_join(targets, shade_groups, by=c("CodeSp")) 
    neighbours <- left_join(neighbours, shade_groups, by=c("CodeSp_NEAR"="CodeSp")) 
    
    targets$shade <- as.factor(targets$shade)
    neighbours$shade <- as.factor(neighbours$shade)
    
    
    # Add classification based on leaf habit
    lhab_groups <- data.frame ("CodeSp"= levels(as.factor(targets$CodeSp)),
                               "lhab"= c("C","D","D","D","D","D","C",
                                         "C","C","C","C","C","C","C",
                                         "C","D","D","C","D"))
    neighbours <- left_join(neighbours, lhab_groups, by=c("CodeSp_NEAR"="CodeSp")) #different name
    
    # Add species identity based on intra vs. interspecific
    neighbours <- mutate(neighbours,intra=ifelse(CodeSp_NEAR==CodeSp_IN,as.character(CodeSp_IN),"zother"))
    neighbours$intra <- as.factor(neighbours$intra)
    
    
# Create files for analyses ---------------------------------------------------------------------

    save(targets,neighbours, file = "./Data/NCI_MAC.Rdata")
    
    
    
    


