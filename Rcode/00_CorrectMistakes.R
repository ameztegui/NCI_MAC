rm(list=ls())

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggvis)

# Load the file with all trees in MAC
all_trees <- read.table('./Data/data2009_2014.txt', header = TRUE, sep="\t", dec=".")
# 
# # Create a dataframe with name, leaf habit, origin and join it to the dataset
namesSp <- data.frame(CodeSp=levels(all_trees$CodeSp),
                      NameSp= c("Abies balsamea","Acer platanoides","Acer rubrum","Acer saccharum",
                                "Betula alleghaniensis","Betula papyrifera","Larix decidua","Larix laricina",
                                "Picea abies","Picea glauca","Picea omorika","Pinus resinosa",
                                "Picea rubens", "Pinus strobus","Pinus sylvestris","Quercus robur",
                                "Quercus rubra","Thuja occidentalis","Tilia cordata"),
                      Lhab=c("C","D","D","D","D","D","C","C",
                             "C","C","C","C","C","C","C","D",
                             "D","C","D"),
                      Origin=c("N","E","N","N","N","N","E","N",
                               "E","N","E","N","N","N","E","E",
                               "N","N","E"))
all_trees <- left_join(all_trees,namesSp,by="CodeSp")
all_trees$Tree_ID <- paste(all_trees$Block, all_trees$Plot, all_trees$Pos, sep="_")

# Correct missing and incorrect values
all_trees <- all_trees %>%
                 mutate(
                        DB2009=replace(DB2009,Tree_ID=="C_M8_H2",6.395),
                        DB2009=replace(DB2009,Tree_ID=="D_M2_D2",3.88),
                        DB2009=replace(DB2009,Tree_ID=="D_LADE_H2",2.94),
                        H2009=replace(H2009,Tree_ID=="D_PIAB_H2",52.07),
                        H2009=replace(H2009,Tree_ID=="D_M1_F7",49.5),
                        H2009=replace(H2009,Tree_ID=="D_ACPL_A1",14.73),
                        H2014=replace(H2014,Tree_ID=="B_ACPL_B2",99.8))


 all_trees$DB2009 [all_trees$Block=="A" & all_trees$Plot=="ACRU" & all_trees$DB2009 >10] <-
         all_trees$DB2009 [all_trees$Block=="A" & all_trees$Plot=="ACRU" & all_trees$DB2009 >10]/2.54
 
 change_values <- function (Tree_ID, height) {
         all_trees$H2014[all_trees$Tree_ID == Tree_ID] <- height
         all_trees <<- all_trees
 }
 
 plot_function <- function (sps) {
         fig <- filter(all_trees,CodeSp==!!sps, Etat2014 == "Alive") %>%
                 ggplot(aes(DB2014,H2014))+
                 geom_point() +
                 geom_text(aes(label=Tree_ID), size = 3)
         fig
 }
 
 
plot_function("ABBA")
 plot_function("ACPL")
         change_values ("B_ACPL_B1", 419)
         change_values ("B_M8_H2", 463)
 
 plot_function("ACRU")
 
 plot_function("ACSA")
         change_values ("A_M7_G1", 455)
         change_values ("B_2NR7A_C2", 389)
 
 
 plot_function("BEAL")
        change_values ("D_M1_D1", 415)
 
 plot_function("BEPA")
         change_values ("D_M7_D5", 310)
         change_values ("B_M6_G6", 513)
         change_values ("B_M8_E8", 558)
         change_values ("D_M7_E2", 570)
         change_values ("B_M7_F4", 592)
         change_values ("B_M8_G7", 654)
         
 plot_function("LADE")
 plot_function("LALA")
        change_values ("A_MM6A_F2", 555)
 
 plot_function("PIAB")
 plot_function("PIGL")
        change_values ("C_12N_A3", 300)
 
 plot_function("PIOM")
 plot_function("PIRE")
 plot_function("PIRU")
 plot_function("PIST")
        change_values ("B_PIST_B6", 370)
 
 plot_function("PISY")
 plot_function("QURO")
        change_values ("A_M6_H5", 200)
        change_values ("B_QURO_E8", 415)
        
 plot_function("QURU")
         change_values ("C_4N2_A1", 320)
 
plot_function("THOC")
        change_values ("D_4N7_G8", 225)
        change_values ("A_M8_D3", 276)

 plot_function("TICO")
 

 # Calculate growth and other variables (volumen in dm3)
all_trees <- mutate(all_trees, 
                    V2009 = (DB2009/100)^2 * H2009/10,
                    V2014 = (DB2014/100)^2 * H2014/10,
                    Dgrowth= DB2014 - DB2009, 
                    Hgrowth = H2014-H2009,
                    Vgrowth=V2014-V2009, 
                    Dmean =(DB2009+DB2014)/2, 
                    Hmean = (H2009+H2014)/2, 
                    Vmean=(V2009+V2014)/2)


write.table(all_trees, "./Data/tidydata_2009_2014.txt", row.names = F, sep="\t")
