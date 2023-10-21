## Utilize phytools::Map.Overlap() to get proportion of edges spent in each FS/Coop state
## Kate Snyder
## 4/4/2022
## Method based on Huelsenbeck et al (2003)
## 
## Edited 6/1/23 - added checkpoint save to CharacterSimmaps
## Edited 9/27/23 - changed "Dummy" simmap generation to use sim.history() with Q rates, ancestral character estimation instead of randomizing tip states; but seems to have gotten totally weird - output values odd
## Edited 9/29/23

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")
library(phytools)
source("subsettreedata.R")
#df <- read.csv("2022-03-10CoopSong_All.csv")
#df <- read.csv("2022-03-10CoopSong_PlusWebbFS.csv")
#newdata = "/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2023-06-01_CoopSongFS_RColumns.csv"

newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data_R.csv"
df = read.csv(newdata)

Hacktree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000.nex")
treelabel <- "Hackett"
#Erictree <-read.nexus("2021-08-31ConsensusPasserineTreeEricson10_1000.nex")
#treelabel <- "Ericson"

subsets = subsettreedata(columns = c("MeanCoopTie2Noncoop", "FemaleSong_Agg01"), newdata = df, newtree = Hacktree)    # decided that for these, will put the pre-subsetted dataset into CharacterSimmaps, to be passed to findQrates, so the Qrates and ancestral character estimation aren't based on any of the non-passerines etc. May want to go back to regular full-subset Q/Anc estimation for other trait combinations
subsetdf = subsets$subsetdf
subsettree = subsets$subsettree

# checking whether the Qs from ace and make.simmap are the same - they are!
tempBinTrait <- subsetdf$MeanCoopTie2Noncoop
names(tempBinTrait) <- subsetdf$species
tempSim = make.simmap(subsettree, tempBinTrait, "ARD")
tempSim$Q
ace(tempBinTrait, subsettree, type = "discrete", model = "ARD")

nsims_real = 5
nsims_dummy = 20

dfout5 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Agg01"), df = subsetdf, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Tie2NonCoop FSAgg_SubsetQs")
#write.csv(dfout5, paste0("Simmap overlap output_Tie2NonCoop FSAgg nsim", nsims_real,"_Hackett REAL ", Sys.Date(),".csv"))
#dfout6 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Agg01"), df = df, tree =  Erictree, dummy = FALSE, nsims = nsims_real, treelabel = "Ericson", datalabel = "Tie2NonCoop FSAgg")

dfDummyMkSimmap <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Agg01"), df = subsetdf, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Tie2NonCoop FSAgg_SubsetQs", dummyMethod = "makeSimmap")
dfDummySimHist <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Agg01"), df = subsetdf, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Tie2NonCoop FSAgg_SubsetQs", dummyMethod = "simHistory")

calcHuel(dfout5,dfDummy5, nsims_real = nsims_real, nsims_dummy = nsims_dummy)
calcHuel(dfout6,dfDummy6, nsims_real = nsims_real, nsims_dummy = nsims_dummy)


# dfout5 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","FemaleSong_Aggregated"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Tie2Coop FSAgg")
# #write.csv(dfout5, "Simmap overlap output_Tie2Coop FSAgg nsim1000_Hackett REAL.csv")
# dfout6 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","FemaleSong_Aggregated"), df = df, tree =  Erictree, dummy = FALSE, nsims = nsims_real, treelabel = "Ericson", datalabel = "Tie2Coop FSAgg")
# #write.csv(dfout6, "Simmap overlap output_Tie2Coop FSAgg nsim1000_Ericson REAL.csv")
# 
# dfDummy5 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","FemaleSong_Aggregated"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Tie2Coop FSAgg")
# write.csv(dfDummy5, "Simmap overlap output_Tie2Coop FSAgg nsim10000_Hackett DUMMY.csv")
# dfDummy6 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","FemaleSong_Aggregated"), df = df, tree =  Erictree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Ericson", datalabel = "Tie2Coop FSAgg")
# write.csv(dfDummy6, "Simmap overlap output_Tie2Coop FSAgg nsim10000_Ericson DUMMY.csv")

# dfout <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","Webb_FemaleSong"), df = df, tree =  Hacktree, dummy = FALSE, nsims = 10, treelabel = "Hackett", datalabel = "CoopTie2Coop FSWebb")
# dfDummy <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","Webb_FemaleSong"), df = df, tree =  Hacktree, dummy = TRUE, nsims = 50, treelabel = "Hackett", datalabel = "CoopTie2Coop FSWebb")

# dfout2 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Aggregated"), df = df, tree =  Erictree, dummy = FALSE, nsims = nsims_real, treelabel = "Ericson", datalabel = "Tie2NonCoop FSAgg")
# dfDummy2 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong_Aggregated"), df = df, tree =  Erictree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Ericson", datalabel = "Tie2NonCoop FSAgg")
# 
# dfout3 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","Webb_FemaleSong"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Tie2NonCoop FSWebb")
# dfDummy3 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","Webb_FemaleSong"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Tie2NonCoop FSWebb")
# 
# dfout4 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","Odom_FemaleSong"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Tie2NonCoop FSOdom")
# dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Tie2NonCoop FSOdom")
# 
# dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","HighConfidence_FemaleSong"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Tie2NonCoop FSHighConf")
# 
# 
# dfout4 <- CharacterSimmaps(columns = c("Final.polygyny","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Polygyny FSAgg")
# dfDummy4 <- CharacterSimmaps(columns = c("Final.polygyny","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Polygyny FSAgg")
# dfout5 <- CharacterSimmaps(columns = c("Final.polygyny","HighConfidence_FemaleSong"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Polygyny FSHighConf")
# dfDummy5 <- CharacterSimmaps(columns = c("Final.polygyny","HighConfidence_FemaleSong"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Polygyny FSHighConf")


df$Kin_NK[which(df$Kin_NK == "Mixed")] <- NA
df$Kin_NK[which(df$Kin_NK == "NonKin")] <- 0
df$Kin_NK[which(df$Kin_NK == "Kin")] <- 1
dfout4 <- CharacterSimmaps(columns = c("Kin_NK","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Kin_NK FSAgg")
dfDummy4 <- CharacterSimmaps(columns = c("Kin_NK","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Kin_NK FSAgg")
KinFSreal = read.csv("Kin_NK FSAgg simmap overlap output nsim 1000 Hackett .csv")
KinFSdummy = read.csv("Kin_NK FSAgg DUMMYResampledCoopFS simmap overlap output nsim 1000 Hackett .csv")
nsims_real = length(KinFSreal$treenum)
nsims_dummy = length(KinFSdummy$treenum)
calcHuel(KinFSreal, KinFSdummy, nsims_real = nsims_real, nsims_dummy = nsims_dummy)

dfout4 <- CharacterSimmaps(columns = c("Griesser2017FamilialLiving","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Familial FSAgg")
dfDummy4 <- CharacterSimmaps(columns = c("Griesser2017FamilialLiving","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Familial FSAgg")
FamFSreal = read.csv("Familial FSAgg simmap overlap output nsim 1000 Hackett .csv")
FamFSdummy = read.csv("Familial FSAgg DUMMYResampledCoopFS simmap overlap output nsim 5000 Hackett .csv")
nsims_real = length(FamFSreal$treenum)
nsims_dummy = length(FamFSdummy$treenum)
calcHuel(FamFSreal, FamFSdummy, nsims_real = nsims_real, nsims_dummy = nsims_dummy)


# dfout <- read.csv("CoopBreed FSWebb simmap overlap output nsim1000 Hackett .csv")
# calcHuel(dfout, dfDummy, nsims_real = 1000, nsims_dummy = 10000)


#### CharacterSimmaps fxn ----
CharacterSimmaps <- function(columns, df, tree, dummy, nsims, treelabel, datalabel, dummyMethod = c("simHistory", "makeSimmap")) {
  
  require(stringr)
  source("findQrates.R")
  cooprates <- findQrates(columns = columns[1], newdata = df, newtree = tree)
  coopQ <- cooprates$qrates
  coopQ01 <- coopQ[3]
  coopQ10 <- coopQ[2]
  coopAnc = str_remove(cooprates$ARDlikanc, "ARDlik.anc ")  # added this for sim.history()
  coopAnc = as.numeric(coopAnc)
  names(coopAnc) <- c("0","1")
  FSrates <- findQrates(columns = columns[2], newdata = df, newtree = tree)
  FSQ <- FSrates$qrates
  FSQAbsPres <- FSQ[3]
  FSQPresAbs <- FSQ[2]
  FSAnc = str_remove(FSrates$ARDlikanc, "ARDlik.anc ")  # added this for sim.history()
  FSAnc = as.numeric(FSAnc)
  names(FSAnc) <- c("0","1")
  subsets <- subsettreedata(columns = columns, newdata = df, newtree = tree)
  subsetdf <- subsets$subsetdf
  subsettree <- subsets$subsettree
  FSvec <- subsetdf[,columns[2]]
  names(FSvec) <- subsetdf$species
  Coopvec <- subsetdf[,columns[1]]
  names(Coopvec) <- subsetdf$species
  
  if (dummy == FALSE) {
    FSsimtrees <- make.simmap(tree = subsettree, x = FSvec, model = "ARD", nsim = nsims, Q = FSQ)
    Coopsimtrees <- make.simmap(tree = subsettree, x = Coopvec, model = "ARD", nsim = nsims, Q = coopQ)
    datalabel = paste(datalabel, "REAL")
  } else {
    # dummy simmaps made from one of two methods
    if (dummyMethod == "simHistory") {
      print(paste("starting Dummy Coop simmaps", Sys.time()))
      Coopsimtrees <- sim.history(tree = subsettree, Q = coopQ, nsim = nsims, anc = coopAnc)
      print(paste("starting Dummy FemSong simmaps", Sys.time()))
      FSsimtrees = sim.history(tree = subsettree, Q = FSQ, nsim = nsims, anc = FSAnc)
      datalabel <- paste(datalabel,"DUMMYSimHist-CoopFS")
    } else if (dummyMethod == "makeSimmap") {
      # Make randomized versions of CoopBreed simmaps / DUMMY data
      CoopsimtreesRand <- list()
      for (j in 1:nsims) { 
        Coopvec <- subsetdf[,columns[1]]
        CoopvecRandom <- sample(Coopvec)
        names(CoopvecRandom) <- subsetdf$species
        Coopsimtree <- make.simmap(tree = subsettree, x = CoopvecRandom, model = "ARD", nsim = 1, Q = coopQ)

        CoopsimtreesRand[[j]] <- Coopsimtree
        print(paste(j, Sys.time()))
      } # end for j in 1:nsims (Coop)
      Coopsimtrees <- CoopsimtreesRand
      
      # Make randomized versions of FemaleSong simmaps / DUMMY data
      FSsimtreesRand <- list()
      print(paste("starting Dummy FemSong simmaps", Sys.time()))
      for (j in 1:nsims) {
        FSvec <- subsetdf[,columns[2]]
        FSvecRandom <- sample(FSvec)
        names(FSvecRandom) <- subsetdf$species
        FSsimtree <- make.simmap(tree = subsettree, x = FSvecRandom, model = "ARD", nsim = 1, Q = FSQ)
        FSsimtreesRand[[j]] <- FSsimtree
        print(j)
      } # end for j in 1:nsims (FS)
      FSsimtrees<- FSsimtreesRand
      
      datalabel <- paste(datalabel,"DUMMYResampledMkSimmap-CoopFS")
    } # end if dummyMethod else
    
    
  } # end if dummy = FALSE else
  
  
  column1 <- columns[1]
  column2 <- columns[2]
  Nspecies <- length(subsetdf$species)
  
  
  # Note: output organization:
  # > Map.Overlap(FSsimtree1, Coopsimtree1)
  #          0          1
  # Absent  0.2182494 0.02578955
  # Present 0.6247627 0.13119829
  # --> Map.Overlap(returned in rows, returned in columns)
  
  dfout <- set.seed(10)
  for (i in 1:nsims) {
    
    FSsimtree1 <- FSsimtrees[[i]]
    Coopsimtree1 <- Coopsimtrees[[i]]
    
    FSdescribed <- describe.simmap(FSsimtree1)
    propFSabsent <- FSdescribed$times[2,1]
    propFSpresent <- FSdescribed$times[2,2]
    Coopdescribed <- describe.simmap(Coopsimtree1)
    propNoncoop <- Coopdescribed$times[2,1]
    propCoop <- Coopdescribed$times[2,2]
    totaltime <- FSdescribed$times[1,3]
    
    OverlapMat <- Map.Overlap(FSsimtree1, Coopsimtree1)
    ObsProp0Absent <- OverlapMat[1,1]
    ObsProp1Absent <- OverlapMat[1,2]
    ObsProp0Present <- OverlapMat[2,1]
    ObsProp1Present <- OverlapMat[2,2]
    # ChiMat <- round(OverlapMat*totaltime)
    # chiOutSim <- chisq.test(ChiMat, simulate.p.value = TRUE)
    # chiStatSim <- chiOutSim$statistic
    # chiPvalSim <- chiOutSim$p.value
    # chiOut <- chisq.test(ChiMat, simulate.p.value = FALSE)
    # chiStat <- chiOut$statistic
    # chiPval <- chiOut$p.value
    
    
    treenum <- i
    
    # temprow <- c(treenum, column1, column2, coopQ01, coopQ10, FSQAbsPres, FSQPresAbs, Nspecies, propFSabsent, propFSpresent, propNoncoop, propCoop, ObsProp0Absent, ObsProp0Present, ObsProp1Absent, ObsProp1Present, totaltime, chiStat, chiPval, chiStatSim, chiPvalSim)
    temprow <- c(treenum, column1, column2, coopQ01, coopQ10, FSQAbsPres, FSQPresAbs, Nspecies, propFSabsent, propFSpresent, propNoncoop, propCoop, ObsProp0Absent, ObsProp0Present, ObsProp1Absent, ObsProp1Present, totaltime)
    dfout <- rbind(dfout, temprow)
    dfout <- as.data.frame(dfout)
    colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
    if (i %in% c(100, 200, 250,500,1000,2000,3000,4000,5000)) {
      write.csv(dfout, file = paste(datalabel, "simmap overlap output nsim", nsims, treelabel,".csv"))
    }
  }
  dfout <- as.data.frame(dfout)
  #colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime", "chiStat", "chiPval", "chiStatSim", "chiPvalSim")
  colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
  write.csv(dfout, file = paste("simmap overlap output nsim", nsims, treelabel, datalabel,".csv"))
  return(dfout)
} # end function


#### Cycle families ----
newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data-wAvoNetFamilies_R.csv"
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
dataIn = read.csv(newdata)
dataIn$X = NULL

CBcolumn = "MeanCoopTie2Coop"
columns = c(CBcolumn, "FemaleSong_Agg01")

subset1 <- subsettreedata(columns = columns, newdata = dataIn, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf
subsetdf1[,columns[1]] <- as.character(subsetdf1[,columns[1]])
subsetdf1[,columns[2]] <- as.character(subsetdf1[,columns[2]])

unique(subsetdf1$Order3_BirdtreeMatchSpecies2)
familyvec = unique(subsetdf1$Family3_BirdtreeMatchSpecies2)
familyvec = na.omit(familyvec)

for (i in 1:length(familyvec)) {
  
  familyToRemove = familyvec[i]
  currentlabel = paste0("remove",familyToRemove)
  print(currentlabel)
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2 != familyToRemove),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  
  
  
} # end cycle through families



#### Rest of procedure from Huelsenbeck et al 2003
## 
#dfout <- read.csv("CoopBreed FSOdom simmap overlap output nsim1000 Hackett .csv")
#dfout <- read.csv("CoopBreed FSWebb simmap overlap output nsim1000 Hackett .csv")
# dfout <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg simmap overlap output nsim 1000 Hackett .csv")
# dfDummy <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg DUMMYResampledCoopFS simmap overlap output nsim 500 Hackett .csv")
# dfdummy2 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg DUMMYResampledCoopFS simmap overlap output nsim 100 Hackett .csv")
# dfdummy400 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg DUMMYResampledCoopFS simmap overlap output nsim 400 Hackett .csv")
# dfDummy500b = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg DUMMYResampledCoopFS simmap overlap output nsim 500 Hackett .csv")
# dfDummy3 = rbind(dfDummy, dfdummy2, dfDummy500b, dfdummy400)
# nsims = 100
# calcHuel(dfout,dfDummy3, nsims_real = 1000, nsims_dummy = 1500)




#### calcHuel fxn ----
calcHuel <- function(dfout, dfDummy, nsims_real, nsims_dummy) {
  ExpPropAbsent0 <- as.numeric(dfout$propFSabsent)*as.numeric(dfout$propNoncoop)
  ExpPropAbsent1 <- as.numeric(dfout$propFSabsent)*as.numeric(dfout$propCoop)
  ExpPropPresent0 <- as.numeric(dfout$propFSpresent)*as.numeric(dfout$propNoncoop)
  ExpPropPresent1 <- as.numeric(dfout$propFSpresent)*as.numeric(dfout$propCoop)
  
  # sum(obs-expected) for each state
  dAbsent0 <- sum(abs(as.numeric(dfout$ObsProp0Absent) - ExpPropAbsent0))
  #as.numeric(dfout$ObsProp0Absent) - ExpPropAbsent0
  dAbsent1 <- sum(abs(as.numeric(dfout$ObsProp1Absent) - ExpPropAbsent1))
  #as.numeric(dfout$ObsProp1Absent) - ExpPropAbsent1
  dPresent0 <- sum(abs(as.numeric(dfout$ObsProp0Present) - ExpPropPresent0))
  #as.numeric(dfout$ObsProp0Present) - ExpPropPresent0
  dPresent1 <- sum(abs(as.numeric(dfout$ObsProp1Present) - ExpPropPresent1))
  #as.numeric(dfout$ObsProp1Present) - ExpPropPresent1
 
   sum(dAbsent0, dAbsent1, dPresent0, dPresent1) 
  D_real <- sum(dAbsent0, dAbsent1, dPresent0, dPresent1)/nsims_real
  print(paste("D_real:",D_real))
  
  # do sum-sum in other order to get distribution of Dsims 
  dAbsent0 <- abs(as.numeric(dfout$ObsProp0Absent) - ExpPropAbsent0)
  dAbsent1 <- abs(as.numeric(dfout$ObsProp1Absent) - ExpPropAbsent1)
  dPresent0 <- abs(as.numeric(dfout$ObsProp0Present) - ExpPropPresent0)
  dPresent1 <- abs(as.numeric(dfout$ObsProp1Present) - ExpPropPresent1)
  Real_dsims <- rowSums(cbind(dAbsent0, dAbsent1, dPresent0, dPresent1))
  
  #dfDummy <- read.csv("CoopBreed FSOdom DUMMYResampledCoopFS simmap overlap output nsim1000 Hackett .csv")
  # dfDummy <- read.csv("CoopBreed FSWebb DUMMYResampledCoopFS simmap overlap output nsim1000 Hackett .csv")
  # dfDummy <- read.csv("CoopTie2NonCoop FSWebb DUMMYResampledCoopFS simmap overlap output nsim 100 Hackett .csv")
  ExpPropAbsent0 <- as.numeric(dfDummy$propFSabsent)*as.numeric(dfDummy$propNoncoop)
  ExpPropAbsent1 <- as.numeric(dfDummy$propFSabsent)*as.numeric(dfDummy$propCoop)
  ExpPropPresent0 <- as.numeric(dfDummy$propFSpresent)*as.numeric(dfDummy$propNoncoop)
  ExpPropPresent1 <- as.numeric(dfDummy$propFSpresent)*as.numeric(dfDummy$propCoop)
  dAbsent0 <- abs(as.numeric(dfDummy$ObsProp0Absent) - ExpPropAbsent0)
  dAbsent1 <- abs(as.numeric(dfDummy$ObsProp1Absent) - ExpPropAbsent1)
  dPresent0 <- abs(as.numeric(dfDummy$ObsProp0Present) - ExpPropPresent0)
  dPresent1 <- abs(as.numeric(dfDummy$ObsProp1Present) - ExpPropPresent1)
  Dummy_dsums <- rowSums(cbind(dAbsent0, dAbsent1, dPresent0, dPresent1))
  
  pval <- sum(Dummy_dsums > D_real)/nsims_dummy
  print(paste("num Dummy dsums > D_real:", sum(Dummy_dsums > D_real)))
  print(paste("pval:",pval))
  
  xmax = max(c(Real_dsims, Dummy_dsums))*1.1
  
  par(mfrow=c(2,1))
  hist(Real_dsims, xlim = c(0,xmax), breaks = 20)
  abline(v = D_real, col = "red")
  hist(Dummy_dsums, xlim = c(0,xmax), breaks = 20)
  
  realDsimDF = cbind(c(rep("real", length(Real_dsims))), Real_dsims)
  dummyDsimDF = cbind(c(rep("dummy", length(Real_dsims))), Real_dsims)
  
}


### Plot hists better
font_size <- 1.25

# Setting layout for 2 plots
par(mfrow=c(2,1), mai=c(0.8,1.0,0.5,0.3), oma=c(2,3,2,3), 
    font.main=1, cex.main=1.5, cex.lab=font_size, cex.axis=1.2)

# Histogram for Real_dsims
hist(Real_dsims, xlim=c(0,xmax), breaks=20, 
     col=rgb(0.2,0.5,0.7,0.5), # semi-transparent blue color
     main="", xlab="D statistic from real data simmaps", ylab="Frequency",
     border="white")
abline(v = D_real, col="red", lwd=2.5) # thicker red line

# Histogram for Dummy_dsums
hist(Dummy_dsums, xlim=c(0,xmax), breaks=20, 
     col=rgb(0.7,0.5,0.2,0.5), # semi-transparent orange color
     main="", xlab="D statistic from simulated independent data simmaps", ylab="Frequency",
     border="white")

dev.off()


### Plot trees and pval hist
pdf(file = paste(datalabel, "simmaps overlap", treelabel, ".pdf"), width = 10, height = 10)
par(mfcol = c(2,3))
par(mar = c(4,1,4,1))
numsig <- sum(as.numeric(dfout$chiPval) < 0.05)
hist(as.numeric(dfout$chiPval), breaks = 40, main = paste("chisqPval distribution, numsig:",numsig,"/1000"))
numsig <- sum(as.numeric(dfout$chiPvalSim) < 0.05)
hist(as.numeric(dfout$chiPvalSim), breaks = 40, main = paste("chisq SimPval distribution, numsig:",numsig,"/1000"))
densityMap(FSsimtrees, fsize = 0.1, lwd = 0.7)
densityMap(Coopsimtrees, fsize = 0.1, lwd = 0.7)
#dfout <- read.csv("OdomFS vs CoopBreed simmap overlap output nsim1000 HackettConsensus.csv")
nonsigvec <- which(dfout$chiPval > 0.05)
sigvec <- which(dfout$chiPval < 0.01)
treevec <- c(sigvec[1:15], nonsigvec[1:15])
for (i in treevec) {
  plotSimmap(FSsimtrees[[i]],fsize = 0.2, lwd = 0.7)
  chiPval <- dfout[which(dfout$treenum == i),"chiPval"]
  if (chiPval < 0.0001) {
    Pvaltext <- "< 0.0001"
  } else if (chiPval < 0.001) {
    Pvaltext <- "< 0.001"
  } else if (chiPval < 0.01) {
    Pvaltext <- "< 0.01"
  } else if (chiPval < 0.05) {
    Pvaltext <- "< 0.05"
  } else {
    Pvaltext <- "> 0.05"
  }
  par(cex.main = 0.9)
  title(main = paste("\ntree num", i, "chiPval", Pvaltext))
  cols <- c("black","orange")
  names(cols) <- c(0,1)
  plotSimmap(Coopsimtrees[[i]], fsize = 0.2, colors = cols, lwd = 0.7)
  title(main = paste("\ntree num", i,"Female song present = red, CoopBreed = orange"))
}
dev.off()


# densityMap to visualize
##densityMap(Coopsimtrees, fsize = 0.1, lwd = 0.7)
# randomize order of coopbreed 1000 times to make 1000 simmaps, make sure pval is nonsig
Coopsimtrees <- list()
for (j in 1:1000) {
  Coopvec <- subsetdf[,columns[1]]
  CoopvecRandom <- sample(Coopvec)
  names(CoopvecRandom) <- subsetdf$species
  Coopsimtree <- make.simmap(tree = subsettree, x = Coopvec, model = "ARD", nsim = 1, Q = coopQ)
  Coopsimtrees[[i]] <- Coopsimtree
}


### Check how state proportions correlate with pvals
dfDummy <- read.csv("CoopBreed FSWebb DUMMYcoop simmap overlap output nsim1000 Hackett .csv")
dfout <- read.csv("CoopBreed FSWebb simmap overlap output nsim1000 Hackett .csv")
logpval <- log(dfDummy$chiPval,10)
hist(logpval,40)
hist(dfDummy$chiPval)

logpval <- log(dfout$chiPval,10)
hist(logpval, 40)
hist(dfout$chiPval)


pvalgroup <- set.seed(10)
pvalgroup[which(dfDummy$chiPval < 0.05)] <- "pSignificant0.05"
pvalgroup[which(dfDummy$chiPval < 0.01)] <- "Significant0.01"
pvalgroup[which(dfDummy$chiPval < 0.001)] <- "VerySignificant0.001"
pvalgroup[which(dfDummy$chiPval > 0.05)] <- "NotSignificant"
dfDummy <- cbind(dfDummy, pvalgroup)

pvalgroup <- set.seed(10)
pvalgroup[which(dfout$chiPval < 0.05)] <- "pSignificant0.05"
pvalgroup[which(dfout$chiPval < 0.01)] <- "Significant0.01"
pvalgroup[which(dfout$chiPval < 0.001)] <- "VerySignificant0.001"
pvalgroup[which(dfout$chiPval > 0.05)] <- "NotSignificant"
dfout <- cbind(dfout, pvalgroup)

library(ggplot2)
library(cowplot)
library(reshape2)
dfDummyMelt <- melt(dfDummy, measure.vars = c("ObsProp0Absent","ObsProp1Absent","ObsProp0Present","ObsProp1Present"), variable.name = "state", value.name = "state.prop")
dfDummyMelt <- melt(dfDummyMelt, measure.vars = c("propFSabsent", "propFSpresent"), variable.name = "FS_state", value.name = "FSstate.prop")
dfDummyMelt <- melt(dfDummyMelt, measure.vars = c("propNoncoop", "propCoop"), variable.name = "CB_state", value.name = "CBstate.prop")

a <- ggplot(dfDummyMelt, aes(x = state, y = state.prop)) + geom_boxplot() + ggtitle("dummy CoopBreed data")
b <- ggplot(dfDummyMelt, aes(x = FS_state, y = FSstate.prop)) + geom_boxplot()+ ggtitle("dummy CoopBreed data")
c <- ggplot(dfDummyMelt, aes(x = CB_state, y = CBstate.prop)) + geom_boxplot()+ ggtitle("dummy CoopBreed data")


dfoutMelt <- melt(dfout, measure.vars = c("ObsProp0Absent","ObsProp1Absent","ObsProp0Present","ObsProp1Present"), variable.name = "state", value.name = "state.prop")
dfoutMelt <- melt(dfoutMelt, measure.vars = c("propFSabsent", "propFSpresent"), variable.name = "FS_state", value.name = "FSstate.prop")
dfoutMelt <- melt(dfoutMelt, measure.vars = c("propNoncoop", "propCoop"), variable.name = "CB_state", value.name = "CBstate.prop")

d <- ggplot(dfoutMelt, aes(x = state, y = state.prop)) + geom_boxplot() + ggtitle("real CoopBreed data")
e <- ggplot(dfoutMelt, aes(x = FS_state, y = FSstate.prop)) + geom_boxplot() + ggtitle("real CoopBreed data")
f <- ggplot(dfoutMelt, aes(x = CB_state, y = CBstate.prop)) + geom_boxplot() + ggtitle("real CoopBreed data")

require(cowplot)
plot_grid(a,b,c,d,e,f, ncol = 3)
