# Adapted from test_trait_overlap_simmaps.R by Kate Snyder
# Created for server by Kate Snyder
# 6/22/2023
# Last edited: 7/28/2023 (pasted onto this doc bc )


library(phytools)
library(ape)
source("subsettreedata.R")

list.files()

newdata = "2023-06-20_CoopBreed-FemaleSong01-Song_Data_R.csv"
newdata = "2023-06-20_CoopBreed-FemaleSong01HighConf-Song_Data_R.csv"
dfAgg = read.csv(newdata)
df = dfAgg

Hacktree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000.nex")
treelabel <- "Hackett"

#df[,"FemaleSong_Aggregated"][which(df[,"FemaleSong_Aggregated"] == "Absent")] <- 0
#df[,"FemaleSong_Aggregated"][which(df[,"FemaleSong_Aggregated"] == "Present")] <- 1

nsims_real = 1000
nsims_dummy = 5000

datalabel = "Tie2Noncoop FSHighConf"
dfout5 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy5 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout5, dfDummy5, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "Tie2Coop FSHighConf"
dfout4 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
dummydf = read.csv("Tie2Coop FSHighConf DUMMYResampledCoopFS simmap overlap output nsim 5000 Hackett .csv")
dummydf = dummydf[,2:18]
dfout = read.csv("Tie2Coop FSHighConf simmap overlap output nsim 1000 Hackett .csv")
dfout = dfout[,2:18]
calcHuel(dfout, dummydf, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "CoopOmitTies FSHighConf"
dfout4 <- CharacterSimmaps(columns = c("MeanCoopOmitTies","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopOmitTies","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout4, dfDummy4, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "CoopOmitTies FSAgg"
dfout4 <- CharacterSimmaps(columns = c("MeanCoopOmitTies","FemaleSong_Agg01"), df = newdata, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopOmitTies","FemaleSong_Agg01"), df = newdata, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout4, dfDummy4, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "AnyCoop FSHighConf"
dfout4 <- CharacterSimmaps(columns = c("AnyCoopEqualsCoop","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("AnyCoopEqualsCoop","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout4, dfDummy4, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "AnyCoop FSAgg"
dfout4 <- CharacterSimmaps(columns = c("AnyCoopEqualsCoop","FemaleSong_Agg01"), df = newdata, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("AnyCoopEqualsCoop","FemaleSong_Agg01"), df = newdata, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout4, dfDummy4, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "KinNK FSAgg"
df$Kin_NK[which(df$Kin_NK == "Mixed")] <- NA
df$Kin_NK[which(df$Kin_NK == "Kin")] <- 1
df$Kin_NK[which(df$Kin_NK == "NonKin")] <- 0
dfout4 <- CharacterSimmaps(columns = c("Kin_NK","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("Kin_NK","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout4, dfDummy4, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "KinNK FSHighConf"
df$Kin_NK[which(df$Kin_NK == "Mixed")] <- NA
df$Kin_NK[which(df$Kin_NK == "Kin")] <- 1
df$Kin_NK[which(df$Kin_NK == "NonKin")] <- 0
dfout4 <- CharacterSimmaps(columns = c("Kin_NK","HighConfidence_FemaleSong"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("Kin_NK","HighConfidence_FemaleSong"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout4, dfDummy4, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "Familial FSHighConf"
dfout4 <- CharacterSimmaps(columns = c("Griesser2017FamilialLiving","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("Griesser2017FamilialLiving","HighConfidence_FemaleSong"), df = newdata, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout4, dfDummy4, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)

datalabel = "Familial FSAgg"
dfout4 <- CharacterSimmaps(columns = c("Griesser2017FamilialLiving","FemaleSong_Agg01"), df = newdata, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = datalabel)
dfDummy4 <- CharacterSimmaps(columns = c("Griesser2017FamilialLiving","FemaleSong_Agg01"), df = newdata, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = datalabel)
calcHuel(dfout, dfDummy, nsims_real = nsims_real, nsims_dummy = nsims_dummy, datalabel = datalabel)



CharacterSimmaps <- function(columns, df, tree, dummy, nsims, treelabel, datalabel) {
  
  source("findQrates.R")
  cooprates <- findQrates(columns = columns[1], newdata = df, newtree = tree)
  coopQ <- cooprates$qrates
  coopQ01 <- coopQ[3]
  coopQ10 <- coopQ[2]
  FSrates <- findQrates(columns = columns[2], newdata = df, newtree = tree)
  FSQ <- FSrates$qrates
  FSQAbsPres <- FSQ[3]
  FSQPresAbs <- FSQ[2]
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
  } else {
    # Make randomized versions of CoopBreed simmaps / DUMMY data
    CoopsimtreesRand <- list()
    print(paste("starting Dummy Coop simmaps", Sys.time()))
    for (j in 1:nsims) {
      Coopvec <- subsetdf[,columns[1]]
      CoopvecRandom <- sample(Coopvec)
      #CoopvecRandom <- sample(c(0,1), length(Coopvec), replace = TRUE)
      names(CoopvecRandom) <- subsetdf$species
      Coopsimtree <- make.simmap(tree = subsettree, x = CoopvecRandom, model = "ARD", nsim = 1, Q = coopQ)
      CoopsimtreesRand[[j]] <- Coopsimtree
      print(j)
    }
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
      print(paste(j, Sys.time()))
    }
    
    FSsimtrees<- FSsimtreesRand
    datalabel <- paste(datalabel,"DUMMYResampled")
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
    
    treenum <- i
    
    temprow <- c(treenum, column1, column2, coopQ01, coopQ10, FSQAbsPres, FSQPresAbs, Nspecies, propFSabsent, propFSpresent, propNoncoop, propCoop, ObsProp0Absent, ObsProp0Present, ObsProp1Absent, ObsProp1Present, totaltime)
    dfout <- rbind(dfout, temprow)
    dfout <- as.data.frame(dfout)
    colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
    if (i %in% c(100,250,500,1000,2000,3000,4000,5000)) {
      write.csv(dfout, file = paste(datalabel, "simmap overlap output nsim", nsims, treelabel,".csv"))
    }
  }
  dfout <- as.data.frame(dfout)
  
  colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
  write.csv(dfout, file = paste(datalabel, "simmap overlap output nsim", nsims, treelabel,".csv"))
  return(dfout)
} # end function



calcHuel <- function(dfout, dfDummy, nsims_real, nsims_dummy, datalabel) {
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
  
  pdf(file = paste0("Simmap overlap ", datalabel, " ", nsims_real," ", nsims_dummy, ".pdf"), height = 8, width = 5)
  par(mfrow=c(2,1))
  hist(Real_dsims, xlim = c(0,0.1), breaks = 20, main = paste("Real_dsims", datalabel, treelabel, "N =", nsims_real))
  abline(v = D_real, col = "red")
  hist(Dummy_dsums, xlim = c(0,0.1), breaks = 20, main = paste("Dummy_dsums pval:", pval, "N =", nsims_dummy))
  dev.off()
  
  realDsimDF = cbind(c(rep("real", length(Real_dsims))), Real_dsims)
  dummyDsimDF = cbind(c(rep("dummy", length(Dummy_dsums))), Dummy_dsums)
  
}
