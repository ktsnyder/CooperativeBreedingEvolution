## Utilize phytools::Map.Overlap() to get proportion of edges spent in each FS/Coop state
## Kate Snyder
## 4/4/2022
## Method based on Huelsenbeck et al (2003)
## 
## Edited 6/1/23 - added checkpoint save to CharacterSimmaps
## Edited 9/27/23 - changed "Dummy" simmap generation to use sim.history() with Q rates, ancestral character estimation instead of randomizing tip states; but seems to have gotten totally weird - output values odd
## Edited 9/29/23
## 10/30/2023 - added boxplot to calcHuel function; 3 ggplots now returned from calcHuel
## 11/22/2023 - added transition counts by state to CharacterSimmaps output
## 11/27/2023 - added boxplots of transition counts (real vs dummy, real vs expected) to calcHuel function, Split/Archived from test_trait_overlap_simmaps.R

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


# df$Kin_NK[which(df$Kin_NK == "Mixed")] <- NA
# df$Kin_NK[which(df$Kin_NK == "NonKin")] <- 0
# df$Kin_NK[which(df$Kin_NK == "Kin")] <- 1
# dfout4 <- CharacterSimmaps(columns = c("Kin_NK","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Kin_NK FSAgg")
# dfDummy4 <- CharacterSimmaps(columns = c("Kin_NK","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Kin_NK FSAgg")
# KinFSreal = read.csv("Kin_NK FSAgg simmap overlap output nsim 1000 Hackett .csv")
# KinFSdummy = read.csv("Kin_NK FSAgg DUMMYResampledCoopFS simmap overlap output nsim 1000 Hackett .csv")
# nsims_real = length(KinFSreal$treenum)
# nsims_dummy = length(KinFSdummy$treenum)
# calcHuel(KinFSreal, KinFSdummy, nsims_real = nsims_real, nsims_dummy = nsims_dummy)

# dfout4 <- CharacterSimmaps(columns = c("Griesser2017FamilialLiving","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = "Familial FSAgg")
# dfDummy4 <- CharacterSimmaps(columns = c("Griesser2017FamilialLiving","FemaleSong_Agg01"), df = df, tree =  Hacktree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = "Familial FSAgg")
# FamFSreal = read.csv("Familial FSAgg simmap overlap output nsim 1000 Hackett .csv")
# FamFSdummy = read.csv("Familial FSAgg DUMMYResampledCoopFS simmap overlap output nsim 5000 Hackett .csv")
# nsims_real = length(FamFSreal$treenum)
# nsims_dummy = length(FamFSdummy$treenum)
# calcHuel(FamFSreal, FamFSdummy, nsims_real = nsims_real, nsims_dummy = nsims_dummy)

# 10/25/2023
nsims_real = 500
nsims_dummy = 2000
OscineTree = read.nexus("/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data-wAvoNetFamilies_R.csv"

dfout4 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","HighConfidence_FemaleSong"), df = newdata, tree =  OscineTree, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL)
dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","HighConfidence_FemaleSong"), df = newdata, tree =  OscineTree, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
calcHuel(dfout4, dfDummy4)

# 10/26/2023
nsims_real = 500
nsims_dummy = 2000
newdata = "2023-10-26_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
OscineTree = read.nexus("/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
otherlabel = "HackettOscine"
dataIn = read.csv(newdata)
dataIn$X = NULL

CBcolumns = c("Griesser2023.Colonial01", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LongSocialBonds", "Griesser2023.MoreThanTwoCaretakers", "Griesser2017FamilialLiving")
FScolumn = "FemaleSong_Agg01"
filelist = list.files("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs", full.names = T)

for (i in 1:5) {
  CBcolumn = CBcolumns[i]
  #RealDF = CharacterSimmaps(columns = c(CBcolumn,FScolumn), df = newdata, tree =  OscineTree, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL)
  #DummyDF = CharacterSimmaps(columns = c(CBcolumn,FScolumn), df = newdata, tree =  OscineTree, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
  
  RealFile <- filelist[which(str_detect(filelist, pattern = CBcolumn) & str_detect(filelist, pattern = FScolumn) & str_detect(filelist, pattern = "REAL"))]
  DummyFile <- filelist[which(str_detect(filelist, pattern = CBcolumn) & str_detect(filelist, pattern = FScolumn) & str_detect(filelist, pattern = "DUMMY"))]
  RealDF = read.csv(RealFile)
  DummyDF = read.csv(DummyFile)

  HuelOut = calcHuel(dfout = RealDF, dfDummy = DummyDF, newplot = FALSE, otherlabel = otherlabel, plot_ggplots_pdf = TRUE)
  
  # pdf(paste("Simmap Overlap Outputs/SimmapOverlap", CBcolumn, FScolumn, "HackettOscine.pdf"), width = 6, height = 14)
  # require(cowplot)
  # print(plot_grid(HuelOut[[1]], HuelOut[[2]], HuelOut[[3]], ncol = 1))
  
  dev.off()
}


# Get Q rates
source(findQrates.R)
for (i in 1:5) {
  columns = c(CBcolumns[i], "FemaleSong_Agg01")
  subset = subsettreedata(columns = columns, newdata = dataIn, newtree = OscineTree)
  subsettree = subset$subsettree
  subsetdf = subset$subsetdf
  #Qout <- findQrates(columns, plot = T, newtree = OscineTree, newdata = dataIn)
  whichnodes = "no"  #"FS" # "Coop" # "no"
  tipsize = 0.1
  filename = paste(columns[1], columns[2], "fan phylo double tips", whichnodes, "nodes.pdf")
  pdf(filename, height = 8, width = 9)
  plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
  py = c("black","orange")
  treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,columns[1]] == 1)]
  tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=tipsize, offset = 1)
  py2 = c("blue","red")
  treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,columns[2]] == 1)]
  tiplabels(pch=21,bg=py2[as.numeric(treetiplabels2)+1], col = py2[as.numeric(treetiplabels2)+1], cex=tipsize, offset = 2)
  if (whichnodes == "Coop") {
    discretetraitvec <- subsetdf[,columns[1]]
    names(discretetraitvec) <- subsetdf[,1]
    pynodes = py
    circles=ace(x=discretetraitvec,phy=subsettree,type="discrete",model="ARD")
    nodelabels(thermo=circles$lik.anc,piecol=pynodes, height = 1.2, width = 1.2, horiz = TRUE, frame = "circle")
  } else if (whichnodes == "FS") {
    discretetraitvec <- subsetdf[,columns[2]]
    names(discretetraitvec) <- subsetdf[,1]
    pynodes = py2
    circles=ace(x=discretetraitvec,phy=subsettree,type="discrete",model="ARD")
    nodelabels(thermo=circles$lik.anc,piecol=pynodes, height = 1.2, width = 1.2, horiz = TRUE, frame = "circle")
  } 
  allpy = c(py, py2)
  alllabs = c(paste("Not", columns[1]), columns[1], "Female Song Absent", "Female Song Present")
  legend("bottomleft", legend = alllabs, cex = 0.9, fill=allpy, bty="n")
  dev.off()
  
  print(columns)
}



#### CharacterSimmaps fxn ----
CharacterSimmaps <- function(columns, df, tree, dummy, nsims, treelabel, datalabel = NULL, dummyMethod = c("simHistory", "makeSimmap")) {
  
  if (is.null(datalabel)) {
    datalabel = paste(columns[1], columns[2])
  }
  
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
    if (!dir.exists("Simmap Overlap Outputs")) {
      dir.create("Simmap Overlap Outputs")
    }
    if (i %in% c(100, 200, 250,500,1000,2000,3000,4000,5000)) {
      write.csv(dfout, file = paste("Simmap Overlap Outputs/", datalabel, "simmap overlap_counts output nsim", nsims, treelabel,".csv"), row.names = FALSE)
    }
  }
  dfout <- as.data.frame(dfout)
  #colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime", "chiStat", "chiPval", "chiStatSim", "chiPvalSim")
  colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
  
  source("find transition counts by state for 2 Discrete traits.R")
  transStateCounts = getTransitionStateCounts(Coopsimtrees = Coopsimtrees, FSsimtrees = FSsimtrees)
  dfout = merge(dfout, transStateCounts, by.x = "treenum", by.y = "TreeNum", )
  
  write.csv(dfout, file = paste("Simmap Overlap Outputs/", datalabel, "simmap overlap_counts output nsim", nsims, treelabel,".csv"), row.names = FALSE)
  return(dfout)
} # end function


#### Cycle families - jackknife ----
newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data-wAvoNetFamilies_R.csv"
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
Hacktree = read.nexus(treefile)
dataIn = read.csv(newdata)
dataIn$X = NULL

CBcolumn = "MeanCoopTie2Coop"
#columns = c(CBcolumn, "FemaleSong_Agg01")
columns = c(CBcolumn, "FemaleSong_Agg01")

subset1 <- subsettreedata(columns = columns, newdata = dataIn, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf
subsetdf1[,columns[1]] <- as.character(subsetdf1[,columns[1]])
subsetdf1[,columns[2]] <- as.character(subsetdf1[,columns[2]])

familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2[which(familycounts$n > 2)]

nsims_real = 50
nsims_dummy = 200


for (j in 1:length(familyvec)) {
  
  familyToRemove = familyvec[j]
  templabel = paste0(columns[1], "_", columns[2], "_", "remove",familyToRemove)
  print(paste(j, templabel))
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2 != familyToRemove),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  dfout4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = templabel, dummyMethod = "makeSimmap")
  dfDummy4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = templabel, dummyMethod = "makeSimmap")
  
  
} # end cycle through families for jackknife


#### Cycle families - single family ----
newdata = "2023-09-14_CoopBreed-FemaleSong-Song_Data-wAvoNetFamilies_R.csv"
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
Hacktree = read.nexus(treefile)
dataIn = read.csv(newdata)
dataIn$X = NULL

CBcolumn = "MeanCoopTie2Noncoop"
columns = c(CBcolumn, "FemaleSong_Agg01")

subset1 <- subsettreedata(columns = columns, newdata = dataIn, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf
subsetdf1[,columns[1]] <- as.character(subsetdf1[,columns[1]])
subsetdf1[,columns[2]] <- as.character(subsetdf1[,columns[2]])

familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2[which(familycounts$n > 4)]

nsims_real = 100
nsims_dummy = 500


for (j in 1:length(familyvec)) {
  
  familyToKeep = familyvec[j]
  templabel = paste0(columns[1], "_", columns[2], "_", "only",familyToKeep)
  print(paste(j, templabel))
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2 == familyToKeep),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  if ( length(unique(subsetdf[,columns[1]])) == 2 & length(unique(subsetdf[,columns[2]])) == 2 ) {
    print(paste(familyToKeep, "has both FS/noFS and CB/nonCB"))
  dfout4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = FALSE, nsims = nsims_real, treelabel = "Hackett", datalabel = templabel, dummyMethod = "makeSimmap")
  dfDummy4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = TRUE, nsims = nsims_dummy, treelabel = "Hackett", datalabel = templabel, dummyMethod = "makeSimmap")
  } # end if family has instances of CB = 0,1 and FS = 0,1
  
} # end cycle through families for single family


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
calcHuel <- function(dfout, dfDummy, nsims_real = NULL, nsims_dummy = NULL, otherlabel = NULL, newplot = TRUE, plot_ggplots_pdf = FALSE) {
  require(tidyverse)
  require(ggplot2)
  if (is.null(nsims_real)) {
    nsims_real = length(dfout[,1])
  }
  if (is.null(nsims_dummy)) {
    nsims_dummy = length(dfDummy[,1])
  }
  
  Nspecies = dfout[1,"Nspecies"]
  trait1 = dfout[1,"column1"]
  trait2 = dfout[1,"column2"]
  
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
  
  plotlabel = paste(trait1, trait2, "\nN species =", Nspecies, otherlabel)
  dummytitle = paste("Nsims =", nsims_dummy, otherlabel, "\nnum Dummy dsums > D_real:", sum(Dummy_dsums > D_real), ", pval =", pval)
  
  xmax = max(c(Real_dsims, Dummy_dsums))*1.1
  
  # par(mfrow=c(2,1))
  # hist(Real_dsims, xlim = c(0,xmax), breaks = 20)
  # abline(v = D_real, col = "red")
  # hist(Dummy_dsums, xlim = c(0,xmax), breaks = 20)
  
  
  ### Plot hists better
  
if (newplot == TRUE) {
  font_size <- 1.0
  # Setting layout for 2 plots
  par(mfrow=c(2,1), mai=c(0.8,1.0,0.5,0.3), oma=c(2,3,2,3), 
      font.main=1, cex.main=1.25, cex.lab=font_size, cex.axis=1)
}
  
  # Histogram for Real_dsims
  hist(Real_dsims, xlim=c(0,xmax), breaks=20, 
       col=rgb(0.2,0.5,0.7,0.5), # semi-transparent blue color
       main=plotlabel, xlab="D statistic from real data simmaps", ylab="Frequency",
       border="white")
  abline(v = D_real, col="red", lwd=2.5) # thicker red line
  
  # Histogram for Dummy_dsums
  hist(Dummy_dsums, xlim=c(0,xmax), breaks=20, 
       col=rgb(0.7,0.5,0.2,0.5), # semi-transparent orange color
       main=dummytitle, xlab="D statistic from simulated independent data simmaps", ylab="Frequency",
       border="white")
  
  real_data <- data.frame(value = Real_dsims)
  dummy_data <- data.frame(value = Dummy_dsums)
  
  # ggplot histograms to return
  p1 <- ggplot(real_data, aes(x = value)) +
    geom_histogram(binwidth = xmax/20, fill = rgb(0.2, 0.5, 0.7, 0.5), color = "white") +
    geom_vline(xintercept = D_real, color = "red") +
    xlim(c(0, xmax)) +
    labs(title = plotlabel, x = "D statistic from real data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)

  # Create the histogram for Dummy_dsums
  p2 <- ggplot(dummy_data, aes(x = value)) +
    geom_histogram(binwidth = xmax/20, fill = rgb(0.7, 0.5, 0.2, 0.5), color = "white") +
    xlim(c(0, xmax)) +
    labs(title = dummytitle, x = "D statistic from simulated independent data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)
  
  # Create mutated dataframes for boxplot
  dfRealMelt <- dfout %>%
    gather("ObservedState", "ObservedState.prop", ObsProp0Absent:ObsProp1Present) %>%
    mutate(Which = "Real")
  
  dfDummyMelt <- dfDummy %>%
    gather("ObservedState", "ObservedState.prop", ObsProp0Absent:ObsProp1Present) %>%
    mutate(Which = "Dummy")
  
  dfCombined <- rbind(dfRealMelt, dfDummyMelt)
  
  dfCombined$ObservedState.prop = as.numeric(dfCombined$ObservedState.prop)
  
  # Create a combined label for x-axis
  dfCombined <- dfCombined %>%
    mutate(Label = paste(ObservedState, Which, sep = "\n"))
  
  # Create the boxplot
  p3 <-  ggplot(dfCombined, aes(x = Label, y = ObservedState.prop, fill = Which)) +
    geom_boxplot(outlier.shape = NA) + # Exclude outliers
    theme_minimal() +
    labs(y = "Observed State Proportion", x = "", fill = "Simulation Data") +
    scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste(trait1, trait2))
  
  if ("Coop0to1inFS0" %in% colnames(dfout) & "Coop0to1inFS0" %in% colnames(dfDummy)) { # compared rates are from Dummy data simmaps
    dfRealMeltCounts <- dfout %>%
      gather("Transition", "Count", Coop0to1inFS0:FS1to0inCoop1) %>%
      mutate(Which = "Real")
    
    dfDummyMeltCounts <- dfDummy %>%
      gather("Transition", "Count", Coop0to1inFS0:FS1to0inCoop1) %>%
      mutate(Which = "Dummy")
    
    dfCombinedCounts <- rbind(dfRealMeltCounts, dfDummyMeltCounts)
    
    dfCombinedCounts$Count = as.numeric(dfCombinedCounts$Count)
    
    # Create a combined label for x-axis
    dfCombinedCounts <- dfCombinedCounts %>%
      mutate(Label = paste(Transition, Which, sep = "\n"))
    
    p4 <-  ggplot(dfCombinedCounts, aes(x = Label, y = Count, fill = Which)) +
      geom_boxplot(outlier.shape = NA) + # Exclude outliers
      theme_minimal() +
      labs(y = "Transition Counts", x = "", fill = "Simulation Data") +
      scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste(trait1, trait2))
    
    if ("Coop0to1inFS0Expected" %in% colnames(dfout)) { # compared rates are calculated "Expected" values from e.g. (Total # transitions trait1=0 to trait1=1) * (total time in trait2 = 0)
      dfMeltCounts <- dfout %>%
        gather("Transition", "Count", c(Coop0to1inFS0:FS1to0inCoop1, Coop0to1inFS0Expected:FS1to0inCoop1Expected))
      dfMeltCounts$ObservedVsExpected = "Observed"
      dfMeltCounts$ObservedVsExpected[which(str_detect(dfMeltCounts$Transition, "Expected"))] = "Expected"
      dfMeltCounts$Count = as.numeric(dfMeltCounts$Count)
      
      dfMeltCounts$Label = dfMeltCounts$Transition
      
      p5 <-  ggplot(dfMeltCounts, aes(x = Label, y = Count, fill = ObservedVsExpected)) +
        geom_boxplot(outlier.shape = NA) + # Exclude outliers
        theme_minimal() +
        labs(y = "Transition Counts", x = "", fill = "Observed/Expected") +
        scale_fill_manual(values = c("Observed" = "blue", "Expected" = "red")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(paste(trait1, trait2))
      
      return(list(p1, p2, p3, p4, p5))
      
    } else {
      return(list(p1, p2, p3, p4))
    }
    
  } else {
    p4 = NULL
    return(list(p1, p2, p3))
  }
  
  
  if (plot_ggplots_pdf == TRUE) {
    require(cowplot)
    pdf(file = paste("Simmap Overlap Outputs/Simmap Overlap",trait1, trait2, nsims_real, nsims_dummy, otherlabel, ".pdf"), height = 14, width = 6)
    print(plot_grid(p1, p2, p3, ncol = 1))
    dev.off()
  }
  #require(cowplot)
  #print(plot_grid(p1, p2, p3, ncol = 1))
  
  return(list(p1, p2, p3, p4))
  
}


#### cycle calcHuel ---- 
# After doing the Cycle families processes above
require(reshape2)
trait1 = "MeanCoopTie2Noncoop"
trait2 = "FemaleSong_Agg01"

## EITHER
filelist = list.files(paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/Jackknife ", trait1, " ", trait2), recursive = T, pattern = ".csv", full.names = T)
tempFiles = jackknifeFiles = filelist[which(str_detect(filelist, "remove") & str_detect(filelist, trait1) & str_detect(filelist, trait2))]
tempUniqueFamilies = unique_remove_strings <- unique(gsub(".*remove([^ ]*) .*", "\\1", jackknifeFiles))
familytreatment = "removed"
pdf(file = paste0("jackknifed Simmap Overlaps ", trait1, " ", trait2, ".pdf"), width = 12, height = 10)

## OR
filelist = list.files(paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/Single Family ", trait1, " ", trait2), recursive = T, pattern = ".csv", full.names = T)
tempFiles = onefamilyFiles = filelist[which(str_detect(filelist, "only") & str_detect(filelist, trait1) & str_detect(filelist, trait2))]
tempUniqueFamilies = unique_only_strings <- unique(gsub(".*only([^ ]*) .*", "\\1", onefamilyFiles))
familytreatment = "only"
pdf(paste0("single family Simmap Overlaps ", trait1, " ", trait2, ".pdf"), width = 12, height = 10)


par(mfcol=c(3,3), mai=c(0.8,1.0,0.5,0.3), oma=c(2,3,2,3), 
    font.main=1, cex.main=1.25, cex.lab=1.1, cex.axis=1.1)

for (i in 1:length(tempUniqueFamilies)) {
  tempFamily = tempUniqueFamilies[i]
  tempRealFile = tempFiles[which(str_detect(tempFiles, tempFamily) & str_detect(tempFiles, "REAL"))]
  tempDummyFile = tempFiles[which(str_detect(tempFiles, tempFamily) & str_detect(tempFiles, "DUMMY"))]
  
  # RealDF = read.csv(paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/",tempRealFile))
  # DummyDF = read.csv(paste0("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/",tempDummyFile))
  RealDF = read.csv(tempRealFile)
  DummyDF = read.csv(tempDummyFile)
  
  templabel = paste0(familytreatment, tempFamily)
  
  calcHuel(dfout = RealDF, dfDummy = DummyDF, newplot = FALSE, otherlabel = templabel)
  
  # Third plot: observed state proportions
  dfDummyMelt <- melt(RealDF, measure.vars = c("ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present"), variable.name = "ObservedState", value.name = "ObservedState.prop")
  
  # Create the basic boxplot without x-axis labels (xaxt = "n") and without outlines
  boxplot(ObservedState.prop ~ ObservedState, data = dfDummyMelt, xlab = "", ylab = "Observed State Proportion", xaxt = "n", outline = FALSE, boxwex = 0.5, col = "lightgray")
  
  # Add points to the plot
  points(jitter(as.numeric(dfDummyMelt$ObservedState), amount = 0.05), dfDummyMelt$ObservedState.prop, pch = 16, col = "darkred", cex = 0.6)
  
  # Add rotated x-axis labels
  axis(1, at = 1:length(unique(dfDummyMelt$ObservedState)), 
       labels = FALSE)
  text(1:length(unique(dfDummyMelt$ObservedState)), 
       par("usr")[3] - 0.001, srt = 45, adj = 1.2, 
       labels = as.character(unique(dfDummyMelt$ObservedState)), xpd = TRUE, cex=0.8)
}
dev.off()


#### Misc ----

### Plot hists better - moved into calcHuel
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
