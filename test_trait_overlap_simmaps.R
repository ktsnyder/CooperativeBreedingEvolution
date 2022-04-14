## Utilize phytools::Map.Overlap() to get proportion of edges spent in each FS/Coop state
## Kate Snyder
## 4/4/2022
## Method based on Huelsenbeck et al (2003)

library(phytools)
source("subsettreedata.R")
df <- read.csv("2022-03-10CoopSong_All.csv")
Hacktree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000.nex")
treelabel <- "Hackett"
Erictree <-read.nexus("2021-08-31ConsensusPasserineTreeEricson10_1000.nex")
treelabel <- "Ericson"

WebbFSdf <- read.csv("Webb et al 2016 Female Song Plumage Data.csv")
dfnew <- merge(df, WebbFSdf, by.x = "species", by.y = "TipLabel")
columns <- c("MeanCoopTie2Noncoop","Female_song_score") # Webb FS
datalabel <- "CoopBreed FSWebb"
 columns <- c("MeanCoopTie2Noncoop","FemaleSong") # Odom FS
 datalabel <- "CoopBreed FSOdom"
# columns <- c("MeanCoopTie2Noncoop", "O.C")
df <- dfnew[which(dfnew[,columns[2]] %in% c("Present","Absent")),]
write.csv(df, "2022-03-10CoopSong_PlusWebbFS.csv")
# df[,columns[2]][which(df[,columns[2]] == "Absent")] <- 0
# df[,columns[2]][which(df[,columns[2]] == "Present")] <- 1



dfout <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","Female_song_score"), df = df, tree =  Hacktree, dummy = FALSE, nsims = 10, treelabel = "Hackett", datalabel = "CoopTie2NonCoop FSWebb")
dfDummy <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","Female_song_score"), df = df, tree =  Hacktree, dummy = TRUE, nsims = 10, treelabel = "Hackett", datalabel = "CoopTie2NonCoop FSWebb")
calcHuel(dfout,dfDummy, nsims = 10)

dfout <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","Female_song_score"), df = df, tree =  Hacktree, dummy = FALSE, nsims = 1000, treelabel = "Hackett", datalabel = "CoopTie2Coop FSWebb")
dfDummy <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","Female_song_score"), df = df, tree =  Hacktree, dummy = TRUE, nsims = 1000, treelabel = "Hackett", datalabel = "CoopTie2Coop FSWebb")

dfout2 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","Female_song_score"), df = df, tree =  Erictree, dummy = FALSE, nsims = 1000, treelabel = "Ericson", datalabel = "CoopTie2NonCoop FSWebb")
dfDummy2 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","Female_song_score"), df = df, tree =  Erictree, dummy = TRUE, nsims = 1000, treelabel = "Ericson", datalabel = "CoopTie2NonCoop FSWebb")

df <- read.csv("2022-03-10CoopSong_All.csv")
dfout3 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong"), df = df, tree =  Hacktree, dummy = FALSE, nsims = 1000, treelabel = "Hackett", datalabel = "CoopTie2NonCoop FSOdom")
dfDummy3 <- CharacterSimmaps(columns = c("MeanCoopTie2Noncoop","FemaleSong"), df = df, tree =  Hacktree, dummy = TRUE, nsims = 1000, treelabel = "Hackett", datalabel = "CoopTie2NonCoop FSOdom")

dfout4 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","FemaleSong"), df = df, tree =  Hacktree, dummy = FALSE, nsims = 1000, treelabel = "Hackett", datalabel = "CoopTie2Coop FSOdom")
dfDummy4 <- CharacterSimmaps(columns = c("MeanCoopTie2Coop","FemaleSong"), df = df, tree =  Hacktree, dummy = TRUE, nsims = 1000, treelabel = "Hackett", datalabel = "CoopTie2Coop FSOdom")

dfout5 <- CharacterSimmaps(columns = c("Final.EPP","Female_song_score"), df = df, tree =  Hacktree, dummy = FALSE, nsims = 100, treelabel = "Hackett", datalabel = "CoopTie2NonCoop FSWebb")
dfDummy5 <- CharacterSimmaps(columns = c("Final.EPP","Female_song_score"), df = df, tree =  Hacktree, dummy = TRUE, nsims = 100, treelabel = "Hackett", datalabel = "CoopTie2NonCoop FSWebb")

calcHuel(dfout,dfDummy, nsims = 1000)
calcHuel(dfout2,dfDummy2, nsims = 1000)
calcHuel(dfout3, dfDummy3, nsims= 1000)
calcHuel(dfout4, dfDummy4, nsims= 1000)
Sys.time()

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
    for (j in 1:nsims) {
      FSvec <- subsetdf[,columns[2]]
      FSvecRandom <- sample(FSvec)
      names(FSvecRandom) <- subsetdf$species
      FSsimtree <- make.simmap(tree = subsettree, x = FSvecRandom, model = "ARD", nsim = 1, Q = FSQ)
      FSsimtreesRand[[j]] <- FSsimtree
      print(j)
    }
    
    FSsimtrees<- FSsimtreesRand
    datalabel <- paste(datalabel,"DUMMYResampledCoopFS")
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
  }
  dfout <- as.data.frame(dfout)
  #colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime", "chiStat", "chiPval", "chiStatSim", "chiPvalSim")
  colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
  write.csv(dfout, file = paste(datalabel, "simmap overlap output nsim", nsims, treelabel,".csv"))
  return(dfout)
} # end function


#### Rest of procedure from Huelsenbeck et al 2003
## 
#dfout <- read.csv("CoopBreed FSOdom simmap overlap output nsim1000 Hackett .csv")
#dfout <- read.csv("CoopBreed FSWebb simmap overlap output nsim1000 Hackett .csv")
dfout <- read.csv("CoopTie2NonCoop FSWebb simmap overlap output nsim 100 Hackett .csv")




calcHuel <- function(dfout, dfDummy, nsims) {
  ExpPropAbsent0 <- as.numeric(dfout$propFSabsent)*as.numeric(dfout$propNoncoop)
  ExpPropAbsent1 <- as.numeric(dfout$propFSabsent)*as.numeric(dfout$propCoop)
  ExpPropPresent0 <- as.numeric(dfout$propFSpresent)*as.numeric(dfout$propNoncoop)
  ExpPropPresent1 <- as.numeric(dfout$propFSpresent)*as.numeric(dfout$propCoop)
  
  # sum(obs-expected) for each state
  dAbsent0 <- sum(abs(as.numeric(dfout$ObsProp0Absent) - ExpPropAbsent0))
  dAbsent1 <- sum(abs(as.numeric(dfout$ObsProp1Absent) - ExpPropAbsent1))
  dPresent0 <- sum(abs(as.numeric(dfout$ObsProp0Present) - ExpPropPresent0))
  dPresent1 <- sum(abs(as.numeric(dfout$ObsProp1Present) - ExpPropPresent1))
  
  D_real <- sum(dAbsent0, dAbsent1, dPresent0, dPresent1)/nsims
  print(paste("D_real:",D_real))
  
  #dfDummy <- read.csv("CoopBreed FSOdom DUMMYResampledCoopFS simmap overlap output nsim1000 Hackett .csv")
  dfDummy <- read.csv("CoopBreed FSWebb DUMMYResampledCoopFS simmap overlap output nsim1000 Hackett .csv")
  dfDummy <- read.csv("CoopTie2NonCoop FSWebb DUMMYResampledCoopFS simmap overlap output nsim 100 Hackett .csv")
  ExpPropAbsent0 <- as.numeric(dfDummy$propFSabsent)*as.numeric(dfDummy$propNoncoop)
  ExpPropAbsent1 <- as.numeric(dfDummy$propFSabsent)*as.numeric(dfDummy$propCoop)
  ExpPropPresent0 <- as.numeric(dfDummy$propFSpresent)*as.numeric(dfDummy$propNoncoop)
  ExpPropPresent1 <- as.numeric(dfDummy$propFSpresent)*as.numeric(dfDummy$propCoop)
  dAbsent0 <- abs(as.numeric(dfDummy$ObsProp0Absent) - ExpPropAbsent0)
  dAbsent1 <- abs(as.numeric(dfDummy$ObsProp1Absent) - ExpPropAbsent1)
  dPresent0 <- abs(as.numeric(dfDummy$ObsProp0Present) - ExpPropPresent0)
  dPresent1 <- abs(as.numeric(dfDummy$ObsProp1Present) - ExpPropPresent1)
  Dummy_dsums <- rowSums(cbind(dAbsent0, dAbsent1, dPresent0, dPresent1))
  
  pval <- sum(Dummy_dsums > D_real)/nsims
  print(paste("pval:",pval))
}







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
dfDummy <- melt(dfDummy, measure.vars = c("ObsProp0Absent","ObsProp1Absent","ObsProp0Present","ObsProp1Present"), variable.name = "state", value.name = "state.prop")
dfDummy <- melt(dfDummy, measure.vars = c("propFSabsent", "propFSpresent"), variable.name = "FS_state", value.name = "FSstate.prop")
dfDummy <- melt(dfDummy, measure.vars = c("propNoncoop", "propCoop"), variable.name = "CB_state", value.name = "CBstate.prop")

a <- ggplot(dfDummy, aes(x = state, y = state.prop, fill = pvalgroup)) + geom_boxplot() + ggtitle("dummy CoopBreed data")
b <- ggplot(dfDummy, aes(x = FS_state, y = FSstate.prop, fill = pvalgroup)) + geom_boxplot()+ ggtitle("dummy CoopBreed data")
c <- ggplot(dfDummy, aes(x = CB_state, y = CBstate.prop, fill = pvalgroup)) + geom_boxplot()+ ggtitle("dummy CoopBreed data")


dfout <- melt(dfout, measure.vars = c("ObsProp0Absent","ObsProp1Absent","ObsProp0Present","ObsProp1Present"), variable.name = "state", value.name = "state.prop")
dfout <- melt(dfout, measure.vars = c("propFSabsent", "propFSpresent"), variable.name = "FS_state", value.name = "FSstate.prop")
dfout <- melt(dfout, measure.vars = c("propNoncoop", "propCoop"), variable.name = "CB_state", value.name = "CBstate.prop")

d <- ggplot(dfout, aes(x = state, y = state.prop, fill = pvalgroup)) + geom_boxplot() + ggtitle("real CoopBreed data")
e <- ggplot(dfout, aes(x = FS_state, y = FSstate.prop, fill = pvalgroup)) + geom_boxplot() + ggtitle("real CoopBreed data")
f <- ggplot(dfout, aes(x = CB_state, y = CBstate.prop, fill = pvalgroup)) + geom_boxplot() + ggtitle("real CoopBreed data")


plot_grid(a,b,c,d,e,f, ncol = 3)