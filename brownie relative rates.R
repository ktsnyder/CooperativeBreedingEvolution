## Process brownie outputs into table
# 1/23/2024
# Kate SNyder
# 
# edited 2/1/2024 to include coop breed tie2noncoop, other song features
# 6/14/2024 - MESSED AROUND ADDING COLUMNS. MIGHT BREAK (I think just the lower section ## summary table by filelist is iffy now)

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")
require(stringr)
require(dplyr)

source("test_trait_overlap_simmaps.R") # for getLabel
filefolder = "OutputFiles" #/Brownie"
filefolder = "/Users/kate/Desktop/CooperativeBreedingEvolution/BrownieJackknifeOutputs"
filelist = list.files(filefolder)
#filelist = list.files()
filelist = filelist[which(str_detect(filelist, "brownie"))]
filelist = filelist[which(str_detect(filelist, ".csv"))]
filelist = filelist[which(str_detect(filelist, "2024-05"))]
#filelist = filelist[which(str_detect(filelist, "Updated"))]

#### summary table by trait vectors ----
#multigrouptraits = c("grouping", "social_bonds", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "MeanCoopTie2Noncoop")
multigrouptraits = c("MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "Griesser2023.Asocial0vsSocial1", "Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny", "Griesser2023.MoreThanTwoCaretakers","BiagoliniCoop", "DowningCoop", "JetzCoop", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop", "RubensteinCoop")
songtraits = "Song.rep.final"
multigrouptraits = c("HighConfidence_Coop")
songtraits = c("Syllable.rep.final", "Song.rep.final", "Syll.song.final", "Duration.final", "Interval.final", "Syllable.rep.min", "Syllable.rep.max", "Song.rep.min", "Song.rep.max", "Syll.song.min", "Syll.song.max")


browniesummary = set.seed(10)
for (i in 1:length(multigrouptraits)) {
  for (j in 1:length(songtraits)) {
  temptrait = multigrouptraits[i]
  tempsongtrait = songtraits[j]
  statelabels = getLabels(temptrait)
  tempfile = filelist[which(str_detect(filelist, temptrait) & str_detect(filelist, tempsongtrait))]  
  if (length(tempfile) == 1) {
    df = read.csv(paste0(filefolder,"/", tempfile))
    
    #subset results for plotting
    browniedf <- df[df$convergence == "Optimization has converged.",]
    browniedf <- browniedf[!is.na(browniedf$convergence),]
    
    #calculate overall mean pval
    ERloglikmean <- mean(browniedf$ERloglik)
    ARDloglikmean <- mean(browniedf$ARDloglik)
    ERARDPval = round(pchisq(2*(ARDloglikmean-ERloglikmean),1,lower.tail=FALSE), digits = 4) #testing whether the two rates of continuous trait evolution are significantly different
    
    nsims = length(df[,1])
    n1greaterthan0 = sum(df$ARDRate0 < df$ARDRate1)
    frac1greaterthan0 = n1greaterthan0/nsims
    if (frac1greaterthan0 > 0.5) {
      fasterRate = statelabels[2]
      slowerRate = statelabels[1]
      fractionGreater = frac1greaterthan0
    } else {
      fasterRate = statelabels[1]
      slowerRate = statelabels[2]
      fractionGreater = 1-frac1greaterthan0
    }
    nPvalUnder0.05 = sum(df$Pval < 0.05)
    medianPval = median(df$Pval)
    Nspecies = NA
    jackedfam = NA
    #temprow = c(tempfile, temptrait, tempsongtrait, nsims, n1greaterthan0, nPvalUnder0.05)
    temprow = c(tempfile, temptrait, tempsongtrait, nsims, fasterRate, slowerRate, fractionGreater, nPvalUnder0.05, medianPval, Nspecies, jackedfam, ERloglikmean, ARDloglikmean, ERARDPval)
    browniesummary = rbind(browniesummary, temprow)
    browniesummary = as.data.frame(browniesummary)
    colnames(browniesummary) = c("File","DiscreteTrait", "ContinuousTrait", "Nsims", "HigherRate", "LowerRate", "FractionOfSims", "num_Significant_Pval", "median_Pval", "numSpeciesInSubset", "removedFamily", "ERloglikmean", "ARDloglikmean", "BrowniePvalFromMeanLogLiks")
  } else if (length(tempfile) == 0) {
    print(paste(temptrait, tempsongtrait, "not found, skipped"))
  } else {
    print(paste(temptrait, tempsongtrait, length(tempfile), "is too many files, gonna do them all tho"))
    tempfiles = tempfile
    for (thisfile in 1:length(tempfiles)) {
      tempfile = tempfiles[thisfile]
      df = read.csv(paste0(filefolder,"/", tempfile))
      
      #subset results for plotting
      browniedf <- df[df$convergence == "Optimization has converged.",]
      browniedf <- browniedf[!is.na(browniedf$convergence),]
      
      #calculate overall mean pval
      ERloglikmean <- mean(browniedf$ERloglik)
      ARDloglikmean <- mean(browniedf$ARDloglik)
      ERARDPval = round(pchisq(2*(ARDloglikmean-ERloglikmean),1,lower.tail=FALSE), digits = 3) #testing whether the two rates of continuous trait evolution are significantly different
      
      nsims = length(df[,1])
      n1greaterthan0 = sum(df$ARDRate0 < df$ARDRate1)
      frac1greaterthan0 = n1greaterthan0/nsims
      if (frac1greaterthan0 > 0.5) {
        fasterRate = statelabels[2]
        slowerRate = statelabels[1]
        fractionGreater = frac1greaterthan0
      } else {
        fasterRate = statelabels[1]
        slowerRate = statelabels[2]
        fractionGreater = 1-frac1greaterthan0
      }
      nPvalUnder0.05 = sum(df$Pval < 0.05)
      medianPval = median(df$Pval)
      if ("numSpecies" %in% colnames(df)) {
        numSpeciesInSubset = df$numSpecies[1]
      } else {
        numSpeciesInSubset = NA
      }
      #temprow = c(tempfile, temptrait, tempsongtrait, nsims, n1greaterthan0, nPvalUnder0.05)
      if (str_detect(tempfile, "jacked")) {
        splitfilename = str_split(tempfile, "jacked")
        jackedfam = str_remove(splitfilename[[1]][2],".csv")
      } else {
        jackedfam = NA
      }
      
      temprow = c(tempfile, temptrait, tempsongtrait, nsims, fasterRate, slowerRate, fractionGreater, nPvalUnder0.05, medianPval, numSpeciesInSubset, jackedfam, ERloglikmean, ARDloglikmean, ERARDPval)
      browniesummary = rbind(browniesummary, temprow)
      browniesummary = as.data.frame(browniesummary)
   #   colnames(browniesummary) = c("File","DiscreteTrait", "ContinuousTrait", "Nsims", "HigherRate", "LowerRate", "FractionOfSims", "num_Significant_Pval", "median_Pval", "numSpeciesInSubset", "removedFamily")
      colnames(browniesummary) = c("File","DiscreteTrait", "ContinuousTrait", "Nsims", "HigherRate", "LowerRate", "FractionOfSims", "num_Significant_Pval", "median_Pval", "numSpeciesInSubset", "removedFamily", "ERloglikmean", "ARDloglikmean", "BrowniePvalFromMeanLogLiks")
    }
  }
  }
}
write.csv(browniesummary, paste0(Sys.Date(),"Brownie summary table for supp_SongRep_jackknifes.csv"), row.names = FALSE)

write.csv(browniesummary, paste0(Sys.Date(),"Brownie summary table logliks.csv"), row.names = FALSE)

#### summary table by filelist ----
source("test_trait_overlap_simmaps.R") # for getLabel
filefolder = "OutputFiles" #/Brownie"
#filefolder = "/Users/kate/Desktop/CooperativeBreedingEvolution/BrownieJackknifeOutputs"
filelist = list.files(filefolder)
filelist = filelist[which(str_detect(filelist, "brownie"))]
filelist = filelist[which(str_detect(filelist, ".csv"))]
filelist = filelist[which(str_detect(filelist, "2024-06"))]


browniesummary = set.seed(10)
for (i in 1:length(filelist)) {
  tempfile = filelist[i]
  tempdf = read.csv(paste0(filefolder,"/",tempfile))
    temptrait = tempdf$DiscreteTrait[1]
    tempsongtrait = tempdf$ContinuousTrait[1]
    statelabels = getLabels(temptrait)
    if (length(tempfile) == 1) {
      df = read.csv(paste0(filefolder,"/", tempfile))
      nsims = length(df[,1])
      n1greaterthan0 = sum(df$ARDRate0 < df$ARDRate1)
      frac1greaterthan0 = n1greaterthan0/nsims
      if (frac1greaterthan0 > 0.5) {
        fasterRate = statelabels[2]
        slowerRate = statelabels[1]
        fractionGreater = frac1greaterthan0
      } else {
        fasterRate = statelabels[1]
        slowerRate = statelabels[2]
        fractionGreater = 1-frac1greaterthan0
      }
      nPvalUnder0.05 = sum(df$Pval < 0.05)
      medianPval = median(df$Pval)
      if ("numSpecies" %in% colnames(df)) {
        numSpeciesInSubset = df$numSpecies[1]
      } else {
        numSpeciesInSubset = NA
      }
      jackedfam =NA
      #temprow = c(tempfile, temptrait, tempsongtrait, nsims, n1greaterthan0, nPvalUnder0.05)
      temprow = c(tempfile, temptrait, tempsongtrait, nsims, fasterRate, slowerRate, fractionGreater, nPvalUnder0.05, medianPval, numSpeciesInSubset, jackedfam, )
      browniesummary = rbind(browniesummary, temprow)
      browniesummary = as.data.frame(browniesummary)
      colnames(browniesummary) = c("File","DiscreteTrait", "ContinuousTrait", "Nsims", "HigherRate", "LowerRate", "FractionOfSims", "num_Significant_Pval", "median_Pval", "numSpeciesInSubset", "removedFamily")
    } else if (length(tempfile) == 0) {
      print(paste(temptrait, tempsongtrait, "not found, skipped"))
    } else {
      print(paste(temptrait, tempsongtrait, length(tempfile), "is too many files, gonna do them all tho"))
      tempfiles = tempfile
      for (thisfile in 1:length(tempfiles)) {
        tempfile = tempfiles[thisfile]
        df = read.csv(paste0(filefolder,"/", tempfile))
        nsims = length(df[,1])
        n1greaterthan0 = sum(df$ARDRate0 < df$ARDRate1)
        frac1greaterthan0 = n1greaterthan0/nsims
        if (frac1greaterthan0 > 0.5) {
          fasterRate = statelabels[2]
          slowerRate = statelabels[1]
          fractionGreater = frac1greaterthan0
        } else {
          fasterRate = statelabels[1]
          slowerRate = statelabels[2]
          fractionGreater = 1-frac1greaterthan0
        }
        nPvalUnder0.05 = sum(df$Pval < 0.05)
        medianPval = median(df$Pval)
        if ("numSpecies" %in% colnames(df)) {
          numSpeciesInSubset = df$numSpecies[1]
        } else {
          numSpeciesInSubset = NA
        }
        #temprow = c(tempfile, temptrait, tempsongtrait, nsims, n1greaterthan0, nPvalUnder0.05)
        if (str_detect(tempfile, "jacked")) {
          splitfilename = str_split(tempfile, "jacked")
          jackedfam = str_remove(splitfilename[[1]][2],".csv")
        } else {
          jackedfam = NA
        }
        
        temprow = c(tempfile, temptrait, tempsongtrait, nsims, fasterRate, slowerRate, fractionGreater, nPvalUnder0.05, medianPval, numSpeciesInSubset, jackedfam)
        browniesummary = rbind(browniesummary, temprow)
        browniesummary = as.data.frame(browniesummary)
        colnames(browniesummary) = c("File","DiscreteTrait", "ContinuousTrait", "Nsims", "HigherRate", "LowerRate", "FractionOfSims", "num_Significant_Pval", "median_Pval", "numSpeciesInSubset", "removedFamily")
      }
    }
}
browniesummary$FractionSignificant = as.numeric(browniesummary$num_Significant_Pval)/as.numeric(browniesummary$Nsims)
write.csv(browniesummary, paste0(Sys.Date()," Brownie summary table for supp_UpdatedSongData_SingleSourceCB.csv"), row.names = FALSE)




# run missing brownie analyses
dfin = read.csv("Brownie summary table for supp_allNonArchived_2024-04-24.csv")
df = dfin[which(dfin$File == "Need"),]

newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
datadf = read.csv(newdata)
datadf$AnyNoncoopEqualsNoncoop = datadf$MeanCoopTie2Noncoop
datadf$AnyNoncoopEqualsNoncoop[which(datadf$SourceDiscrepancy == 1)] = 0

source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")

currentlabel <- "_HackettOscine_"
nsim = 500
for (p in c(7,10)) {
    DiscreteTrait = df$DiscreteTrait[p]
    ContinuousTrait = df$ContinuousTrait[p]
    print(Sys.time())
    newdata = datadf
    treefile = treefile
    currentlabel <- currentlabel
    browniefunction(columns = c(DiscreteTrait, ContinuousTrait), newdata = newdata, newtree = treefile, nsim = nsim, islog = ContinuousTrait, plotsimmaps = FALSE, otherlabel = currentlabel)
}


#### run other coop w updates ----
newdata = "2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

datadf = read.csv(newdata)
tree = read.nexus(treefile)

currentlabel = "_HackettOscine_UpdatedSongData_"

ContinuousTraitvec = c(rep("Song.rep.final",5), rep("Syllable.rep.final",5))
DiscreteTraitvec = rep(c("MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "MeanCoopTie2Noncoop", "AnyNoncoopEqualsNoncoop"), 2)

nsim = 500
for (p in 6:10) {
  DiscreteTrait = DiscreteTraitvec[p]
  ContinuousTrait = ContinuousTraitvec[p]
  print(paste(Sys.time(), DiscreteTrait, ContinuousTrait, p))
  browniefunction(columns = c(DiscreteTrait, ContinuousTrait), newdata = newdata, newtree = treefile, nsim = nsim, islog = ContinuousTrait, plotsimmaps = FALSE, otherlabel = currentlabel)
}
