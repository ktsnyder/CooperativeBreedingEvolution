## Process brownie outputs into table
# 1/23/2024
# Kate SNyder
# 
# edited 2/1/2024 to include coop breed tie2noncoop, other song features

source("test_trait_overlap_simmaps.R") # for getLabel
filefolder = "OutputFiles/Brownie"
filefolder = "/Users/kate/Desktop/CooperativeBreedingEvolution/BrownieJackknifeOutputs"
filelist = list.files(filefolder)
filelist = filelist[which(str_detect(filelist, "brownie"))]
filelist = filelist[which(str_detect(filelist, ".csv"))]

#multigrouptraits = c("grouping", "social_bonds", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "MeanCoopTie2Noncoop")
multigrouptraits = c("Griesser2023.Asocial0vsSocial1", "Griesser2023.Asocial0VsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny", "Griesser2023.MoreThanTwoCaretakers", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "MeanCoopTie2Noncoop", "HighConfidence_Coop", "AnyNoncoopEqualsNoncoop")
songtraits = c("Song.rep.final", "Song.rep.min", "Song.rep.max","Syllable.rep.final", "Syllable.rep.min", "Syllable.rep.max", "Syll.song.final", "Duration.final", "Interval.final")

browniesummary = set.seed(10)
for (i in 1:length(multigrouptraits)) {
  for (j in 1:length(songtraits)) {
  temptrait = multigrouptraits[i]
  tempsongtrait = songtraits[j]
  statelabels = getLabels(temptrait)
  tempfile = filelist[which(str_detect(filelist, temptrait) & str_detect(filelist, tempsongtrait))]  
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
    #temprow = c(tempfile, temptrait, tempsongtrait, nsims, n1greaterthan0, nPvalUnder0.05)
    temprow = c(tempfile, temptrait, tempsongtrait, nsims, fasterRate, slowerRate, fractionGreater, nPvalUnder0.05)
    browniesummary = rbind(browniesummary, temprow)
    browniesummary = as.data.frame(browniesummary)
    colnames(browniesummary) = c("File","DiscreteTrait", "ContinuousTrait", "Nsims", "HigherRate", "LowerRate", "FractionOfSims", "num_Significant_Pval")
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
}
write.csv(browniesummary, "Brownie summary table for supp_jackknifes.csv", row.names = FALSE)


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

currentlabel <- "_HackettOscine"
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

