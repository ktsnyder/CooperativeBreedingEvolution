## Process brownie outputs into table
# 1/23/2024
# Kate SNyder
# 
# edited 2/1/2024 to include coop breed tie2noncoop, other song features

source("test_trait_overlap_simmap.R") # for getLabel
filelist = list.files("OutputFiles")
filelist = filelist[which(str_detect(filelist, "brownie"))]
filelist = filelist[which(str_detect(filelist, ".csv"))]

#multigrouptraits = c("grouping", "social_bonds", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "MeanCoopTie2Noncoop")
multigrouptraits = c("Griesser2023.Asocial0vsSocial1", "Griesser2023.GroupsLargerThanPair", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.LongSocialBonds", "Griesser2017FamilialLiving", "Griesser2023.Colonial01", "Final.polygyny", "Griesser2023.MoreThanTwoCaretakers", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "MeanCoopTie2Noncoop")
songtraits = c("Song.rep.final","Syllable.rep.final", "Syll.song.final", "Duration.final", "Interval.final")

browniesummary = set.seed(10)
for (i in 1:length(multigrouptraits)) {
  for (j in 1:length(songtraits)) {
  temptrait = multigrouptraits[i]
  tempsongtrait = songtraits[j]
  statelabels = getLabels(temptrait)
  tempfile = filelist[which(str_detect(filelist, temptrait) & str_detect(filelist, tempsongtrait))]  
  if (length(tempfile) == 1) {
    df = read.csv(paste0("OutputFiles/", tempfile))
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
    print(paste(temptrait, tempsongtrait, length(tempfile), "is too many files, skipped"))
  }
  }
}
#write.csv(browniesummary, "Brownie summary table for supp.csv", row.names = FALSE)
