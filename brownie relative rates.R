## Process brownie outputs into table
# 1/23/2024
# Kate SNyder
# 
# 


filelist = list.files("OutputFiles")
filelist = filelist[which(str_detect(filelist, "brownie"))]
filelist = filelist[which(str_detect(filelist, ".csv"))]

multigrouptraits = c("grouping", "social_bonds", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0vsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop")
#songtraits = c("Song.rep.final","Syllable.rep.final", "Syll.song.final", "Duration.final", "Interval.final")

browniesummary = set.seed(10)
for (i in 1:length(multigrouptraits)) {
  temptrait = multigrouptraits[i]
  tempfile = filelist[which(str_detect(filelist, temptrait))]  
  if (length(tempfile) == 1) {
  df = read.csv(paste0("OutputFiles/", tempfile))
  nsims = length(df[,1])
  n1greaterthan0 = sum(df$ARDRate0 < df$ARDRate1)
  nPvalUnder0.05 = sum(df$Pval < 0.05)
  temprow = c(tempfile, temptrait, nsims, n1greaterthan0, nPvalUnder0.05)
  browniesummary = rbind(browniesummary, temprow)
  browniesummary = as.data.frame(browniesummary)
  colnames(browniesummary) = c("File","DiscreteTrait", "Nsims", "num_State1Rate_greaterthan_State0Rate", "num_Significant_Pval")
  }
}
write.csv(browniesummary, "Song.rep.final Brownie summary table for supp.csv", row.names = FALSE)
