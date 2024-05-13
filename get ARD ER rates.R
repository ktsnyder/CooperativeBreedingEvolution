# ARD/ER table

source("findQrates.R")

newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
dfIn = read.csv(newdata)
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

SocialColumns <- c("Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.Asocial0VsSocial1", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving", "Final.polygyny", "HighConfidence_Coop", "FemaleSong_Agg01", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop")

SocialColumns %in% colnames(dfIn)

dfIn$AnyNoncoopEqualsNoncoop = dfIn$MeanCoopTie2Noncoop
dfIn$AnyNoncoopEqualsNoncoop[which(dfIn$SourceDiscrepancy == 1)] = 0
df = dfIn

allQout = set.seed(10)

for ( i in 1:length(SocialColumns)) {
  temptrait = SocialColumns[i]
  
  Qoutput = findQrates(columns = temptrait, newtree = treefile, newdata = df, plot = F)
  Qoutput
  
  browniedata = set.seed(10)
  
  browniedata$trait = temptrait
  browniedata$ERsimmapQ = gsub("ERrates ", "", Qoutput$ERrates)
  browniedata$ARDsimmapQ0to1 = gsub("ARDrates ", "", Qoutput$ARDrates[2])
  browniedata$ARDsimmapQ1to0 = gsub("ARDrates ", "", Qoutput$ARDrates[1])
  browniedata$ERsimmapQ.LogLik = Qoutput$anovaERARD$`Log lik.`[1]
  browniedata$ARDsimmapQ.LogLik = Qoutput$anovaERARD$`Log lik.`[2]
  browniedata$ARDvERsimmapQ.LRtestPval = Qoutput$anovaERARD$`Pr(>|Chi|)`[2]
  
  browniedata = as.data.frame(browniedata)
  allQout = rbind(allQout, browniedata)
}
allQout$ARDvERsimmapQ.LRtestPval.abbr = as.numeric(allQout$ARDvERsimmapQ.LRtestPval)
allQout$ARDvERsimmapQ.LRtestPval.abbr = round(allQout$ARDvERsimmapQ.LRtestPval.abbr, digits = 3)
allQout$ARDvERsimmapQ.LRtestPval.abbr[which(allQout$ARDvERsimmapQ.LRtestPval.abbr < 0.001)] <- "<0.001"
write.csv(allQout,"binary trait Qrates for supp_rounded.csv")
