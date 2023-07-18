# Test versions of btw and bayestraits
# 
# Did repeated tests of data on BTv2, v3, v4. Independent, dependent, independent via restriction

source("subsettreedata.R")
source("~/Desktop/CooperativeBreedingEvolution/btwDiscreteKTS.R")
newdata = "~/Desktop/CooperativeBreedingEvolution/2023-06-20_CoopBreed-FemaleSong01-Song_Data_R.csv"
dataNoSongless = read.csv(newdata)
columns = c(CBcolumn, "FemaleSong_Agg01")
treefile = "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
currentlabel <- "Hackett-TieNoncoop-FSAgg mlt100_noRes"
subsetbtw <- subsettreedata(columns = columns, newdata = dataNoSongless, newtree = treefile, skinnydata = TRUE)
subsetdf <- subsetbtw$subsetdf
subsetdf[,CBcolumn] <- as.character(subsetdf[,CBcolumn])
subsettree <- subsetbtw$subsettree



# btwV1 - BayesTraitsV2 - having trouble getting BayesTraits to read the .nex file (says tree has branch lengths of 0)
.BayesTraitsPath = "~/Documents/BayesTraitsV2"
nocorrDv2 <- DiscreteKTS(subsettree, subsetdf, KeepBTInputFiles = TRUE)
corrDv2 <- DiscreteKTS(subsettree, subsetdf, dependent=TRUE, KeepBTInputFiles = TRUE, silent = FALSE)
lrtestresultsV2 <- lrtestV1(corrDv2, nocorrDv2)

# btwV1 - BayesTraitsV3
.BayesTraitsPath = "~/Documents/BayesTraitsV3"
nocorrDv3 <- DiscreteKTS(subsettree, subsetdf, KeepBTInputFiles = TRUE)
corrDv3 <- DiscreteKTS(subsettree, subsetdf, dependent=TRUE, KeepBTInputFiles = TRUE, silent = FALSE)
lrtestresultsV3 <- lrtestV1(corrDv3, nocorrDv3)

# btwV1 - BayesTraitsV4
.BayesTraitsPath = "~/Documents/BayesTraitsV4"
nocorrDv4 <- DiscreteKTS(subsettree, subsetdf)
corrDv4 <- DiscreteKTS(subsettree, subsetdf, dependent=TRUE)
lrtestresultsV4 <- lrtestV1(corrDv4, nocorrDv4)


# btwV2 (BayesTraitsV3)
set.seed(10)
commandIndML <- c("2","1", "mlt 10")#, "res q31 q42 15") # replace nocorrD
IndMLout <- bayestraits(df,tree,commandIndML, remove_files = TRUE)
nocorrDbtwV2 = IndMLout$Log$results

set.seed(10)
commandDepML <- c("3","1", "mlt 10", "Se 10")#, "res q31 q42 15")  # replace corrD
DepMLout <- bayestraits(df,tree,commandDepML, remove_files = FALSE) # , silent = FALSE
corrDbtwV2 = DepMLout$Log$results
lrtestresultsBTWv2 <- lrtest(corrDbtwV2, nocorrDbtwV2)
lrtestresultsBTWv2reverse <- lrtest(nocorrDbtwV2, corrDbtwV2)

# bayestraitsKTS
commandDepML <- c("3","1", "mlt 10", "Se 10")
DepMLout <- bayestraitsKTS(subsetdf,subsettree,commandDepML, remove_files = FALSE, version = "V3", BTdirpath = "~/Documents")


# Major differences between corrDv3 and corrDbtwV2/corrDv4: q13, q31, q34, q43

# IN terminal, try BTv3 both with "Independent" selection and with BTW's version of Independence --> 
#     input = c(input, "res q12 q34")
#input = c(input, "res q21 q43")
#input = c(input, "res q13 q24")
#input = c(input, "res q31 q42"
#"res q12 q34", "res q21 q43", "res q13 q24", "res q31 q42"
#
#
#
# alpha1 = q13, q24
# alpha2 = q12, q34
# beta1 = q42, q31
# beta2 = q43, q21


# Test BTv2 in Terminal with smaller trees (e.g. from NatComms paper) using both the Independent and Dependent+restrictions methods
# --> see Test BayesTraits Discrete Methods.xslx


# Test bayestraitsKTS many times in each version
# FINDINGS: Seeds not necessarily equal across versions - but effectvely equal if mlt is 100?
# 

newdata = "/Users/kate/Desktop/CooperativeBreedingEvolution/2023-06-20_CoopBreed-FemaleSong01HighConf-Song_Data_R.csv"
dataNoSongless = read.csv(newdata)
CBcolumn = "MeanCoopOmitTies"
columns = c(CBcolumn, "HighConfidence_FemaleSong")
columns = c(CBcolumn, "Final.polygyny")
currentlabel <- "Hackett-OmitTies-Polygyny"
treefile = "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
subsetbtw <- subsettreedata(columns = columns, newdata = dataNoSongless, newtree = treefile, skinnydata = TRUE)
subsetdf <- subsetbtw$subsetdf
subsetdf[,CBcolumn] <- as.character(subsetdf[,CBcolumn])
subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
subsettree <- subsetbtw$subsettree

seeds = 101:500

# Do Independent first
outputdf = set.seed(10)

for (i in seeds[253:300]) {
  print(paste("Independent, Seed:", i))
  tempdf = set.seed(i)
  Seed = i
  Model = "Independent"
  Method = "ML"
  commandVector = c("2", "1", "mlt 100", paste("Se", i)) # mlt number of tries
  
  outV2Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V2", silent = T, remove_files = T)
  resultsV2 = outV2Ind$Log$results
  Version = "V2"
  temprow = cbind(Seed, Version, Model, Method, resultsV2)
  tempdf = rbind(tempdf, temprow)
  
  outV3Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V3", remove_files = T)
  resultsV3 = outV3Ind$Log$results
  Version = "V3"
  temprow = cbind(Seed, Version, Model, Method, resultsV3)
  tempdf = rbind(tempdf, temprow)
  
  outV4Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V4")
  resultsV4 = outV4Ind$Log$results
  Version = "V4"
  temprow = cbind(Seed, Version, Model, Method, resultsV4)
  tempdf = rbind(tempdf, temprow)
  
  outdf = tempdf[,c("Seed","Version","Model","Method","Tree.No", "Lh")]
  outdf$q12 = tempdf$alpha2
  outdf$q13 = tempdf$alpha1
  outdf$q21 = tempdf$beta2
  outdf$q24 = tempdf$alpha1
  outdf$q31 = tempdf$beta1
  outdf$q34 = tempdf$alpha2
  outdf$q42 = tempdf$beta1
  outdf$q43 = tempdf$beta2
  outdf = cbind(outdf, tempdf[,c("Root...P.0.0.", "Root...P.0.1.", "Root...P.1.0.", "Root...P.1.1.")])
  
  outputdf = rbind(outputdf,outdf)
}
write.csv(outputdf, paste0("TestBayesTraits_", currentlabel, ".csv"))


# Next, Dependent with restricted rates to make it Independent and Dependent
for (i in seeds[257:300]) {
  print(paste("Dependent, Seed:", i))
  tempdf = set.seed(i)
  Seed = i
  Method = "ML"
  
  Model = "Independent-DependentRestrictedRates"
  commandVector = c("3", "1", "res q12 q34", "res q21 q43", "res q13 q24", "res q31 q42", "mlt 100", paste("Se", i))
  
  outV2Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V2", silent = T, remove_files = T)
  resultsV2 = outV2Ind$Log$results
  Version = "V2"
  temprow = cbind(Seed, Version, Model, Method, resultsV2)
  tempdf = rbind(tempdf, temprow)
  
  outV3Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V3", remove_files = T)
  resultsV3 = outV3Ind$Log$results
  Version = "V3"
  temprow = cbind(Seed, Version, Model, Method, resultsV3)
  tempdf = rbind(tempdf, temprow)
  
  outV4Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V4")
  resultsV4 = outV4Ind$Log$results
  Version = "V4"
  temprow = cbind(Seed, Version, Model, Method, resultsV4)
  tempdf = rbind(tempdf, temprow)
  
  
  Model = "Dependent"
  commandVector = c("3", "1", "mlt 100",paste("Se", i))
  
  outV2Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V2", silent = T, remove_files = T)
  resultsV2 = outV2Ind$Log$results
  Version = "V2"
  temprow = cbind(Seed, Version, Model, Method, resultsV2)
  tempdf = rbind(tempdf, temprow)
  
  outV3Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V3", remove_files = T)
  resultsV3 = outV3Ind$Log$results
  Version = "V3"
  temprow = cbind(Seed, Version, Model, Method, resultsV3)
  tempdf = rbind(tempdf, temprow)
  
  outV4Ind <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = "V4", remove_files = T)
  resultsV4 = outV4Ind$Log$results
  Version = "V4"
  temprow = cbind(Seed, Version, Model, Method, resultsV4)
  tempdf = rbind(tempdf, temprow)
  
  outputdf = rbind(outputdf, tempdf)
  
  if (Seed %in% c(110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 240, 260, 280, 300, 350, 400, 450, 500)) {
  write.csv(outputdf, paste0("TestBayesTraits_", currentlabel, ".csv"))
  }
}
write.csv(outputdf, paste0("TestBayesTraits_", currentlabel, ".csv"))



# outputdf100 = read.csv("/Users/kate/Documents/TestBayesTraits_Hackett-OmitTies-Polygyny_mlt100.csv")
# outputdf10 = read.csv("/Users/kate/Documents/TestBayesTraits_Hackett-OmitTies-Polygyny.csv")
# currentlabel = "Hackett-OmitTies-Polygyny_mlt10"
# outputdf10$X = NULL
# summ100 = outputdf100 %>% group_by(Version, Model) %>% summarize(n=n(), meanLh = mean(Lh))
# summ10 = outputdf10 %>% group_by(Version, Model) %>% summarize(n=n(), meanLh = mean(Lh))
# merge(summ100, summ10, by = c("Version","Model"), suffixes = c("_mlt100", "_mlt10"))
# hist(outputdf100$Lh[which(outputdf$Version=="V3" & outputdf$Model == "Independent")])
# hist(outputdf10$Lh[which(outputdf$Version=="V3" & outputdf$Model == "Independent")])
# hist(outputdf100$Lh[which(outputdf$Version=="V4" & outputdf$Model == "Independent")])
# hist(outputdf10$Lh[which(outputdf$Version=="V4" & outputdf$Model == "Independent")])

outputdf = read.csv("/Users/kate/Documents/TestBayesTraits_Hackett-OmitTies-Polygyny_mlt10.csv")
outputdf$X = NULL
currentlabel = "Hackett-OmitTies-Polygyny_mlt10_400reps"
Models = c("Independent", "Independent-DependentRestrictedRates", "Dependent")
booldf = set.seed(10)
seeds = unique(outputdf$Seed)
for (Model in Models) {
for (Seed in seeds) {
  seeddf = outputdf[which(outputdf$Seed == Seed & outputdf$Model == Model),]
  bool = seeddf[which(seeddf$Version == "V3"), 6:14] == seeddf[which(seeddf$Version == "V4"), 6:14]
  boolsum = sum(bool)
  temprow = cbind(Seed, Model, bool, boolsum)
  print(temprow)
  booldf = as.data.frame(booldf)
  booldf = rbind(booldf, temprow)
}
}

booldf %>% group_by(Model, boolsum) %>% summarize(n=n())

SeedsBoolsum9 = booldf$Seed[which(booldf$boolsum == "9" & booldf$Model == "Independent-DependentRestrictedRates")]
SeedsBoolsum0 = booldf$Seed[which(booldf$boolsum == "0" & booldf$Model == "Independent-DependentRestrictedRates")]
SeedsBoolsum9Ind = booldf$Seed[which(booldf$boolsum == "9" & booldf$Model == "Independent")]
SeedsBoolsum0Ind = booldf$Seed[which(booldf$boolsum == "0" & booldf$Model == "Independent")]
SeedsBoolsum9Dep = booldf$Seed[which(booldf$boolsum == "9" & booldf$Model == "Dependent")]
SeedsBoolsum0Dep = booldf$Seed[which(booldf$boolsum == "0" & booldf$Model == "Dependent")]
SeedsBoolsumList = list(SeedsBoolsum9, SeedsBoolsum9Ind, SeedsBoolsum9Dep)
BoolSumsOrder = c(9,9,9)
ModelOrder = c("Independent-DependentRestrictedRates", "Independent", "Dependent")

#V3Ind = outputdf[which(outputdf$Version == "V3" & outputdf$Model == "Independent"),]
#V4Ind = outputdf[which(outputdf$Version == "V4" & outputdf$Model == "Independent"),]
rates = colnames(outputdf)[c(6:9,11)]
pdf(paste("TestBayesTraitsV3vV4-Rates", currentlabel, ".pdf"), height = 9, width = 8)
par(mfrow = c(5,2))
par(mar = c(3,2,2,1))
for (j in 1:3) {
  tempseeds = SeedsBoolsumList[[j]]
  Boolsum = BoolSumsOrder[j]
  Model = ModelOrder[j]
  V3Ind = outputdf[which(outputdf$Version == "V3" & outputdf$Model == Model & outputdf$Seed %in% tempseeds),]
  V4Ind = outputdf[which(outputdf$Version == "V4" & outputdf$Model == Model & outputdf$Seed %in% tempseeds),]
  for (i in rates) {
    hist(V3Ind[,i], main = paste(Model, "V3\n",i, "NumValuesSame =", Boolsum))
    hist(V4Ind[,i], main = paste("V4\n",i, "Num Seeds =", length(tempseeds)))
  }
}
dev.off()

V3df = outputdf[which(outputdf$Version == "V3"),]
V4df = outputdf[which(outputdf$Version == "V4"),]

ratecounts = outputdf[which(outputdf$Version %in% c("V3","V4")),rates] %>% group_by_all() %>% count
unique(ratecounts$n)


outputdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/BayesTraitsDiscrete_Hackett-TieNoncoop-FSAgg mlt100_noRes.csv")
outputdf$X = NULL
library(tidyverse)
Modelsummary = outputdf %>% group_by(Version, Model) %>% summarize(MeanLh = mean(Lh), n=n(), MinLh = min(Lh), MaxLh = max(Lh), Meanq12 = mean(q12), Minq12 = min(q12), Maxq12 = max(q12), Meanq13 = mean(q13), Minq13 = min(q13), Maxq13 = max(q13), Meanq21 = mean(q21), Minq21 = min(q21), Maxq21 = max(q21), Meanq24 = mean(q24), Minq24 = min(q24), Maxq24 = max(q24), Meanq31 = mean(q31), Minq31 = min(q31), Maxq31 = max(q31),Meanq34 = mean(q34), Minq34 = min(q34), Maxq34 = max(q34), Meanq42 = mean(q42), Minq42 = min(q42), Maxq42 = max(q42), Meanq43 = mean(q43), Minq43 = min(q43), Maxq43 = max(q43) )


Lh2 = LhDependent = Modelsummary$MeanLh[which(Modelsummary$Model == "Dependent")]  
Lh1 = LhIndependent = Modelsummary$MeanLh[which(Modelsummary$Model == "Independent")]  
#Lh1 = LhIndependent = Modelsummary$MeanLh[which(Modelsummary$Model == "Independent-DependentRestrictedRates")]  
  # LRstat from btwV1

version = c("V4")
LRstatdf = c()
for (IndModel in c("Independent")) { # ,"Independent-DependentRestrictedRates"
  Lh1 = LhIndependent = Modelsummary$MeanLh[which(Modelsummary$Model == IndModel)]  
  Lh2 = LhDependent = Modelsummary$MeanLh[which(Modelsummary$Model == "Dependent")]
  version = Modelsummary$Version[which(Modelsummary$Model == IndModel)]
  LRstat = c()
  pval = c()
  Model = c()
  for (n in 1:length(Lh1)) {
    lrs = 2*(Lh1[n] - Lh2[n])
  
    if (Lh1[n] < Lh2[n]) {lrs = -lrs}
      pv = pchisq(lrs, df=1, lower.tail=F) 
      LRstat = c(LRstat, lrs)
      pval = c(pval, pv)
      Model = c(Model, IndModel)
  }
  tempdf = cbind(version, IndModel, LhIndependent, LhDependent, LRstat, pval)
  LRstatdf = rbind(LRstatdf, tempdf)
}
LRstatdf = as.data.frame(LRstatdf)


df_long <- outputdf %>%
  pivot_longer(cols = c(q12,	q13,	q21,	q24,	q31,	q34,	q42,	q43), 
               names_to = "Rate", 
               values_to = "Value")

ggplot(df_long, aes(x = Rate, y = log(Value), color = Model)) +
  geom_point() +
  facet_grid(~Version) +
  labs(color = "Model") +
  theme_minimal()

hist(outputdf$Lh[which(outputdf$Version=="V4" & outputdf$Model == "Independent")])
dev.off()


df_summary <- df_long %>%
  group_by(Version, Model, Rate) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE),
            se = sd(Value, na.rm = TRUE) / sqrt(n()), min_value = min(Value, na.rm = TRUE), max_value = max(Value, na.rm = TRUE), .groups = "drop")

pdf(paste("BayesTraits-MeanRates", currentlabel, ".pdf"), width = 11, height = 5)
ggplot(df_summary, aes(x = Rate, y = log(mean_value), color = Model)) +
  geom_point() +
  geom_errorbar(aes(ymin = log(mean_value - se), ymax = log(mean_value + se)), width = 0.2) +
  facet_wrap(~Version) +
  labs(color = "Model") +
  theme_minimal()
ggplot(outputdf, aes(x = Lh, color = Model, fill = Model)) +
  geom_histogram() +
  facet_wrap(~Version)
dev.off()

# ggplot(df_summary, aes(x = Rate, y = log(mean_value), color = Model)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = log(min_value), ymax = log(max_value)), width = 0.2) +
#   facet_wrap(~Version) +
#   labs(color = "Model") +
#   theme_minimal()
#   

# Plot rates
outputdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/BayesTraitsDiscrete_Hackett-Tie2Noncoop-FSHighConf mlt100_noRes.csv")
outputdf$X = NULL
currentlabel = "Hackett-Tie2Noncoop-FSHighConf mlt100_noRes"
source("~/Desktop/CooperativeBreedingEvolution/btwDiscreteKTS.R")
Models = c("Dependent",  "Independent") #"Independent-DependentRestrictedRates",
Versions = c("V2","V3","V4")
pdf(file = paste("BayesTraits-TransitionPlots", currentlabel,".pdf"), height = 7, width = 5)
par(mfrow = c(2,1))
par(mar = c(1,1,4,1))
for (i in 3) { 
  Version = Versions[i]
  print(Version)
  for (j in 1:2) { 
    Model = Models[j]
    print(Model)
    subdf = outputdf[which(outputdf$Model == Model & outputdf$Version == Version),]
    nsims = length(subdf$Seed)
    if (Model == "Dependent") {
      LhDependent = mean(subdf$Lh)
      LhLabel = paste("Mean Lh =", round(LhDependent,4))
      subdfDependent = subdf
    } else {
      LhIndependent = mean(subdf$Lh)
      subdfIndependent = subdf
      Lhs = c(LhDependent, LhIndependent)
      ind = sort(Lhs, index.return=TRUE)$ix
      Lh1 = Lhs[ind[1]]
      Lh2 = Lhs[ind[2]]
      LRstat = c()
      pval = c()
        lrs = 2*(Lh1 - Lh2)
        if (Lh1 < Lh2) {lrs = -lrs}
        pv = pchisq(lrs, df=1, lower.tail=F) 
        LRstat = c(LRstat, lrs)
        pval = c(pval, pv)
        print(LRstat)
        print(pval)
        LhLabel = paste("Mean Lh =", round(LhIndependent,4), "\nLRstat =", round(LRstat,4) , "Pval =", round(pval,6))
        # LRtestV1out = lrtestV1(subdfDependent, subdfIndependent)
        # LRtestStat = round(mean(LRtestV1out$LRstat),4)
        # LRtestPval = round(mean(LRtestV1out$pval),6)
        # # mockmodelDep = list()
        # # mockmodelDep$Log$results = subdfDependent[,5:18]
        # # mockmodelInd = list()
        # # mockmodelInd$Log$results = subdfIndependent[,5:18]
        # # LRtestV2out = lrtest(mockmodelDep, mockmodelInd)
        # # LRtestStat = round(mean(LRtestV2out$LRstat),4)
        # # LRtestPval = round(mean(LRtestV2out$pval),6)
        
        # LhLabel = paste(LhLabel, "\nLRtestStat =", LRtestStat, "LRtestPval =", LRtestPval)
    }
    plotdiscrete(subdf, main = paste(Version, Model, "N =", nsims, LhLabel))
  }
}
ggplot(outputdf, aes(x = Lh, color = Model, fill = Model)) +
  geom_histogram() +
  facet_wrap(~Version, dir = "v")
dev.off()


#### MCMC ----
library(psych)
commandVector = c("3", "2", "Prior q21 exp 10", "Prior q43 exp 10", "burnin 220000", "Stones 100 1000") 
outDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, version = Version, remove_files = F, BTdirpath = "~/Documents", silent = FALSE)
outDepLog = outDep$Log$results
outDepStonesLh = outDep$Stones$logMarLH
harmonic.mean(outDepLog$Lh)
mean(outDepLog$Lh)

commandVector = c("2", "2", "Prior alpha?? exp 10", "burnin 220000", "Stones 100 1000") 
outInd <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = FALSE)

# PriorAll - individual runs
commandVector = c("3", "2", "PriorAll exp 10", "burnin 220000", "Stones 100 1000") 
outPriorAllDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = FALSE, OutputFolderPath = "PriorAll-Exp-10_3")
commandVector = c("2", "2", "PriorAll exp 10", "burnin 220000", "Stones 100 1000") 
outPriorAllInd <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = FALSE, OutputFolderPath = "PriorAll-Exp-10_3")

# PriorAll - loop
columns = c("MeanCoopTie2Noncoop", "HighConfidence_FemaleSong")
#columns = c("MeanCoopTie2Noncoop", "FemaleSong_Agg01")
#columns = c("MeanCoopTie2Coop", "FemaleSong_Agg01")
newdata = "2023-06-20_CoopBreed-FemaleSong01HighConf-Song_Data_R.csv"
dataNoSongless = read.csv(newdata)
dataNoSongless$X = NULL
treefile <- "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex"
subsetbtw = subsettreedata(columns, newdata = dataNoSongless, newtree = treefile, skinnydata = T)
subsettree = subsetbtw$subsettree
subsetdf = subsetbtw$subsetdf
nsims = 40

AdditionalCommandsDep = AdditionalCommandsInd = c("PriorAll exp 10", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "PriorAll-Exp-10"
#AdditionalCommands = NULL
#OutputFolderPath = "MaxLikelihood-Defaults"
AdditionalCommandsDep = c("Prior q12 exp 11", "Prior q13 exp 4", "Prior q21 exp 7", "Prior q24 exp 11", "Prior q31 exp 43", "Prior q34 exp 20", "Prior q42 exp 60", "Prior q43 exp 7", "burnin 220000", "Stones 100 1000")
AdditionalCommandsInd = c("Prior alpha1 exp 8", "Prior beta1 exp 51", "Prior alpha2 exp 6", "Prior beta2 exp 4", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "Priors-exp-MaxLikValues20230712"

# alpha1 = q13, q24
# alpha2 = q12, q34
# beta1 = q42, q31
# beta2 = q43, q21
AdditionalCommandsDep = c("Prior q12 exp 6", "Prior q13 exp 8", "Prior q21 exp 4", "Prior q24 exp 8", "Prior q31 exp 51", "Prior q34 exp 6", "Prior q42 exp 51", "Prior q43 exp 4", "burnin 220000", "Stones 100 1000")
AdditionalCommandsInd = c("Prior alpha1 exp 8", "Prior beta1 exp 51", "Prior alpha2 exp 6", "Prior beta2 exp 4", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "Priors-exp-MaxLikIndependentValues20230712"

# qrates from ace
column1qrates= findQrates(columns = columns[1], plot=FALSE, newtree = treefile, newdata = dataNoSongless, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL)
column2qrates= findQrates(columns = columns[2], plot=FALSE, newtree = treefile, newdata = dataNoSongless, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL)
AdditionalCommandsDep = c("Prior q12 exp 0.062", "Prior q13 exp 0.009", "Prior q21 exp 0.037", "Prior q24 exp 0.009", "Prior q31 exp 0.059", "Prior q34 exp 0.062", "Prior q42 exp 0.059", "Prior q43 exp 0.037", "burnin 220000", "Stones 100 1000")
AdditionalCommandsInd = c("Prior alpha1 exp 0.009", "Prior beta1 exp 0.059", "Prior alpha2 exp 0.062", "Prior beta2 exp 0.037", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "Priors-exp-AceQrates"

# 
AdditionalCommandsDep = AdditionalCommandsInd = c("PriorAll exp 0.01", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "PriorAll-Exp-0.01"

AdditionalCommandsDep = AdditionalCommandsInd = c("burnin 220000", "Stones 100 1000")
OutputFolderPath = "NoPriors_Stones100-1000"

AdditionalCommandsDep = AdditionalCommandsInd = c("burnin 220000", "Stones 100 10000")
OutputFolderPath = "NoPriors_Stones100-10000"

nsims = 40
AdditionalCommandsDep = AdditionalCommandsInd = c("PriorAll uniform 0 50", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "PriorAll-Uniform-0-50_2"

AdditionalCommandsDep = AdditionalCommandsInd = c("PriorAll uniform 0 10", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "PriorAll-Uniform-0-10_2"

AdditionalCommandsDep = AdditionalCommandsInd = c("PriorAll uniform 0 5", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "PriorAll-Uniform-0-5"

AdditionalCommandsDep = AdditionalCommandsInd = c("PriorAll gamma 0 5", "burnin 220000", "Stones 100 1000")
OutputFolderPath = "PriorAll-Gamma-0-5"

TestPrior = TRUE
if (TestPrior) {
  AdditionalCommandsDep <- c(AdditionalCommandsDep, "TestPrior q12 1000")
  AdditionalCommandsInd <- c(AdditionalCommandsInd, "TestPrior alpha1 1000")
}

outDepdf = set.seed(10)
outInddf = set.seed(10)

for (i in 1:nsims) { 
  print(i)
  # Dependent
  commandVector = c("3", "2", AdditionalCommandsDep) 
  outPriorAllDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = T, OutputFolderPath = OutputFolderPath, TestPrior = TestPrior)
  
  outPriorAllDepOptions <- outPriorAllDep$Log$options
  outPriorAllDepResults <- outPriorAllDep$Log$results
  outPriorAllDepStonesLh <- outPriorAllDep$Stones$logMarLH
  Model = str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Model")]), "Model: ")
  if (sum(str_detect(outPriorAllDepOptions, "Iterations")) == 1) {
    Iterations = str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Iterations")]), "Iterations: ")
  } else {Iterations = NA}
  if (sum(str_detect(outPriorAllDepOptions, "Burn in")) == 1) {
    BurnIn = str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Burn in")]), "Burn in: ")
  } else {BurnIn = NA}
  Seed = str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Seed")]), "Seed: ")
  ScheduleFile = NA #str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Schedule File:")]), "Schedule File: ")
  if (!is.null(outPriorAllDepStonesLh)) {
    StonesLh = outPriorAllDepStonesLh
  } else {StonesLh = NA}
  meanResults = apply(outPriorAllDepResults[,c("Lh", "q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")], 2, mean)
  names(meanResults)[which(names(meanResults) == "Lh")] <- "SamplingMeanLh"
  outDeprow = c(i, Model, Iterations, BurnIn, Seed, ScheduleFile, StonesLh, meanResults)
  names(outDeprow) = c("Sim", "Model", "Iterations", "BurnIn", "Seed", "ScheduleFile", "StonesLh", "SamplingMeanLh", "q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
  outDepdf = rbind(outDepdf, outDeprow)
  outDepdf = as.data.frame(outDepdf)
  
  # Independent
  commandVector = c("2", "2", AdditionalCommandsInd) 
  outPriorAllInd <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = T, OutputFolderPath = OutputFolderPath, TestPrior = TestPrior)
  
  outPriorAllIndOptions <- outPriorAllInd$Log$options
  outPriorAllIndResults <- outPriorAllInd$Log$results
  outPriorAllIndStonesLh <- outPriorAllInd$Stones$logMarLH
  Model = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Model")]), "Model: ")
  #Iterations = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Iterations")]), "Iterations: ")
  #BurnIn = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Burn in")]), "Burn in: ")
  if (sum(str_detect(outPriorAllIndOptions, "Iterations")) == 1) {
    Iterations = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Iterations")]), "Iterations: ")
  } else {Iterations = NA}
  if (sum(str_detect(outPriorAllIndOptions, "Burn in")) == 1) {
    BurnIn = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Burn in")]), "Burn in: ")
  } else {BurnIn = NA}
  Seed = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Seed")]), "Seed: ")
  ScheduleFile = NA #str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Schedule File:")]), "Schedule File: ")
  #StonesLh = outPriorAllIndStonesLh
  if (!is.null(outPriorAllIndStonesLh)) {
    StonesLh = outPriorAllIndStonesLh
  } else {StonesLh = NA}
  meanResultsInd = apply(outPriorAllIndResults[,c("Lh", "alpha1", "beta1", "alpha2", "beta2")], 2, mean)
  names(meanResultsInd)[which(names(meanResultsInd) == "Lh")] <- "SamplingMeanLh"
  outIndrow = c(i, Model, Iterations, BurnIn, Seed, ScheduleFile, StonesLh, meanResultsInd)
  names(outIndrow) = c("Sim", "Model", "Iterations", "BurnIn", "Seed", "ScheduleFile", "StonesLh", "SamplingMeanLh", "alpha1", "beta1", "alpha2", "beta2")
  outInddf = rbind(outInddf, outIndrow)
  outInddf = as.data.frame(outInddf)
  
  if (i %in% seq(1,100, by = 9)) {
    filenameDep = paste0(OutputFolderPath, "/",columns[1], "-", columns[2], "_", "Dependent_", Sys.Date(),".csv")
    write.csv(outDepdf, file = filenameDep, row.names = FALSE)
    filenameInd = paste0(OutputFolderPath, "/", columns[1], "-", columns[2], "_", "Independent_", Sys.Date(),".csv")
    write.csv(outInddf, file = filenameInd, row.names = FALSE)
    print(paste("Finished loop ", i, "with columns", columns[1], columns[2], "at", Sys.time()))
  }
}
filenameDep = paste0(OutputFolderPath, "/",columns[1], "-", columns[2], "_", "Dependent_", Sys.Date(),".csv")
write.csv(outDepdf, file = filenameDep, row.names = FALSE)
filenameInd = paste0(OutputFolderPath, "/", columns[1], "-", columns[2], "_", "Independent_", Sys.Date(),".csv")
write.csv(outInddf, file = filenameInd, row.names = FALSE)
print(paste("Finished loop ", i, "with columns", columns[1], columns[2], "at", Sys.time(), "in", OutputFolderPath))




# Hyperpriors
commandVector = c("3", "2", "HyperPriorAll exp 0 10", "burnin 250000", "Stones 100 1000") 
outHyperDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = FALSE, OutputFolderPath = "HyperPriorAll-Exp-0-10_3")
commandVector = c("2", "2", "HyperPriorAll exp 0 10", "burnin 250000", "Stones 100 1000") 
outHyperInd <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = FALSE, OutputFolderPath = "HyperPriorAll-Exp-0-10_3")
print(paste(Sys.time(), "both hyper jumps ended"))

# Reverse Jump
commandVector = c("3", "2", "RevJump exp 10", "burnin 250000", "Stones 100 1000") 
outHyperDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = FALSE, OutputFolderPath = "RevJump-exp-10_2")

commandVector = c("3", "2", "RevJumpHP exp 0 10", "burnin 250000", "Stones 100 1000") 
outHyperRevJump <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = FALSE, OutputFolderPath = "RevJumpHP-exp-0-100_2")
print(paste(Sys.time(), "both rev jumps ended"))




Results <- outHyperRevJump$Log$results
Options <- outHyperDep$Log$options
logMarLH <- outHyperDep$Stones$logMarLH

write.csv()

commandVector = c("TestPrior", "gamma 0 10") 
TestPriorOff <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, Model)


AdditionalCommandsDep = AdditionalCommandsInd = c("PriorAll uniform 0 5", "TestPrior q12 1000", "burnin 220000", "Stones 100 1000")
commandVector = c("3", "2", AdditionalCommandsDep) 
outPriorAllDep <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = F, OutputFolderPath = OutputFolderPath)

commandVector = c("2", "2", AdditionalCommandsInd) 
outPriorAllInd <- bayestraitsKTS(data = subsetdf, tree = subsettree, commands = commandVector, remove_files = F, BTdirpath = "~/Documents", silent = F, OutputFolderPath = OutputFolderPath, TestPrior = T)


my_wd = "/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Test_BayesTraits_MCMC_settings/PriorAll-Uniform-0-5/Discrete-Independent_MCMC"
system(paste(paste0(BTdir, "/BayesTraits", BTversionNum), paste0(my_wd, "/tree.nex"), paste0(my_wd, "/data.txt"), paste0("< ", my_wd, "/inputfile.txt"), paste0(">", my_wd, "/LogAllOutputs.txt")), ignore.stdout = silent)
