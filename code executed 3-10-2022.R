## Temp - copied from Run_CoopBreed_Analyses for running 3/10/2022-3/10/2022
# 12:30pm - re-ran CoopBreed_species_summary and 
# 
setwd("~/Desktop/CooperativeBreedingEvolution/")
setwd("~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB")
require(dplyr)

source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")

songfeatures <- c("Syll.song.final", "Song.rep.final", "Syll.rep.final")
classmethods = c("MeanCoopOmitTies", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "AnyCoopEqualsCoop")
newdata = "2022-03-10CoopSong_All.csv"
olddata = "2022-03-09CoopSong__All_preJetzCorrection.csv"
treefile = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" 
treelabel = "PasserineTreeEricson-"

dfnew <- read.csv(newdata)
dfnew$species[which(dfnew$SourceDiscrepancy == TRUE)]
dfnew$species[which(dfnew$numSourcesNonCoop == dfnew$numSourcesCoop & dfnew$numSourcesCoop > 0)]
sum(dfnew$MeanCoopTie2Noncoop == 1, na.rm = TRUE)

dfold <- read.csv(olddata)
dfold$species[which(dfold$SourceDiscrepancy == TRUE)]
dfold$species[which(dfold$numSourcesNonCoop == dfold$numSourcesCoop & dfold$numSourcesCoop > 0)]
sum(dfold$MeanCoopTie2Noncoop == 1, na.rm = TRUE)

dfnew$species[which(dfnew$Jetz != dfnew$Cockburn)]
unique(dfnew$Jetz)
dfnew %>% group_by(Jetz, Cockburn) %>% summarize(n=n())

jetzview <- dfnew %>% group_by(Jetz, Cockburn, JetzSource) %>% summarize(n=n())

dfnew$species[which(dfnew$numSourcesNonCoop == dfnew$numSourcesCoop & dfnew$numSourcesNonCoop != 0)]

dfnew %>% group_by(MeanCoopTie2Noncoop ,Griesser2017FamilialLiving) %>% summarize(n = n())

source("simplebtwDiscrete.R")
classmethods <- c("MeanCoopOmitTies",  "MeanCoopTie2Noncoop", "MeanCoopTie2Coop",  "AnyCoopEqualsCoop" )
nsim = 100
for (i in classmethods) {
  currentclassmethod <- i
  columns <- c(currentclassmethod, "FemaleSong")
  print(columns)
  print(Sys.time())
  simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, treelabel = treelabel, nsim = nsim)
  plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = NULL)
}


columns <- c("Griesser2017FamilialLiving", "FemaleSong")
print(columns)
print(Sys.time())
simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, treelabel = treelabel, nsim = nsim)
plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = " 'COOPERATIVE' ACTUALLY 'FAMILIAL LIVING' ")

columns <- c("Kin_NK", "FemaleSong")
print(columns)
print(Sys.time())
#simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, treelabel = treelabel, nsim = nsim)  - this seems to have broken BT, be careful
plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = " 'COOPERATIVE' ACTUALLY 'KIN COOP v NONKIN COOP' ")

columns <- c("MeanCoopTie2Noncoop", "O.C")
print(columns)
print(Sys.time())
simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, treelabel = treelabel, nsim = nsim)
plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = "")

columns <- c("Final.polygyny", "FemaleSong")
print(columns)
print(Sys.time())
simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, treelabel = treelabel, nsim = nsim)
plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = " 'COOPERATIVE' ACTUALLY 'POLYGYNY v MONOGAMY' ")

columns <- c("Final.EPP", "FemaleSong") 
print(columns)
print(Sys.time())
simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, treelabel = treelabel, nsim = nsim)
plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = " 'COOPERATIVE' ACTUALLY 'EPP' ")

columns <- c("MeanCoopTie2Noncoop" ,"Griesser2017FamilialLiving")
print(columns)
print(Sys.time())
simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, treelabel = treelabel, nsim = nsim)
plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = "")




# require("btw")
# currentclassmethod = classmethods[2]   # MeanCoopTie2Noncoop
# subsetbtw <- subsettreedata(columns = c(currentclassmethod,"FemaleSong"), newdata = "2022-03-09CoopSong__All.csv", newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", skinnydata = TRUE)
# currentlabel <- paste0("PasserTreeEric-",currentclassmethod)
# subsetdf <- subsetbtw$subsetdf
# colnames(subsetdf)[which(colnames(subsetdf) == currentclassmethod)] <- "CoopBreed"
# subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Present")] <- "1"
# subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Absent")] <- "0"
# subsetdf$CoopBreed <- as.character(subsetdf$CoopBreed)
# subsettree <- subsetbtw$subsettree
# subsetdf %>% group_by(CoopBreed,FemaleSong) %>% summarise(n=n())
# simplebtwOut <- set.seed(10)
# nsim = 2000
# for (n in 1:nsim) {
#   nocorrD <- Discrete(subsettree, subsetdf)
#   corrD <- Discrete(subsettree, subsetdf, dependent=TRUE)
#   lrtestresults <- lrtest(corrD, nocorrD)
#   tempRow <- cbind(corrD, lrtestresults)
#   simplebtwOut <- rbind(simplebtwOut, tempRow)
# }
# simplebtwOut <- as.data.frame(simplebtwOut)
# means <- apply(X = simplebtwOut,MARGIN = 2,FUN = mean)
# meansdf <- as.data.frame(rbind(means,means))
# #meansdf <- as.data.frame(as.matrix(means))
# pvalMed <- median(simplebtwOut$pval)
# pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs FemaleSong ", nsim, "sims.pdf"))
# plotdiscrete(meansdf[1,1:14], main = paste(currentlabel, "vs FemSong, \nnsims =",nsim, "median pval =", pvalMed))
# dev.off()
# 
# # BayesTraits female song - MeanCoopTie2Coop
# currentclassmethod = classmethods[3]
# subsetbtw <- subsettreedata(columns = c(currentclassmethod,"FemaleSong"), newdata = "2022-03-09CoopSong__All.csv", newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex", skinnydata = TRUE)
# currentlabel <- paste0("PasserTreeEric-",currentclassmethod)
# subsetdf <- subsetbtw$subsetdf
# colnames(subsetdf)[which(colnames(subsetdf) == currentclassmethod)] <- "CoopBreed"
# subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Present")] <- "1"
# subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Absent")] <- "0"
# subsetdf$CoopBreed <- as.character(subsetdf$CoopBreed)
# subsettree <- subsetbtw$subsettree
# subsetdf %>% group_by(CoopBreed,FemaleSong) %>% summarise(n=n())
# simplebtwOut <- set.seed(10)
# nsim = 2000
# for (n in 1:nsim) {
#   nocorrD <- Discrete(subsettree, subsetdf)
#   corrD <- Discrete(subsettree, subsetdf, dependent=TRUE)
#   lrtestresults <- lrtest(corrD, nocorrD)
#   tempRow <- cbind(corrD, lrtestresults)
#   simplebtwOut <- rbind(simplebtwOut, tempRow)
# }
# simplebtwOut <- as.data.frame(simplebtwOut)
# means <- apply(X = simplebtwOut,MARGIN = 2,FUN = mean)
# meansdf <- as.data.frame(rbind(means,means))
# #meansdf <- as.data.frame(as.matrix(means))
# pvalMed <- median(simplebtwOut$pval)
# pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs FemaleSong ", nsim, "sims.pdf"))
# plotdiscrete(meansdf[1,1:14], main = paste(currentlabel, "vs FemSong, \nnsims =",nsim, "median pval =", pvalMed))
# dev.off()
# 
# 
# # BayesTraits song features
# source("btwfunction.R")
# source("BayesPlots_choosebin.R")
# songfeatures <- c("Syllable.rep.final", "Syll.song.final", "Song.rep.final", "Duration.final", "Interval.final", "Song.rate")
# newdata = "2022-03-09CoopSong__All.csv"
# #newdata = "2022-03-08CoopSong_MeanCoop_All.csv"
# treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
# currentclassmethod = classmethods[2]
# currentlabel <- paste0("PasserTreeEric-",currentclassmethod)
# nsim = 50
# for (k in 3) {  # only did SongRep, "MeanCoopTie2Noncoop"
#   feature <- songfeatures[k]
#   print(feature)
#   print(currentlabel)
#   btwfunction(MateParam = currentclassmethod,SongParam = feature, plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = nsim, newtreefile = treefile, newdata = newdata)
#   feature <- songfeatures[k]
#   filename <- paste0(Sys.Date(),"Bayes", currentclassmethod,feature, nsim, "reps.csv") 
#   BTdf <- read.csv(filename)
#   #BTdf <- BTdf[,which(colnames(BTdf) != "X")]
#   colnames(BTdf)[16:18] <- c("LRstat", "LRpval", "songcontvec")
#   transitionBinplots(MateParam = currentclassmethod,SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3)
# } 
# 
# 
# ## Brownie
nsim = 100
songfeatures <- c("Syllable.rep.final","Syll.song.final", "Song.rep.final")
for (k in c(2,3,1)) {
  currentclassmethod = classmethods[2]
  currentlabel <- paste0("PasserTreeEric-",currentclassmethod)
  feature <- songfeatures[k]
  print(feature)
  print(currentlabel)
  browniefunction(columns = c(currentclassmethod, feature), newdata = newdata, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)

  plotbrownie(data = paste0(Sys.Date(),currentclassmethod,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(currentclassmethod,feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)


  # currentclassmethod = classmethods[3]
  # currentlabel <- paste0("PasserTreeEric-",currentclassmethod)
  # print(currentlabel)
  # browniefunction(columns = c(currentclassmethod, feature), newdata = newdata, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)
  # 
  # plotbrownie(data = paste0(Sys.Date(),currentclassmethod,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(currentclassmethod,feature), discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
}

# 
# ## ACE tree
# newdata = "2022-03-09CoopSong__All.csv"
# treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex"
# treefile <- "2020-10-11ConsensusPasserineTreeHack100.nex"
# source("plotACEtree.R")
# for (feature in songfeatures) {  # repeated Syllrep, Songrep, Interval with "ER"
#   currentclassmethod <- classmethods[2]
#   currentlabel <- paste0("PasserTreeHack-ARD-",currentclassmethod)
#   print(feature)
#   print(currentlabel)
# plotACEtree(columns = c(currentclassmethod,feature), newdata = newdata, newtree = treefile, islog = feature, discretelabels = c("NonCoop","Coop"), discretemodel = "ARD", otherlabel = currentlabel)
# 
#   currentclassmethod <- classmethods[3]
#   currentlabel <- paste0("PasserTreeHack-ARD-",currentclassmethod)
#   print(feature)
#   print(currentlabel)
#   plotACEtree(columns = c(currentclassmethod,feature), newdata = newdata, newtree = treefile, islog = feature, discretelabels = c("NonCoop","Coop"), discretemodel = "ARD", otherlabel = currentlabel)
# }
# 
# ### Did not get below here
# newdata = "2022-03-09CoopSong__All.csv"
# treefile <- "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" # 3/8/2022
# currentclassmethod = classmethods[3]
# currentlabel <- paste0("PasserTreeEric-",currentclassmethod)
# print(currentlabel)
# nsim = 150
# for (k in 1:6) { 
#   feature <- songfeatures[k]
#   print(feature)
#   print(currentlabel)
#   btwfunction(MateParam = currentclassmethod,SongParam = feature, plot=FALSE, jackknife = FALSE, csvsout = TRUE, nsim = nsim, newtreefile = treefile, newdata = newdata)
#   feature <- songfeatures[k]
#   filename <- paste0(Sys.Date(),"Bayes", currentclassmethod,feature, nsim, "reps.csv") 
#   BTdf <- read.csv(filename)
#   #BTdf <- BTdf[,which(colnames(BTdf) != "X")]
#   colnames(BTdf)[16:18] <- c("LRstat", "LRpval", "songcontvec")
#   transitionBinplots(MateParam = currentclassmethod,SongParam = feature, df = BTdf,newpdf = TRUE, nsim = nsim, binnum = 3)
# } 
# 
# 
# ## scratch corHMM etc
# require(corHMM)
# source("subsettreedata.R")
# classmethods = c("MeanCoopOmitTies", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "AnyCoopEqualsCoop")
# newdata = "2022-03-09CoopSong__All.csv"
# treefile = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex" 
# subsetout <- subsettreedata(columns = c(classmethods[2],"FemaleSong"), newdata = newdata, newtree = treefile, skinnydata = TRUE )
# subsetdf <- subsetout$subsetdf
# subsettree <- subsetout$subsettree
# subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Present")] <- 1
# subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Absent")] <- 0
# subsetdf$FemaleSong <- as.integer(subsetdf$FemaleSong)
# 
# corHMM_results<-corDISC(subsettree,subsetdf,ntraits=2,model="ARD", node.states="marginal", diagn=FALSE)
# 
# coopvec <- subsetdf$MeanCoopTie2Noncoop
# names(coopvec) <- subsetdf$species
# femvec <- subsetdf$FemaleSong
# names(femvec) <- subsetdf$species
# fitPagel_results <- fitPagel(subsettree,x=coopvec,y=femvec)
# plot(fitPagel_results)
# 
# library(phylolm)
# fit0 <- phyloglm(coopvec ~ 1, phy = subsettree)
# fit2 = phyloglm(coopvec ~ femvec, phy=subsettree)
# data.frame(model=c("Coop ~ FemSong","Null model"),
#            log_lik=c(logLik(fit2)$logLik,logLik(fit0)$logLik),
#            df=c(logLik(fit2)$df,logLik(fit0)$df),
#            AIC=c(AIC(fit2),AIC(fit0)))
# fit0 <- phyloglm(femvec ~ 1, phy = subsettree)
# fit2 = phyloglm(femvec ~ coopvec, phy=subsettree)
# data.frame(model=c("FemSong ~ Coop","Null model"),
#            log_lik=c(logLik(fit2)$logLik,logLik(fit0)$logLik),
#            df=c(logLik(fit2)$df,logLik(fit0)$df),
#            AIC=c(AIC(fit2),AIC(fit0)))
# 
# 
# 
# ape::compar.gee(coopvec ~ femvec, phy = subsettree)  #not working
# 
# # this works
# pic.coop <- pic(coopvec, subsettree)
# pic.fem <- pic(femvec, subsettree)
# cor.test(pic.coop,pic.fem)
# cor.test(pic.fem,pic.coop)
# lm(pic.fem ~ pic.coop-1)
# lm(pic.coop ~ pic.fem-1)
# 
# require("nlme")
# tips <- subsettree$tip.label
# tree.corr <- corBrownian(phy=subsettree, form = ~ tips)
# # The PGLS ANOVA
# test.anova <- gls(coopvec ~ femvec - 1, correlation = tree.corr)
# summary(test.anova)
# test.anova <- gls(femvec ~ coopvec - 1, correlation = tree.corr)
# summary(test.anova)
# 
# 
# phylANOVA(subsettree, x = coopvec, y = femvec, nsim = 50000)
# phylANOVA(subsettree, x = femvec, y = coopvec, nsim = 50000)
# 
# # felsenstein's threshold method, via phytools
# df <- cbind(coopvec, femvec)
# row.names(df) <- names(coopvec)
# thresh.results <- threshBayes(subsettree, df, ngen = 20000)
# plot(density(thresh.results))
# 
# fakedata <- sample(c(0,1), size = length(coopvec), replace = TRUE)
# df <- cbind(coopvec, fakedata)
# row.names(df) <- names(coopvec)
# thresh.results <- threshBayes(subsettree, df, ngen = 20000)
# plot(density(thresh.results))
# 
# 
# subset <- subsettreedata(columns = c(classmethods[2],classmethods[4]), newdata = newdata, newtree = treefile, skinnydata = TRUE )
# df <- subset$subsetdf[,2:3]
# row.names(df) <- subset$subsetdf[,1]
# tree <- subset$subsettree
# thresh.results <- threshBayes(subsettree, df, ngen = 20000)
# plot(density(thresh.results))
