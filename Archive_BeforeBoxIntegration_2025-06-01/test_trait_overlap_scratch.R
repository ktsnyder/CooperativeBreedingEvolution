# test trait overlap scratch

#dfReal50 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg simmap overlap output nsim 50 Hackett .csv")
#dfReal50$X = NULL
dfReal50 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg_SubsetQs simmap overlap output nsim 50 Hackett .csv")
dfReal50$X = NULL
dfReal500 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg simmap overlap output nsim 500 Hackett .csv")
dfReal500$X = NULL
dfReal2000 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg_SubsetQs simmap overlap output nsim 500 Hackett .csv")
dfReal2000$X = NULL
dfDummy200 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg DUMMYSimHist-CoopFS simmap overlap output nsim 200 Hackett .csv") # from data with dummy using sim.history()
dfDummy200$X = NULL
dfDummyMkSimmap = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg DUMMYResampledCoopFS simmap overlap output nsim 2500 Hackett .csv") # from data with dummy using make.simmap()
dfDummyMkSimmap$X = NULL

dfDummy200 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg_SubsetQs DUMMYSimHist-CoopFS simmap overlap output nsim 200 Hackett .csv")
dfDummy200$X = NULL
dfDummy2000 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg_SubsetQs DUMMYSimHist-CoopFS simmap overlap output nsim 2000 Hackett .csv")
dfDummy2000$X = NULL

#dfoutReal[,4:17] = apply(dfoutReal[,4:17], MARGIN = 2, FUN = as.numeric)

dfReal50mod = melt(dfReal50, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfReal50mod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1)) + ggtitle("Real simmaps using make.simmap(), subsettedQs, n=50")

dfReal2000mod = melt(dfReal2000, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfReal2000mod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1)) + ggtitle("Real simmaps using make.simmap(), subsettedQs, n=2000")

dfReal500mod = melt(dfReal500, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfReal500mod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1)) + ggtitle("Real simmaps using make.simmap(), FullSetQs, n=500")

calcHuel(dfReal500, dfDummy2000, nsims_real = nsims_real, nsims_dummy = nsims_dummy)



dfDummy200mod = melt(dfDummy200, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfDummy200mod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1)) + ggtitle("Dummy simmaps using sim.history(), FullSetQs, n=200")

dfDummyMkSimmapmod = melt(dfDummyMkSimmap, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfDummyMkSimmapmod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1))+ ggtitle("Dummy simmaps using make.simmap(), FullSetQs, n=2500")

dfDummy2000mod = melt(dfDummy2000, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfDummy2000mod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1)) + ggtitle("Dummy simmaps using sim.history(), subsettedQs, n=2000")


#FSsimtrees
pdf(file = paste0(Sys.Date(), " FSAgg egDummySimmaps w SimHistory.pdf"), height=12,width=8)
layout(matrix(1:8,nrow = 2,ncol=4))
for (i in 1:8) {
  simmap <- FSsimtrees[[i]]
  py = c("black","red")
  pynamed <- py
  names(pynamed) <- c(0,1)
  
  # Plot simmap 
  plotSimmap(simmap,fsize=0.2, lwd = 0.8, colors = pynamed)
  title(paste(" ", "\nFemaleSong_Agg01 Sim.history()", "Treenum:", i),cex.main = 0.5)
}
dev.off()

#Coopsimtrees
pdf(file = paste0(Sys.Date(), "MeanCoopTie2NonCoop egDummySimmaps w SimHistory.pdf"),height=12,width=8)
layout(matrix(1:8,nrow = 2,ncol=4))
for (i in 1:8) {
  simmap <- Coopsimtrees[[i]]
  py = c("black","red")
  pynamed <- py
  names(pynamed) <- c(0,1)
  
  # Plot simmap 
  plotSimmap(simmap,fsize=0.2, lwd = 0.8, colors = pynamed)
  title(paste(" ", "\nMeanCoopTie2NonCoop Sim.history()", "Treenum:", i),cex.main = 0.5)
}
dev.off()
  

# Get tip states of Coopsimtrees
tipdf = data.frame()
for (i in 1:length(Coopsimtrees)) {
  temptree = Coopsimtrees[[i]]  
  tempmappededge = temptree$mapped.edge
  tempmappededge = as.data.frame(tempmappededge)
  tempmappededge$edges = rownames(tempmappededge)
  
  tempedges = str_split(tempmappededge$edges, "\\,", simplify = TRUE)
  tempedges = as.data.frame(tempedges)
  tempedges$TipwardsNode = as.numeric(tempedges$V2)
  tempedges$RootwardsNode = as.numeric(tempedges$V1)
  tempmappededge = cbind(tempedges, tempmappededge)
  tempmappededge$TipState = NA
  colnames(tempmappededge)[which(colnames(tempmappededge) == 0)] <- "state0"
  colnames(tempmappededge)[which(colnames(tempmappededge) == 1)] <- "state1"
  tempmappededge[which(tempmappededge$TipwardsNode %in% 1:1080 & tempmappededge$state1 > tempmappededge$state0),"TipState"] <- 1
  tempmappededge[which(tempmappededge$TipwardsNode %in% 1:1080 & tempmappededge$state1 < tempmappededge$state0),"TipState"] <- 0
  tempmappededge[which(tempmappededge$TipwardsNode %in% 1:1080 & tempmappededge$state1 == tempmappededge$state0),"TipState"] <- 0.5
  tipcount = tempmappededge %>% group_by(TipState) %>% count
  num0 = tipcount$n[which(tipcount$TipState == 0)]
  num1 = tipcount$n[which(tipcount$TipState == 1)]
  temprow = c(i, num0, num1)
  tipdf = rbind(tipdf, temprow)
}
colnames(tipdf) <- c("treenum","num0", "num1")
tipdf$Prop0 = tipdf$num0/(tipdf$num1+tipdf$num0)
CoopTipDF = tipdf
hist(CoopTipDF$Prop0, breaks = 20) # actual 0.8713
abline(v=0.8713, col = "red")
hist(CoopTipDF$num1, breaks = 20) # actual 139
abline(v=139, col = "red")


# Get tip states of FSsimtrees
tipdf = data.frame()
for (i in 1:length(FSsimtrees)) {
  temptree = FSsimtrees[[i]]  
  tempmappededge = temptree$mapped.edge
  tempmappededge = as.data.frame(tempmappededge)
  tempmappededge$edges = rownames(tempmappededge)
  
  tempedges = str_split(tempmappededge$edges, "\\,", simplify = TRUE)
  tempedges = as.data.frame(tempedges)
  tempedges$TipwardsNode = as.numeric(tempedges$V2)
  tempedges$RootwardsNode = as.numeric(tempedges$V1)
  tempmappededge = cbind(tempedges, tempmappededge)
  tempmappededge$TipState = NA
  colnames(tempmappededge)[which(colnames(tempmappededge) == 0)] <- "state0"
  colnames(tempmappededge)[which(colnames(tempmappededge) == 1)] <- "state1"
  tempmappededge[which(tempmappededge$TipwardsNode %in% 1:1080 & tempmappededge$state1 > tempmappededge$state0),"TipState"] <- 1
  tempmappededge[which(tempmappededge$TipwardsNode %in% 1:1080 & tempmappededge$state1 < tempmappededge$state0),"TipState"] <- 0
  tempmappededge[which(tempmappededge$TipwardsNode %in% 1:1080 & tempmappededge$state1 == tempmappededge$state0),"TipState"] <- 0.5
  tipcount = tempmappededge %>% group_by(TipState) %>% count
  num0 = tipcount$n[which(tipcount$TipState == 0)]
  num1 = tipcount$n[which(tipcount$TipState == 1)]
  temprow = c(i, num0, num1)
  tipdf = rbind(tipdf, temprow)
}
colnames(tipdf) <- c("treenum","num0", "num1")
tipdf$Prop0 = tipdf$num0/(tipdf$num1+tipdf$num0)
FSTipDF = tipdf
hist(FSTipDF$Prop0, breaks = 20) # actual 0.3676
abline(v=0.3676, col = "red")
hist(FSTipDF$num1, breaks = 20) # actual 683
abline(v=683, col = "red")


hist(dfoutReal$propFSabsent)
hist(dfDummy$propFSabsent)

hist(dfoutReal$propNoncoop)
hist(dfDummy$propNoncoop)

hist(dfoutReal$ObsProp0Absent)
hist(dfDummy$ObsProp0Absent)

hist(dfoutReal$ObsProp1Absent)
hist(dfDummy$ObsProp1Absent)

hist(dfoutReal$ObsProp0Present)
hist(dfDummy$ObsProp0Present)

hist(dfoutReal$ObsProp1Present)
hist(dfDummy$ObsProp1Present)


