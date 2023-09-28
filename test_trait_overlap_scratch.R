# test trait overlap scratch

#dfReal50 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg simmap overlap output nsim 50 Hackett .csv")
#dfReal50$X = NULL
dfReal50 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg_SubsetQs simmap overlap output nsim 50 Hackett .csv")
dfReal50$X = NULL
dfReal500 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Tie2NonCoop FSAgg simmap overlap output nsim 500 Hackett .csv")
dfReal500$X = NULL
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

dfReal500mod = melt(dfReal500, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfReal500mod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1)) + ggtitle("Real simmaps using make.simmap(), FullSetQs, n=500")

dfDummy200mod = melt(dfDummy200, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfDummy200mod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1)) + ggtitle("Dummy simmaps using sim.history(), n=200")

dfDummyMkSimmapmod = melt(dfDummyMkSimmap, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfDummyMkSimmapmod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1))+ ggtitle("Dummy simmaps using make.simmap(), n=2500")

dfDummy2000mod = melt(dfDummy2000, id.vars = 'treenum', measure.vars = c('propFSabsent', 'propFSpresent', 'propNoncoop', 'propCoop', 'ObsProp0Absent', 'ObsProp0Present', 'ObsProp1Absent', 'ObsProp1Present'))
ggplot(dfDummy2000mod) + geom_boxplot(aes(x = variable, y = value)) + scale_y_continuous(limits=c(0,1)) + ggtitle("Dummy simmaps using sim.history(), subsettedQs, n=2000")


#FSsimtrees
pdf(file = paste0(Sys.Date(), " FSAgg egDummySimmaps w SimHistory.pdf"),height=12,width=8)
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


