## BayesTraits with multiPhylo
## Kate T Snyder
## Created: 3/24/2022
## Last edited: 3/30/2022

# Testing BayesTraitsv4, MCMC

library("devtools")
install_github("rgriff23/btw")
library("btw")

# library(devtools)
# install_github("rgriff23/btw", ref="v1")
# library(btw)

require(phytools)
require(btw)

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution")

EricConsensus <- read.nexus("2021-08-31ConsensusPasserineTreeEricson10_1000.nex")
HackConsensus <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000.nex")
HackMulti <- read.nexus("PasserineMultiphy1000Hackett4_nondicho.nex")

BTMulti <- read.nexus("CoopVfemsongWebb_BT_multitree.nex") # 997 tips; Webb FS classification
BT1 <- BTMulti[[1]]


newdata = "2022-03-10CoopSong_All.csv"
df <- read.csv(newdata)
WebbFSdf <- read.csv("Webb et al 2016 Female Song Plumage Data.csv")
dfnew <- merge(df, WebbFSdf, by.x = "species", by.y = "TipLabel")
columns <- c("MeanCoopTie2Noncoop","Female_song_score") # Webb FS
columns <- c("MeanCoopTie2Noncoop","FemaleSong") # Odom FS
columns <- c("MeanCoopTie2Noncoop", "O.C")
subsetdf <- dfnew[which(dfnew[,columns[2]] %in% c("Present","Absent")),]

source("subsettreedata.R")
output <- subsettreedata(columns = columns, newdata = subsetdf, newtree = HackConsensus, skinnydata = TRUE)  # newtree = BT1
multitree <- output$subsettree
#write.nexus(multitree, file = "CoopVfemsongOdom_BT_multitree.nex")
df <- output$subsetdf
df[,columns[2]][which(df[,columns[2]] == "Absent")] <- 0
df[,columns[2]][which(df[,columns[2]] == "Present")] <- 1
df[,columns[1]] <- as.character(df[,columns[1]])
df[,columns[2]] <- as.character(df[,columns[2]])

source("btwDiscreteKTS.R")
.BayesTraitsPath <- "~/Documents/BayesTraitsV4.0.0-OSX/BayesTraitsV4"
.BayesTraitsPath <- "~/Documents/BayesTraitsV3"

nocorrD <- DiscreteKTS(tree = multitree, data = df, mode = "MCMC", dependent = FALSE, silent = FALSE, pr = "PriorAll exp 10")
corrD <- DiscreteKTS(tree = multitree, data = df, mode = "MCMC", dependent = TRUE, silent = FALSE, pr = "PriorAll exp 10")


# testing new btw function
require(btw)
onetree <- multitree[[1]]
attributes(multitree)$TipLabel <- onetree$tip.label  # check -  is this ok??  - it might not be, but see solution below
multitreeTipLabel1 <- multitree
attributes(multitree)$TipLabel <- c(onetree$tip.label[200:315],onetree$tip.label[1:199])
multitreeTipLabelScramble <- multitree

#one tree test - good
# commandIndMCMC <- c("2","2")
# IndMCMCout <- bayestraits(df,onetree,commandIndMCMC)
# IndMCMCout$Log

commandIndML <- c("2","1", "Seed 10")
IndMLTipLabel1<- bayestraits(df,multitreeTipLabel1,commandIndML)
IndMLTipLabelScramble <- bayestraits(df,multitreeTipLabelScramble,commandIndML)
IndML1Log <- IndMLTipLabel1$Log$results
IndMLScrambleLog <- IndMLTipLabelScramble$Log$results

Da1 <- density(IndML1Log$alpha1)
Da2 <- density(IndML1Log$alpha2)
Db1 <- density(IndML1Log$beta1)
Db2 <- density(IndML1Log$beta2)
plot(Da1, col="blue",
     xlim=c(0,100),
     ylim=c(0,4), main = "TipLabel1")
lines(Da2, col = "red")
lines(Db1, col = "green")
lines(Db2, col = "purple")
legend("topright",legend = c("alpha1", "alpha2","beta1","beta2"), lwd=1,col=c("blue","red", "green","purple"))

Da1 <- density(IndMLScrambleLog$alpha1,)
Da2 <- density(IndMLScrambleLog$alpha2)
Db1 <- density(IndMLScrambleLog$beta1)
Db2 <- density(IndMLScrambleLog$beta2)
plot(Da1, col="blue",
     xlim=c(0,100),
     ylim=c(0,1), main = "TipLabelScramble")
lines(Da2, col = "red")
lines(Db1, col = "green")
lines(Db2, col = "purple")
legend("topright",legend = c("alpha1", "alpha2","beta1","beta2"), lwd=1,col=c("blue","red", "green","purple"))

mammaltrees <- read.nexus("Mammal.trees")


# SOLUTION: So, assigning TipLabel by just pulling the vector out of one of the trees may not be ok. BUT writing the multiPhylo as nexus and then reading it back in puts it in the right format
#OdomMultitreeNex <- read.nexus("CoopVfemsongOdom_BT_multitree.nex")
columns <- c("MeanCoopTie2Noncoop","FemaleSong") # Odom FS
columns <-  c("MeanCoopTie2Noncoop","Female_song_score") # Webb FS 
HackConsensus <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000.nex")

# setwd("~/Documents")
.BayesTraitsPath <- "~/Documents/BayesTraitsV4"

commandIndMCMC <- c("2","2", "Seed 10", "Stones 10 10000")
IndMCMCout <- bayestraits(df,BTMulti,commandIndMCMC)
commandDepMCMC <- c("3","2", "Seed 10", "Stones 10 10000")
DepMCMCout <- bayestraits(df,BTMulti,commandDepMCMC)
logMarLH_ind <- IndMCMCout$Stones$logMarLH
IndLog <- IndMCMCout$Log$results
a1 <- IndLog$alpha1
a2 <- IndLog$alpha2
b1 <- IndLog$beta1
b2 <- IndLog$beta2
qIndRates <- cbind(a2, a1, b2, a1, b1, a2, b1, b2)
colnames(qIndRates) <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
IndLogQ <- cbind(IndLog, qIndRates)
logMarLH_dep <- DepMCMCout$Stones$logMarLH
DepLog <- DepMCMCout$Log$results

2*(logMarLH_dep - logMarLH_ind)

plotDiscreteBayes(columns = columns, simplebtwOut = DepLog, nocorrDdf = IndLogQ, newpdf = TRUE, treelabel = "NoPrior_MCMC_Webb")


# ML
output <- subsettreedata(columns = columns, newdata = subsetdf, newtree = HackConsensus, skinnydata = TRUE)
tree <- output$subsettree

df <- output$subsetdf
df[,columns[2]][which(df[,columns[2]] == "Absent")] <- 0
df[,columns[2]][which(df[,columns[2]] == "Present")] <- 1
df[,columns[1]] <- as.character(df[,columns[1]])
df[,columns[2]] <- as.character(df[,columns[2]])

setwd("~/Documents")
commandIndML <- c("2","1", "Seed 10")
IndMLout <- bayestraits(df,tree,commandIndML)
commandDepML <- c("3","1", "Seed 10")#, "res q31 q42 1.5")
DepMLout <- bayestraits(df,tree,commandDepML, silent = FALSE)
logMarLH_ind <- IndMLout$Stones$logMarLH
IndLog <- IndMCMCout$Log$results
a1 <- IndLog$alpha1
a2 <- IndLog$alpha2
b1 <- IndLog$beta1
b2 <- IndLog$beta2
qIndRates <- cbind(a2, a1, b2, a1, b1, a2, b1, b2)
colnames(qIndRates) <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
IndLogQ <- cbind(IndLog, qIndRates)
logMarLH_dep <- DepMCMCout$Stones$logMarLH
DepLog <- DepMCMCout$Log$results

commandIndML <- c("2","1", "Seed 10")
IndMLout <- bayestraits(df,tree,commandIndML)



# try with Webb data - doneish
# compare ML and MCMC outputs
# try setting rates in ML - can only do each of 8 q rates, not 4 alpha/beta 

# 
newdata = "2022-03-10CoopSong_All.csv"
df <- read.csv(newdata)
columns <- c("MeanCoopTie2Noncoop", "Syll.song.final")
subsetout <- subsettreedata(columns = columns, newdata = df, newtree = HackConsensus, islog = columns[2])
subsetdf <- subsetout$subsetdf
subsettree <- subsetout$subsettree
x <- subsetdf[,columns[2]]
names(x) <- subsetdf$species
ancOut<- fastAnc(subsettree, x)
xDisc <- subsetdf[,columns[1]]
names(xDisc) <- subsetdf$species
simmap <- make.simmap(tree = subsettree, x = xDisc, model = "ARD")

