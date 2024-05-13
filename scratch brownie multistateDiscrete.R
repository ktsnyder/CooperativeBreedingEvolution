# multistate brownie scratch
# 
# Building function in 2nd section - status: update plot to allow nGroups != 4

newdata = "2024-01-08_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
newdata = "/Users/kate/Desktop/CooperativeBreedingEvolution/2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

df = read.csv(newdata)
tree = read.nexus(treefile)

# subsetout = subsettreedata(columns = c("Griesser2023.GroupsLargerThanPair", "Song.rep.final"), newdata = newdata, newtree = tree, islog = "Song.rep.final")
# subsetdf = subsetout$subsetdf
# subsettree = subsetout$subsettree

#### Groups larger than pair ----
subsetout = subsettreedata(columns = "Griesser2023.GroupsLargerThanPair", newdata = newdata, newtree = tree)
subsetdf = subsetout$subsetdf
subsettree = subsetout$subsettree

unique(subsetdf$grouping)
discretetraitvec = subsetdf$grouping
names(discretetraitvec) = subsetdf$species

simmapER <- make.simmap(subsettree,discretetraitvec,nsim=1,model = "ER") 

ERmodel <- ace(discretetraitvec,subsettree, type="discrete",model = "ER")
ARDmodel <- ace(discretetraitvec,subsettree, type="discrete",model = "ARD")
anovaERARD <- anova(ERmodel,ARDmodel)
anovaERARD
SYMmodel <- ace(discretetraitvec,subsettree, type="discrete",model = "SYM")
# this is the SYM model rate matrix: matrix(c(0, 1, 2, 3, 1, 0, 4, 5, 2, 4, 0, 6, 3,5,6,0), 4); may want to try changing some of the rates to 0
anova(ERmodel,SYMmodel)
anovaSYMARD <- anova(SYMmodel,ARDmodel)
anovaSYMARD


# get rates from ace() output
aceARDrates = ARDmodel$rates
aceARDrates = cbind(1:12, aceARDrates)
aceARDrates = as.data.frame(aceARDrates)
colnames(aceARDrates) <- c("rate_index", "rates")
rate_index_matrix = ARDmodel$index.matrix
groupnames = colnames(ARDmodel$lik.anc)

rate_matrix = matrix(rep(NA,16), nrow = 4)
rownames(rate_matrix) = colnames(rate_matrix) = groupnames

for (i in 1:4) {
  for (j in 1:4) {
    index = rate_index_matrix[i,j]
    if (!is.na(index)) {
      rate_matrix[i,j] = aceARDrates$rates[which(aceARDrates$rate_index == index)]
    }
  }
}
diagvals = rowSums(rate_matrix, na.rm = T)*-1
diag(rate_matrix) <- diagvals
rate_matrix

# get rates from make.simmap
simmapARD = make.simmap(subsettree,discretetraitvec, nsim = 1, model = "ARD")
ARDoscineQ = simmapARD$Q


## use the above ARDoscineQ rates to make simmaps on subset 
##  Note: found that those rates are too high, kind of don't follow trend from ace rates
subsetout = subsettreedata(columns = c("Griesser2023.GroupsLargerThanPair", "Song.rep.final"), newdata = newdata, newtree = tree, islog = "Song.rep.final")
subsetdf = subsetout$subsetdf
subsettree = subsetout$subsettree

unique(subsetdf$grouping)
discretetraitvec = subsetdf$grouping
names(discretetraitvec) = subsetdf$species

continuoustraitvec = subsetdf$Song.rep.final
names(continuoustraitvec) = subsetdf$species

# use ARDoscineQ
simmappyARDq = make.simmap(subsettree, discretetraitvec, nsim = 100, Q= ARDoscineQ, type = "discrete") # rates too fast but trying brownie anyway

simmappy = make.simmap(subsettree, discretetraitvec, nsim = 100, Q= rate_matrix, type = "discrete") 

plotSimmap(simmappy[[2]])

simmapfor = simmappy[[1]]
brownieliteresults <- brownie.lite(simmapfor,continuoustraitvec,maxit=75000)
statenames = names(brownieliteresults$sig2.multiple)


browniedata <- data.frame(DiscreteTrait=character(nsim),ContinuousTrait=character(nsim),Pval=numeric(nsim),ERRate=numeric(nsim),ERloglik=numeric(nsim),ERace=numeric(nsim),ARDRateAsocial=numeric(nsim),ARDRateLargeGroups=numeric(nsim),ARDRatePair=numeric(nsim),ARDRateSmallGroups=numeric(nsim),ARDloglik=numeric(nsim),ARDace=numeric(nsim),k2=numeric(nsim),convergence=character(nsim),simmapnumber=integer(nsim),phylanovaP=numeric(nsim),stringsAsFactors = FALSE)
brownied <- set.seed(10)
browniedf <- set.seed(10)

browniedata$DiscreteTrait = "grouping"
browniedata$ContinuousTrait = "logSong.rep.final"

for (i in 1:nsim) {
  simmapfor <- simmappy[[i]]
  brownieliteresults <- set.seed(10)
  
 # tryCatch(
 #   expr = {
 #     withTimeout(expr={
      
      brownieliteresults <- brownie.lite(simmapfor,continuoustraitvec,maxit=75000)
      browniedata[i,3] <- brownieliteresults$P.chisq
      browniedata[i,4] <- brownieliteresults$sig2.single
      browniedata[i,5] <- brownieliteresults$logL1
      browniedata[i,6] <- brownieliteresults$a.single
      browniedata[i,7] <- brownieliteresults$sig2.multiple[1]
      browniedata[i,8] <- brownieliteresults$sig2.multiple[2]
      browniedata[i,9] <- brownieliteresults$sig2.multiple[3]
      browniedata[i,10] <- brownieliteresults$sig2.multiple[4]
      browniedata[i,11] <- brownieliteresults$logL.multiple
      browniedata[i,12] <- brownieliteresults$a.multiple
      browniedata[i,13] <- brownieliteresults$k2
      browniedata[i,14] <- as.character(brownieliteresults$convergence)
      browniedata[i,15] <- i #}, timeout = 16, cpu=Inf, onTimeout = "error")
  # },
  # TimeoutException = function(ex) {browniedata[i,2:15]<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"timeout",i);
  # print(paste("timeout",i));
  # },
  # error = function(e) {browniedata[i,2:15]<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"error",i);
  # print(paste("compute error",i));
  # })
}


write.csv(browniedata, paste(Sys.Date(), "test brownie multistate grouping_ARD rates from ace.csv"), row.names = FALSE)


D0 <- density(browniedata$ARDRateAsocial)
D1 <- density(browniedata$ARDRateLargeGroups)
D2 <- density(browniedata$ARDRatePair)
D3 <- density(browniedata$ARDRateSmallGroups)

par(mar = c(3.8,3.5,4,1))
plot(D0,col="blue",
     xlim=c(min(c(D0$x,D1$x, D2$x, D3$x)),
            max(c(D0$x,D1$x, D2$x, D3$x))),
     ylim=c(min(c(D0$y,D1$y, D2$y, D3$y)),
            max(c(D0$y,D1$y, D2$y, D3$y))),
     main="", xlab="",ylab="", cex.axis=1.5)     
title(main="", cex.main = 2, line = 1)
title(ylab = "Frequency",line=2.5, cex.lab=1.15)
#   axis(1, cex.axis=1.2)
#    axis(2, cex.axis=1.2)
lines(D1, col="red")
lines(D2, col="green")
lines(D3, col="purple")
abline(v=browniedata$ERRate[1], lty = 2)

legend("topright",legend = c("Asocial","Pair", "Small Groups", "Large Groups","Equal Rates"), lwd=1,col=c("blue","green", "purple", "red", "black"), lty = c(1,1,1,1,2), cex=1.2)
pval = round(P.chisqAll,4)
#text(x = min(c(D0$x,D1$x)) + (max(c(D0$x,D1$x))-min(c(D0$x,D1$x)))*0.15, y=max(c(D0$y,D1$y))*0.65, labels = bquote(italic(p) == .(pval)), cex=2)
title(xlab=paste0("Rate of evolution of ", "log Song rep"),line = 2.5, cex.lab = 1.8)

# post-hoc idea: # of sims where rate in e.g. pair is the highest rate


#### Griesser 2017 social_system (coop x fam) ----
#### Now being reworked for flex/functionality
#df =  read.csv("2024-01-08_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv")
df = read.csv("2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv")
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
DiscreteTrait = "social_system_incl_nk_coop_Griesser2017"
DiscreteTrait = "grouping"
ContinuousTrait = "Song.rep.final"
ContinuousTrait = "Syllable.rep.final"
nsim = 500
plotsimmaps = TRUE
unique(df$social_system_incl_nk_coop_Griesser2017)
df$social_system_incl_nk_coop_Griesser2017[which(df$social_system_incl_nk_coop_Griesser2017 == "nk-coop")] <- "nk.coop"


subsetout = subsettreedata(columns = DiscreteTrait, newdata = df, newtree = treefile)
subsetDisctree = subsetout$subsettree
subsetDiscdf = subsetout$subsetdf
#subsetDiscdf %>% group_by(social_system_incl_nk_coop_Griesser2017) %>% count

discretetraitvecDisc = subsetDiscdf[,DiscreteTrait]
names(discretetraitvecDisc) = subsetDiscdf$species

#simmapER <- make.simmap(subsetDisctree,discretetraitvecDisc,nsim=1,model = "ER") 

ERmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ER")
ARDmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ARD")
anovaERARD <- anova(ERmodel,ARDmodel)
SYMmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "SYM")
anovaERSYM = anova(ERmodel,SYMmodel)
anovaSYMARD <- anova(SYMmodel,ARDmodel)

print(anovaERSYM)
print(anovaSYMARD)


# get rates from ace() output
aceARDratesVec = ARDmodel$rates
aceARDrates = cbind(1:length(aceARDratesVec), aceARDratesVec)
aceARDrates = as.data.frame(aceARDrates)
colnames(aceARDrates) <- c("rate_index", "rates")
rate_index_matrix = ARDmodel$index.matrix
groupnames = colnames(ARDmodel$lik.anc)
nGroups = length(groupnames)
rate_matrix = matrix(rep(NA,nGroups^2), nrow = nGroups)
rownames(rate_matrix) = colnames(rate_matrix) = groupnames

for (i in 1:nGroups) {
  for (j in 1:nGroups) {
    index = rate_index_matrix[i,j]
    if (!is.na(index)) {
      rate_matrix[i,j] = aceARDrates$rates[which(aceARDrates$rate_index == index)]
    }
  }
}
diagvals = rowSums(rate_matrix, na.rm = T)*-1
diag(rate_matrix) <- diagvals
print(rate_matrix)

# get rates from make.simmap
#simmapARD = make.simmap(subsettree,discretetraitvec, nsim = 1, model = "ARD")
#makesimmapARDrates = simmapARD$Q

# Make simmaps from data subsetted to those with song reps
subsetSong = subsettreedata(columns = c(DiscreteTrait, ContinuousTrait), newdata = df, newtree = treefile, islog = ContinuousTrait)
subsetdf = subsetSong$subsetdf
subsettree = subsetSong$subsettree

discretetraitvec = subsetdf[,DiscreteTrait]
names(discretetraitvec) = subsetdf$species

continuoustraitvec = subsetdf[,ContinuousTrait]
names(continuoustraitvec) = subsetdf$species

simmappy = make.simmap(subsettree, discretetraitvec, nsim = nsim, Q= rate_matrix, type = "discrete") 

plotSimmap(simmappy[[1]])

if (plotsimmaps) { 
  simmapFileName = paste0(DiscreteTrait, " multistate simmap plots ", ContinuousTrait, " subset", Sys.Date(),".pdf")
  pdf(simmapFileName, height = 9, width = 12)
  par(mfrow = c(2,3))
  par(mar = c(3.8,3.8,3,1))
  for (simmapNum in 1:5) {
    simmap = simmappy[[simmapNum]]
    simmapQ = simmap$Q
    SimmapStates = colnames(simmap$Q)
    py = c("orange", "#009E73", "blue", "#CC79A7")
    pynamed <- py[1:length(SimmapStates)]
    names(pynamed) <- SimmapStates
    tipcols = rep(NA, length(simmap$tip.label))
    for (stateNum in 1:length(SimmapStates)) { # get color vector of tips
      tempstate = SimmapStates[stateNum]
      tipcols[which(simmap$tip.label %in% names(which(discretetraitvec==tempstate)))] <- pynamed[tempstate]
    }
    # Plot simmap 
    plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)
    tiplabels(pch=21,bg=tipcols, col = tipcols, cex=0.3)
    legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
  }
  # plot to add rate matrix
  # Create an empty plot
  par(mar = c(5,5,4,1))
  plot(1, type = "n", xlim = c(0, ncol(simmapQ)+1), ylim = c(0, nrow(simmapQ)+1), 
       xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  
  # Add column and row names
  axis(1, at = 1:ncol(simmapQ), labels = colnames(simmapQ), las = 2, tick = FALSE)
  axis(2, at = 1:nrow(simmapQ), labels = rev(rownames(simmapQ)), las = 2, tick = FALSE)
  # Add the matrix values
  for (i in 1:nrow(simmapQ)) {
    for (j in 1:ncol(simmapQ)) {
      text(j, nrow(simmapQ) - i + 1, round(simmapQ[i, j], 6))
    }
  }
  # Add label re transitions
  axis(1, at = 0.1, labels = "To:", las = 1, tick = FALSE, font = 2)
  axis(2, at = nrow(simmapQ)+0.9, labels = "From:", las = 2, tick = FALSE, font = 2)
  title("Transition Rates")
  dev.off()
}


egSimmap = simmappy[[1]]
brownieliteresults <- brownie.lite(egSimmap,continuoustraitvec,maxit=75000)
statenames = names(brownieliteresults$sig2.multiple)

ARDRateColNames = paste0("ARDRate_", statenames)

browniedata <- data.frame(DiscreteTrait=character(nsim),ContinuousTrait=character(nsim),Pval=numeric(nsim),ERRate=numeric(nsim),ERloglik=numeric(nsim),ERace=numeric(nsim),ARDloglik=numeric(nsim),ARDace=numeric(nsim),k2=numeric(nsim),convergence=character(nsim),simmapnumber=integer(nsim),phylanovaP=numeric(nsim),stringsAsFactors = FALSE)
browniedata[,ARDRateColNames] = NA
ncolsBrownieData = length(colnames(browniedata))
brownied <- set.seed(10)
browniedf <- set.seed(10)

browniedata$DiscreteTrait = DiscreteTrait
browniedata$ContinuousTrait = paste0("log",ContinuousTrait)


for (i in 1:nsim) {
  if (i %in% seq(0,2500,by=50)) {
    print(i)
  }
  simmapfor <- simmappy[[i]]
  brownieliteresults <- set.seed(10)
  
  # tryCatch(
  #   expr = {
  #     withTimeout(expr={
  
  brownieliteresults <- brownie.lite(simmapfor,continuoustraitvec,maxit=75000)
  browniedata[i,3] <- brownieliteresults$P.chisq
  browniedata[i,4] <- brownieliteresults$sig2.single
  browniedata[i,5] <- brownieliteresults$logL1
  browniedata[i,6] <- brownieliteresults$a.single
  browniedata[i,7] <- brownieliteresults$logL.multiple
  browniedata[i,8] <- brownieliteresults$a.multiple
  browniedata[i,9] <- brownieliteresults$k2
  browniedata[i,10] <- as.character(brownieliteresults$convergence)
  browniedata[i,11] <- i 
  browniedata[i,12] <- NA 
  for (ARDcolumn in 1:length(ARDRateColNames)) {
    browniedata[i,ARDRateColNames[ARDcolumn]] <- brownieliteresults$sig2.multiple[ARDcolumn]
  }
  
      #}, timeout = 16, cpu=Inf, onTimeout = "error")
  # },
  # TimeoutException = function(ex) {browniedata[i,3:ncolsBrownieData]<-c(NA,NA,NA,NA,NA,NA,NA,NA,i,"timeout", rep(NA,times = length(ARDRateColNames)));
  # print(paste("timeout",i));
  # },
  # error = function(e) {browniedata[i,3:ncolsBrownieData]<-c(NA,NA,NA,NA,NA,NA,NA,NA,i,"error", rep(NA,times = length(ARDRateColNames)));
  # print(paste("compute error",i));
  # })
}

# add columns that say which rates are higher in each sim - added 1/10/2024
for (ARDcolumn1Num in 1:length(ARDRateColNames)) {
  ARDcolumn = ARDRateColNames[ARDcolumn1Num]
  if (ARDcolumn1Num != length(ARDRateColNames)) {
    for (ARDcolumn2Num in (ARDcolumn1Num+1):length(ARDRateColNames)) {
      ARDcolumn2 = ARDRateColNames[ARDcolumn2Num] 
      newcolname = paste0(ARDcolumn,"_greater_than_", ARDcolumn2)
      browniedata[,newcolname] = browniedata[,ARDcolumn] > browniedata[,ARDcolumn2]
    }
  }
}

nsim = length(browniedata[,1])
CompareColNames = colnames(browniedata)[grep("_greater_than_", colnames(browniedata))]
CompareColSums = data.frame(DiscreteTrait = rep(DiscreteTrait, length(CompareColNames)), ContinuousTrait = rep(ContinuousTrait, length(CompareColNames)), nsims = rep(nsim, length(CompareColNames)))
CompareColSums$RateComparison = CompareColNames
CompareColSums$Sums = colSums(browniedata[,CompareColNames])
CompareColSums$Fraction1 = CompareColSums$Sums/CompareColSums$nsims
SplitRates = str_remove_all(CompareColSums$RateComparison, "ARDRate_")
SplitRates = str_split(SplitRates, "_greater_than_")
CompareColSums$HigherRate = NA
CompareColSums$LowerRate = NA
CompareColSums$Fraction = NA
for (i in 1:length(SplitRates)) { # added 1/24/2024
  if (CompareColSums$Fraction1[i] >= 0.5) {
    CompareColSums$HigherRate[i] = SplitRates[[i]][1]
    CompareColSums$LowerRate[i] = SplitRates[[i]][2]
    CompareColSums$Fraction[i] = CompareColSums$Fraction1[i]
  } else {
    CompareColSums$HigherRate[i] = SplitRates[[i]][2]
    CompareColSums$LowerRate[i] = SplitRates[[i]][1]
    CompareColSums$Fraction[i] = 1-CompareColSums$Fraction1[i]
  }
}
compareCSVname = paste(Sys.Date(), DiscreteTrait, ContinuousTrait, "multistate aceARD Brownie", nsim, "sims COMPARE RATES.csv")
write.csv(CompareColSums, compareCSVname, row.names = F)
# end added 1/10/2024


csvname = paste(Sys.Date(), DiscreteTrait, ContinuousTrait, "multistate aceARD Brownie", nsim, "sims.csv")
pdfname = paste0(DiscreteTrait, " ", ContinuousTrait, " multistate aceARD Brownie ", nsim, " sims ", Sys.Date(), ".pdf")

#write.csv(browniedata, csvname, row.names = FALSE)


browniedata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/MultistateBrownie/2024-01-08 social_system_incl_nk_coop_Griesser2017 Song.rep.final multistate aceARD Brownie 500 sims.csv")
browniedata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/MultistateBrownie/2024-01-09 grouping Song.rep.final multistate aceARD Brownie 500 sims.csv")
browniedata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/MultistateBrownie/2024-04-02 social_system_incl_nk_coop_Griesser2017 Syllable.rep.final multistate aceARD Brownie 500 sims.csv")
browniedata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/MultistateBrownie/2024-04-02 grouping_Griesser2023 Syllable.rep.final multistate aceARD Brownie 500 sims.csv")
sum(browniedata$Pval < 0.05)/500

ARDRateColNames = colnames(browniedata)[which(str_detect(colnames(browniedata), "ARDRate") & str_detect(colnames(browniedata), "greater", negate = T))]
DiscreteTrait = browniedata$DiscreteTrait[1]
ContinuousTrait = browniedata$ContinuousTrait[1]
nsim = length(browniedata$DiscreteTrait)
if (browniedata$DiscreteTrait[1] %in% c("grouping", "grouping_Griesser2023")) {
  ARDRateColNames = c("ARDRate_asocial", "ARDRate_pair", "ARDRate_small_groups", "ARDRate_large_groups")
}

pdfname = paste0(DiscreteTrait, " ", ContinuousTrait, " multistate aceARD Brownie ", nsim, " sims ", Sys.Date(), ".pdf")

# make rate distribution plots
for (ARDcolumn in 1:length(ARDRateColNames)) {
  assign(x = paste0("D",ARDcolumn-1), value = density(browniedata[,ARDRateColNames[ARDcolumn]]))
}

pdf(pdfname, width = 8, height = 7)
par(mar = c(3.8,3.8,3,1))
par(mfrow = c(2,1))
if (length(ARDRateColNames) == 4) {
  xmax = max(c(D0$x,D1$x, D2$x, D3$x))
  if (xmax > 1.5) {xmax = .75}
  plot(D0,col="orange",
       xlim=c(min(c(D0$x,D1$x, D2$x, D3$x)),
              xmax),
       ylim=c(min(c(D0$y,D1$y, D2$y, D3$y)),
              max(c(D0$y,D1$y, D2$y, D3$y))),
       main="", xlab="",ylab="", cex.axis=1)     
  title(main="", cex.main = 2, line = 1)
  title(ylab = "Frequency",line=2.5, cex.lab=1.15)
  #   axis(1, cex.axis=1.2)
  #    axis(2, cex.axis=1.2)
  lines(D1, col="#009E73")
  lines(D2, col="blue")
  lines(D3, col="#CC79A7")
  abline(v=browniedata$ERRate[1], lty = 2)
  
  legend("topright",legend = c(ARDRateColNames,"Equal Rates"), lwd=1,col=c("orange", "#009E73", "blue", "#CC79A7", "black"), lty = c(1,1,1,1,2), cex=1) # check order
  title(xlab=paste0("Rate of evolution of ", ContinuousTrait),line = 2.5, cex.lab = 1)
  
} else if (length(ARDRateColNames) == 3) {
  plot(D0,col="orange",
       xlim=c(min(c(D0$x,D1$x, D2$x)),
              max(c(D0$x,D1$x, D2$x))),
       ylim=c(min(c(D0$y,D1$y, D2$y)),
              max(c(D0$y,D1$y, D2$y))),
       main="", xlab="",ylab="", cex.axis=1)     
  title(main="", cex.main = 2, line = 1)
  title(ylab = "Frequency",line=2.5, cex.lab=1.15)
  #   axis(1, cex.axis=1.2)
  #    axis(2, cex.axis=1.2)
  lines(D1, col="#009E73")
  lines(D2, col="blue")
  abline(v=browniedata$ERRate[1], lty = 2)
  
  legend("topright",legend = c(ARDRateColNames,"Equal Rates"), lwd=1,col=c("orange", "#009E73", "blue", "black"), lty = c(1,1,1,2), cex=1) # check order
  title(xlab=paste0("Rate of evolution of ", ContinuousTrait),line = 2.5, cex.lab = 1)
  
}
pdens = density(browniedata$Pval)
plot(pdens, main = paste("N =", nsim), ylab = "", xlab = "")
title(ylab = "Frequency",line=2.5, cex.lab=1.15)
title(xlab = "p-value", line = 2.5, cex.lab = 1)
abline(v=0.05, col = "gray")
dev.off()
