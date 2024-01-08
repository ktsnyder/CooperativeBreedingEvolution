# multistate brownie scratch
# 

newdata = "2024-01-08_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv"
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
df =  read.csv("2024-01-08_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_R.csv")
dfOs = df[which(df$species %in% tree$tip.label),]
dfOs %>% group_by(social_system_incl_nk_coop_Griesser2017) %>% count

subsetout = subsettreedata(columns = "social_system_incl_nk_coop_Griesser2017", newdata = dfOs, newtree = tree)
subsetDisctree = subsetout$subsettree
subsetDiscdf = subsetout$subsetdf
subsetDiscdf %>% group_by(social_system_incl_nk_coop_Griesser2017) %>% count
# songrepsubset = subsetdf[which(!is.na(subsetdf$Song.rep.final)), ]
# songrepsubset %>% group_by(social_system_incl_nk_coop_Griesser2017) %>% count
# songrepsubset %>% group_by(social_system_Griesser2017) %>% count

unique(subsetDiscdf$social_system_incl_nk_coop_Griesser2017)
discretetraitvecDisc = subsetDiscdf$social_system_incl_nk_coop_Griesser2017
names(discretetraitvecDisc) = subsetDiscdf$species

simmapER <- make.simmap(subsetDisctree,discretetraitvecDisc,nsim=1,model = "ER") 

ERmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ER")
ARDmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ARD")
anovaERARD <- anova(ERmodel,ARDmodel)
anovaERARD
SYMmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "SYM")
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
makesimmapARDrates = simmapARD$Q

# Make simmaps from data subsetted to those with song reps
subsetSong = subsettreedata(columns = c("social_system_incl_nk_coop_Griesser2017", "Song.rep.final"), newdata = df, newtree = tree, islog = "Song.rep.final")
subsetdf = subsetSong$subsetdf
subsettree = subsetSong$subsettree

unique(subsetdf$social_system_incl_nk_coop_Griesser2017)
discretetraitvec = subsetdf$social_system_incl_nk_coop_Griesser2017
names(discretetraitvec) = subsetdf$species

continuoustraitvec = subsetdf$Song.rep.final
names(continuoustraitvec) = subsetdf$species

simmappy = make.simmap(subsettree, discretetraitvec, nsim = 100, Q= rate_matrix, type = "discrete") 

plotSimmap(simmappy[[2]])

simmapfor = simmappy[[1]]
brownieliteresults <- brownie.lite(simmapfor,continuoustraitvec,maxit=75000)
statenames = names(brownieliteresults$sig2.multiple)

browniedata <- data.frame(DiscreteTrait=character(nsim),ContinuousTrait=character(nsim),Pval=numeric(nsim),ERRate=numeric(nsim),ERloglik=numeric(nsim),ERace=numeric(nsim),ARDRateCoopFam=numeric(nsim),ARDRateNoncoopFam=numeric(nsim),ARDRateCoopNonkin=numeric(nsim),ARDRateNoncoopNonkin=numeric(nsim),ARDloglik=numeric(nsim),ARDace=numeric(nsim),k2=numeric(nsim),convergence=character(nsim),simmapnumber=integer(nsim),phylanovaP=numeric(nsim),stringsAsFactors = FALSE)
brownied <- set.seed(10)
browniedf <- set.seed(10)

browniedata$DiscreteTrait = "social_system_incl_nk_coop_Griesser2017"
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

write.csv(browniedata, paste(Sys.Date(), "test brownie multistate social_system_with_nk_Griesser2017 _ARD rates from ace2.csv"), row.names = FALSE)



D0 <- density(browniedata$ARDRateCoopFam)
D1 <- density(browniedata$ARDRateNoncoopFam)
D2 <- density(browniedata$ARDRateCoopNonkin)
D3 <- density(browniedata$ARDRateNoncoopNonkin)

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

legend("topright",legend = c("Coop/Familial","Noncoop/Familial", "Coop/Nonkin", "Noncoop/Nonfamilial","Equal Rates"), lwd=1,col=c("blue","red", "green", "purple", "black"), lty = c(1,1,1,1,2), cex=1.2) # check order
pval = round(P.chisqAll,4)
#text(x = min(c(D0$x,D1$x)) + (max(c(D0$x,D1$x))-min(c(D0$x,D1$x)))*0.15, y=max(c(D0$y,D1$y))*0.65, labels = bquote(italic(p) == .(pval)), cex=2)
title(xlab=paste0("Rate of evolution of ", "log Song rep"),line = 2.5, cex.lab = 1.8)

pdens = density(browniedata$Pval)
plot(pdens)
