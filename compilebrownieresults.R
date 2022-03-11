########
#Coded by Kate T. Snyder
#Last Modified 3-11-2022
#Built using RStudio Version 1.1.453
#R Version 3.4.2
#
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#BayesTraitsV2
########


#Compile brownie jackknife plots and data

compilebrownieresults <- function(matensim, allcsvs = FALSE) {

require(ape)
require(phytools)
source(file = "Supplement_jackknifingbrownie.R")
source(file = "Supplement_plotbrowniejacks.R")
  
MateParams = c("Polygyny","EPP")
SongParams = c("Song","Syllsong","Syllrep","Interval","Duration","Rate","Continuity")

for (i in 1:2) {
  MateParam = MateParams[i]
  for (j in 1:7) {
    SongParam = SongParams[j]
browniejackout <- jackbrowniefunction(MateParam = MateParam, SongParam = SongParam, matensim = matensim, allcsvs = allcsvs)

browniejacks = browniejackout$brownielist
familyvec = browniejackout$familyvec

plotbrowniejacks(MateParam, SongParam, browniejacks = browniejacks)
  }
}
outputbrowniejacks(browniejacks = browniejacks)

}

outputbrowniejacks <- function(browniejacks = NA) {
  allbrownieresults <- data.frame()
  for (i in 1:length(browniejacks)) {
    onejackbrowniedata <- data.frame()
    brownied  <- browniejacks[[1]]
    tempfamily = brownied[1,14]
    browniedf <- brownied[brownied$convergence == "Optimization has converged.",]
    browniedf <- browniedf[!is.na(browniedf$convergence),]
    ARDratio0to1 <- browniedf$ARDRate0/browniedf$ARDRate1
    browniedf <- cbind(browniedf,ARDratio0to1)
    all0over1df <- browniedf[which(browniedf$ARDratio0to1 > 1),]
    all1over0df <- browniedf[which(browniedf$ARDratio0to1 < 1),]
    
    ERloglikmean <- mean(browniedf$ERloglik)
    ARDloglikmean <- mean(browniedf$ARDloglik)
    P.chisqAll=pchisq(2*(ARDloglikmean-as.numeric(ERloglikmean)),1,lower.tail=FALSE) #testing whether the two rates of song evolution (i.e. within polygyny or monogamy) are significantly different
    ERloglikmean0over1 <- mean(all0over1df$ERloglik)
    ARDloglikmean0over1 <- mean(all0over1df$ARDloglik)
    P.chisq0over1=pchisq(2*(ARDloglikmean0over1-as.numeric(ERloglikmean0over1)),1,lower.tail=FALSE)
    ERloglikmean1over0 <- mean(all1over0df$ERloglik)
    ARDloglikmean1over0 <- mean(all1over0df$ARDloglik)
    P.chisq1over0=pchisq(2*(ARDloglikmean1over0-as.numeric(ERloglikmean1over0)),1,lower.tail=FALSE)
    
    onejackbrowniedata[i,1] <- tempfamily
    onejackbrowniedata[i,2] <- MateParam
    onejackbrowniedata[i,3] <- SongParam
    onejackbrowniedata[i,4] <- length(browniedf$MatePar)
    onejackbrowniedata[i,5] <- ERloglikmean
    onejackbrowniedata[i,6] <- ARDloglikmean
    onejackbrowniedata[i,7] <- P.chisqAll
    onejackbrowniedata[i,8] <- ERloglikmean0over1
    onejackbrowniedata[i,9] <- ARDloglikmean0over1
    onejackbrowniedata[i,10] <- P.chisq0over1
    onejackbrowniedata[i,11] <- ERloglikmean1over0
    onejackbrowniedata[i,12] <- ARDloglikmean1over0
    onejackbrowniedata[i,13] <- P.chisq1over0
    
    
    allbrownieresults <- rbind(allbrownieresults,onejackbrowniedata)
  }
  colnames(allbrownieresults) <- c("JackedFamily","MateParam","SongParam","TotalConverged","ERmeanloglik","ARDmeanloglik","pvalOverall","ERmeanloglik0over1","ARDmeanloglik0over1","pval0over1","ERmeanloglik1over0","ARDmeanloglik1over0","pval1over0")
  write.csv(allbrownieresults,file = "AllBrownieChiSqResults.csv")
}