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
# 3/11/2022 allcsvs used as another string to search for in filelist for csvs (must be start of file - e.g. the date)


#Brownie Jack Plot


plotbrowniejacks <- function(columns, browniejacks = NA, allcsvs = FALSE, nsim = NULL, otherlabel = NULL) {   #browniejacks is a list of dataframes, all the same Mate/Song combo, but each dataframe the test with one family removed
require(ape)
require(phytools)
  
  MateParam <- columns[1]
  SongParam <- columns[2]
  if (allcsvs != FALSE) {
    require(stringr)
    filelist <- list.files(getwd())
    # 2022-03-11 BrownieJack MeanCoopTie2Noncoop Syllable.rep.final Wilsonia_canadensis .csv
    csvPatterns <- c(allcsvs, "BrownieJack", MateParam, SongParam, ".csv")
    #filelist[which(str_detect(filelist, paste(csvPatterns, collapse = "&") ) )]
    browniejacksfiles <- filelist[which(str_detect(filelist, paste(allcsvs,"BrownieJack",MateParam,SongParam)))]
    browniejacks = list()
    for (l in 1:length(browniejacksfiles)) {
      tempfile <- browniejacksfiles[l]
      browniejacks[[l]] <- read.csv(tempfile)
    }
  }
  
    pdf(paste(Sys.Date(),"PlotBrownieJacks",otherlabel,MateParam, SongParam,nsim,"sim.pdf"), height = 10, width = 8)
  par(mfrow = c(6,4), cex = 0.3, mar=rep(2,4))

  if (MateParam == "Polygyny") {
    state0 <- "Monogamy"
    state1 <- "Polygyny"
  } else if (MateParam == "EPP") {
    state0 <- "Low EPP"
    state1 <- "High EPP"
  } else {
    state0 <- "NonCooperative"
    state1 <- "Cooperative"
  }

  for (k in 1:length(browniejacks)) {
    dftempNotConv <- browniejacks[[k]] 
    tempfamily = dftempNotConv[1,"jackedfam"]
    dftemp <- dftempNotConv[dftempNotConv$convergence == "Optimization has converged.",]
    dftemp <- dftemp[!is.na(dftemp$convergence),]

    ARDratio0to1 <- dftemp$ARDRate0/dftemp$ARDRate1
    dftemp <- cbind(dftemp,ARDratio0to1)
    all0over1df <- dftemp[which(dftemp$ARDratio0to1 > 1),]
    all1over0df <- dftemp[which(dftemp$ARDratio0to1 < 1),]
    ERloglikmean <- mean(dftemp$ERloglik)
    ARDloglikmean <- mean(dftemp$ARDloglik)
    P.chisqAll=pchisq(2*(ARDloglikmean-as.numeric(ERloglikmean)),1,lower.tail=FALSE)  #testing whether the two rates of song evolution (i.e. within polygyny or monogamy) are significantly different
    ERloglikmean0over1 <- mean(all0over1df$ERloglik)
    ARDloglikmean0over1 <- mean(all0over1df$ARDloglik)
    P.chisq0over1=pchisq(2*(ARDloglikmean0over1-as.numeric(ERloglikmean0over1)),1,lower.tail=FALSE)
    ERloglikmean1over0 <- mean(all1over0df$ERloglik)
    ARDloglikmean1over0 <- mean(all1over0df$ARDloglik)
    P.chisq1over0=pchisq(2*(ARDloglikmean1over0-as.numeric(ERloglikmean1over0)),1,lower.tail=FALSE)
    

    familylabel = paste(tempfamily)

    D0 <- density(dftemp$ARDRate0)
    D1 <- density(dftemp$ARDRate1)
    
    par(mar = c(3.8,3.5,4,1))
    plot(D0,col="blue",
         xlim=c(min(c(D0$x,D1$x)),
                max(c(D0$x,D1$x))),
         ylim=c(min(c(D0$y,D1$y)),
                max(c(D0$y,D1$y))),
         main="", xlab="",ylab="", cex.axis=1.5)     
    title(main=paste("Removed:", familylabel), cex.main = 2, line = 1)
    title(ylab = "Frequency",line=2.5, cex.lab=1.15)
 #   axis(1, cex.axis=1.2)
#    axis(2, cex.axis=1.2)
    lines(D1, col="red")
    abline(v=dftemp$ERRate[1], lty = 2)
    
    if (k==1) {
      #if (SongParam == "Syllrep") {
      legend("topright",legend = c(paste(state0),paste(state1),"Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2), cex=1.9)
      pval = round(P.chisqAll,4)
      text(x = min(c(D0$x,D1$x)) + (max(c(D0$x,D1$x))-min(c(D0$x,D1$x)))*0.15, y=max(c(D0$y,D1$y))*0.65, labels = bquote(italic(p) == .(pval)), cex=2)
      title(xlab=paste("Rate of evolution of syllable repertoire"),line = 2.5, cex.lab = 1.8)

    } else {  # end if k == 1 (i.e. the "None" removed test)
      # if (SongParam == "Syllrep") {
        pval = round(P.chisqAll,4)
        text(x = min(c(D0$x,D1$x)) + (max(c(D0$x,D1$x))-min(c(D0$x,D1$x)))*0.25, y=max(c(D0$y,D1$y))*0.65,labels = bquote(italic(p) == .(pval)), cex=2)
        title(xlab=paste("Rate of evolution of", SongParam),line = 2.5, cex.lab = 1.8)
        

    } # end else (i.e. k =/= 1)
    
    #plot brownie pvalue 
    D0 <- density(dftemp$Pval)
    sdev <- sd(dftemp$Pval)
    meanphy <- mean(dftemp$Pval)
    plot(D0,col="black",
         xlim=c(min(D0$x),
                max(D0$x)),
         ylim=c(min(D0$y),
                max(D0$y)),
         main=paste(columns[1], columns[2], "Brownie pvals", ", # sims =", nsim, " \nMean =", round(meanphy,4), "/ StdDev =", round(sdev,4), otherlabel), cex.main = 0.85, xlab="Pval" ,ylab="Frequency") 
    abline(v=0.05, col = "gray")
    
    rm(dftemp)
  } # end for k in 1:length(browniejacks) - going through each dataframe of results from each removed family
  dev.off()
  }

