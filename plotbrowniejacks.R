########
#Coded by Kate T. Snyder
#Last Modified 6-22-2018
#Built using RStudio Version 1.1.453
#R Version 3.4.2
#
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#BayesTraitsV2
########

#Brownie Jack Plot


plotbrowniejacks <- function(MateParam, SongParam, browniejacks = NA, allcsvs = FALSE) {   #browniejacks is a list of dataframes, all the same Mate/Song combo, but each dataframe the test with one family removed
require(ape)
require(phytools)
  
    pdf(paste(Sys.Date(),"BrownieJacks",MateParam, SongParam,".pdf"), height = 10, width = 8)
  par(mfrow = c(6,4), cex = 0.3, mar=rep(2,4))

  if (SongParam == "Song") {
    SongParam = "SongRep"
  }
  if (MateParam == "Polygyny") {
    state0 <- "Monogamy"
    state1 <- "Polygyny"
  } else if (MateParam == "EPP") {
    state0 <- "Low EPP"
    state1 <- "High EPP"
  } else {
    state0 <- "Binary State 0"
    state1 <- "Binary State 1"
  }

  for (k in 1:length(browniejacks)) {
    dftempNotConv <- browniejacks[[k]] 
    tempfamily = dftempNotConv[1,14]
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
    
    par(mar = c(3.5,3.5,4,1))
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
      text(x = 0.05,y=100,labels = bquote(italic(p) == .(pval)), cex=2)
      title(xlab=paste("Rate of evolution of syllable repertoire"),line = 2.5, cex.lab = 1.8)
    #  } else if (SongParam == "Duration") {
      #   legend("topleft",legend = c(paste(state0),paste(state1),"Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2), cex = 1.9)
      #   pval = round(P.chisqAll,4)
      # text(x = 0.005,y=300,labels = bquote(italic(p) == .(pval)), cex=2)
      # title(xlab=paste("Rate of evolution of song duration"),line = 2.5, cex.lab = 1.8)
      # }
    } else {  # end if k == 1 (i.e. the "None" removed test)
      # if (SongParam == "Syllrep") {
        pval = round(P.chisqAll,4)
        text(x = 0.05,y=max(c(D0$y,D1$y))*0.65,labels = bquote(italic(p) == .(pval)), cex=2)
        title(xlab=paste("Rate of evolution of", SongParam),line = 2.5, cex.lab = 1.8)
        
      # } else if (SongParam == "Duration") {
      #   pval = round(P.chisqAll,4)
      #   text(x = 0.005,y=max(c(D0$y,D1$y))*0.65,labels = bquote(italic(p) == .(pval)), cex=2)
      #   title(xlab=paste("Rate of evolution of song duration"),line = 2.5, cex.lab = 1.8)
      # } # end if SongParam == Duration
    } # end else (i.e. k =/= 1)
    
    rm(dftemp)
  } # end for k in 1:length(browniejacks) - going through each dataframe of results from each removed family
  dev.off()
  }

