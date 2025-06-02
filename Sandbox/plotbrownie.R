#Coded by Kate T. Snyder
#Last Modified 8-18-2021
#Built using RStudio Version 1.1.453
#R Version 3.5.2?
#
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#
#
#5/12/2020- copied bulk of below from browniefunction.R. Some plotting issues, not sure what to do about Pvals
#5/14/2020 - add mean log likelihood-based pval calculation, added to title. Next: fix axes
#8/20/2020 - figured out
#10/9/2020 - reordered loglabel 
#8/18/2021 - commented out example, added otherlabel arg
#8/25/2021 - add otherlabel to y axis

#plotbrownie(data = "2020-10-10CoopBreedSyllable.rep.finalbrownie500sim.csv", columns = c("Final.polygyny","Syllable.rep.final"), discreteCategoryLabels = c("Monogamy","Polygyny"), newpdf = FALSE)  #discrete category labels will be e.g. c("Monogamy","Polygyny") # default is to make a new PDF


plotbrownie <- function(data, columns, discreteCategoryLabels = c("state0","state1"), cladesubsetvalue = NULL, nsim = 500, islog = FALSE, newpdf = TRUE, otherlabel = NULL) {
  
if (is.data.frame(data)) {
  brownied <- data
  browniefile <- NULL
} else {
  datafile <- paste0(getwd(), "/OutputFiles/", data)
  brownied <- as.data.frame(read.csv(datafile, stringsAsFactors = FALSE))
  browniefile <- data
}
  
  #subset results for plotting
  browniedf <- brownied[brownied$convergence == "Optimization has converged.",]
  browniedf <- browniedf[!is.na(browniedf$convergence),]
  
  #calculate overall mean pval
  ERloglikmean <- mean(browniedf$ERloglik)
  ARDloglikmean <- mean(browniedf$ARDloglik)
  ERARDPval = round(pchisq(2*(ARDloglikmean-ERloglikmean),1,lower.tail=FALSE), digits = 3) #testing whether the two rates of continuous trait evolution are significantly different
  
  if (islog == FALSE) {   #make text for title/axis if log
    loglabel = NULL
  } else if (islog == TRUE) {
    loglabel = "log"
  } else {
    loglabel = islog
  }
  
##### plot brownie distribution
  
  if (newpdf == TRUE) {
    pdf(file = paste0(getwd(),"/OutputFiles/", Sys.Date(),columns[1], loglabel, columns[2], otherlabel, "brownie.pdf"), width = 10, height = 5)
    par(mar = c(4,4,2,1))
    par(mfrow = c(1,2)) 
  } else {
    par(mar = c(2,3,2,1))
  }
  
  
  titlelabel <- paste(columns[1], loglabel, columns[2], "all converged runs from \n", browniefile)
  
  D0 <- density(browniedf$ARDRate0)
  D1 <- density(browniedf$ARDRate1)
  
  plot(D0,col="blue",
       xlim=c(min(c(D0$x,D1$x)),
              max(c(D0$x,D1$x))),
       ylim=c(min(c(D0$y,D1$y)),
              max(c(D0$y,D1$y))),
       main=titlelabel, 
       xlab = "", 
       ylab = "",
       cex.main = 0.6)
  lines(D1, col="red")
  abline(v=browniedf$ERRate[1], lty = 2)
  title(xlab=paste("Rate of", loglabel, columns[2],"evolution"),
        ylab= paste("Number of Observations", otherlabel), line = 2)
  
  state0 <- discreteCategoryLabels[1]
  state1 <- discreteCategoryLabels[2]
  
  #legend
  if (c(D0$x,D1$x)[which(c(D0$y,D1$y) == max(D0$y,D1$y))] > browniedf$ERRate[1]) {
    legend("topleft",legend = c(paste(state0),paste(state1),"Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2))
  } else if (c(D0$x,D1$x)[which(c(D0$y,D1$y) == max(D0$y,D1$y))] < browniedf$ERRate[1]) {
    legend("topright",legend = c(paste(state0),paste(state1), "Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2))
  }

  
  #plot brownie pvalue 
  D0 <- density(browniedf$Pval)
  sdev <- sd(browniedf$Pval)
  meanphy <- mean(browniedf$Pval)
  plot(D0,col="black",
       xlim=c(min(D0$x),
              max(D0$x)),
       ylim=c(min(D0$y),
              max(D0$y)),
       main=paste(columns[1], columns[2], "Brownie pvals", ", # sims =", nsim, " \nMean =", round(meanphy,4), "/ StdDev =", round(sdev,4), otherlabel), cex.main = 0.75, xlab="Pval" ,ylab="Frequency") 
  abline(v=0.05, col = "gray")
  
  if (newpdf == TRUE) {
    dev.off()
  }
  
  # writeLines(paste0("OneRate=", round(ERloglikmean, digits = 4)))
  # writeLines(paste0("TwoRates=", round(ARDloglikmean, digits = 4)))
  # writeLines(paste0("pVal", ERARDPval))
  # writeLines(paste("",sep="\n\n"))
  # writeLines(paste("",sep="\n\n"))  

  
}  #end plotbrownie function


####Cristina's
##

# browniedata <- data.frame(Pval=numeric(nsim), ERRate=numeric(nsim), ERloglik=numeric(nsim),
#                           ERace=numeric(nsim), ARDloglik=numeric(nsim),ARDace=numeric(nsim),
#                           stringsAsFactors = FALSE)
# 
# 
# 
# BrowniePlotRates <- function(dataset, title="BrowniePlot", col=c("blue","red"), Groups=c("Rate0", "Rate1"),
#                              Xlim =c(min(XMins),
#                                      max(XMaxes)) ){
#   #data cleaning
#   dataset<- dataset[!is.na(dataset$Pval),]
#   dataset<- dataset[dataset$convergence == "Optimization has converged.",]
#   
#   #plots
#   ARDs <- length(grep("ARDRate", colnames(dataset)))
#   D <- list()
#   YMaxes <- numeric(length=ARDs)
#   YMins <- numeric(length=ARDs)
#   XMaxes <- numeric(length=ARDs)
#   XMins <- numeric(length=ARDs)
#   for(i in 1:ARDs){
#     D[[i]] <- density(dataset[,paste0('ARDRate', i-1)])
#     YMaxes[i] <- max(D[[i]]$y)
#     YMins[i] <- min(D[[i]]$y)
#     XMaxes[i] <- max(D[[i]]$x)
#     XMins[i] <- min(D[[i]]$x)
#   }
#   plot(D[[1]],col=col[1],
#        xlim=Xlim,
#        ylim=c(min(YMins),
#               max(YMaxes)),
#        main=title, xlab="Rates", font.lab=2,
#        ylab="Number of Observations",
#        cex.main=1, lwd=2)
#   for(i in 2:ARDs){
#     lines(D[[i]], col=col[i], lwd=2)
#   }
#   abline(v=dataset$ERRate[1])
#   legend("topright", legend=Groups, col=col, lty=1, lwd=2)
#   #stats
#   MeanArd <- mean(dataset$ARDloglik)
#   MeanER <- mean(dataset$ERloglik)
#   pval <- round(pchisq(2*(MeanArd - MeanER),1,lower.tail=FALSE),digits=3)
#   
# 
#   
#   ifelse(pval == 0,pval <- "<0.001", pval <- paste0("=",pval))
#   writeLines(title)
#   writeLines(paste0("OneRate=", round(MeanER, digits = 4)))
#   writeLines(paste0("TwoRates=", round(MeanArd, digits = 4)))
#   writeLines(paste0("pVal", pval))
#   writeLines(paste("",sep="\n\n"))
#   writeLines(paste("",sep="\n\n"))
# }