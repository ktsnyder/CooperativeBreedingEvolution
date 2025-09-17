#Coded by Kate T. Snyder
#Last Modified 8-18-2021
#Built using RStudio Version 1.1.453
#R Version 4.3.1
#
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#


#plotbrownie(data = "2020-10-10CoopBreedSyllable.rep.finalbrownie500sim.csv", columns = c("Final.polygyny","Syllable.rep.final"), discreteCategoryLabels = c("Monogamy","Polygyny"), newpdf = FALSE)  #discrete category labels will be e.g. c("Monogamy","Polygyny") # default is to make a new PDF


plotbrownie <- function(data, columns, discreteCategoryLabels = c("state0","state1"), cladesubsetvalue = NULL, nsim = 500, islog = FALSE, newpdf = TRUE, otherlabel = NULL) {
  
if (is.data.frame(data)) {
  brownied <- data
  browniefile <- NULL
} else {
  datafile <- file.path("Outputs", "Brownie_outputs", data)
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
    # Generate filename base for both PDF and PNG
    filename_base <- paste0(Sys.Date(), "_", columns[1], "_", loglabel, columns[2], "_brownie", otherlabel)
    
    dir.create(file.path("Outputs", "Brownie_outputs"), recursive = T)
    
    # Create PDF output
    pdf(file = file.path("Outputs", "Brownie_outputs", paste0(filename_base, ".pdf")), width = 10, height = 5)
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
  D0_pval <- density(browniedf$Pval)
  sdev <- sd(browniedf$Pval)
  meanphy <- mean(browniedf$Pval)
  plot(D0_pval,col="black",
       xlim=c(min(D0_pval$x),
              max(D0_pval$x)),
       ylim=c(min(D0_pval$y),
              max(D0_pval$y)),
       main=paste(columns[1], columns[2], "Brownie pvals", ", # sims =", nsim, " \nMean =", round(meanphy,4), "/ StdDev =", round(sdev,4), otherlabel), cex.main = 0.75, xlab="Pval" ,ylab="Frequency") 
  abline(v=0.05, col = "gray")
  
  if (newpdf == TRUE) {
    dev.off()
    
    # Now create PNG output with same content
    png(file = file.path("Outputs", "Brownie_outputs", paste0(filename_base, ".png")), 
        width = 10*150, height = 5*150, res = 150)
    par(mar = c(4,4,2,1))
    par(mfrow = c(1,2))
    
    # Recreate the first plot
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
    
    # Add legend
    if (c(D0$x,D1$x)[which(c(D0$y,D1$y) == max(D0$y,D1$y))] > browniedf$ERRate[1]) {
      legend("topleft",legend = c(paste(state0),paste(state1),"Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2))
    } else if (c(D0$x,D1$x)[which(c(D0$y,D1$y) == max(D0$y,D1$y))] < browniedf$ERRate[1]) {
      legend("topright",legend = c(paste(state0),paste(state1), "Equal Rates"), lwd=1,col=c("blue","red", "black"), lty = c(1,1,2))
    }
    
    # Recreate the second plot
    plot(D0_pval,col="black",
         xlim=c(min(D0_pval$x),
                max(D0_pval$x)),
         ylim=c(min(D0_pval$y),
                max(D0_pval$y)),
         main=paste(columns[1], columns[2], "Brownie pvals", ", # sims =", nsim, " \nMean =", round(meanphy,4), "/ StdDev =", round(sdev,4), otherlabel), cex.main = 0.75, xlab="Pval" ,ylab="Frequency") 
    abline(v=0.05, col = "gray")
    
    dev.off()
  }
  
  
}  #end plotbrownie function
