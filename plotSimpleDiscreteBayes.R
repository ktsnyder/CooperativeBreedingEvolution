## Plot Discrete BayesTraits
## From simplebtwDiscrete.R
## Kate Snyder
## 7/13/2023
## Last edited: 7/21/2023 - add parameter roundDigits (labeling arrow rates)


#bayesIndFiles <- list.files(pattern = "FemaleSong_Independent", recursive = T)
#bayesDepFiles <- list.files(pattern = "FemaleSong_Dependent", recursive = T)
bayesIndFiles <- list.files(pattern = "Independent_2023-", recursive = T)
bayesDepFiles <- list.files(pattern = "Dependent_2023-", recursive = T)

bayesIndSplit = as.data.frame(str_split(bayesIndFiles, "/", simplify = T))
bayesDepSplit = as.data.frame(str_split(bayesDepFiles, "/", simplify = T))
colnames(bayesIndSplit) <- colnames(bayesDepSplit) <- c("Folder", "Filename")

loglikType <- "StonesLh"
loglikType <- "SamplingMeanLh"


pdf(file = paste0(Sys.Date(), " bayestraits AceQrates etc ", loglikType,"_scaledArrows.pdf"), width = 15, height = 5)
par(mfrow=c(1,3))
par(mar = c(4,3,3,1))
for (j in 1:length(bayesIndFiles)) {
  tempfolder = bayesIndSplit$Folder[j]
  outInddf = read.csv(bayesIndFiles[j])
  if (length(bayesDepFiles[which(bayesDepSplit$Folder == tempfolder)]) == 1) {
    outDepdf = read.csv(bayesDepFiles[which(bayesDepSplit$Folder == tempfolder)])
  } else {
    outDepdf = read.csv(bayesDepFiles[which(bayesDepSplit$Folder == tempfolder)][2])
  }
  if (length(outDepdf$Sim) > 1 & length(outInddf$Sim) > 1 &   sum(!is.na(outDepdf[,loglikType])) > 2  ) {
    plotSimpleDiscreteBayes(columns = c("MeanCoopTie2Noncoop","HighConfidence_FemaleSong"), df = outDepdf, nocorrDdf = outInddf, LhCol = loglikType, nsim = NULL, treelabel = NULL, newpdf = FALSE, otherlabel = tempfolder, ylabel = tempfolder, arrowmod = 1, roundDigits = 3)
  }
}
dev.off()

#### Top ----
plotSimpleDiscreteBayes <- function(columns, df, nocorrDdf = NULL, LhCol = NULL, nsim = NULL, treelabel = NULL, newpdf = TRUE, cladesubsetvalue = NULL, ylabel = NULL, arrowmod = 1, otherlabel = NULL, roundDigits = 2) {
  
  traitColsLabel <- paste(columns[1], columns[2])
  currentlabel <- paste(treelabel,traitColsLabel)
  
  if (newpdf == TRUE) {
    pdf(file = paste0(getwd(),"/OutputFiles/BayesTraitsDiscrete/",Sys.Date(),"BayesTraits_",currentlabel, "_", nsim, "sims", cladesubsetvalue, otherlabel, ".pdf"), width = 14, height = 5) 
  }
  
  tempdfDep <- df
  trait1 <- columns[1]
  trait2 <- columns[2]
  ## replace plotdiscrete - from BayesPlots_choosebin
  if ("Model" %in% colnames(df)) {
    if ("Dependent" %in% df$Model | "Discrete: Dependent" %in% df$Model) {
      tempdfDep = df[which(df$Model %in% c("Dependent", "Discrete: Dependent")),]
      qDepColumns <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
    } else {
      tempdfDep = df
      print("Warning: no 'Dependent' in Model column")
    }
  } else if ("q12" %in% colnames(df)) {
    tempdfDep = df
    qDepColumns <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
  } else if (is.null(nocorrDdf)) {
    qDepColumns <- 3:10
  } else {
    qDepColumns <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
    tempdfDep = df
    print("Warning: columns might not match")
  }
  if (is.null(nsim)) {
    nsim = length(tempdfDep[,1])
  }
  tempdfDep[,qDepColumns] <- apply(X = tempdfDep[,qDepColumns],MARGIN = 2,FUN = as.numeric)
  means <- apply(X = tempdfDep[,qDepColumns],MARGIN = 2,FUN = mean)
  ttests <- apply(X = tempdfDep[,qDepColumns],MARGIN = 2,FUN = t.test)
  confInts = list()
  segmentmeans <- set.seed(10)
  segmentmins <- set.seed(10)
  segmentmaxs <- set.seed(10)
  for (k in 1:8) {
    confInts[[k]] <- ttests[[k]]$conf.int[1:2]
  }
  confIntsdf <- as.data.frame(confInts)
  names(means) <- colnames(df[,qDepColumns])
  mins <- confIntsdf[1,]
  names(mins) <- colnames(df[,qDepColumns])
  maxs <- confIntsdf[2,]
  names(maxs) <- colnames(df[,qDepColumns])
  segmentmeans <- rbind(segmentmeans,means) #stores each(all) segment's mean rates
  segmentmins <- rbind(segmentmins,mins) #stores each(all) segment's lower 95CI rates
  segmentmaxs <- rbind(segmentmaxs,maxs) #stores each(all) segment's upper 95CI rates
  #} # end for (j in 1:length(minvec))
  segmentmeansdf <- as.data.frame(segmentmeans)
  segmentminsdf <- as.data.frame(segmentmins)
  segmentmaxsdf <- as.data.frame(segmentmaxs)
  arrowcols = c("red","blue","green","purple")
  
  m=1 
  model <- segmentmeansdf[m,]
  # CREATE TRANSITION RATE MATRIces - from btw function "plotdiscrete"
  meanratesmat = matrix(0, 4, 4, dimnames=list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))
  for (i in 1:4) {
    for (j in 1:4) {
      if (i != j  && sum(i,j) != 5) {
        meanratesmat[i,j] = mean(model[,grep(tail(paste("q", i, j, sep=""), 1), names(model))])
      } else meanratesmat[i,j] = NaN
    } # end for j in 1:4 (make transition matrix) mean
  } #end for i in 1:4 (make transition matrix) mean
  
  model <- segmentminsdf[m,] 
  minratesmat = matrix(0, 4, 4, dimnames=list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))
  for (i in 1:4) {
    for (j in 1:4) {
      if (i != j  && sum(i,j) != 5) {
        minratesmat[i,j] = mean(model[,grep(tail(paste("q", i, j, sep=""), 1), names(model))])
      } else minratesmat[i,j] = NaN
    } # end for j in 1:4 (make transition matrix) min
  } #end for i in 1:4 (make transition matrix) min
  
  model <- segmentmaxsdf[m,] 
  maxratesmat = matrix(0, 4, 4, dimnames=list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))
  for (i in 1:4) {
    for (j in 1:4) {
      if (i != j  && sum(i,j) != 5) {
        maxratesmat[i,j] = mean(model[,grep(tail(paste("q", i, j, sep=""), 1), names(model))])
      } else maxratesmat[i,j] = NaN
    } # end for j in 1:4 (make transition matrix) max
  } #end for i in 1:4 (make transition matrix) max
  
  
  ##### repeat rate processing for nocorrD, if present ----
  if ("alpha1" %in% colnames(df)) {
    tempdfInd = df
    qColumns <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
    tempdfInd$q12 = tempdfInd$alpha2
    tempdfInd$q13 = tempdfInd$alpha1
    tempdfInd$q21 = tempdfInd$beta2
    tempdfInd$q24 = tempdfInd$alpha1
    tempdfInd$q31 = tempdfInd$beta1
    tempdfInd$q34 = tempdfInd$alpha2
    tempdfInd$q42 = tempdfInd$beta1
    tempdfInd$q43 = tempdfInd$beta2
  } else if (!is.null(nocorrDdf)) {
    tempdfInd = nocorrDdf
    if ("alpha1" %in% colnames(tempdfInd)) {
      tempdfInd$q12 = tempdfInd$alpha2
      tempdfInd$q13 = tempdfInd$alpha1
      tempdfInd$q21 = tempdfInd$beta2
      tempdfInd$q24 = tempdfInd$alpha1
      tempdfInd$q31 = tempdfInd$beta1
      tempdfInd$q34 = tempdfInd$alpha2
      tempdfInd$q42 = tempdfInd$beta1
      tempdfInd$q43 = tempdfInd$beta2
      qColumns <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
    } else if ("q12" %in% colnames(nocorrDf)) {
      qColumns <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
    }
  }  else if ("Model" %in% colnames(df)) {
    if ("Independent" %in% df$Model | "Discrete: Independent" %in% df$Model) {
      tempdfInd = df[which(df$Model %in% c("Independent", "Discrete: Independent")),]
      qColumns <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
    }
  } else if (is.null(nocorrDdf)) {
    tempdfInd = df
    qColumns = 21:28
  } else {
    qColumns = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
    tempdfInd = nocorrDdf
  }
  tempdfInd[,qColumns] <- apply(X = tempdfInd[,qColumns], MARGIN = 2, FUN = as.numeric)
  means_nocorr <- apply(X = tempdfInd[,qColumns],MARGIN = 2,FUN = mean)
  ttests_nocorr <- apply(X = tempdfInd[,qColumns],MARGIN = 2,FUN = t.test)
  confInts_nocorr = list()
  segmentmeans_nocorr <- set.seed(10)
  segmentmins_nocorr <- set.seed(10)
  segmentmaxs_nocorr <- set.seed(10)
  for (k in 1:length(ttests_nocorr)) {
    confInts_nocorr[[k]] <- ttests_nocorr[[k]]$conf.int[1:2]
  }
  confIntsdf_nocorr <- as.data.frame(confInts_nocorr)
  names(means_nocorr) <- colnames(df[,qColumns])
  mins_nocorr <- confIntsdf_nocorr[1,]
  names(mins_nocorr) <- colnames(df[,qColumns])
  maxs_nocorr <- confIntsdf_nocorr[2,]
  names(maxs_nocorr) <- colnames(df[,qColumns])
  segmentmeans_nocorr <- rbind(segmentmeans_nocorr,means_nocorr) #stores each(all) segment's mean rates  # this isn't necessary for discrete
  segmentmins_nocorr <- rbind(segmentmins_nocorr,mins_nocorr) #stores each(all) segment's lower 95CI rates
  segmentmaxs_nocorr <- rbind(segmentmaxs_nocorr,maxs_nocorr) #stores each(all) segment's upper 95CI rates
  #} # end for (j in 1:length(minvec))
  segmentmeansdf_nocorr <- as.data.frame(segmentmeans_nocorr)
  segmentminsdf_nocorr <- as.data.frame(segmentmins_nocorr)
  segmentmaxsdf_nocorr <- as.data.frame(segmentmaxs_nocorr)
  
  
  ## get arrow mod based on max rates
  allrates = c(means, means_nocorr)
  if (max(allrates) < 1) {
    arrowmod = arrowmod*50
  } else if (max(allrates) < 5) {
    arrowmod = 20
  } else if (max(allrates) < 20) {
    arrowmod = 5
  }
  
  
  arrowcols = c("red","blue","green","purple")
  #  for (m in 1:length(segmentmeansdf$q12)) {  # goes til end
  model <- segmentmeansdf_nocorr[m,]
  # CREATE TRANSITION RATE MATRIces - from btw function "plotdiscrete"
  meanratesmat_nocorr = matrix(0, 4, 4, dimnames=list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))
  for (i in 1:4) {
    for (j in 1:4) {
      if (i != j  && sum(i,j) != 5) {
        meanratesmat_nocorr[i,j] = mean(model[,grep(tail(paste("q", i, j, sep=""), 1), names(model))])
      } else meanratesmat_nocorr[i,j] = NaN
    } # end for j in 1:4 (make transition matrix) mean
  } #end for i in 1:4 (make transition matrix) mean
  
  model <- segmentminsdf_nocorr[m,] 
  minratesmat_nocorr = matrix(0, 4, 4, dimnames=list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))
  for (i in 1:4) {
    for (j in 1:4) {
      if (i != j  && sum(i,j) != 5) {
        minratesmat_nocorr[i,j] = mean(model[,grep(tail(paste("q", i, j, sep=""), 1), names(model))])
      } else minratesmat_nocorr[i,j] = NaN
    } # end for j in 1:4 (make transition matrix) min
  } #end for i in 1:4 (make transition matrix) min
  
  model <- segmentmaxsdf_nocorr[m,] 
  maxratesmat_nocorr = matrix(0, 4, 4, dimnames=list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))
  for (i in 1:4) {
    for (j in 1:4) {
      if (i != j  && sum(i,j) != 5) {
        maxratesmat_nocorr[i,j] = mean(model[,grep(tail(paste("q", i, j, sep=""), 1), names(model))])
      } else maxratesmat_nocorr[i,j] = NaN
    } # end for j in 1:4 (make transition matrix) max
  } #end for i in 1:4 (make transition matrix) max
  
  
  # calculate number runs significant
  if ("pval" %in% colnames(df)) {
    dfrangesig <- df[which(df$pval < 0.05),]
    meannumbersig <- length(dfrangesig$Tree.No)
    LRplot = FALSE
  } else if (!is.null(LhCol)) {
    LRplot = TRUE
  }
  
  if (newpdf == TRUE) {  #used when not plotting jackknifes
    par(mar = rep(2, 4))
    par(mfrow = c(1, 3))
    runsperthresh <- paste("/",nsim,sep="")
    arrowmod <- arrowmod
    yfamilylabel <- ylabel
  } else if (newpdf == FALSE) {
    par(mar = c(1.9,1.9,2.4,1.9))
    runsperthresh <- paste("/",nsim,sep="")
    arrowmod <- arrowmod*0.6
    yfamilylabel = paste(ylabel, cladesubsetvalue)
  }
  
  labx0 = paste(trait2,"Absent")
  labx1 = paste(trait2,"Present")
  
  if (trait1 == "Final.polygyny") {
    lab0x = paste("Monogamy")
    lab1x = paste("Polygyny")
  } else if (grepl("coop", trait1, ignore.case = T)) {  #(str_detect(trait1, "coop")) {
    lab0x = paste("Non-Cooperative")
    lab1x = paste("Cooperative")
  }
  
  ## nocorrD rates
  mat <- meanratesmat_nocorr
  rates <- meanlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),roundDigits)  
  
  ## nocorrD plot  (plot 1)
  plot(c(0,100), c(0,100), type = "n", xaxt = "n", yaxt = "n", xlab = "", main=paste(currentlabel,"Independent Model"), cex.main=1, ylab = "")
  title(ylab = yfamilylabel, line = 1)
  text(x=c(15, 85, 15, 85), y=c(80, 80, 20, 20), labels=c(paste(lab0x,"\n",labx0,sep=""), paste(lab0x,"\n",labx1,sep=""), paste(lab1x,"\n",labx0,sep=""), paste(lab1x,"\n",labx1,sep="")), cex=0.8)
  if (trait1 == "EPP") {
    rates <- rates/2
  }
  arrowcolvec <- rep(arrowcols[2],times=8)
  arrowcolvec[which(rates == 0)] <- "gray"
    arrows(x0=c(35, 65, 85, 75, 65, 35, 15, 25), y0=c(85, 75, 65, 35, 15, 25, 35, 65), x1=c(65, 35, 85, 75, 35, 65, 15, 25), y1=c(85, 75, 35, 65, 15, 25, 65, 35), lwd=rates*arrowmod, col=arrowcolvec, length = 0.7) 
    
    mat <- maxratesmat_nocorr
    maxlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),roundDigits) 
    mat <- minratesmat_nocorr
    minlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),roundDigits) 
    
    labs = set.seed(10)
    for (y in 1:8) {
      labs[y] <- paste(meanlabs[y],"\n(",minlabs[y],", ",maxlabs[y],")", sep = "")} #end for y in 1:8
    
    text(x=c(50, 50, 94, 66, 50, 50, 6, 34), y=c(93, 67, 50, 50, 7,33, 50, 50), labels=labs, cex=0.85)
    #end plot nocorrD
    
    
    ## corrD rates
    mat <- meanratesmat
    rates <- meanlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),roundDigits)  
    
    ## corrD plot
    plot(c(0,100), c(0,100), type = "n", xaxt = "n", yaxt = "n", xlab = "", main = paste(columns[1], columns[2], "Dependent Model"), cex.main=1, ylab = "")
         #main=paste(currentlabel,"\n# Runs significant: ", round(meannumbersig, digits = 1), runsperthresh, sep = ""), cex.main=1, ylab = "")
    title(ylab = yfamilylabel, line = 1)
    text(x=c(15, 85, 15, 85), y=c(80, 80, 20, 20), labels=c(paste(lab0x,"\n",labx0,sep=""), paste(lab0x,"\n",labx1,sep=""), paste(lab1x,"\n",labx0,sep=""), paste(lab1x,"\n",labx1,sep="")), cex=0.8)
    if (trait1 == "EPP") {
      rates <- rates/2
    }
    arrowcolvec <- rep(arrowcols[m],times=8)
    arrowcolvec[which(rates == 0)] <- "gray"
      arrows(x0=c(35, 65, 85, 75, 65, 35, 15, 25), y0=c(85, 75, 65, 35, 15, 25, 35, 65), x1=c(65, 35, 85, 75, 35, 65, 15, 25), y1=c(85, 75, 35, 65, 15, 25, 65, 35), lwd=rates*arrowmod, col=arrowcolvec, length = 0.7) 
      
      mat <- maxratesmat
      maxlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),roundDigits) 
      mat <- minratesmat
      minlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),roundDigits) 
      
      labs = set.seed(10)
      for (y in 1:8) {
        labs[y] <- paste(meanlabs[y],"\n(",minlabs[y],", ",maxlabs[y],")", sep = "")} #end for y in 1:8
      
      text(x=c(50, 50, 94, 66, 50, 50, 6, 34), y=c(93, 67, 50, 50, 7,33, 50, 50), labels=labs, cex=0.85)
      #end plot corrD
      
      if (!is.null(df$pval)) {
        D0 <- density(df$pval)
        
        plot(D0,col="black",
             xlim=c(min(D0$x),
                    max(D0$x)),
             ylim=c(min(D0$y),
                    max(D0$y)),
             main=paste(columns[1], columns[2], "BayesTraits pvals", ", # sims =", nsim, " \n", otherlabel), cex.main = 0.85, xlab="Pval" ,ylab="Frequency") 
        abline(v=0.05, col = "gray")
      } else if (LRplot == TRUE) {
        D0 <- density(as.numeric(tempdfDep[,LhCol]), na.rm = T)
        D1 <- density(as.numeric(tempdfInd[,LhCol]), na.rm = T)
        MeanLhDep = round(mean(as.numeric(tempdfDep[,LhCol]), na.rm = T),2)
        MeanLhInd = round(mean(as.numeric(tempdfInd[,LhCol]), na.rm = T),2)
        BayesFactor = round(2*(MeanLhDep-MeanLhInd),2)
        plot(D0,col=arrowcols[m],
             xlim=c(min(c(D0$x,D1$x)),
                    max(c(D0$x,D1$x))),
             ylim=c(min(c(D0$y,D1$y)),
                    max(c(D0$y,D1$y))),
             main=paste(columns[1], columns[2], "BayesTraits", LhCol, ", # sims =", nsim, " \nMean Independent:", MeanLhInd, "Mean Dependent:", MeanLhDep, "BayesFactor:", BayesFactor, otherlabel), cex.main = 0.85, xlab=LhCol ,ylab="Frequency") 
        lines(D1, col = arrowcols[2])
      }
      
      if (newpdf == TRUE) {
        dev.off()
      }
} # end function
