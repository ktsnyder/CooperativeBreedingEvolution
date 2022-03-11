


plotBTjacks <- function(MateParam, SongParam, d, sigonly = FALSE, familysplit = NULL, nsim = 10) {
  if (sigonly == TRUE) {
    sigtext = "Sig"
  } else if (sigonly == FALSE) {
    sigtext <- "ALL"
  }
  #  pdf(paste(Sys.Date(),MateParam,SongParam,sigtext,"BayesJacks4plots.pdf",sep=""), height = 10, width = 8)
  #  par(mfrow = c(6,4), cex = 0.3)
  #  layout(matrix(1:(length(csvfilevec)*4), nrow = length(csvfilevec), ncol = 4, byrow=TRUE))
  dfforminmax <- d
  if (SongParam == "Interval") {
    dfforminmax$songcontvec[which(dfforminmax$songcontvec == 0.01)] <- 0.25
  }
  dfforminmax <- dfforminmax[order(dfforminmax$songcontvec),]
  alldfvals <- unique(dfforminmax$songcontvec)
  #print(alldfvals)
  
  minvec = set.seed(10)
  maxvec = set.seed(10)
  minvec[1] = round(min(alldfvals),2)
  minvec[2] <- alldfvals[ceiling(length(alldfvals)/3)]
  maxvec[1] <- alldfvals[ceiling(length(alldfvals)/3)]
  minvec[3] <- alldfvals[ceiling(length(alldfvals)*2/3)]
  maxvec[2] <- alldfvals[ceiling(length(alldfvals)*2/3)]
  maxvec[3] <- max(alldfvals)
  print(minvec)
  print(maxvec)
  # for (k in 1:length(csvfilevec)) {
  #   csvfiletemp <- csvfilevec[k]
  transitionplots(MateParam, SongParam, df = dfforminmax,newpdf = FALSE, sigonly = sigonly, minvec = minvec, maxvec = maxvec, nsim=nsim)
  # }
  # dev.off()
}



BTbinplotsfunc <- function(MateParam, SongParam, csvfile) { 
thisdf <- as.data.frame(read.csv(csvfile))
pdf(file = paste0(MateParam,SongParam,"BTplotsBins.pdf"), height = 35, width = 7)
layout(matrix(c(1,2,3,3,0,0,4,5,6,7,0,0,8,9,10,11,12,12,0,0,13,14,15,16,17,18),13,2, byrow = TRUE),widths= c(3,3),heights=c(3,3.5,0.5,3,3,0.5,3,3,3.5,0.5,3,3,3))
  for (i in c(2,3,4,5)) {
    transitionBinplots(MateParam = MateParam, SongParam = SongParam, df = thisdf, newpdf = FALSE, nsim = 100, binnum = i)
  }
dev.off()
} # end function



############### arrow plots AND pval vs threshold plot
transitionBinplots <- function(MateParam,SongParam,df = df, newpdf = TRUE, sigonly = FALSE, minvec = NULL, maxvec = NULL, onetransplot = FALSE, familysplit = NULL, nsim = 100, binnum = NULL) {

  
if (newpdf == TRUE) {
  pdf(file = paste(Sys.Date(),MateParam,SongParam,nsim,"BayesTraits.pdf"), height = 10, width = 10)
  par(mfrow = c(2,2))
}

if (SongParam == "Interval") {
  df$songcontvec[which(df$songcontvec == 0.01)] <- 0.25
}

df <- df[order(df$songcontvec),]
alldfvals <- unique(df$songcontvec)

if (sigonly == FALSE) {
dfallsig <- df 
} else if (sigonly == TRUE) {
  dfallsig <- df[df$LRpval < 0.05,]
}

songcontvec <- dfallsig$songcontvec
LRpval <- dfallsig$LRpval
if (onetransplot == TRUE) {
  means <- apply(X = dfallsig[,4:11],MARGIN = 2,FUN = mean)
  ttests <- apply(X = dfallsig[,4:11],MARGIN = 2,FUN = t.test)
  pdf(file=paste(Sys.Date(),MateParam,SongParam,"Bayes1plot.pdf",sep=""), width = 5, height = 5)
  labx0 = paste("Low",SongParam)
  labx1 = paste("High",SongParam)
  if (MateParam == "OC") {
    lab0x = "Song-stable"
    lab1x = "Song-plastic"
  } #end if MateParam == "OC"
  else if (MateParam == "Polygyny") {
    lab0x = "Monogamous"
    lab1x = "Polygynous"
  } else if (MateParam == "EPP") {#end if MateParam == "Polygyny"
    lab0x = "Low EPP"
    lab1x = "High EPP"
  } else  { #if (MateParam == "CoopBreed") { #end if MateParam == "EPP"
    lab0x = "Noncooperative"
    lab1x = "Cooperative"
    } 
}

labx0 = "Low"
labx1 = "High"
if (MateParam == "OC") {
  lab0x = "Song-stable"
  lab1x = "Song-plastic"
} #end if MateParam == "OC"
else if (MateParam == "Polygyny") {
  lab0x = "Monogamous"
  lab1x = "Polygynous"
} else if (MateParam == "EPP") {#end if MateParam == "Polygyny"
  lab0x = "Low EPP"
  lab1x = "High EPP"
} else { #if (MateParam == "CoopBreed") { #end if MateParam == "EPP"
lab0x = "Noncooperative"
lab1x = "Cooperative"
} 


songvals <- unique(dfallsig$songcontvec)

  minvec = set.seed(10)
  maxvec = set.seed(10)
  minvec[1] = round(min(alldfvals),2)
  minvec[2] <- alldfvals[ceiling(length(alldfvals)/binnum)]
  maxvec[1] <- alldfvals[ceiling(length(alldfvals)/binnum)]
  if (binnum==2) {
    maxvec[2] <- max(alldfvals)
  } 
  if (binnum > 2) {
  minvec[3] <- alldfvals[ceiling(length(alldfvals)*2/binnum)]
  maxvec[2] <- alldfvals[ceiling(length(alldfvals)*2/binnum)]
    if (binnum == 3) {
      maxvec[3] <- max(alldfvals)
    }
  } 
  if (binnum > 3) {
    minvec[4] <- alldfvals[ceiling(length(alldfvals)*3/binnum)]
    maxvec[3] <- alldfvals[ceiling(length(alldfvals)*3/binnum)]
    if (binnum==4) {
      maxvec[4] <- max(alldfvals)
    }
  }
  if (binnum > 4) {
    minvec[5] <- alldfvals[ceiling(length(alldfvals)*4/binnum)]
    maxvec[4] <- alldfvals[ceiling(length(alldfvals)*4/binnum)]
    if (binnum==5) {
      maxvec[5] <- max(alldfvals)
    }
  }  
  # maxvec[3] <- max(alldfvals)  #this "threshold" value is not included in calculations because it is not actually a threshold separating two groups. Since the original bayestraits function set the threshold and all lower numbers as "0" and only higher numbers as "1", the highest value of threshold in the plots is not actually comparing two groups.

segmentmeans <- set.seed(10)
segmentmins <- set.seed(10)
segmentmaxs <- set.seed(10)

#FIRST for loop, make the transition matrices for each segment
#Within loop, makes 3 transition matrices - lower 95CI, upper 95CI, mean
for (j in 1:length(minvec)) { #when doing 3rds, this is 1:3
  dfrangetemp <- dfallsig[which(dfallsig$songcontvec >= minvec[j] & dfallsig$songcontvec < maxvec[j]),]
  
  numberofthreshs <- length(unique(dfrangetemp$songcontvec))
  meannumbersig <- length(dfrangetemp$Tree.No)/numberofthreshs
  means <- apply(X = dfrangetemp[,4:11],MARGIN = 2,FUN = mean)
  ttests <- apply(X = dfrangetemp[,4:11],MARGIN = 2,FUN = t.test)
  confInts = list()
  for (k in 1:8) {
    confInts[[k]] <- ttests[[k]]$conf.int[1:2]
  }
  confIntsdf <- as.data.frame(confInts)
  names(means) <- colnames(df[,4:11])
  mins <- confIntsdf[1,]
  names(mins) <- colnames(df[,4:11])
  maxs <- confIntsdf[2,]
  names(maxs) <- colnames(df[,4:11])
  segmentmeans <- rbind(segmentmeans,means) #stores each(all) segment's mean rates
  segmentmins <- rbind(segmentmins,mins) #stores each(all) segment's lower 95CI rates
  segmentmaxs <- rbind(segmentmaxs,maxs) #stores each(all) segment's upper 95CI rates
} # end for (j in 1:length(minvec))
segmentmeansdf <- as.data.frame(segmentmeans)
segmentminsdf <- as.data.frame(segmentmins)
segmentmaxsdf <- as.data.frame(segmentmaxs)
arrowcols = c("red","orange","blue","green","purple")
for (m in 1:length(segmentmeansdf$q12)) {
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
  
  #Now you have 3 transition matrices, mean, min (lower 95% CI), max (upper 95 CI) within one segment
  
  ##This part repeated from earlier for loop so the plot titles have the right values
  dfrangetemp <- dfallsig[which(dfallsig$songcontvec >= minvec[m] & dfallsig$songcontvec < maxvec[m]),]
  numberofthreshs <- length(unique(dfrangetemp$songcontvec))
  dfrangesig <- dfrangetemp[which(dfrangetemp$LRpval < 0.05),]
  meannumbersig <- length(dfrangesig$LRpval)/numberofthreshs
  ##end repeated code
  if (SongParam == "Song") {
    SongParam = "SongRep"
  }
  
  yfamilylabel <- ""
  if (newpdf == TRUE) {  #used when not plotting jackknifes
    par(mar = rep(2, 4))
    runsperthresh <- paste("/",nsim,sep="")
    arrowmod <- 1
  } else if (newpdf == FALSE) {
    par(mar = c(1.9,1.9,2.4,1.9))
    runsperthresh <- paste("/",nsim,sep="")
    arrowmod <- 0.5
  }

    labx0 = paste("Lower",SongParam)
    labx1 = paste("Higher",SongParam)

  runsperthresh <- paste("/",nsim,sep="")
  
  #PLOT (arrows taken from btw function plotdiscrete)
  mat <- meanratesmat
  rates <- meanlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2)  
 
  plot(c(0,100), c(0,100), type = "n", xaxt = "n", yaxt = "n", xlab = "", main=paste("Thresholds >=", round(minvec[m], digits = 3)," - <",round(maxvec[m], digits = 3), " (", numberofthreshs,")", "\nMean # Runs significant: ", round(meannumbersig, digits = 1), runsperthresh, sep = ""), cex.main=1, ylab = "")
  title(ylab = yfamilylabel, line = 1)
  text(x=c(15, 85, 15, 85), y=c(80, 80, 20, 20), labels=c(paste(lab0x,"\n",labx0,sep=""), paste(lab0x,"\n",labx1,sep=""), paste(lab1x,"\n",labx0,sep=""), paste(lab1x,"\n",labx1,sep="")), cex=0.8)
  if (MateParam == "EPP") {
    rates <- rates/2
  }
  arrowcolvec <- rep(arrowcols[m],times=8)
  arrowcolvec[which(rates == 0)] <- "gray"
  arrows(x0=c(35, 65, 85, 75, 65, 35, 15, 25), y0=c(85, 75, 65, 35, 15, 25, 35, 65), x1=c(65, 35, 85, 75, 35, 65, 15, 25), y1=c(85, 75, 35, 65, 15, 25, 65, 35), lwd=rates*arrowmod/2*15, col=arrowcolvec, length = arrowmod/5) 
  
  mat <- maxratesmat
  maxlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2) 


  mat <- minratesmat
  minlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2) 
  
  labs = set.seed(10)
  for (y in 1:8) {
    labs[y] <- paste(meanlabs[y],"\n(",minlabs[y],", ",maxlabs[y],")", sep = "")} #end for y in 1:8
  
  text(x=c(50, 50, 94, 66, 50, 50, 6, 34), y=c(93, 67, 50, 50, 7,33, 50, 50), labels=labs, cex=0.85)
  #end plot
  
} #end for (m in 1:length(segmentmeans$q12)), i.e. end of going through each segment

MateSongParams = NULL
par(mar=c(4,3,2,0.5),cex.axis = 1, cex.lab = 1)
if (newpdf == TRUE) {
  par(mar=c(4,3,2,2), cex.axis = 1, cex.lab = 1)
  MateSongParams <- paste(MateParam,SongParam)
} #end if newpdf = TRUE

#plot pvals vs thresholds
with(df, plot(songcontvec, LRpval, pch=20, 
              main = paste(MateSongParams, "Original threshold range: ", round(min(alldfvals),4)," - ",round(max(alldfvals),3)," (",length(alldfvals)-1,")", sep = ""), xlim = c(min(songvals)-0.1*min(songvals),max(songvals)+0.1*max(songvals)), log = "x", xlab = "", ylab="", cex.main = 0.7
))
title(ylab = "p-value", line = 2)
title(xlab = paste("High/Low Threshold:",SongParam), line = 2.2)

segments(x0 = minvec, x1 = maxvec, y0 = rep(-0.03,times=length(minvec)), y1=rep(-0.03,times=length(minvec)), lwd = 5, col = c("red","orange","blue","green","purple"))
abline(h=0.05,col="blue")
if (newpdf == TRUE) {
dev.off()
  } # end if newpdf == TRUE
}