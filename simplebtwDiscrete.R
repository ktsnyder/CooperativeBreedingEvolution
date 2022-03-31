# Simple BayesTraits Discrete run
# Created 3/9/2022
# Kate T Snyder
# Last Edited: 3/30/2022
# Added   res = c("q10 q00 1.2", "q11 q01 1.2") # remove
# Still hard coded to do Female Song 
# Saves nocorrD output to csv and plots this
# Added optional input nocorrDdf to plotDiscreteBayes() for separate df corrD vs nocorrD; flexible calling of columns for q12, q13, etc; doesn't automatically plot pval density

# Purposes: 
#   Make transition plots between two discrete characters without having to use btw::plotdiscrete()
#   Use multitree input
#   Put into other functions to e.g. jackknife

# setwd("~/Desktop/CooperativeBreedingEvolution")


source("subsettreedata.R")
require(phytools)
require(btw)
require(dplyr)

# newdata = "2022-03-09CoopSong__All.csv"
# newtree = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex"
.BayesTraitsPath <- "~/Documents/BayesTraitsV4.0.0-OSX/BayesTraitsV4"
.BayesTraitsPath <- "~/Documents/BayesTraitsV2"

# 
# ## Switch CoopBreed classifications of disputed birds
# datain <- read.csv(newdata)
# #switchspecies <- datain$species[which(is.na(datain$MeanCoopOmitTies) & !is.na(datain$FemaleSong) & !is.na(datain$AnyCoopEqualsCoop))]
# switchspecies <- datain$species[which(datain$SourceDiscrepancy == "1" & !is.na(datain$FemaleSong))]
# switchspeciesdf <- datain[which(datain$species %in% switchspecies),]
# orderedspecies <- switchspeciesdf[order(switchspeciesdf[,currentclassmethod]),]
# switchedClass <- set.seed(10)
# switchedClass[which(orderedspecies[,currentclassmethod] == 0)] <- 1
# switchedClass[which(orderedspecies[,currentclassmethod] == 1)] <- 0
# switchedClassesdf <- cbind(orderedspecies$species, orderedspecies[,currentclassmethod], switchedClass)
# colnames(switchedClassesdf) <- c("species", "originalClass","switchedClass")
# switchedClassesdf <- as.data.frame(switchedClassesdf)
# 
# 
source("btwDiscreteKTS.R")
classmethods <- c("MeanCoopOmitTies",  "MeanCoopTie2Noncoop", "MeanCoopTie2Coop") #,  "AnyCoopEqualsCoop" )
currentclassmethod <- classmethods[2]
nsim = 100
nsim = 14
treefile = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex"
treelabel = "PasserineTreeEricson-"
newdata = "2022-03-10CoopSong_All.csv"
df <- read.csv(newdata)
WebbFSdf <- read.csv("Webb et al 2016 Female Song Plumage Data.csv")
dfnew <- merge(df, WebbFSdf, by.x = "species", by.y = "TipLabel")
dfnew <- dfnew[which(dfnew$Female_song_score %in% c("Present","Absent")),]
secondcol <- "Female_song_score"
secondcol <- "FemaleSong"
for (i in classmethods) {
  currentclassmethod <- i
  columns <- c(currentclassmethod, secondcol)
  print(columns)
  print(Sys.time())
  simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = dfnew, newtree = treefile, treelabel = treelabel, nsim = nsim, savecsvs = TRUE, KeepBTInputFiles = TRUE)
  plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = "Old FS Data, BayesTraitsV4", arrowmod = 1)
}



simplebtwDiscrete <- function(columns, newdata, newtree, treelabel, nsim, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, savecsvs = FALSE, KeepBTInputFiles = FALSE) {
  require(btw)
  currentclassmethod <- columns[1]
  currentlabel <- paste0(treelabel,currentclassmethod)
  subsetbtw <- subsettreedata(columns = columns, newdata = newdata, newtree = newtree, skinnydata = TRUE, cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue)
  currentlabel <- paste0(treelabel,currentclassmethod)
  subsetdf <- subsetbtw$subsetdf
  colnames(subsetdf)[which(colnames(subsetdf) == currentclassmethod)] <- "CoopBreed"
  if ("Present" %in% subsetdf[,columns[2]]) {
  subsetdf[,columns[2]][which(subsetdf[,columns[2]] == "Present")] <- "1"
  subsetdf[,columns[2]][which(subsetdf[,columns[2]] == "Absent")] <- "0"
  }
  subsetdf$CoopBreed <- as.character(subsetdf$CoopBreed)
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  # tipsToDrop <- thousandtrees[[1]]$tip.label[which(!thousandtrees[[1]]$tip.label %in% subsetdf$species)]
  # subsettree <- drop.tip.multiPhylo(phy = thousandtrees, tip=tipsToDrop)
  subsettree <- subsetbtw$subsettree
  #subsettree5 <- as.multiPhylo(subsettree5)
  
  simplebtwOut <- set.seed(10)
  nsim = nsim
  for (n in 1:nsim) {
    if (n %in% c(1, 2, 5,10,25,50,100,150,200, 500, 1000, 1500, 2000)) {
      print(n)
      print(Sys.time())
    }
    nocorrD <- DiscreteKTS(subsettree, subsetdf, KeepBTInputFiles = KeepBTInputFiles)
    corrD <- DiscreteKTS(subsettree, subsetdf, dependent=TRUE, KeepBTInputFiles = KeepBTInputFiles)
    lrtestresults <- lrtest(corrD, nocorrD)
    tempRow <- cbind(corrD, lrtestresults, nocorrD)
    simplebtwOut <- rbind(simplebtwOut, tempRow)
  }
  simplebtwOut <- as.data.frame(simplebtwOut)
  if (savecsvs == TRUE) {
    write.csv(simplebtwOut, file = paste0(getwd(),"/OutputFiles/BayesTraitsDiscrete/",Sys.Date(), "BayesTraitsDiscrete_", currentlabel, " ",columns[2],nsim, "sims", cladesubsetvalue,".csv"))
  }
  
  Output <- simplebtwOut
  return(Output)
} # end function

# means <- apply(X = simplebtwOut,MARGIN = 2,FUN = mean)
# meansdf <- as.data.frame(rbind(means,means))
# #meansdf <- as.data.frame(as.matrix(means))
# pvalMed <- median(simplebtwOut$pval)
# pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs FemaleSong ", nsim, "sims.pdf"))
# plotdiscrete(meansdf[1,1:14], main = paste(currentlabel, "vs FemSong, \nnsims =",nsim, "median pval =", pvalMed))
# dev.off()




plotDiscreteBayes <- function(columns, simplebtwOut, nocorrDdf = NULL, nsim = NULL, treelabel = NULL, newpdf, cladesubsetvalue = NULL, ylabel = NULL, arrowmod = 1) {
  
  currentclassmethod <- columns[1]
  currentlabel <- paste0(treelabel,currentclassmethod)
  
  if (newpdf == TRUE) {
    pdf(file = paste0(getwd(),"/OutputFiles/BayesTraitsDiscrete/",Sys.Date(),"BayesTraits_",currentlabel," vs ", columns[2], nsim, "sims", cladesubsetvalue,".pdf"), width = 14, height = 5) 
  }
  df <- dfrangetemp <- simplebtwOut
  MateParam <- columns[1]
  SongParam <- columns[2]
  ## replace plotdiscrete - from BayesPlots_choosebin
  if (is.null(nocorrDdf)) {
    qDepColumns <- 3:10
  } else {
    qDepColumns <- c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
  }
  
  means <- apply(X = dfrangetemp[,qDepColumns],MARGIN = 2,FUN = mean)
  ttests <- apply(X = dfrangetemp[,qDepColumns],MARGIN = 2,FUN = t.test)
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
#  for (m in 1:length(segmentmeansdf$q12)) {  # goes til end
    m=1  # insted of above line
    
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
    
    
    ##### repeat rate processing for nocorrD, if present
    if (is.null(nocorrDdf)) {
        qColumns = 21:28
    } else {
      qColumns = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
      dfrangetemp = nocorrDdf
    }
    means_nocorr <- apply(X = dfrangetemp[,qColumns],MARGIN = 2,FUN = mean)
    ttests_nocorr <- apply(X = dfrangetemp[,qColumns],MARGIN = 2,FUN = t.test)
    confInts_nocorr = list()
    segmentmeans_nocorr <- set.seed(10)
    segmentmins_nocorr <- set.seed(10)
    segmentmaxs_nocorr <- set.seed(10)
    for (k in 1:8) {
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
    dfrangesig <- dfrangetemp[which(dfrangetemp$pval < 0.05),]
    meannumbersig <- length(dfrangesig$Tree.No)
    
    
    if (newpdf == TRUE) {  #used when not plotting jackknifes
      par(mar = rep(2, 4))
      par(mfrow = c(1, 3))
      runsperthresh <- paste("/",nsim,sep="")
      arrowmod <- arrowmod
      yfamilylabel <- ylabel
    } else if (newpdf == FALSE) {
      par(mar = c(1.9,1.9,2.4,1.9))
      runsperthresh <- paste("/",nsim,sep="")
      arrowmod <- 0.6
      yfamilylabel = paste(ylabel, cladesubsetvalue)
    }
    
    labx0 = paste(SongParam,"Absent")
    labx1 = paste(SongParam,"Present")
    
    if (MateParam == "Final.polygyny") {
      lab0x = paste("Monogamy")
      lab1x = paste("Polygyny")
    } else {
      lab0x = paste("Non-Cooperative")
      lab1x = paste("Cooperative")
    }
    
    ## nocorrD rates
    mat <- meanratesmat_nocorr
    rates <- meanlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2)  
    
    ## nocorrD plot  (plot 1)
    plot(c(0,100), c(0,100), type = "n", xaxt = "n", yaxt = "n", xlab = "", main=paste(currentlabel,"no-correlation model"), cex.main=1, ylab = "")
    title(ylab = yfamilylabel, line = 1)
    text(x=c(15, 85, 15, 85), y=c(80, 80, 20, 20), labels=c(paste(lab0x,"\n",labx0,sep=""), paste(lab0x,"\n",labx1,sep=""), paste(lab1x,"\n",labx0,sep=""), paste(lab1x,"\n",labx1,sep="")), cex=0.8)
    if (MateParam == "EPP") {
      rates <- rates/2
    }
    arrowcolvec <- rep(arrowcols[2],times=8)
    arrowcolvec[which(rates == 0)] <- "gray"
    arrows(x0=c(35, 65, 85, 75, 65, 35, 15, 25), y0=c(85, 75, 65, 35, 15, 25, 35, 65), x1=c(65, 35, 85, 75, 35, 65, 15, 25), y1=c(85, 75, 35, 65, 15, 25, 65, 35), lwd=rates*arrowmod/2*15, col=arrowcolvec, length = arrowmod/4) 
    
    mat <- maxratesmat_nocorr
    maxlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2) 
    mat <- minratesmat_nocorr
    minlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2) 
    
    labs = set.seed(10)
    for (y in 1:8) {
      labs[y] <- paste(meanlabs[y],"\n(",minlabs[y],", ",maxlabs[y],")", sep = "")} #end for y in 1:8
    
    text(x=c(50, 50, 94, 66, 50, 50, 6, 34), y=c(93, 67, 50, 50, 7,33, 50, 50), labels=labs, cex=0.85)
    #end plot nocorrD
    

    ## corrD rates
    mat <- meanratesmat
    rates <- meanlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2)  
    
    ## corrD plot
    plot(c(0,100), c(0,100), type = "n", xaxt = "n", yaxt = "n", xlab = "", main=paste(currentlabel,"\n# Runs significant: ", round(meannumbersig, digits = 1), runsperthresh, sep = ""), cex.main=1, ylab = "")
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
    #end plot corrD
    
    if (!is.null(dfrangetemp$pval)) {
    D0 <- density(dfrangetemp$pval)
    
    plot(D0,col="black",
         xlim=c(min(D0$x),
                max(D0$x)),
         ylim=c(min(D0$y),
                max(D0$y)),
         main=paste(columns[1], columns[2], "BayesTraits pvals", ", # sims =", nsim, " \n", otherlabel), cex.main = 0.85, xlab="Pval" ,ylab="Frequency") 
    abline(v=0.05, col = "gray")
    }

  if (newpdf == TRUE) {
  dev.off()
  }
} # end function
