# Simple BayesTraits Discrete run
# Created 3/9/2022
# Kate T Snyder
# Last Edited: 3/9/2022

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
.BayesTraitsPath <- "~/Documents/BayesTraitsV2"

# # multitree
# startReadingTrees <- Sys.time()
# startReadingTrees
# #thousandtrees = read.tree(file="BirdzillaEricson10.tre")
# #thousandtrees = read.tree(file="BirdzillaHackett3_Stage2_1000trees.tre")
# thousandtrees = read.tree(file = "/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/Birdsong - Life History Evolution/BirdzillaHackett4_Stage2_1000trees.tre")
# endReadingTrees <- Sys.time()
# endReadingTrees
# timeReadingTrees <- endReadingTrees-startReadingTrees
# timeReadingTrees #Now only 1 minute! Nice.
# 
# drop.tip.multiPhylo<-function(phy, tip, ...){
#   if(!inherits(phy,"multiPhylo"))
#     stop("phy is not an object of class \"multiPhylo\".")
#   else {
#     trees<-lapply(phy,drop.tip,tip=tip,...)
#     class(trees)<-"multiPhylo"
#   }
#   trees
# }
# 
# # running overnight 3/9/22-3/10/22
# newdata = "2022-03-09CoopSong__All.csv"
# currentclassmethod = classmethods[2]
# columns = c(currentclassmethod, "FemaleSong")
# nsim = 50
# 
# 
# pdf(file = paste0(Sys.Date(),"BayesTraitsDiscrete_", currentclassmethod, "_",columns[2], "_50trees.pdf"), width = 8, height = 11)
# par(mfrow = c(3,2))
# par(oma = c(2,3,2,1))
# for (i in 1:50) {
#   temptree <- thousandtrees[[i]]
#   treelabel = paste0("Hackett4-",i,"_")
#   subset <- subsettreedata(columns = columns, newdata = newdata, newtree = temptree)
#   subsetHack <- subset$subsettree
#   simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = temptree, treelabel = treelabel, nsim = nsim)
#   plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = FALSE, ylabel = NULL)
# }
# dev.off()
# 
# treefile = "2021-08-31ConsensusPasserineTreeEricson10_1000.nex"
# nsim = 200
# treelabel = "PasserineTreeEricson-"
# currentclassmethod = classmethods[4]
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
# pdf(file = paste0(Sys.Date(),"BayesTraitsDiscrete", treelabel, columns[2], nsim, "sim_switchedTies.pdf"), width = 8, height = 11)
# par(mfrow = c(3,2))
# par(oma = c(2,3,2,1))
# for (i in 1:length(switchedClassesdf$species)) {
#   switchBird <- switchedClassesdf$species[i]
#   print(switchBird)
#   oldValue <- switchedClassesdf$originalClass[i]
#   switchedValue <- switchedClassesdf$switchedClass[i]
#   datain[which(datain$species == switchBird),currentclassmethod] <- switchedValue
#   tempY <- paste(switchBird,"CoopBreed was", oldValue, "now", switchedValue)
#   
#   simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = datain, newtree = treefile, treelabel = treelabel, nsim = nsim)
#   plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = FALSE, ylabel = tempY)
# }
# dev.off()
# 
# 
# classmethods <- c("MeanCoopOmitTies",  "MeanCoopTie2Noncoop", "MeanCoopTie2Coop",  "AnyCoopEqualsCoop" )
# nsim = 250
# for (i in classmethods) {
#   currentclassmethod <- i
#   columns <- c(currentclassmethod, "FemaleSong")
#   print(columns)
#   print(Sys.time())
#   simplebtwOutput <- simplebtwDiscrete(columns = columns, newdata = newdata, newtree = treefile, treelabel = treelabel, nsim = nsim)
#   plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOutput, nsim = nsim, treelabel = treelabel, newpdf = TRUE, ylabel = NULL)
# }


simplebtwDiscrete <- function(columns, newdata, newtree, treelabel, nsim, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, savecsvs = FALSE) {
  
  currentclassmethod <- columns[1]
  currentlabel <- paste0(treelabel,currentclassmethod)
  subsetbtw <- subsettreedata(columns = columns, newdata = newdata, newtree = newtree, skinnydata = TRUE, cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue)
  currentlabel <- paste0(treelabel,currentclassmethod)
  subsetdf <- subsetbtw$subsetdf
  colnames(subsetdf)[which(colnames(subsetdf) == currentclassmethod)] <- "CoopBreed"
  if (columns[2] == "FemaleSong") {
  subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Present")] <- "1"
  subsetdf$FemaleSong[which(subsetdf$FemaleSong == "Absent")] <- "0"
  }
  subsetdf$CoopBreed <- as.character(subsetdf$CoopBreed)
  # tipsToDrop <- thousandtrees[[1]]$tip.label[which(!thousandtrees[[1]]$tip.label %in% subsetdf$species)]
  # subsettree <- drop.tip.multiPhylo(phy = thousandtrees, tip=tipsToDrop)
  subsettree <- subsetbtw$subsettree
  #subsettree5 <- as.multiPhylo(subsettree5)
  
  simplebtwOut <- set.seed(10)
  nsim = nsim
  for (n in 1:nsim) {
    nocorrD <- Discrete(subsettree, subsetdf)
    corrD <- Discrete(subsettree, subsetdf, dependent=TRUE)
    lrtestresults <- lrtest(corrD, nocorrD)
    tempRow <- cbind(corrD, lrtestresults)
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

plotDiscreteBayes <- function(columns, simplebtwOut, nsim, treelabel, newpdf, cladesubsetvalue = NULL, ylabel = NULL) {
  
  currentclassmethod <- columns[1]
  currentlabel <- paste0(treelabel,currentclassmethod)
  
  if (newpdf == TRUE) {
    pdf(file = paste0(Sys.Date(),"BayesTraits_",currentlabel," vs ", columns[2], nsim, "sims", cladesubsetvalue,".pdf")) 
  }
  df <- dfrangetemp <- simplebtwOut
  MateParam <- columns[1]
  SongParam <- columns[2]
  ## replace plotdiscrete - from BayesPlots_choosebin
  means <- apply(X = dfrangetemp[,3:10],MARGIN = 2,FUN = mean)
  ttests <- apply(X = dfrangetemp[,3:10],MARGIN = 2,FUN = t.test)
  confInts = list()
  segmentmeans <- set.seed(10)
  segmentmins <- set.seed(10)
  segmentmaxs <- set.seed(10)
  for (k in 1:8) {
    confInts[[k]] <- ttests[[k]]$conf.int[1:2]
  }
  confIntsdf <- as.data.frame(confInts)
  names(means) <- colnames(df[,3:10])
  mins <- confIntsdf[1,]
  names(mins) <- colnames(df[,3:10])
  maxs <- confIntsdf[2,]
  names(maxs) <- colnames(df[,3:10])
  segmentmeans <- rbind(segmentmeans,means) #stores each(all) segment's mean rates
  segmentmins <- rbind(segmentmins,mins) #stores each(all) segment's lower 95CI rates
  segmentmaxs <- rbind(segmentmaxs,maxs) #stores each(all) segment's upper 95CI rates
  #} # end for (j in 1:length(minvec))
  segmentmeansdf <- as.data.frame(segmentmeans)
  segmentminsdf <- as.data.frame(segmentmins)
  segmentmaxsdf <- as.data.frame(segmentmaxs)
  arrowcols = c("red","orange","blue","green","purple")
  for (m in 1:length(segmentmeansdf$q12)) {  # goes til end
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
    
    dfrangesig <- dfrangetemp[which(dfrangetemp$pval < 0.05),]
    meannumbersig <- length(dfrangesig$Tree.No)
    
    
    if (newpdf == TRUE) {  #used when not plotting jackknifes
      par(mar = rep(2, 4))
      runsperthresh <- paste("/",nsim,sep="")
      arrowmod <- 1
      yfamilylabel <- ylabel
    } else if (newpdf == FALSE) {
      par(mar = c(1.9,1.9,2.4,1.9))
      runsperthresh <- paste("/",nsim,sep="")
      arrowmod <- 0.6
      yfamilylabel = paste(ylabel, cladesubsetvalue)
    }
    
    labx0 = paste(SongParam,"Absent")
    labx1 = paste(SongParam,"Present")
    lab0x = paste("Non-Cooperative")
    lab1x = paste("Cooperative")
    
    mat <- meanratesmat
    rates <- meanlabs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2)  
    
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
    #end plot
    
  } #end for (m in 1:length(segmentmeans$q12)), i.e. end of going through each segment
  if (newpdf == TRUE) {
  dev.off()
  }
} # end function
