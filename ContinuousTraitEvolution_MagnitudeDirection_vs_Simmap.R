## Get direction and magnitude of continuous character change
## For confirming Brownie results

# Next: Take amount change/branch length for each SEGMENT of the simmap
#   Each simmmap will have a distribution of chnge/segment values (roughly around 0) - take the mode/mean/median/location of the max density() for each tree; sum distributions?
# Maybe next: test significance by comparing to random tip trees. Or ask whether two distributions for each tree are different

setwd("~/Desktop/CooperativeBreedingEvolution")
nsims = 1000
state0 <- "Noncooperative"
state1 <- "Cooperative"

newdata = "2022-03-10CoopSong_All.csv"
df <- read.csv(newdata)
columns <- c("MeanCoopTie2Noncoop", "Syll.song.final")
tree <- HackConsensus <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000.nex")

songfeatures <- c("Syllable.rep.final", "Syll.song.final","Song.rep.final","Duration.final","Interval.final")

for (j in songfeatures) {
columns <- c("MeanCoopTie2Noncoop", j)
confirmbrownie(columns = columns, df = df, tree = tree, nsims = 1001, state0 = "Noncooperative", state1 = "Cooperative")
}



confirmbrownie <- function(columns, df, tree, nsims, state0=NULL, state1=NULL) {
  print(columns)
  source("findQrates.R")
  cooprates <- findQrates(columns = columns[1], newdata = df, newtree = tree)
  coopQ <- cooprates$qrates
  print(coopQ)
  
  subsetout <- subsettreedata(columns = columns, newdata = df, newtree = HackConsensus, islog = columns[2])
  subsetdf <- subsetout$subsetdf
  subsettree <- subsetout$subsettree
  contvec <- subsetdf[,columns[2]]
  names(contvec) <- subsetdf$species
  ancOut<- fastAnc(subsettree, contvec)
  ancOutdf <- as.data.frame(cbind(names(ancOut), ancOut))
  colnames(ancOutdf) <- c("label", "value")
  
  
  tipnum <- 1:length(subsettree$tip.label)
  tipsdf <- as.data.frame(cbind(tipnum, subsettree$tip.label))
  subsetdftips <- merge(subsetdf, tipsdf, by.x = "species", by.y = "V2")
  subsetdftips <- subsetdftips[,c("tipnum",columns[2])]
  colnames(subsetdftips) <- c("label", "value")
  
  contdf <- rbind(ancOutdf, subsetdftips)
  
  xDisc <- subsetdf[,columns[1]]
  names(xDisc) <- subsetdf$species
  
  simmaps <- make.simmap(tree = subsettree, x = xDisc, model = "ARD", nsim = nsims, Q = coopQ)
 # plotSimmap(simmap, fsize = 0.8)
  # tip label numbers are in order of simmap$tip.labels
  # tip.label is the same (and in same order) for both simmap and subsettree
  # in simmap$mapped.edge, the second number corresponds to tip.label number 
  # ancOut provides values at nodes only, would need to add in tip values
#  head(simmap$mapped.edge)
 # class(rownames(simmap$mapped.edge))
  
  ## Next: make below into a loop through each simmap in simmap
  
  
  require(stringr)
  require(tidyr)
  
  outdf <- set.seed(10)
  for (i in 1:length(simmaps)) {
    simmap <- simmaps[[i]]
    nodes <- str_split(rownames(simmap$mapped.edge), ",")
    simmap_edges <- as.data.frame(simmap$mapped.edge)
    simmap_edges <- cbind(rownames(simmap_edges),simmap_edges)
    colnames(simmap_edges) <- c("nodes","timeIn0","timeIn1")
    simmap_edgesdf <- separate(simmap_edges, "nodes", into = c("node1", "node2"), sep = ",")
    
    simmap_edgesdf <- merge(simmap_edgesdf, contdf, by.x = "node1", by.y = "label")
    simmap_edgesdf <- merge(simmap_edgesdf, contdf, by.x = "node2", by.y = "label", suffixes = c("_node1","_node2"), sort = FALSE)
    
    delta_contval <- as.numeric(simmap_edgesdf$value_node1) - as.numeric(simmap_edgesdf$value_node2)
    edgesdf <- cbind(simmap_edgesdf, delta_contval)
    
    changeDuring0 <- edgesdf$delta_contval*edgesdf$timeIn0
    changeDuring1 <- edgesdf$delta_contval*edgesdf$timeIn1
    edgesdf <- cbind(edgesdf, changeDuring0, changeDuring1)
    
    netChangeDuring0 <- sum(changeDuring0)
    totalChangeDuring0 <- sum(abs(changeDuring0))
    totalTimeIn0 <- sum(edgesdf$timeIn0)
    MeanChangeDuring0 <- netChangeDuring0/totalTimeIn0
    netChangeDuring1 <- sum(changeDuring1)
    totalChangeDuring1 <- sum(abs(changeDuring1))
    totalTimeIn1 <- sum(edgesdf$timeIn1)
    MeanChangeDuring1 <- netChangeDuring1/totalTimeIn1
    netChangeOverall <- netChangeDuring0+netChangeDuring1
    totalChangeOverall <- sum(abs(c(changeDuring0,changeDuring1)))
    totalTime <- totalTimeIn0+totalTimeIn1
    MeanChangeOverall <- netChangeOverall/totalTime
    
    totalPositiveChangeDuring0 <- sum(changeDuring0[which(changeDuring0 > 0)])
    totalNegativeChangeDuring0 <- sum(changeDuring0[which(changeDuring0 < 0)])
    totalPositiveChangeDuring1 <- sum(changeDuring1[which(changeDuring1 > 0)])
    totalNegativeChangeDuring1 <- sum(changeDuring1[which(changeDuring1 < 0)])
    
    row <- c(i, totalPositiveChangeDuring0, totalNegativeChangeDuring0, netChangeDuring0, totalChangeDuring0 , totalTimeIn0, MeanChangeDuring0, totalPositiveChangeDuring1, totalNegativeChangeDuring1, netChangeDuring1, totalChangeDuring1, totalTimeIn1, MeanChangeDuring1, netChangeOverall, totalChangeOverall, totalTime, MeanChangeOverall)
    
    outdf <- rbind(outdf,row)
    
  } # end for loop
  
  colnames(outdf) <- c("treenum", "totalPositiveChangeDuring0", "totalNegativeChangeDuring0", "netChangeDuring0", "totalChangeDuring0","totalTimeIn0", "MeanChangeDuring0", "totalPositiveChangeDuring1", "totalNegativeChangeDuring1", "netChangeDuring1", "totalChangeDuring1","totalTimeIn1", "MeanChangeDuring1","netChangeOverall", "totalChangeOverall","totalTime", "MeanChangeOverall")
  outdf <- as.data.frame(outdf)
  write.csv(outdf, file = paste0(columns[1],columns[2],"_RatePerState_",nsims, "nsims.csv"))
  
  
  # ttestout <- t.test(outdf$MeanChangeDuring0, outdf$MeanChangeDuring1)
  # pval <- ttestout$p.value
  
  D0 <- density(outdf$totalChangeDuring0/outdf$totalTimeIn0)
  D1 <- density(outdf$totalChangeDuring1/outdf$totalTimeIn1)
  overallRate <- outdf$totalChangeOverall[1]/outdf$totalTime[1]
  
  pdf(file = paste0(columns[1],columns[2],"_RatePerState_",nsims, "nsims.pdf"), height = 8, width = 5)
  par(mfrow=c(3,1))
  
  # plot net change in 0 and 1 per branch length
  plot(D0,col="blue",
       xlim=c(min(c(D0$x,D1$x)),
              max(c(D0$x,D1$x))),
       ylim=c(min(c(D0$y,D1$y)),
              max(c(D0$y,D1$y))),
       main=paste(columns[1], columns[2], "nsims:", nsims), 
       xlab = paste("Total change in log", columns[2],"per branch time per state"), 
       ylab = "Density",
       cex.main = 0.6)
  lines(D1, col="red")
  abline(v = overallRate, col = "black")
  print(overallRate)
  
  legend("topleft",legend = c(paste(state0),paste(state1),"Overall Rate"), lwd=1,col=c("blue","red","black"), lty = c(1,1,2))
  
  
  # plot distributions of total positive and negative change in 0 and 1
  Dpos0 <- density(outdf$totalPositiveChangeDuring0/outdf$totalTimeIn0)
  Dneg0 <- density(outdf$totalNegativeChangeDuring0/outdf$totalTimeIn0)
  Dpos1 <- density(outdf$totalPositiveChangeDuring1/outdf$totalTimeIn1)
  Dneg1 <- density(outdf$totalNegativeChangeDuring1/outdf$totalTimeIn1)
  
  plot(Dpos0,col="blue",
       xlim=c(min(c(Dpos0$x,Dpos1$x,Dneg0$x,Dneg1$x)),
              max(c(Dpos0$x,Dpos1$x,Dneg0$x,Dneg1$x))),
       ylim=c(min(c(Dpos0$y,Dpos1$y,Dneg0$y,Dneg1$y)),
              max(c(Dpos0$y,Dpos1$y,Dneg0$y,Dneg1$y))),
       #main=paste(columns[1], columns[2], "nsims:", nsims), 
       xlab = paste("Change in log", columns[2], "in pos or neg direction per branch length"), 
       ylab = "Density",
       cex.main = 0.6)
  lines(Dpos1, col="red")
  lines(Dneg0, col="blue", lty = 2)
  lines(Dneg1, col="red", lty = 2)
  
  legend("topleft",legend = c(paste(state0,"Pos"),paste(state0, "Neg"), paste(state1,"Pos"), paste(state1,"Neg")), lwd=1,col=c("blue","blue","red", "red"), lty = c(1,2,1,2))
  
  # plot total time in 0 and 1 
  
  Dtotal0 <- density(outdf$totalTimeIn0)
  Dtotal1 <- density(outdf$totalTimeIn1)
  
  plot(Dtotal0,col="blue",
       xlim=c(min(c(Dtotal0$x,Dtotal1$x)),
              max(c(Dtotal0$x,Dtotal1$x))),
       ylim=c(min(c(Dtotal0$y,Dtotal1$y)),
              max(c(Dtotal0$y,Dtotal1$y))),
       #main=paste(columns[1], columns[2], "nsims:", nsims), 
       xlab = "Total time in state", 
       ylab = "Density",
       cex.main = 0.6)
  lines(Dtotal1, col="red")
  
  legend("topleft",legend = c(paste(state0),paste(state1)), lwd=1,col=c("blue","red"), lty = c(1,1,2))
  
  dev.off()
} # end function

