## Utilize phytools::Map.Overlap() to get proportion of edges spent in each FS/Coop state
## Kate Snyder
## 4/4/2022
## Method based on Huelsenbeck et al (2003)
## 
## Edited 6/1/23 - added checkpoint save to CharacterSimmaps
## Edited 9/27/23 - changed "Dummy" simmap generation to use sim.history() with Q rates, ancestral character estimation instead of randomizing tip states; but seems to have gotten totally weird - output values odd
## Edited 9/29/23
## 10/30/2023 - added boxplot to calcHuel function; 3 ggplots now returned from calcHuel
## 11/22/2023 - added transition counts by state to CharacterSimmaps output
## 11/27/2023 - added boxplots of transition counts (real vs dummy, real vs expected) to calcHuel function; removed all non-function sections (Cycle families - jackknife, Cycle families - single family, cycle calcHuel, Misc, as well as many executions at top) - can be found in test_trait_overlap_simmaps_2023-11-27ArchiveInclNonfunctions.R
## 11/28/2023  - calcHuel: stats on real vs expected transition counts, removed "Expected" from x-axis labels in plots 5-6, named elements in returned list; CharacterSimmaps: added plotSampleSimmaps to args
## 1/9/2024 - added calculation of Nsims Actual greater than Expected to calcHuel - transition counts

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")
library(phytools)
source("subsettreedata.R")

# Get Q rates
source("findQrates.R")
# for (i in 1:5) {
#   columns = c(CBcolumns[i], "FemaleSong_Agg01")
#   subset = subsettreedata(columns = columns, newdata = dataIn, newtree = OscineTree)
#   subsettree = subset$subsettree
#   subsetdf = subset$subsetdf
#   #Qout <- findQrates(columns, plot = T, newtree = OscineTree, newdata = dataIn)
#   whichnodes = "no"  #"FS" # "Coop" # "no"
#   tipsize = 0.1
#   filename = paste(columns[1], columns[2], "fan phylo double tips", whichnodes, "nodes.pdf")
#   pdf(filename, height = 8, width = 9)
#   plot.phylo(subsettree, type = "f", show.tip.label = FALSE, align.tip.label = TRUE, cex = 0.01)
#   py = c("black","orange")
#   treetiplabels = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,columns[1]] == 1)]
#   tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=tipsize, offset = 1)
#   py2 = c("blue","red")
#   treetiplabels2 = subsettree$tip.label %in% subsetdf$species[which(subsetdf[,columns[2]] == 1)]
#   tiplabels(pch=21,bg=py2[as.numeric(treetiplabels2)+1], col = py2[as.numeric(treetiplabels2)+1], cex=tipsize, offset = 2)
#   if (whichnodes == "Coop") {
#     discretetraitvec <- subsetdf[,columns[1]]
#     names(discretetraitvec) <- subsetdf[,1]
#     pynodes = py
#     circles=ace(x=discretetraitvec,phy=subsettree,type="discrete",model="ARD")
#     nodelabels(thermo=circles$lik.anc,piecol=pynodes, height = 1.2, width = 1.2, horiz = TRUE, frame = "circle")
#   } else if (whichnodes == "FS") {
#     discretetraitvec <- subsetdf[,columns[2]]
#     names(discretetraitvec) <- subsetdf[,1]
#     pynodes = py2
#     circles=ace(x=discretetraitvec,phy=subsettree,type="discrete",model="ARD")
#     nodelabels(thermo=circles$lik.anc,piecol=pynodes, height = 1.2, width = 1.2, horiz = TRUE, frame = "circle")
#   } 
#   allpy = c(py, py2)
#   alllabs = c(paste("Not", columns[1]), columns[1], "Female Song Absent", "Female Song Present")
#   legend("bottomleft", legend = alllabs, cex = 0.9, fill=allpy, bty="n")
#   dev.off()
#   
#   print(columns)
# }


#### CharacterSimmaps fxn ----
CharacterSimmaps <- function(columns, df, tree, dummy, nsims, treelabel, datalabel = NULL, dummyMethod = c("simHistory", "makeSimmap"), plotSampleSimmaps = FALSE) {
  
  if (is.null(datalabel)) {
    datalabel = paste(columns[1], columns[2])
  }
  
  require(stringr)
  source("findQrates.R")
  cooprates <- findQrates(columns = columns[1], newdata = df, newtree = tree)
  coopQ <- cooprates$qrates
  coopQ01 <- coopQ[3]
  coopQ10 <- coopQ[2]
  coopAnc = str_remove(cooprates$ARDlikanc, "ARDlik.anc ")  # added this for sim.history()
  coopAnc = as.numeric(coopAnc)
  names(coopAnc) <- c("0","1")
  FSrates <- findQrates(columns = columns[2], newdata = df, newtree = tree)
  FSQ <- FSrates$qrates
  FSQAbsPres <- FSQ[3]
  FSQPresAbs <- FSQ[2]
  FSAnc = str_remove(FSrates$ARDlikanc, "ARDlik.anc ")  # added this for sim.history()
  FSAnc = as.numeric(FSAnc)
  names(FSAnc) <- c("0","1")
  subsets <- subsettreedata(columns = columns, newdata = df, newtree = tree)
  subsetdf <- subsets$subsetdf
  subsettree <- subsets$subsettree
  FSvec <- subsetdf[,columns[2]]
  names(FSvec) <- subsetdf$species
  Coopvec <- subsetdf[,columns[1]]
  names(Coopvec) <- subsetdf$species
  
  if (dummy == FALSE) {
    FSsimtrees <- make.simmap(tree = subsettree, x = FSvec, model = "ARD", nsim = nsims, Q = FSQ)
    Coopsimtrees <- make.simmap(tree = subsettree, x = Coopvec, model = "ARD", nsim = nsims, Q = coopQ)
    datalabel = paste(datalabel, "REAL")
    if (plotSampleSimmaps) {
      pdf(file = paste0("Simmap Overlap Outputs/",columns[1], " ", columns[2]," ",cladesubsetvalue, treelabel, " egSimmaps.pdf"),height=12,width=6)
      layout(matrix(1:6,nrow = 2,ncol=3))
      for (i in 1:3) {
        simmap1 <- Coopsimtrees[[i]]
        py = c("black","red")
        pynamed <- py
        names(pynamed) <- c(0,1)
        
        # Plot columns[1] data simmap 
        plotSimmap(simmap1,fsize=0.2, lwd = 0.8, colors = pynamed)
        #add tips
        treetiplabels <- simmap1$tip.label %in% names(Coopvec[Coopvec==1]) 
        tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
        numrates1 <- lapply(coopQ,round,digits=6)
        title(paste(" ", "\nmake.simmap","Qrates:", numrates1[2],numrates1[3], columns[1]),cex.main = 0.5)
        
        simmapQset <- FSsimtrees[[i]]
        plotSimmap(simmapQset,fsize=0.2,lwd=0.8, colors = pynamed)
        #add tips
        treetiplabels2 <- simmapQset$tip.label %in% names(FSvec[FSvec==1]) 
        tiplabels(pch=21,bg=py[as.numeric(treetiplabels2)+1], col = py[as.numeric(treetiplabels2)+1], cex=0.1)
        numrates2 <- lapply(FSQ,round,digits=6)
        title(main=paste(" ","\nARDmodel","Qrates:",numrates2[2],numrates2[3], columns[2]),cex.main = 0.5)
      }
      dev.off()
    }
  } else {
    # dummy simmaps made from one of two methods
    if (dummyMethod == "simHistory") {
      print(paste("starting Dummy Coop simmaps", Sys.time()))
      Coopsimtrees <- sim.history(tree = subsettree, Q = coopQ, nsim = nsims, anc = coopAnc)
      print(paste("starting Dummy FemSong simmaps", Sys.time()))
      FSsimtrees = sim.history(tree = subsettree, Q = FSQ, nsim = nsims, anc = FSAnc)
      datalabel <- paste(datalabel,"DUMMYSimHist-CoopFS")
    } else if (dummyMethod == "makeSimmap") {
      # Make randomized versions of CoopBreed simmaps / DUMMY data
      CoopsimtreesRand <- list()
      for (j in 1:nsims) { 
        Coopvec <- subsetdf[,columns[1]]
        CoopvecRandom <- sample(Coopvec)
        names(CoopvecRandom) <- subsetdf$species
        Coopsimtree <- make.simmap(tree = subsettree, x = CoopvecRandom, model = "ARD", nsim = 1, Q = coopQ)

        CoopsimtreesRand[[j]] <- Coopsimtree
        print(paste(j, Sys.time()))
      } # end for j in 1:nsims (Coop)
      Coopsimtrees <- CoopsimtreesRand
      
      # Make randomized versions of FemaleSong simmaps / DUMMY data
      FSsimtreesRand <- list()
      print(paste("starting Dummy FemSong simmaps", Sys.time()))
      for (j in 1:nsims) {
        FSvec <- subsetdf[,columns[2]]
        FSvecRandom <- sample(FSvec)
        names(FSvecRandom) <- subsetdf$species
        FSsimtree <- make.simmap(tree = subsettree, x = FSvecRandom, model = "ARD", nsim = 1, Q = FSQ)
        FSsimtreesRand[[j]] <- FSsimtree
        print(j)
      } # end for j in 1:nsims (FS)
      FSsimtrees<- FSsimtreesRand
      
      datalabel <- paste(datalabel,"DUMMYResampledMkSimmap-CoopFS")
    } # end if dummyMethod else
    
    
  } # end if dummy = FALSE else
  
  
  column1 <- columns[1]
  column2 <- columns[2]
  Nspecies <- length(subsetdf$species)
  
  
  # Note: output organization:
  # > Map.Overlap(FSsimtree1, Coopsimtree1)
  #          0          1
  # Absent  0.2182494 0.02578955
  # Present 0.6247627 0.13119829
  # --> Map.Overlap(returned in rows, returned in columns)
  
  dfout <- set.seed(10)
  for (i in 1:nsims) {
    
    FSsimtree1 <- FSsimtrees[[i]]
    Coopsimtree1 <- Coopsimtrees[[i]]
    
    FSdescribed <- describe.simmap(FSsimtree1)
    propFSabsent <- FSdescribed$times[2,1]
    propFSpresent <- FSdescribed$times[2,2]
    Coopdescribed <- describe.simmap(Coopsimtree1)
    propNoncoop <- Coopdescribed$times[2,1]
    propCoop <- Coopdescribed$times[2,2]
    totaltime <- FSdescribed$times[1,3]
    
    OverlapMat <- Map.Overlap(FSsimtree1, Coopsimtree1)
    ObsProp0Absent <- OverlapMat[1,1]
    ObsProp1Absent <- OverlapMat[1,2]
    ObsProp0Present <- OverlapMat[2,1]
    ObsProp1Present <- OverlapMat[2,2]
    # ChiMat <- round(OverlapMat*totaltime)
    # chiOutSim <- chisq.test(ChiMat, simulate.p.value = TRUE)
    # chiStatSim <- chiOutSim$statistic
    # chiPvalSim <- chiOutSim$p.value
    # chiOut <- chisq.test(ChiMat, simulate.p.value = FALSE)
    # chiStat <- chiOut$statistic
    # chiPval <- chiOut$p.value
    
    
    treenum <- i
    
    # temprow <- c(treenum, column1, column2, coopQ01, coopQ10, FSQAbsPres, FSQPresAbs, Nspecies, propFSabsent, propFSpresent, propNoncoop, propCoop, ObsProp0Absent, ObsProp0Present, ObsProp1Absent, ObsProp1Present, totaltime, chiStat, chiPval, chiStatSim, chiPvalSim)
    temprow <- c(treenum, column1, column2, coopQ01, coopQ10, FSQAbsPres, FSQPresAbs, Nspecies, propFSabsent, propFSpresent, propNoncoop, propCoop, ObsProp0Absent, ObsProp0Present, ObsProp1Absent, ObsProp1Present, totaltime)
    dfout <- rbind(dfout, temprow)
    dfout <- as.data.frame(dfout)
    colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
    if (!dir.exists("Simmap Overlap Outputs")) {
      dir.create("Simmap Overlap Outputs")
    }
    if (i %in% c(100, 200, 250,500,1000,2000,3000,4000,5000)) {
      write.csv(dfout, file = paste("Simmap Overlap Outputs/", datalabel, "simmap overlap_counts output nsim", nsims, treelabel,".csv"), row.names = FALSE)
    }
  }
  dfout <- as.data.frame(dfout)
  #colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime", "chiStat", "chiPval", "chiStatSim", "chiPvalSim")
  colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
  
  source("find transition counts by state for 2 Discrete traits.R")
  transStateCounts = getTransitionStateCounts(Coopsimtrees = Coopsimtrees, FSsimtrees = FSsimtrees)
  dfout = merge(dfout, transStateCounts, by.x = "treenum", by.y = "TreeNum")
  
  write.csv(dfout, file = paste("Simmap Overlap Outputs/", datalabel, "simmap overlap_counts output nsim", nsims, treelabel,".csv"), row.names = FALSE)
  return(dfout)
} # end function

#### calcHuel fxn ----
calcHuel <- function(dfout, dfDummy, nsims_real = NULL, nsims_dummy = NULL, otherlabel = NULL, newplot = TRUE, plot_ggplots_pdf = FALSE) {
  require(tidyverse)
  require(ggplot2)

  if (is.null(nsims_real)) {
    nsims_real = length(dfout[,1])
  }
  if (is.null(nsims_dummy)) {
    nsims_dummy = length(dfDummy[,1])
  }
  
  Nspecies = dfout[1,"Nspecies"]
  trait1 = dfout[1,"column1"]
  trait2 = dfout[1,"column2"]
  
  PDFname = paste("Simmap Overlap Counts",trait1, trait2, nsims_real, nsims_dummy, otherlabel)
  
  ExpPropAbsent0 <- as.numeric(dfout$propFSabsent)*as.numeric(dfout$propNoncoop)
  ExpPropAbsent1 <- as.numeric(dfout$propFSabsent)*as.numeric(dfout$propCoop)
  ExpPropPresent0 <- as.numeric(dfout$propFSpresent)*as.numeric(dfout$propNoncoop)
  ExpPropPresent1 <- as.numeric(dfout$propFSpresent)*as.numeric(dfout$propCoop)
  
  # sum(obs-expected) for each state
  dAbsent0 <- sum(abs(as.numeric(dfout$ObsProp0Absent) - ExpPropAbsent0))
  #as.numeric(dfout$ObsProp0Absent) - ExpPropAbsent0
  dAbsent1 <- sum(abs(as.numeric(dfout$ObsProp1Absent) - ExpPropAbsent1))
  #as.numeric(dfout$ObsProp1Absent) - ExpPropAbsent1
  dPresent0 <- sum(abs(as.numeric(dfout$ObsProp0Present) - ExpPropPresent0))
  #as.numeric(dfout$ObsProp0Present) - ExpPropPresent0
  dPresent1 <- sum(abs(as.numeric(dfout$ObsProp1Present) - ExpPropPresent1))
  #as.numeric(dfout$ObsProp1Present) - ExpPropPresent1
 
   sum(dAbsent0, dAbsent1, dPresent0, dPresent1) 
  D_real <- sum(dAbsent0, dAbsent1, dPresent0, dPresent1)/nsims_real
  print(paste("D_real:",D_real))
  
  # do sum-sum in other order to get distribution of Dsims 
  dAbsent0 <- abs(as.numeric(dfout$ObsProp0Absent) - ExpPropAbsent0)
  dAbsent1 <- abs(as.numeric(dfout$ObsProp1Absent) - ExpPropAbsent1)
  dPresent0 <- abs(as.numeric(dfout$ObsProp0Present) - ExpPropPresent0)
  dPresent1 <- abs(as.numeric(dfout$ObsProp1Present) - ExpPropPresent1)
  Real_dsims <- rowSums(cbind(dAbsent0, dAbsent1, dPresent0, dPresent1))
  
  # Dummy data
  ExpPropAbsent0 <- as.numeric(dfDummy$propFSabsent)*as.numeric(dfDummy$propNoncoop)
  ExpPropAbsent1 <- as.numeric(dfDummy$propFSabsent)*as.numeric(dfDummy$propCoop)
  ExpPropPresent0 <- as.numeric(dfDummy$propFSpresent)*as.numeric(dfDummy$propNoncoop)
  ExpPropPresent1 <- as.numeric(dfDummy$propFSpresent)*as.numeric(dfDummy$propCoop)
  dAbsent0 <- abs(as.numeric(dfDummy$ObsProp0Absent) - ExpPropAbsent0)
  dAbsent1 <- abs(as.numeric(dfDummy$ObsProp1Absent) - ExpPropAbsent1)
  dPresent0 <- abs(as.numeric(dfDummy$ObsProp0Present) - ExpPropPresent0)
  dPresent1 <- abs(as.numeric(dfDummy$ObsProp1Present) - ExpPropPresent1)
  Dummy_dsums <- rowSums(cbind(dAbsent0, dAbsent1, dPresent0, dPresent1))
  
  pval <- sum(Dummy_dsums > D_real)/nsims_dummy
  print(paste("num Dummy dsums > D_real:", sum(Dummy_dsums > D_real)))
  print(paste("pval:",pval))
  
  plotlabel = paste(trait1, trait2, "\nN species =", Nspecies, otherlabel)
  dummytitle = paste("Nsims =", nsims_dummy, otherlabel, "\nnum Dummy dsums > D_real:", sum(Dummy_dsums > D_real), ", pval =", pval)
  
  xmax = max(c(Real_dsims, Dummy_dsums))*1.1
  
  ### Plot hists 
  
if (newplot == TRUE) {
  font_size <- 1.0
  # Setting layout for 2 plots
  par(mfrow=c(2,1), mai=c(0.8,1.0,0.5,0.3), oma=c(2,3,2,3), 
      font.main=1, cex.main=1.25, cex.lab=font_size, cex.axis=1)
}
  
  # Histogram for Real_dsims
  hist(Real_dsims, xlim=c(0,xmax), breaks=20, 
       col=rgb(0.2,0.5,0.7,0.5), # semi-transparent blue color
       main=plotlabel, xlab="D statistic from real data simmaps", ylab="Frequency",
       border="white")
  abline(v = D_real, col="red", lwd=2.5) # thicker red line
  
  # Histogram for Dummy_dsums
  hist(Dummy_dsums, xlim=c(0,xmax), breaks=20, 
       col=rgb(0.7,0.5,0.2,0.5), # semi-transparent orange color
       main=dummytitle, xlab="D statistic from simulated independent data simmaps", ylab="Frequency",
       border="white")
  
  real_data <- data.frame(value = Real_dsims)
  dummy_data <- data.frame(value = Dummy_dsums)
  
  # ggplot histograms to return
  p1 <- ggplot(real_data, aes(x = value)) +
    geom_histogram(binwidth = xmax/20, fill = rgb(0.2, 0.5, 0.7, 0.5), color = "white") +
    geom_vline(xintercept = D_real, color = "red") +
    xlim(c(0, xmax)) +
    labs(title = plotlabel, x = "D statistic from real data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)

  # Create the histogram for Dummy_dsums
  p2 <- ggplot(dummy_data, aes(x = value)) +
    geom_histogram(binwidth = xmax/20, fill = rgb(0.7, 0.5, 0.2, 0.5), color = "white") +
    xlim(c(0, xmax)) +
    labs(title = dummytitle, x = "D statistic from simulated independent data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)
  
  # Create mutated dataframes for boxplot
  dfRealMelt <- dfout[,which(colnames(dfout) %in% colnames(dfDummy))] %>%
    gather("ObservedState", "ObservedState.prop", ObsProp0Absent:ObsProp1Present) %>%
    mutate(Which = "Real")
  
  dfDummyMelt <- dfDummy %>%
    gather("ObservedState", "ObservedState.prop", ObsProp0Absent:ObsProp1Present) %>%
    mutate(Which = "Dummy")
  
  dfCombined <- rbind(dfRealMelt, dfDummyMelt)
  
  dfCombined$ObservedState.prop = as.numeric(dfCombined$ObservedState.prop)
  
  # Create a combined label for x-axis
  dfCombined <- dfCombined %>%
    mutate(Label = paste(ObservedState, Which, sep = "\n"))
  
  # Create the boxplot
  p3 <-  ggplot(dfCombined, aes(x = Label, y = ObservedState.prop, fill = Which)) +
    geom_boxplot(outlier.shape = NA) + # Exclude outliers
    theme_minimal() +
    labs(y = "Observed State Proportion", x = "", fill = "Simulation Data") +
    scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste(trait1, trait2))
  
  if ("Coop0to1inFS0" %in% colnames(dfout) & "Coop0to1inFS0" %in% colnames(dfDummy)) { # compared rates are from Dummy data simmaps
    dfRealMeltCounts <- dfout[,which(colnames(dfout) %in% colnames(dfDummy))] %>%
      gather("Transition", "Count", Coop0to1inFS0:FS1to0inCoop1) %>%
      mutate(Which = "Real")
    
    dfDummyMeltCounts <- dfDummy %>%
      gather("Transition", "Count", Coop0to1inFS0:FS1to0inCoop1) %>%
      mutate(Which = "Dummy")
    
    dfCombinedCounts <- rbind(dfRealMeltCounts, dfDummyMeltCounts)
    
    dfCombinedCounts$Count = as.numeric(dfCombinedCounts$Count)
    
    # Create a combined label for x-axis
    dfCombinedCounts <- dfCombinedCounts %>%
      mutate(Label = paste(Transition, Which, sep = "\n"))
    
    dfCombinedCounts$Which = factor(dfCombinedCounts$Which, levels = c("Real", "Dummy"))
    dfCombinedCounts$Label[which(str_detect(dfCombinedCounts$Label, "Real"))] <- gsub("\n", "\n ", dfCombinedCounts$Label[which(str_detect(dfCombinedCounts$Label, "Real"))])
    
    p4 <-  ggplot(dfCombinedCounts, aes(x = Label, y = Count, fill = Which)) +
      geom_boxplot(outlier.shape = NA) + # Exclude outliers
      theme_minimal() +
      labs(y = "Transition Counts", x = "", fill = "Simulation Data") +
      scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5), title = element_text(size = 8)) +
      #ggtitle(paste(trait1, trait2))
      ggtitle(paste(trait1, trait2, "\nnSimsReal =", nsims_real, "    nSimsDummy =", nsims_dummy))
    
    if ("Coop0to1inFS0Expected" %in% colnames(dfout)) { # compared rates are calculated "Expected" values from e.g. (Total # transitions trait1=0 to trait1=1) * (total time in trait2 = 0)
      dfMeltCounts <- dfout %>%
        gather("Transition", "Count", c(Coop0to1inFS0:FS1to0inCoop1, Coop0to1inFS0Expected:FS1to0inCoop1Expected))
      dfMeltCounts$ObservedVsExpected = "Observed"
      dfMeltCounts$ObservedVsExpected[which(str_detect(dfMeltCounts$Transition, "Expected"))] = "Expected"
      dfMeltCounts$Count = as.numeric(dfMeltCounts$Count)
      dfMeltCounts$Count[which(dfMeltCounts$Count == 0)] = 0.001
      dfMeltCounts$Label = dfMeltCounts$Transition
      dfMeltCounts$Transition = gsub("Expected", "", dfMeltCounts$Transition)
      dfMeltCounts$logCount = log(dfMeltCounts$Count)
      
      # Perform LM analysis on log-transformed Counts
      lmLogMult = lm(formula = logCount ~ ObservedVsExpected*Transition, data = dfMeltCounts)
      logLmANOVA= anova(lmLogMult)
      pvalInteractionLog = logLmANOVA$`Pr(>F)`[3]
      print(pvalInteractionLog)
      if (pvalInteractionLog < 0.0001) {
        pvalInteractionLogLabel = "p < 0.0001"
      } else {
        pvalInteractionLogLabel = paste("p =", round(pvalInteractionLog, digits = 5))
      }
      require(emmeans)
      emmLog <- emmeans(lmLogMult, pairwise ~ ObservedVsExpected | Transition)
      contrastLog <- emmLog$contrasts
      logPairwisePostHoc = summary(contrastLog, adjust = "tukey")
      
      # Perform LM analysis on non-transformed Counts
      lmMult = lm(formula = Count ~ ObservedVsExpected*Transition, data = dfMeltCounts)
      LmANOVA = anova(lmMult)
      pvalInteraction = LmANOVA$`Pr(>F)`[3]
      print(pvalInteraction)
      if (pvalInteraction < 0.0001) {
        pvalInteractionLabel = "p < 0.0001"
      } else {
        pvalInteractionLabel = paste("p =", round(pvalInteraction, digits = 5))
      }
      require(emmeans)
      require(ggpubr)
      emm <- emmeans(lmMult, pairwise ~ ObservedVsExpected | Transition)
      contrast <- emm$contrasts
      PairwisePostHoc = summary(contrast, adjust = "tukey")
      
      # added 1/9/24
      timesActualGreaterThanExpected = dfout[,transitioncols] > dfout[,ExpectedCols]
      colnames(timesActualGreaterThanExpected) <- paste0(transitioncols)
      timesActualGreaterThanExpected = as.data.frame(timesActualGreaterThanExpected)
      NtimesActualGreaterThanExpected = colSums(timesActualGreaterThanExpected)
      # end added 1/9/2024
      
      TransitionStats = list(logLmANOVA = logLmANOVA, logPairwisePostHoc = logPairwisePostHoc, LmANOVA = LmANOVA, PairwisePostHoc = PairwisePostHoc, NtimesActualGreaterThanExpected = NtimesActualGreaterThanExpected)
      
      #all.equal(names(NtimesActualGreaterThanExpected), c("Coop0to1inFS0", "Coop1to0inFS0", "Coop0to1inFS1", "Coop1to0inFS1", "FS0to1inCoop0", "FS1to0inCoop0", "FS0to1inCoop1", "FS1to0inCoop1"))
      dfMeltCounts$Transition = factor(dfMeltCounts$Transition, levels = c("Coop0to1inFS0", "Coop1to0inFS0", "Coop0to1inFS1", "Coop1to0inFS1", "FS0to1inCoop0", "FS1to0inCoop0", "FS0to1inCoop1", "FS1to0inCoop1"))
      groupLabels = paste0(levels(dfMeltCounts$Transition), "\nNsims Actual greater than\nExpected: ", NtimesActualGreaterThanExpected, "/", nsims_real)
      
      p5 <-  ggplot(dfMeltCounts, aes(x = Transition, y = Count, fill = ObservedVsExpected)) +
        geom_boxplot(outlier.shape = NA) + # Exclude outliers
        theme_minimal() +
        labs(y = "Transition Counts", x = "", fill = "Observed/Expected") +
        scale_fill_manual(values = c("Observed" = "blue", "Expected" = "red")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), title = element_text(size = 8)) +
        ggtitle(paste(trait1, trait2, "\nnSimsObserved =", nsims_real, "    nSimsExpected =", nsims_real, "    Obs/Exp:TransCounts", pvalInteractionLabel)) +
        stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
        scale_x_discrete(labels = groupLabels)
      #ns: p > 0.05
      # *: p <= 0.05
      # **: p <= 0.01
      # ***: p <= 0.001
      # ****: p <= 0.0001
        
      p6 <-  ggplot(dfMeltCounts, aes(x = Transition, y = logCount, fill = ObservedVsExpected)) +
        geom_boxplot(outlier.shape = NA) + # Exclude outliers
        theme_minimal() +
        labs(y = "log(Transition Counts)", x = "", fill = "Observed/Expected") +
        scale_fill_manual(values = c("Observed" = "blue", "Expected" = "red")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), title = element_text(size = 8)) +
        ggtitle(paste(trait1, trait2, "\nnSimsObserved =", nsims_real, "    nSimsExpected =", nsims_real, "    Obs/Exp:TransCounts", pvalInteractionLogLabel)) +
        stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
        scale_x_discrete(labels = groupLabels)
      
      return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, filename = PDFname, TransitionStats = TransitionStats))
      
    } else {
      return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, filename = PDFname))
    }
    
  } else {
    p4 = NULL
    return(list(p1 = p1, p2 = p2, p3 = p3, filename = PDFname))
  }
  
  
  if (plot_ggplots_pdf == TRUE) {
    require(cowplot)
    pdf(file = paste("Simmap Overlap Outputs/Simmap Overlap",trait1, trait2, nsims_real, nsims_dummy, otherlabel, ".pdf"), height = 14, width = 6)
    print(plot_grid(p1, p2, p3, ncol = 1))
    dev.off()
  }
  #require(cowplot)
  #print(plot_grid(p1, p2, p3, ncol = 1))
  
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, filename = PDFname))
  
}
