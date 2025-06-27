## Utilize phytools::Map.Overlap() to get proportion of edges spent in each FS/Coop state
## Kate Snyder
## 4/4/2022
## Method based on Huelsenbeck et al (2003)
## Last modified 6/19/2024
## Last edited 12/11/2024 - calcHuelflex - changed columns made as.numeric to start with numericColumnStart based on the position of the last column containing "trait"
## Last edited 2/3/2025 - add "setQratesData" arg to CharacterSimmaps, to work with setQratesTree to allow Q rates to be calculated for all characters using a specific tree and data subset. NOT to be used with columnForGlobalQ and columnGlobalQrates


library(phytools)
source("subsettreedata.R")
source("findQrates.R")
source("find transition counts by state for 2 Discrete traits.R")


#### CharacterSimmaps fxn ----
CharacterSimmaps <- function(columns, df, tree, dummy, nsims, treelabel, datalabel = NULL, dummyMethod = c("simHistory", "makeSimmap"), plotSampleSimmaps = FALSE, cladesubsetvalue = NULL, setQratesTree = NULL, setQratesData = NULL, columnForGlobalQ = NULL, columnGlobalQrates = NULL) {
  
  if (is.null(datalabel)) {
    datalabel = paste(columns[1], columns[2])
  }
  
  if (is.null(setQratesTree)) {
    Qtree = tree
  } else if (is.character(setQratesTree)) {
    Qtree = read.nexus(setQratesTree)
  } else {
    Qtree = setQratesTree
  }
  
  if (is.null(setQratesData)) {
    Qdata = df
  } else if (is.character(setQratesData)) {
    Qdata = read.csv(setQratesData)
  } else {
    Qdata = setQratesData
  }
  
  if (!dir.exists("Simmap Overlap Outputs")) {
    dir.create("Simmap Overlap Outputs")
  }
  
  require(stringr)
  source("findQrates.R")
  cooprates <- findQrates(columns = columns[1], newdata = Qdata, newtree = Qtree)
  coopQ <- cooprates$qrates
  coopQ01 <- coopQ[3]
  coopQ10 <- coopQ[2]
  coopAnc = str_remove(cooprates$ARDlikanc, "ARDlik.anc ")  # added this for sim.history()
  coopAnc = as.numeric(coopAnc)
  names(coopAnc) <- c("0","1")
  FSrates <- findQrates(columns = columns[2], newdata = Qdata, newtree = Qtree)
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
  
  if (is.matrix(columnGlobalQrates)) {
    if (columnForGlobalQ == 1) {
      coopQ <- columnGlobalQrates
      coopQ01 <- coopQ[3]
      coopQ10 <- coopQ[2]
    } else {
      FSQ <- columnGlobalQrates
      FSQAbsPres <- FSQ[3]
      FSQPresAbs <- FSQ[2]
    }
  }
  
  
  if (dummy == FALSE) {
    FSsimtrees <- make.simmap(tree = subsettree, x = FSvec, model = "ARD", nsim = nsims, Q = FSQ)
    Coopsimtrees <- make.simmap(tree = subsettree, x = Coopvec, model = "ARD", nsim = nsims, Q = coopQ)
    datalabel = paste(datalabel, "REAL")
    if (plotSampleSimmaps) {
      
      # Generate filename base for both PDF and PNG
      filename_base <- paste(columns[1], columns[2], cladesubsetvalue, treelabel, "egSimmaps")
      
      pdf(file = file.path("Simmap Overlap Outputs", paste0(filename_base, ".pdf")),height=12,width=6)
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
      
      # Now create PNG output with same content
      png(file = file.path("Outputs/Figures/PNG", paste0(filename_base, ".png")),
          width = 6*150, height = 12*150, res = 150)
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
  
  set.seed(6000)
  dfout <- data.frame()
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
    
    
    treenum <- i
    
    temprow <- c(treenum, column1, column2, coopQ01, coopQ10, FSQAbsPres, FSQPresAbs, Nspecies, propFSabsent, propFSpresent, propNoncoop, propCoop, ObsProp0Absent, ObsProp0Present, ObsProp1Absent, ObsProp1Present, totaltime)
    dfout <- rbind(dfout, temprow)
    dfout <- as.data.frame(dfout)
    colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
    if (!dir.exists("Simmap Overlap Outputs")) {
      dir.create("Simmap Overlap Outputs")
    }
    if (i %in% c(100, 200, 250,500,1000,2000,3000,4000,5000)) {
      write.csv(dfout, file = file.path("Simmap Overlap Outputs", paste(datalabel, "simmap overlap_counts output nsim", nsims, treelabel,".csv")), row.names = FALSE)
    }
  }
  dfout <- as.data.frame(dfout)

  colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
  
  source("find transition counts by state for 2 Discrete traits.R")
  transStateCounts = getTransitionStateCounts(Coopsimtrees = Coopsimtrees, FSsimtrees = FSsimtrees)
  dfout = merge(dfout, transStateCounts, by.x = "treenum", by.y = "TreeNum")
  
  write.csv(dfout, file = file.path("Simmap Overlap Outputs", paste(datalabel, "simmap overlap_counts output nsim", nsims, treelabel,".csv")), row.names = FALSE)
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
  
  ObsProp0Absent_FractionDummyLessThanMedianReal = sum(dfDummy$ObsProp0Absent <= median(dfout$ObsProp0Absent))/length(dfDummy$treenum)
  ObsProp0Present_FractionDummyLessThanMedianReal = sum(dfDummy$ObsProp0Present <= median(dfout$ObsProp0Present))/length(dfDummy$treenum)
  ObsProp1Absent_FractionDummyLessThanMedianReal = sum(dfDummy$ObsProp1Absent <= median(dfout$ObsProp1Absent))/length(dfDummy$treenum)
  ObsProp1Present_FractionDummyLessThanMedianReal = sum(dfDummy$ObsProp1Present <= median(dfout$ObsProp1Present))/length(dfDummy$treenum)
  
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
  trait1StateLabels = getLabels(trait1)
  trait2StateLabels = getLabels(trait2)
  
  dfCombined$trait1StateLabel = NA
  dfCombined$trait2StateLabel = NA
  
  dfCombined$trait1StateLabel[which(dfCombined$ObservedState %in% c("ObsProp0Absent", "ObsProp0Present"))] = trait1StateLabels[1]
  dfCombined$trait1StateLabel[which(dfCombined$ObservedState %in% c("ObsProp1Absent", "ObsProp1Present"))] = trait1StateLabels[2]
  dfCombined$trait2StateLabel[which(dfCombined$ObservedState %in% c("ObsProp0Absent", "ObsProp1Absent"))] = trait2StateLabels[1]
  dfCombined$trait2StateLabel[which(dfCombined$ObservedState %in% c("ObsProp0Present", "ObsProp1Present"))] = trait2StateLabels[2]
  
 # dfCombined <- dfCombined %>%
 #  mutate(Label = paste(ObservedState, Which, sep = "\n"))
  dfCombined <- dfCombined %>%
    mutate(Label = paste(trait1StateLabel, trait2StateLabel, sep = "\n"))
  #dfCombined$Label[which(dfCombined$Which == "Real")] <- paste0(dfCombined$Label[which(dfCombined$Which == "Real")], " ") 
  levels_ordered = c(
    paste(trait1StateLabels[1], trait2StateLabels[1], sep = "\n"),
    paste(trait1StateLabels[1], trait2StateLabels[2], sep = "\n"),
    paste(trait1StateLabels[2], trait2StateLabels[1], sep = "\n"),
    paste(trait1StateLabels[2], trait2StateLabels[2], sep = "\n")
  )
  dfCombined$Label <- factor(dfCombined$Label, levels = levels_ordered, ordered = TRUE)
  
  
  # Create outputs for overall multi-comparison table
  dfout <- dfout %>%
    mutate(across(coopQ01:ObsProp1Present, as.numeric))
  medians = dfout %>% select(coopQ01:ObsProp1Present) %>% summarise_all(median, na.rm = TRUE)
  names(medians)[6:13] <- paste0(names(medians)[6:13],"_MedianReal")
  
  # quartiles <- dfout %>%
  #   select(ObsProp0Absent:ObsProp1Present) %>%
  #   summarise_all(function(x) list(quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)))
  # quartiles <- unnest(quartiles, cols = everything())
  
  
  mediansReal <- cbind(trait1, trait2, Nspecies, nsims_real, nsims_dummy, pval, medians)
  dfDummy <- dfDummy %>%
    mutate(across(coopQ01:ObsProp1Present, as.numeric))
  mediansDummy = dfDummy %>% select(propFSabsent:ObsProp1Present) %>% summarise_all(median, na.rm = TRUE) 
  names(mediansDummy) <- paste0(names(mediansDummy), "_MedianDummy")
  mediansRow = cbind(mediansReal, mediansDummy)
  
  require(emmeans)
  #require(ggpubr)
  lmStates = lm(ObservedState.prop ~ Which*ObservedState, data = dfCombined)
  emm <- emmeans(lmStates, pairwise ~ Which | ObservedState)
  contrast <- emm$contrasts
  PairwisePostHoc = summary(contrast, adjust = "tukey")
  PairwisePvals = PairwisePostHoc$p.value
  PairwisePvals = as.data.frame(t(as.data.frame(PairwisePvals)))
  colnames(PairwisePvals) <- paste0(PairwisePostHoc$ObservedState, "_DummyVsRealPval")
  mediansRow = cbind(mediansRow, PairwisePvals)
  mediansRow = cbind(mediansRow, ObsProp0Absent_FractionDummyLessThanMedianReal, ObsProp0Present_FractionDummyLessThanMedianReal, ObsProp1Absent_FractionDummyLessThanMedianReal, ObsProp1Present_FractionDummyLessThanMedianReal)
  
  
  # Create the boxplot
  p3 <-  ggplot(dfCombined, aes(x = Label, y = ObservedState.prop, fill = Which)) +
    geom_boxplot(outlier.shape = NA) + # Exclude outliers
    # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.05, position=position_dodge(width=0.75), color="black", binwidth = 0.003, stackratio = 0.5) + # consider adding dots to plot - will probably need arg tweaking based on number of dots, however.
    theme_minimal() +
    labs(y = "Observed State Proportion", x = "", fill = "Simulation Data") +
    scale_fill_manual(values = c("Real" = "#762a83", "Dummy" = "#1b7837")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste(trait1, trait2, "p =", pval))
  
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
      scale_fill_manual(values = c("Real" = "#762a83", "Dummy" = "#1b7837")) +
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
      
      
      transitioncols = c("Coop0to1inFS0", "Coop1to0inFS0", "Coop0to1inFS1", "Coop1to0inFS1", "FS0to1inCoop0", "FS1to0inCoop0", "FS0to1inCoop1", "FS1to0inCoop1")  # added 2/3/2024... weren't these supposed to be made already? What happened here? Also ExpectedCols...
      ExpectedCols = paste0(transitioncols, "Expected")
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
        scale_fill_manual(values = c("Observed" = "#762a83", "Expected" = "#1b7837")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), title = element_text(size = 8)) +
        ggtitle(paste(trait1, trait2, "\nnSimsObserved =", nsims_real, "    nSimsExpected =", nsims_real, "    Obs/Exp:TransCounts", pvalInteractionLabel)) +
        stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") 
      #ns: p > 0.05
      # *: p <= 0.05
      # **: p <= 0.01
      # ***: p <= 0.001
      # ****: p <= 0.0001
        
      p6 <-  ggplot(dfMeltCounts, aes(x = Transition, y = logCount, fill = ObservedVsExpected)) +
        geom_boxplot(outlier.shape = NA) + # Exclude outliers
        theme_minimal() +
        labs(y = "log(Transition Counts)", x = "", fill = "Observed/Expected") +
        scale_fill_manual(values = c("Observed" = "#762a83", "Expected" = "#1b7837")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), title = element_text(size = 8)) +
        ggtitle(paste(trait1, trait2, "\nnSimsObserved =", nsims_real, "    nSimsExpected =", nsims_real, "    Obs/Exp:TransCounts", pvalInteractionLogLabel)) +
        stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
        scale_x_discrete(labels = groupLabels)
      
    } else {
      p5 = NULL
      p6 = NULL
    }
    
  } else {
    p4 = NULL
    p5 = NULL
    p6 = NULL
  }
  
  # Generate PNG and PDF outputs if requested
  if (plot_ggplots_pdf == TRUE) {
    require(cowplot)
    # Generate filename base for both PDF and PNG
    filename_base <- paste0("Simmap_Overlap_", trait1, "_", trait2, "_", nsims_real, "sims_real_", nsims_dummy, "sims_dummy", otherlabel)
    
    pdf(file = file.path("Simmap Overlap Outputs", paste0(filename_base, ".pdf")), height = 14, width = 6)
    print(plot_grid(p1, p2, p3, ncol = 1))
    dev.off()
    
    # Now create PNG output with same content
    png(file = file.path("Outputs/Figures/PNG", paste0(filename_base, ".png")),
        width = 6*150, height = 14*150, res = 150)
    print(plot_grid(p1, p2, p3, ncol = 1))
    dev.off()
  }
  
  # Return appropriate list based on what plots were created
  if (!is.null(p6)) {
    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, filename = PDFname, TransitionStats = TransitionStats, dfMeltCounts = dfMeltCounts, mediansRow = mediansRow))
  } else if (!is.null(p4)) {
    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, filename = PDFname, mediansRow = mediansRow))
  } else {
    return(list(p1 = p1, p2 = p2, p3 = p3, filename = PDFname, mediansRow = mediansRow))
  }
  
}


#### calcHuelflex fxn for multistate categorical traits ----
# only works with output csvs that have the state combinations with DUMMY and REAL in the column names
calcHuelflex = function(overlapdf) { # 
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(ggplot2)
  overlapdf$X=NULL
  colnames(overlapdf)
  
  nsims = length(overlapdf[,1])
  numericColumnStart = max(which(str_detect(colnames(overlapdf), "trait")))+1
  overlapdf[,numericColumnStart:length(colnames(overlapdf))] = apply(overlapdf[,numericColumnStart:length(colnames(overlapdf))], MARGIN = 2, FUN = as.numeric) 
  overlaplonger = overlapdf %>%   pivot_longer(
    cols = !c(tree, trait1, trait2),  # may need to be c(tree, trait1, trait2) for some old files?
    names_to = "state",  
    values_to = "proportion"          # The name of the new column for the values
  )
  overlaplonger$Which = NA
  overlaplonger$Which[which(str_detect(overlaplonger$state, "DUMMY"))] = "Dummy"
  overlaplonger$Which[which(str_detect(overlaplonger$state, "REAL"))] = "Real"
  overlaplonger$state <- gsub( "_DUMMY", "", overlaplonger$state)
  overlaplonger$state <- gsub( "_REAL", "", overlaplonger$state)
  
  df_long <- overlaplonger %>%
    separate(state, into = c("trait_state", "FS"), sep = "_FS") 
  
  total_times_trait <- df_long %>%
    group_by(tree, trait1, trait2, trait_state, Which) %>%
    summarize(total_trait = sum(proportion), .groups = "drop")
  
  total_times_FS <- df_long %>%
    group_by(tree, trait1, trait2, FS, Which) %>%
    summarize(total_FS = sum(proportion), .groups = "drop")
  
  df_long <- df_long %>%
    left_join(total_times_trait, by = c("tree", "trait1", "trait2", "trait_state", "Which")) %>%
    left_join(total_times_FS, by = c("tree", "trait1", "trait2", "FS", "Which"))
  
  df_long <- df_long %>%
    mutate(expected_proportion = total_trait * total_FS)
  
  df_real <- df_long %>%
    filter(Which == "Real")
  df_real <- df_real %>%
    mutate(abs_diff = abs(proportion - expected_proportion))
  
  D_real <- df_real %>%
    summarize(total_abs_diff = sum(abs_diff)) %>%
    pull(total_abs_diff) / nsims
  
  Real_dsims <- df_real %>%
    select(tree, trait1, trait2, trait_state, FS, abs_diff) %>%
    group_by(tree, trait1, trait2) %>%
    summarize(Real_dsim = sum(abs_diff), .groups = "drop")
  
  # Dummy data
  # Step 1: Filter for "Dummy" data
  df_dummy <- df_long %>%
    filter(Which == "Dummy")
  
  # Step 2: Calculate the absolute differences
  df_dummy <- df_dummy %>%
    mutate(abs_diff = abs(proportion - expected_proportion))
  
  # Step 3: Calculate Dummy_dsums (row-wise sums of the absolute differences)
  Dummy_dsums <- df_dummy %>%
    group_by(tree, trait1, trait2) %>%
    summarize(Dummy_dsum = sum(abs_diff), .groups = "drop")
  
  #hist(Real_dsims$Real_dsim)
  #hist(Dummy_dsums$Dummy_dsum)
  #abline(v = D_real)
  numGreater = sum(Dummy_dsums$Dummy_dsum > D_real)
  pval = sum(Dummy_dsums$Dummy_dsum > D_real)/nsims
  
  ## Get num Dummy greater than median real
  # Calculate medians for "Real" data
  medians_real <- df_long %>%
    filter(Which == "Real") %>%
    group_by(trait_state, FS) %>%
    summarize(median_real = median(proportion), .groups = "drop")
  
  # Filter for "Dummy" data
  df_dummy <- df_long %>%
    filter(Which == "Dummy")
  
  # Join the medians back to the "Dummy" data
  df_dummy <- df_dummy %>%
    left_join(medians_real, by = c("trait_state", "FS"))
  
  # Calculate the fraction for each state in "Dummy" data
  fraction_dummy_less_than_median_real <- df_dummy %>%
    group_by(trait_state, FS) %>%
    summarize(fraction = sum(proportion <= median_real) / n(), .groups = "drop")
  
  
  trait1 = overlaplonger$trait1[1]
  trait2 = overlaplonger$trait2[1]
  
  plotlabel = paste(trait1, trait2)
  dummytitle = paste("Nsims =", nsims, "\nnum Dummy dsums > D_real:", numGreater, ", pval =", pval)
  
  xmax = max(c(Real_dsims$Real_dsim, Dummy_dsums$Dummy_dsum))*1.1
  
  # ggplot histograms to return
  p1 <- ggplot(Real_dsims, aes(x = Real_dsim)) +
    geom_histogram(binwidth = xmax/20, fill = rgb(0.2, 0.5, 0.7, 0.5), color = "white") +
    geom_vline(xintercept = D_real, color = "red") +
    xlim(c(0, xmax)) +
    labs(title = plotlabel, x = "D statistic from real data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)
  
  # Create the histogram for Dummy_dsums
  p2 <- ggplot(Dummy_dsums, aes(x = Dummy_dsum)) +
    geom_histogram(binwidth = xmax/20, fill = rgb(0.7, 0.5, 0.2, 0.5), color = "white") +
    xlim(c(0, xmax)) +
    labs(title = dummytitle, x = "D statistic from simulated independent data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)
  
  #pdf(paste("simmap overlap Real Dummy multiState Proportions Boxplot -", trait1, trait2, nsims, "sims.pdf"))
  boxplotStates <- ggplot(overlaplonger, aes(x = state, y = proportion, fill = Which)) +
    geom_boxplot(outlier.shape = NA) + # Exclude outliers
    theme_minimal() +
    labs(y = "Observed State Proportion", x = "", fill = "Simulation Data") +
    scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste(trait1, trait2, "p =", pval))
  #dev.off()
  
  return(list(p1=p1, p2=p2, boxplotStates = boxplotStates, fraction_dummy_less_than_median_real = fraction_dummy_less_than_median_real))
  
}

#### getLabels ----
getLabels <- function(trait) {
  require(stringr)
  if (trait == "Final.polygyny") {
    return(c("Monogamy", "Polygyny"))
  } else if (grepl("coop", trait, ignore.case = TRUE)) {
    return(c("Non-Cooperative", "Cooperative"))
  } else if (str_detect(trait, "Kin")) {
    return(c("Non-kin", "Kin"))
  } else if (str_detect(trait, "Familial")) {
    return(c("Non-Familial Living", "Familial Living"))
  } else if (str_detect(trait, "Colonial")) {
    return(c("Non-Colonial", "Colonial"))
  } else if (str_detect(trait, "GroupsLargerThanPair")) {
    return(c("Asocial or pair", "Small or large groups"))
  } else if (str_detect(trait, "LongSocialBonds")) {
    return(c("Bonds last one season or less", "Multi-year bonds"))
  } else if (str_detect(trait, "MoreThanTwoCaretakers")) {
    return(c("Two or fewer caretakers", "More than two caretakers"))
  } else if (str_detect(trait, "TwoOrMoreCaretakers")) {
    return(c("Fewer than two caretakers", "Two or more caretakers"))
  } else if (str_detect(trait, "Asocial")) {
    return(c("Asocial", "Pair or group sociality"))
  } else if (str_detect(trait, "SeasonOrLonger")) {
    return(c("Bonds last less than one season", "Season or longer social bonds"))
  } else if (str_detect(trait, "LargestGroupSizes")) {
    return(c("Asocial, pair, or small groups", "Large groups"))
  } else if (str_detect(trait, "FemaleSong")) {
    return(c("Female Song Absent", "Female Song Present"))  
  } else if (str_detect(trait, "HighConfidence_Coop")) {
    return(c("Non-Cooperative", "Cooperative"))
  } else {
    return(c(paste(trait, "0"), paste(trait, "1")))  # Return NA if no condition matches
  }
}
