# Brownie multistate discrete
# Kate Snyder

# returns several csv files, a pdf plot, and a dataframe with the result statistics
# added requirements 2/8/2025
# 3/18/2025 -  changed plotting to be more flexible for different numbers of states of the discrete trait
# 6/15/2025 - added merge_simmaps() and getSimmapSegments() from merge_simmaps.R
# 7/1/2025 - added calculation of overall p-value from median logliks to summary output
require(stringr)
require(phytools)
source("subsettreedata.R")

# Example: BrownieMultistate(DiscreteTrait = "TerrCoop", ContinuousTrait = "Syllable.rep.final", newdata = newdata, treefile = treefile, nsim = 500, importSimmaps = TerrCoopsims)

BrownieMultistate <- function(DiscreteTrait, ContinuousTrait, newdata, treefile, nsim, plotsimmaps = F, plotResults = T, filename = NULL, otherlabel = "", importSimmaps = NULL)
{
  
  if (is.null(importSimmaps[1])) {
    
    subsetout = subsettreedata(columns = DiscreteTrait, newdata = newdata, newtree = treefile)
    subsetDisctree = subsetout$subsettree
    subsetDiscdf = subsetout$subsetdf
    
    discretetraitvecDisc = subsetDiscdf[,DiscreteTrait]
    names(discretetraitvecDisc) = subsetDiscdf$species
    
    
    ERmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ER")
    ARDmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ARD")
    anovaERARD <- anova(ERmodel,ARDmodel)
    SYMmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "SYM")
    anovaERSYM = anova(ERmodel,SYMmodel)
    anovaSYMARD <- anova(SYMmodel,ARDmodel)
    
    print(anovaERSYM)
    print(anovaSYMARD)
    
    
    # get rates from ace() output
    aceARDratesVec = ARDmodel$rates
    aceARDrates = cbind(1:length(aceARDratesVec), aceARDratesVec)
    aceARDrates = as.data.frame(aceARDrates)
    colnames(aceARDrates) <- c("rate_index", "rates")
    rate_index_matrix = ARDmodel$index.matrix
    groupnames = colnames(ARDmodel$lik.anc)
    nGroups = length(groupnames)
    rate_matrix = matrix(rep(NA,nGroups^2), nrow = nGroups)
    rownames(rate_matrix) = colnames(rate_matrix) = groupnames
    
    for (i in 1:nGroups) {
      for (j in 1:nGroups) {
        index = rate_index_matrix[i,j]
        if (!is.na(index)) {
          rate_matrix[i,j] = aceARDrates$rates[which(aceARDrates$rate_index == index)]
        }
      }
    }
    diagvals = rowSums(rate_matrix, na.rm = T)*-1
    diag(rate_matrix) <- diagvals
    print(rate_matrix)
    
    
    # Make simmaps from data subsetted to those with song data
    subsetSong = subsettreedata(columns = c(DiscreteTrait, ContinuousTrait), newdata = newdata, newtree = treefile, islog = ContinuousTrait)
    subsetdf = subsetSong$subsetdf
    subsettree = subsetSong$subsettree
    
    discretetraitvec = subsetdf[,DiscreteTrait]
    names(discretetraitvec) = subsetdf$species
    
    continuoustraitvec = subsetdf[,ContinuousTrait]
    names(continuoustraitvec) = subsetdf$species
    
    simmappy = make.simmap(subsettree, discretetraitvec, nsim = nsim, Q= rate_matrix, type = "discrete") 
    
    plotSimmap(simmappy[[1]])
    
    if (plotsimmaps) { 
      simmapFileName = paste0(DiscreteTrait, " multistate simmap plots ", ContinuousTrait, " subset ", otherlabel,  Sys.Date(),".pdf")
      pdf(simmapFileName, height = 9, width = 12)
      par(mfrow = c(2,3))
      par(mar = c(3.8,3.8,3,1))
      for (simmapNum in 1:5) {
        simmap = simmappy[[simmapNum]]
        simmapQ = simmap$Q
        SimmapStates = colnames(simmap$Q)
        py = c("orange", "#009E73", "blue", "#CC79A7")
        pynamed <- py[1:length(SimmapStates)]
        names(pynamed) <- SimmapStates
        tipcols = rep(NA, length(simmap$tip.label))
        for (stateNum in 1:length(SimmapStates)) { # get color vector of tips
          tempstate = SimmapStates[stateNum]
          tipcols[which(simmap$tip.label %in% names(which(discretetraitvec==tempstate)))] <- pynamed[tempstate]
        }
        # Plot simmap 
        plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)
        tiplabels(pch=21,bg=tipcols, col = tipcols, cex=0.3)
        legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
      }
      # plot to add rate matrix
      # Create an empty plot
      par(mar = c(5,5,4,1))
      plot(1, type = "n", xlim = c(0, ncol(simmapQ)+1), ylim = c(0, nrow(simmapQ)+1), 
           xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", xaxs = "i", yaxs = "i")
      
      # Add column and row names
      axis(1, at = 1:ncol(simmapQ), labels = colnames(simmapQ), las = 2, tick = FALSE)
      axis(2, at = 1:nrow(simmapQ), labels = rev(rownames(simmapQ)), las = 2, tick = FALSE)
      # Add the matrix values
      for (i in 1:nrow(simmapQ)) {
        for (j in 1:ncol(simmapQ)) {
          text(j, nrow(simmapQ) - i + 1, round(simmapQ[i, j], 6))
        }
      }
      # Add label re transitions
      axis(1, at = 0.1, labels = "To:", las = 1, tick = FALSE, font = 2)
      axis(2, at = nrow(simmapQ)+0.9, labels = "From:", las = 2, tick = FALSE, font = 2)
      title("Transition Rates")
      dev.off()
    }
    
  } else if ("multiSimmap" %in% class(importSimmaps)) {
    
    subsetSong = subsettreedata(columns = c(ContinuousTrait), newdata = newdata, newtree = treefile, islog = ContinuousTrait)
    subsetdf = subsetSong$subsetdf
    subsettree = subsetSong$subsettree
    
    
    pruned_SimsAndData <- prune_multiSimmap(multiSimmap = importSimmaps, df = subsetdf)
    
    subsetdf = pruned_SimsAndData$pruned_df
    simmappy = pruned_SimsAndData$pruned_simmaps
    
    continuoustraitvec = subsetdf[,ContinuousTrait]
    names(continuoustraitvec) = subsetdf$species
    
    
  } else {
    stop("importSimmaps must be of class c('multiSimmap', 'multiPhylo')")
  }
  
  egSimmap = simmappy[[1]]
  brownieliteresults <- brownie.lite(egSimmap,continuoustraitvec,maxit=75000)
  statenames = names(brownieliteresults$sig2.multiple)
  
  ARDRateColNames = paste0("ARDRate_", statenames)
  
  browniedata <- data.frame(DiscreteTrait=character(nsim),ContinuousTrait=character(nsim),Pval=numeric(nsim),ERRate=numeric(nsim),ERloglik=numeric(nsim),ERace=numeric(nsim),ARDloglik=numeric(nsim),ARDace=numeric(nsim),k2=numeric(nsim),convergence=character(nsim),simmapnumber=integer(nsim),phylanovaP=numeric(nsim),stringsAsFactors = FALSE)
  browniedata[,ARDRateColNames] = NA
  ncolsBrownieData = length(colnames(browniedata))
  brownied <- set.seed(10)
  browniedf <- set.seed(10)
  
  browniedata$DiscreteTrait = DiscreteTrait
  browniedata$ContinuousTrait = paste0("log",ContinuousTrait)
  
  
  for (i in 1:nsim) {
    if (i %in% seq(0,2500,by=50)) {
      print(i)
    }
    simmapfor <- simmappy[[i]]
    brownieliteresults <- set.seed(10)
    
    # tryCatch(
    #   expr = {
    #     withTimeout(expr={
    
    brownieliteresults <- brownie.lite(simmapfor,continuoustraitvec,maxit=75000)
    browniedata[i,3] <- brownieliteresults$P.chisq
    browniedata[i,4] <- brownieliteresults$sig2.single
    browniedata[i,5] <- brownieliteresults$logL1
    browniedata[i,6] <- brownieliteresults$a.single
    browniedata[i,7] <- brownieliteresults$logL.multiple
    browniedata[i,8] <- brownieliteresults$a.multiple
    browniedata[i,9] <- brownieliteresults$k2
    browniedata[i,10] <- as.character(brownieliteresults$convergence)
    browniedata[i,11] <- i 
    browniedata[i,12] <- NA 
    for (ARDcolumn in 1:length(ARDRateColNames)) {
      browniedata[i,ARDRateColNames[ARDcolumn]] <- brownieliteresults$sig2.multiple[ARDcolumn]
    }
    
    #}, timeout = 16, cpu=Inf, onTimeout = "error")
    # },
    # TimeoutException = function(ex) {browniedata[i,3:ncolsBrownieData]<-c(NA,NA,NA,NA,NA,NA,NA,NA,i,"timeout", rep(NA,times = length(ARDRateColNames)));
    # print(paste("timeout",i));
    # },
    # error = function(e) {browniedata[i,3:ncolsBrownieData]<-c(NA,NA,NA,NA,NA,NA,NA,NA,i,"error", rep(NA,times = length(ARDRateColNames)));
    # print(paste("compute error",i));
    # })
  }
  
  # add columns that say which rates are higher in each sim
  for (ARDcolumn1Num in 1:length(ARDRateColNames)) {
    ARDcolumn = ARDRateColNames[ARDcolumn1Num]
    if (ARDcolumn1Num != length(ARDRateColNames)) {
      for (ARDcolumn2Num in (ARDcolumn1Num+1):length(ARDRateColNames)) {
        ARDcolumn2 = ARDRateColNames[ARDcolumn2Num] 
        newcolname = paste0(ARDcolumn,"_greater_than_", ARDcolumn2)
        browniedata[,newcolname] = browniedata[,ARDcolumn] > browniedata[,ARDcolumn2]
      }
    }
  }
  brownieCSVname = paste(Sys.Date(), DiscreteTrait, ContinuousTrait, "multistate aceARD Brownie", otherlabel, nsim, "sims.csv")
  write.csv(browniedata, brownieCSVname, row.names = F)
  
  # make compare rates csv
  CompareColNames = colnames(browniedata)[grep("_greater_than_", colnames(browniedata))]
  CompareColSums = data.frame(DiscreteTrait = rep(DiscreteTrait, length(CompareColNames)), ContinuousTrait = rep(ContinuousTrait, length(CompareColNames)), nsims = rep(nsim, length(CompareColNames)))
  CompareColSums$RateComparison = CompareColNames
  CompareColSums$Sums = colSums(browniedata[,CompareColNames])
  CompareColSums$Fraction1 = CompareColSums$Sums/CompareColSums$nsims
  SplitRates = str_remove_all(CompareColSums$RateComparison, "ARDRate_")
  SplitRates = str_split(SplitRates, "_greater_than_")
  CompareColSums$HigherRate = NA
  CompareColSums$LowerRate = NA
  CompareColSums$Fraction = NA
  for (i in 1:length(SplitRates)) { # added 1/24/2024
    if (CompareColSums$Fraction1[i] >= 0.5) {
      CompareColSums$HigherRate[i] = SplitRates[[i]][1]
      CompareColSums$LowerRate[i] = SplitRates[[i]][2]
      CompareColSums$Fraction[i] = CompareColSums$Fraction1[i]
    } else {
      CompareColSums$HigherRate[i] = SplitRates[[i]][2]
      CompareColSums$LowerRate[i] = SplitRates[[i]][1]
      CompareColSums$Fraction[i] = 1-CompareColSums$Fraction1[i]
    }
  }
  compareCSVname = paste(Sys.Date(), DiscreteTrait, ContinuousTrait, "multistate aceARD Brownie", otherlabel, nsim, "sims COMPARE RATES.csv")
  write.csv(CompareColSums, compareCSVname, row.names = F)
  
  
  ARDRateColNames = colnames(browniedata)[which(str_detect(colnames(browniedata), "ARDRate") & str_detect(colnames(browniedata), "greater", negate = T))]
  DiscreteTrait = browniedata$DiscreteTrait[1]
  ContinuousTrait = browniedata$ContinuousTrait[1]
  nsim = length(browniedata$DiscreteTrait)
  if (browniedata$DiscreteTrait[1] %in% c("grouping", "grouping_Griesser2023")) {
    ARDRateColNames = c("ARDRate_asocial", "ARDRate_pair", "ARDRate_small_groups", "ARDRate_large_groups")
  }
  
  ## make rate distribution plots ----
  if (plotResults == TRUE) {
    pdfname = paste0(DiscreteTrait, " ", ContinuousTrait, " multistate aceARD Brownie ", otherlabel," ", nsim, " sims ", Sys.Date(), ".pdf")
    pdf(pdfname, width = 8, height = 7)
    par(mar = c(3.8,3.8,3,1))
    par(mfrow = c(2,1))
    
    #### old plot method - inflexible ----
    # for (ARDcolumn in 1:length(ARDRateColNames)) {
    #   assign(x = paste0("D",ARDcolumn-1), value = density(browniedata[,ARDRateColNames[ARDcolumn]]))
    # }
    # 
    # if (length(ARDRateColNames) == 4) {
    #   xmax = max(c(D0$x,D1$x, D2$x, D3$x))
    #   if (xmax > 1.5) {xmax = .75}
    #   plot(D0,col="orange",
    #        xlim=c(min(c(D0$x,D1$x, D2$x, D3$x)),
    #               xmax),
    #        ylim=c(min(c(D0$y,D1$y, D2$y, D3$y)),
    #               max(c(D0$y,D1$y, D2$y, D3$y))),
    #        main="", xlab="",ylab="", cex.axis=1)     
    #   title(main="", cex.main = 2, line = 1)
    #   title(ylab = "Frequency",line=2.5, cex.lab=1.15)
    #   #   axis(1, cex.axis=1.2)
    #   #    axis(2, cex.axis=1.2)
    #   lines(D1, col="#009E73")
    #   lines(D2, col="blue")
    #   lines(D3, col="#CC79A7")
    #   abline(v=browniedata$ERRate[1], lty = 2)
    #   
    #   legend("topright",legend = c(ARDRateColNames,"Equal Rates"), lwd=1,col=c("orange", "#009E73", "blue", "#CC79A7", "black"), lty = c(1,1,1,1,2), cex=1) # check order
    #   title(xlab=paste0("Rate of evolution of ", ContinuousTrait),line = 2.5, cex.lab = 1)
    #   
    # } else if (length(ARDRateColNames) == 3) {
    #   plot(D0,col="orange",
    #        xlim=c(min(c(D0$x,D1$x, D2$x)),
    #               max(c(D0$x,D1$x, D2$x))),
    #        ylim=c(min(c(D0$y,D1$y, D2$y)),
    #               max(c(D0$y,D1$y, D2$y))),
    #        main="", xlab="",ylab="", cex.axis=1)     
    #   title(main="", cex.main = 2, line = 1)
    #   title(ylab = "Frequency",line=2.5, cex.lab=1.15)
    #   lines(D1, col="#009E73")
    #   lines(D2, col="blue")
    #   abline(v=browniedata$ERRate[1], lty = 2)
    #   
    #   legend("topright",legend = c(ARDRateColNames,"Equal Rates"), lwd=1,col=c("orange", "#009E73", "blue", "black"), lty = c(1,1,1,2), cex=1) # check order
    #   title(xlab=paste0("Rate of evolution of ", ContinuousTrait),line = 2.5, cex.lab = 1)
    
    #### new plot ----
    # Generate density objects dynamically
    ARDRateColNames = sort(ARDRateColNames)
    densities <- lapply(ARDRateColNames, function(colname) {
      density(browniedata[[colname]])
    })
    
    # Determine x and y axis limits dynamically
    xmax <- max(sapply(densities, function(d) max(d$x)))
    if (xmax > 1.5) xmax <- 0.75
    
    xlim_vals <- range(unlist(lapply(densities, function(d) d$x)))
    ylim_vals <- range(unlist(lapply(densities, function(d) d$y)))
    
    # Define colors dynamically (expand as needed)
    color_palette <- c("orange", "#009E73", "blue", "#CC79A7", "purple", "red", "cyan", "magenta")
    
    # Plot first density
    plot(densities[[1]], col=color_palette[1], xlim=xlim_vals, ylim=ylim_vals, 
         main="", xlab="", ylab="", cex.axis=1)
    
    # Overlay remaining densities
    for (i in 2:length(densities)) {
      lines(densities[[i]], col=color_palette[i])
    }
    
    # Add title, labels, and legend dynamically
    title(main="", cex.main=2, line=1)
    title(ylab="Frequency", line=2.5, cex.lab=1.15)
    title(xlab=paste0("Rate of evolution of ", ContinuousTrait), line=2.5, cex.lab=1)
    
    # Add vertical line for equal rates
    abline(v=browniedata$ERRate[1], lty=2)
    
    # Create legend dynamically
    legend("topright", legend=c(ARDRateColNames, "Equal Rates"), lwd=1, 
           col=c(color_palette[1:length(ARDRateColNames)], "black"), 
           lty=c(rep(1, length(ARDRateColNames)), 2), cex=1)
    
    
    #### pval plot ----
    pdens = density(browniedata$Pval)
    plot(pdens, main = paste("N =", nsim), ylab = "", xlab = "")
    title(ylab = "Frequency",line=2.5, cex.lab=1.15)
    title(xlab = "p-value", line = 2.5, cex.lab = 1)
    abline(v=0.05, col = "gray")
    dev.off()
  } # end if plotResults == TRUE
  
  ## summary table ----
  #subset results 
  summarydf = set.seed(10)
  browniedf <- browniedata[browniedata$convergence == "Optimization has converged.",]
  browniedf <- browniedf[!is.na(browniedf$convergence),]
  
  #calculate overall mean pval
  ERloglikmean <- mean(browniedf$ERloglik, na.rm = T)
  ARDloglikmean <- mean(browniedf$ARDloglik, na.rm = T)
  ERARDPval = round(pchisq(2*(ARDloglikmean-ERloglikmean),1,lower.tail=FALSE), digits = 3) #testing whether the two rates of continuous trait evolution are significantly different
  ERARDPval_from_Mean_logliks = round(pchisq(2*(ARDloglikmean-ERloglikmean),1,lower.tail=FALSE), digits = 3) #testing whether the two rates of continuous trait evolution are significantly different
  
  ERloglikmedian <- median(browniedf$ERloglik, na.rm = T)
  ARDloglikmedian <- median(browniedf$ARDloglik, na.rm = T)
  ERARDPval_from_Median_logliks = round(pchisq(2*(ARDloglikmedian-ERloglikmedian),1,lower.tail=FALSE), digits = 3) #testing whether the two rates of continuous trait evolution are significantly different
  
  nSig = sum(browniedf$Pval < 0.05, na.rm = T)
  medianPval = median(browniedf$Pval)
  fractionSig = nSig/nsim
  temprow = c(brownieCSVname, DiscreteTrait, ContinuousTrait, nsim, nSig, fractionSig, medianPval, ERloglikmean, ARDloglikmean, ERARDPval, ERARDPval_from_Mean_logliks, ERloglikmedian, ARDloglikmedian, ERARDPval_from_Median_logliks)
  as.data.frame(temprow)
  summarydf = as.data.frame(rbind(summarydf, temprow))
  colnames(summarydf) = c("filename", "DiscreteTrait", "ContinuousTrait", "nsims", "nSig", "fractionSig", "medianPval", "ERloglikmean", "ARDloglikmean", "ERARDPval", "ERARDPval_from_Mean_logliks", "ERloglikmedian", "ARDloglikmedian", "ERARDPval_from_Median_logliks")
  return(summarydf)
}


prune_multiSimmap <- function(multiSimmap, df) {
  require(phytools)
  
  # Get the list of species present in the dataframe
  species_in_df <- as.character(df$species)
  
  # Initialize a new list for pruned simmaps
  pruned_simmaps <- list()
  
  for (i in seq_along(multiSimmap)) {
    simmap <- multiSimmap[[i]]
    
    # Find tips to drop (those not in the dataframe)
    tips_to_drop <- setdiff(simmap$tip.label, species_in_df)
    
    # Drop tips using phytools function
    pruned_simmap <- drop.tip.simmap(simmap, tips_to_drop)
    
    # Store in new list
    pruned_simmaps[[i]] <- pruned_simmap
  }
  
  # Convert to multiSimmap format
  class(pruned_simmaps) <- c("multiSimmap", "multiPhylo")
  
  # Prune the dataframe to match species in the new tree
  pruned_species <- pruned_simmaps[[1]]$tip.label  # Use first tree as reference
  pruned_df <- df[df$species %in% pruned_species, ]
  
  return(list(pruned_simmaps = pruned_simmaps, pruned_df = pruned_df))
}




merge_simmaps <- function(simmap1, simmap2) {
  require(phytools)
  
  # Ensure both simmaps use the same phylogeny
  if (!all(simmap1$edge == simmap2$edge)) {
    stop("Simmap trees do not match in edge structure!")
  }
  
  # Initialize the merged simmap object
  merged_simmap <- simmap1
  merged_simmap$maps <- vector("list", length(simmap1$edge.length))
  
  # Loop through each edge and merge state maps
  for (edge in seq_along(simmap1$maps)) {
    map1 <- simmap1$maps[[edge]]
    map2 <- simmap2$maps[[edge]]
    
    states1 <- names(map1)
    states2 <- names(map2)
    
    seg1 <- cumsum(map1)
    seg2 <- cumsum(map2)
    
    new_map <- list()
    
    i <- 1
    j <- 1
    last_pos <- 0
    
    while (i <= length(map1) & j <= length(map2)) {
      # Get segment start & end positions
      start_pos <- max(ifelse(i == 1, 0, seg1[i - 1]), ifelse(j == 1, 0, seg2[j - 1]))
      end_pos <- min(seg1[i], seg2[j])
      
      if (start_pos < end_pos) {
        combined_state <- paste0(states1[i], states2[j])
        new_map[[length(new_map) + 1]] <- end_pos - start_pos
        names(new_map)[length(new_map)] <- combined_state
      }
      
      # Move to the next segment
      if (seg1[i] == end_pos) i <- i + 1
      if (seg2[j] == end_pos) j <- j + 1
    }
    
    # Assign the new map to the merged simmap
    merged_simmap$maps[[edge]] <- unlist(new_map)
  }
  
  # Update mapped.edge matrix to reflect new state labels
  unique_states <- unique(unlist(lapply(merged_simmap$maps, names)))
  merged_simmap$mapped.edge <- matrix(0, nrow = nrow(simmap1$mapped.edge), ncol = length(unique_states))
  colnames(merged_simmap$mapped.edge) <- unique_states
  
  for (edge in seq_along(merged_simmap$maps)) {
    for (state in names(merged_simmap$maps[[edge]])) {
      merged_simmap$mapped.edge[edge, state] <- merged_simmap$maps[[edge]][state]
    }
  }
  
  return(merged_simmap)
}


getSimmapSegments = function(tree) {
  # Initialize an empty dataframe
  
  outlist = list()
  
  transitions_df <- data.frame(EdgeNumber = integer(),
                               Transition01or10 = character(),
                               DistanceFromRootwardNode = numeric(),
                               stringsAsFactors = FALSE)
  
  segments_df <- data.frame(EdgeNumber = integer(),
                            State = character(),
                            Start = numeric(), End = numeric(),
                            stringsAsFactors = FALSE)
  
  # Loop through each edge in maps
  for (edge in 1:length(tree$maps)) {
    edge_data <- tree$maps[[edge]]
    states <- names(edge_data)
    
    # Initialize distance from rootward node
    distance = 0
    # Check for transitions
    if (length(states) == 1) {
      # No transition
      transitions_df <- rbind(transitions_df, data.frame(EdgeNumber=edge, Transition01or10=NA, DistanceFromRootwardNode=NA))
    } else {
      # Transition exists, calculate and record
      for (i in seq_along(edge_data)) {
        distance <- distance + edge_data[i]
        if (i < length(edge_data)) {
          transition <- paste0(states[i], "to", states[i+1])
          transitions_df <- rbind(transitions_df, data.frame(EdgeNumber=edge, Transition01or10=transition, DistanceFromRootwardNode=distance))
        }
      }
    }
    
    # Initialize start position
    start_position = 0
    
    # Loop through each segment in the edge
    for (i in seq_along(edge_data)) {
      end_position <- start_position + edge_data[i]
      segments_df <- rbind(segments_df, data.frame(EdgeNumber=edge, State=as.numeric(states[i]), Start=start_position, End=end_position))
      start_position <- end_position
    }
    
  }
  
  outlist$transitions_df = transitions_df
  outlist$segments_df = segments_df
  # View the resulting dataframe
  return(outlist)
}
