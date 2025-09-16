## Make multiple-trait simmaps from two simmaps
## 

allRealList <- readRDS("Simmap Overlap Outputs/500 multisimmaps Territory FemaleSong_Agg01 HighConfidence_Coop subsettree_AllRealData 2025-03-04 .rds")
FSsims = allRealList$multisimmapFSAgg01Real
Terrsims = allRealList$multisimmapMultistateTerritoryReal
Coopsims = allRealList$multisimmapHCCoopReal

TerrFS1 = merge_simmaps(simmap1 = Terrsims[[1]], simmap2 = FSsims[[1]])
TerrFS2 = merge_simmaps(simmap1 = Terrsims[[2]], simmap2 = FSsims[[2]])
TerrCB1 = merge_simmaps(simmap1 = Terrsims[[1]], simmap2 = Coopsims[[1]])
TerrCB2 = merge_simmaps(simmap1 = Terrsims[[2]], simmap2 = Coopsims[[2]])

nsims = 500
TerrCoopsims = list()
for (i in 1:nsims) {
  TerrCoopsims[[i]] <- merge_simmaps(simmap1 = Terrsims[[i]], simmap2 = Coopsims[[i]])
}
class(TerrCoopsims) <- c("multiSimmap", "multiPhylo")


## plot simmaps
simmapFileName = paste0("Territory CoopBreed intersection", " multistate simmap plots ", Sys.Date(),".pdf")
pdf(simmapFileName, height = 30, width = 12)
par(mfrow = c(2,3))
par(mar = c(3.8,3.8,3,1))
  # plot combo 1
  simmap = TerrCB1
  simmapQ = simmap$Q
  SimmapStates = sort(colnames(simmap$mapped.edge))
  py = c("black", "darkgrey", "#CD0BBC", "orange", "#2297E6", "#28E2E5")
  pynamed <- py[1:length(SimmapStates)]
  names(pynamed) <- SimmapStates
  plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)
  legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
  
  # plot Terr 1
  simmap = Terrsims[[1]]
  simmapQ = simmap$Q
  SimmapStates = sort(colnames(simmap$mapped.edge))
  py = c("black", "orange", "#2297E6")
  pynamed <- py[1:length(SimmapStates)]
  names(pynamed) <- SimmapStates
  plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)
  legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
  
  # plot FS 1
  simmap = Coopsims[[1]]
  simmapQ = simmap$Q
  SimmapStates = sort(colnames(simmap$mapped.edge))
  py = c("#CD0BBC", "#28E2E5")
  pynamed <- py[1:length(SimmapStates)]
  names(pynamed) <- SimmapStates
  plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)
  legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
  

  # plot combo 2
  simmap = TerrCB2
  simmapQ = simmap$Q
  SimmapStates = sort(colnames(simmap$mapped.edge))
  py = c("black", "darkgrey", "#CD0BBC", "orange", "#2297E6", "#28E2E5")
  pynamed <- py[1:length(SimmapStates)]
  names(pynamed) <- SimmapStates
  plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)
  legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
  
  # plot Terr 2
  simmap = Terrsims[[2]]
  simmapQ = simmap$Q
  SimmapStates = sort(colnames(simmap$mapped.edge))
  py = c("black", "orange", "#2297E6")
  pynamed <- py[1:length(SimmapStates)]
  names(pynamed) <- SimmapStates
  plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)
  legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
  
  # plot FS 2
  simmap = Coopsims[[2]]
  simmapQ = simmap$Q
  SimmapStates = sort(colnames(simmap$mapped.edge))
  py = c("#CD0BBC", "#28E2E5")
  pynamed <- py[1:length(SimmapStates)]
  names(pynamed) <- SimmapStates
  plotSimmap(simmap, fsize=0.2, lwd = 0.8, colors = pynamed)  
  legend("bottomleft",legend = names(pynamed), lwd=1,col=pynamed, lty = c(rep(1,length(SimmapStates)),2), cex=1) 
dev.off()




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
