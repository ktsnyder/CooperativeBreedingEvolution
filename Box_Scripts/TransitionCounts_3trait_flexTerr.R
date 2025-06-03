# This update to TransitionCounts_3trait.R runs when using the binary Territory (or another binary trait) as the third trait multisimmap 


library(phytools)
library(dplyr)
library(tidyr)
library(stringr)
library(emmeans)
library(ggpubr)


Qdata = df = read.csv('/Users/kate/Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/scratchData_TerrWeakStrong_wJiayingDuet2025-03-04.csv')
Qdata = df = read.csv('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-05-07.csv')

Qdata = df = read.csv("Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-04-29.csv")
Qtree = tree = read.nexus('ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')
colnames(df)
columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrong")
columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrongHighConf")
columns = c("HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrongHighConf")

setwd('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/')


### get rates ----


require(stringr)
source("findQrates.R")
cooprates <- findQrates(columns = columns[1], newdata = Qdata, newtree = Qtree)
coopQ <- cooprates$qrates
coopQ01 <- coopQ[3]
coopQ10 <- coopQ[2]
coopAnc = str_remove(cooprates$ARDlikanc, "ARDlik.anc ")  # added this for sim.history()
coopAnc = as.numeric(coopAnc)
names(coopAnc) <- c("0","1")
trait1StateLabels = getLabels(columns[1])

FSrates <- findQrates(columns = columns[2], newdata = Qdata, newtree = Qtree)
FSQ <- FSrates$qrates
FSQAbsPres <- FSQ[3]
FSQPresAbs <- FSQ[2]
FSAnc = str_remove(FSrates$ARDlikanc, "ARDlik.anc ")  # added this for sim.history()
FSAnc = as.numeric(FSAnc)
names(FSAnc) <- c("0","1")
trait2StateLabels = getLabels(columns[2])

multistateTrait = columns[3]
subsetTerr = subsettreedata(columns = multistateTrait, newdata = df, newtree = tree)
subsetTerrdf = subsetTerr$subsetdf
subsetTerrtree = subsetTerr$subsettree
discretetraitvecDisc = subsetTerrdf[,multistateTrait]
names(discretetraitvecDisc) = subsetTerrdf$species
ARDmodel <- ace(discretetraitvecDisc,subsetTerrtree, type="discrete",model = "ARD")
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

for (k in 1:nGroups) {
  for (j in 1:nGroups) {
    index = rate_index_matrix[k,j]
    if (!is.na(index)) {
      rate_matrix[k,j] = aceARDrates$rates[which(aceARDrates$rate_index == index)]
    }
  }
}
diagvals = rowSums(rate_matrix, na.rm = T)*-1
diag(rate_matrix) <- diagvals
print(rate_matrix)

#### subset data ----
subsets <- subsettreedata(columns = columns, newdata = df, newtree = tree)
subsetdf <- subsets$subsetdf
subsettree <- subsets$subsettree
FSvec <- subsetdf[,columns[2]]
names(FSvec) <- subsetdf$species
Coopvec <- subsetdf[,columns[1]]
names(Coopvec) <- subsetdf$species
Terrvec <- subsetdf[,columns[3]]
names(Terrvec) <- subsetdf$species

nsims = 500
FSsimtrees <- make.simmap(tree = subsettree, x = FSvec, model = "ARD", nsim = nsims, Q = FSQ)
Coopsimtrees <- make.simmap(tree = subsettree, x = Coopvec, model = "ARD", nsim = nsims, Q = coopQ)
Terrsimtrees <- make.simmap(tree = subsettree, x = Terrvec, model = "ARD", nsim = nsims, Q = rate_matrix)

dfout = getTransitionStateCounts3(Coopsimtrees = Coopsimtrees, FSsimtrees = FSsimtrees, Terrsimtrees = Terrsimtrees)
write.csv(dfout, paste0("Simmap Overlap Outputs/", Sys.Date(), "_transition-counts", columns[1], " ", columns[2], " ", columns[3], "_", nsims, "sims.csv"), row.names = F)
#dfout = read.csv('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Simmap Overlap Outputs/transition-counts_Coop_FS_Territory_1000sims.csv')
outlist = processTransitionStats3(dfout = dfout)

TransCountBoxplots = outlist$TransCountBoxplots
logTransCountBoxplots = outlist$logTransCountBoxplots
plotTerr1 = outlist$terr1Transitions
plotTerr2 = outlist$terr2Transitions
plotTerr3 = outlist$terr3Transitions

outlist$Terr0Transitions$transitionplot_GrayNS
outlist$Terr0Transitions$transition_plot

outlist$Terr1Transitions$transitionplot_GrayNS
outlist$Terr1Transitions$transition_plot

plotTerr0 = outlist$Terr0Transitions
plotTerr1 = outlist$Terr1Transitions


#pdf("test Coop FS Territory trans counts boxplot 500sims.pdf", width = 24, height = 6)
#print(TransCountBoxplots)
#dev.off()

# Each of these (plotTerr1, plotTerr2, plotTerr3) is a list where:
#  - $transition_df = processed data
#  - $transition_plot = your ggplot object
#  - $transitionplot_GrayNS = ggplot with non-sig arrows in gray (if p-value is used)

# To display the main color-coded plots:
print(plotTerr1$transitionplot_GrayNS)
print(plotTerr2$transitionplot_GrayNS)
print(plotTerr3$transitionplot_GrayNS)

plotlist = list(TransCountBoxplots, logTransCountBoxplots, plotTerr1$transitionplot_GrayNS, plotTerr1$transition_plot, plotTerr2$transitionplot_GrayNS, plotTerr2$transition_plot, plotTerr3$transitionplot_GrayNS, plotTerr3$transition_plot)

plotlist = list(TransCountBoxplots, logTransCountBoxplots, plotTerr0$transitionplot_GrayNS, plotTerr0$transition_plot, plotTerr1$transitionplot_GrayNS, plotTerr1$transition_plot)

require(gridExtra)
require(grid)
nTransPlots= length(plotlist)
nTransPlots

layout_matrix <- rbind(
  c(1, 1),      # First plot spans both columns
  c(2, 2),      # Second plot spans both columns
  c(3, 4),      # Third and fourth plots side by side
  c(5, 6),      # Fifth plot in first column of last row
  c(7, 8)      # Fifth plot in first column of last row
)

grobs_with_margins <- lapply(plotlist, function(plot) {
  plot_with_margin <- plot + 
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")) # Adjust margins as needed
  ggplotGrob(plot_with_margin)
})

# Arrange the grobs using the layout matrix
single_page_plot <- arrangeGrob(
  grobs = grobs_with_margins,
  layout_matrix = layout_matrix,
  heights = unit(c(1, 1, 1, 1,1), "null")  # Adjust relative heights as needed
)

# Save the arranged plot to a file
ggsave(paste0(Sys.Date(), " ", nsims, "sim ", columns[1], " ", columns[2], " ", columns[3], " transition grobs.pdf"), single_page_plot, width = 15, height = 22, units = "in")


#### Helper functions below here ----

# Currently, Territory can be either 2 state or 3 state, but both are hard-coded. It would be better to have more flexible coding to accommodate any number of states. 
# It assumes that each trait’s simmaps have been created (with make.simmap())
# and that getSimmapSegments() is available 
getTransitionStateCounts3 <- function(Coopsimtrees, FSsimtrees, Terrsimtrees) {
  # use the minimum number of trees available across all three traits
  nTrees <- min(length(Coopsimtrees), length(FSsimtrees), length(Terrsimtrees))
  
  # First determine if territory is 2-state or 3-state
  terr_states <- sort(unique(unlist(sapply(Terrsimtrees[[1]]$maps, names))))
  is_terr_2state <- length(terr_states) == 2
  
  ## Predefine all possible keys for the observed transitions 
  # For Coop transitions: the focal event is a change in Coop (two types: 0to1 and 1to0)
  # and we record the state of FS (0 or 1) and Terr (1,2,3 or 0,1) at that moment.
  if (is_terr_2state) {
    terr_states_for_labels <- 0:1
    coop_keys <- c(paste0("Coop0to1inFS", rep(c(0,1), each=2), "Terr", rep(terr_states_for_labels, times=2)),
                   paste0("Coop1to0inFS", rep(c(0,1), each=2), "Terr", rep(terr_states_for_labels, times=2)))
    
    # For FS transitions (focal event FS change) we record the state of Coop and Terr:
    fs_keys <- c(paste0("FS0to1inCoop", rep(c(0,1), each=2), "Terr", rep(terr_states_for_labels, times=2)),
                 paste0("FS1to0inCoop", rep(c(0,1), each=2), "Terr", rep(terr_states_for_labels, times=2)))
    
    # For Terr transitions (focal event: change among territorial states) - only one transition type for 2-state
    terr_transition_types <- c("0to1", "1to0")
  } else {
    terr_states_for_labels <- 1:3  # Original 3-state territory
    coop_keys <- c(paste0("Coop0to1inFS", rep(c(0,1), each=3), "Terr", rep(terr_states_for_labels, times=2)),
                   paste0("Coop1to0inFS", rep(c(0,1), each=3), "Terr", rep(terr_states_for_labels, times=2)))
    
    # For FS transitions (focal event FS change) we record the state of Coop and Terr:
    fs_keys <- c(paste0("FS0to1inCoop", rep(c(0,1), each=3), "Terr", rep(terr_states_for_labels, times=2)),
                 paste0("FS1to0inCoop", rep(c(0,1), each=3), "Terr", rep(terr_states_for_labels, times=2)))
    
    # For Terr transitions (focal event: change among territorial states) we use the six possible 
    # pairs (e.g. 1to2, 2to1, 1to3, 3to1, 2to3, 3to2)
    terr_transition_types <- c("1to2", "2to1", "1to3", "3to1", "2to3", "3to2")
  }
  
  # For Terr transitions we record the state of Coop and FS
  terr_keys <- c()
  for (tt in terr_transition_types) {
    for (c in 0:1) {
      for (f in 0:1) {
        terr_keys <- c(terr_keys, paste0("Terr", tt, "inCoop", c, "FS", f))
      }
    }
  }
  
  source("MapOverlapThree.R")
  
  # Prepare to store one row of output per simmap triple:
  results_list <- list()
  
  # Loop over each simmap triple (each "simulation")
  for (t in 1:nTrees) {
    
    FSsim = FSsimtrees[[t]]
    Coopsim = Coopsimtrees[[t]]
    Terrsim = Terrsimtrees[[t]]
    OverlapMatCoopFS <- Map.Overlap(FSsim, Coopsim)
    ObsPropsCoopFS <- c(
      ObsPropCoop0FS0 = OverlapMatCoopFS[1,1],
      ObsPropCoop1FS0 = OverlapMatCoopFS[1,2],
      ObsPropCoop0FS1 = OverlapMatCoopFS[2,1],
      ObsPropCoop1FS1 = OverlapMatCoopFS[2,2]
    )
    
    OverlapMatCoopTerr <- Map.Overlap(Terrsim, Coopsim)
    
    # Handle matrix based on number of territory states
    if (is_terr_2state) {
      ObsPropsCoopTerr <- c(
        ObsPropCoop0Terr0 = OverlapMatCoopTerr[1,1],
        ObsPropCoop0Terr1 = OverlapMatCoopTerr[2,1],
        ObsPropCoop1Terr0 = OverlapMatCoopTerr[1,2],
        ObsPropCoop1Terr1 = OverlapMatCoopTerr[2,2]
      )
    } else {
      ObsPropsCoopTerr <- c(
        ObsPropCoop0Terr1 = OverlapMatCoopTerr[1,1],
        ObsPropCoop0Terr2 = OverlapMatCoopTerr[2,1],
        ObsPropCoop0Terr3 = OverlapMatCoopTerr[3,1],
        ObsPropCoop1Terr1 = OverlapMatCoopTerr[1,2],
        ObsPropCoop1Terr2 = OverlapMatCoopTerr[2,2],
        ObsPropCoop1Terr3 = OverlapMatCoopTerr[3,2]
      )
    }
    
    OverlapMatFSTerr <- Map.Overlap(Terrsim, FSsim)
    
    # Handle matrix based on number of territory states
    if (is_terr_2state) {
      ObsPropsFSTerr <- c(
        ObsPropFS0Terr0 = OverlapMatFSTerr[1,1],
        ObsPropFS0Terr1 = OverlapMatFSTerr[2,1],
        ObsPropFS1Terr0 = OverlapMatFSTerr[1,2],
        ObsPropFS1Terr1 = OverlapMatFSTerr[2,2]
      )
    } else {
      ObsPropsFSTerr <- c(
        ObsPropFS0Terr1 = OverlapMatFSTerr[1,1],
        ObsPropFS0Terr2 = OverlapMatFSTerr[2,1],
        ObsPropFS0Terr3 = OverlapMatFSTerr[3,1],
        ObsPropFS1Terr1 = OverlapMatFSTerr[1,2],
        ObsPropFS1Terr2 = OverlapMatFSTerr[2,2],
        ObsPropFS1Terr3 = OverlapMatFSTerr[3,2]
      )
    }
    
    OverlapCoopFSTerr = Map.Overlap.Three(FSsim, Coopsim, Terrsim)
    
    # Handle 3D matrix based on number of territory states
    if (is_terr_2state) {
      ObsPropsAll <- c(
        # Terr0 (array slice 1)
        ObsPropFS0Coop0Terr0 = OverlapCoopFSTerr[1,1,1],  # [FS=0, Coop=0, Terr=0]
        ObsPropFS0Coop1Terr0 = OverlapCoopFSTerr[1,2,1],  # [FS=0, Coop=1, Terr=0]
        ObsPropFS1Coop0Terr0 = OverlapCoopFSTerr[2,1,1],  # [FS=1, Coop=0, Terr=0]
        ObsPropFS1Coop1Terr0 = OverlapCoopFSTerr[2,2,1],  # [FS=1, Coop=1, Terr=0]
        
        # Terr1 (array slice 2)
        ObsPropFS0Coop0Terr1 = OverlapCoopFSTerr[1,1,2],
        ObsPropFS0Coop1Terr1 = OverlapCoopFSTerr[1,2,2],
        ObsPropFS1Coop0Terr1 = OverlapCoopFSTerr[2,1,2],
        ObsPropFS1Coop1Terr1 = OverlapCoopFSTerr[2,2,2]
      )
    } else {
      ObsPropsAll <- c(
        # Terr1 (array slice 1)
        ObsPropFS0Coop0Terr1 = OverlapCoopFSTerr[1,1,1],  # [FS=0, Coop=0, Terr=1]
        ObsPropFS0Coop1Terr1 = OverlapCoopFSTerr[1,2,1],  # [FS=0, Coop=1, Terr=1]
        ObsPropFS1Coop0Terr1 = OverlapCoopFSTerr[2,1,1],  # [FS=1, Coop=0, Terr=1]
        ObsPropFS1Coop1Terr1 = OverlapCoopFSTerr[2,2,1],  # [FS=1, Coop=1, Terr=1]
        
        # Terr2 (array slice 2)
        ObsPropFS0Coop0Terr2 = OverlapCoopFSTerr[1,1,2],
        ObsPropFS0Coop1Terr2 = OverlapCoopFSTerr[1,2,2],
        ObsPropFS1Coop0Terr2 = OverlapCoopFSTerr[2,1,2],
        ObsPropFS1Coop1Terr2 = OverlapCoopFSTerr[2,2,2],
        
        # Terr3 (array slice 3)
        ObsPropFS0Coop0Terr3 = OverlapCoopFSTerr[1,1,3],
        ObsPropFS0Coop1Terr3 = OverlapCoopFSTerr[1,2,3],
        ObsPropFS1Coop0Terr3 = OverlapCoopFSTerr[2,1,3],
        ObsPropFS1Coop1Terr3 = OverlapCoopFSTerr[2,2,3]
      )
    }
    
    # Get the overall proportion of branch lengths in each state for each trait
    FS_desc <- describe.simmap(FSsim)
    ObsPropFS0 <- FS_desc$times["prop", "0"]
    ObsPropFS1 <- FS_desc$times["prop", "1"]
    
    Coop_desc <- describe.simmap(Coopsim)
    ObsPropCoop0 <- Coop_desc$times["prop", "0"]
    ObsPropCoop1 <- Coop_desc$times["prop", "1"]
    
    Terr_desc <- describe.simmap(Terrsim)
    
    # Create territory proportions based on available states
    if (is_terr_2state) {
      ObsPropTerr0 <- Terr_desc$times["prop", "0"]
      ObsPropTerr1 <- Terr_desc$times["prop", "1"]
      ObsPropTerr <- c(ObsPropTerr0 = ObsPropTerr0, ObsPropTerr1 = ObsPropTerr1)
    } else {
      ObsPropTerr1 <- Terr_desc$times["prop", "1"]
      ObsPropTerr2 <- Terr_desc$times["prop", "2"]
      ObsPropTerr3 <- Terr_desc$times["prop", "3"]
      ObsPropTerr <- c(ObsPropTerr1 = ObsPropTerr1, ObsPropTerr2 = ObsPropTerr2, ObsPropTerr3 = ObsPropTerr3)
    }
    
    # Get segment and transition dataframes (using your helper function)
    FS_segs   <- getSimmapSegments(FSsim)
    FS_trans  <- FS_segs$transitions_df
    FS_seg    <- FS_segs$segments_df
    
    Coop_segs <- getSimmapSegments(Coopsim)
    Coop_trans<- Coop_segs$transitions_df
    Coop_seg  <- Coop_segs$segments_df
    
    Terr_segs <- getSimmapSegments(Terrsim)
    Terr_trans<- Terr_segs$transitions_df
    Terr_seg  <- Terr_segs$segments_df
    
    ## Initialize counters (all possible keys are preset to 0)
    coop_counts <- setNames(rep(0, length(coop_keys)), coop_keys)
    fs_counts   <- setNames(rep(0, length(fs_keys)), fs_keys)
    terr_counts <- setNames(rep(0, length(terr_keys)), terr_keys)
    
    ## (A) Count Coop transitions
    if(nrow(Coop_trans) > 0) {
      for (i in 1:nrow(Coop_trans)) {
        row_i <- Coop_trans[i, ]
        if (is.na(row_i$Transition01or10)) next  # skip if no change
        # Get the FS state at the precise location on the branch:
        fs_state <- FS_seg$State[ FS_seg$EdgeNumber == row_i$EdgeNumber & 
                                    FS_seg$Start < row_i$DistanceFromRootwardNode & 
                                    FS_seg$End   >= row_i$DistanceFromRootwardNode ]
        if(length(fs_state) == 0) next
        # Get the Terr state at that same point:
        terr_state <- Terr_seg$State[ Terr_seg$EdgeNumber == row_i$EdgeNumber & 
                                        Terr_seg$Start < row_i$DistanceFromRootwardNode & 
                                        Terr_seg$End   >= row_i$DistanceFromRootwardNode ]
        if(length(terr_state) == 0) next
        key <- paste0("Coop", row_i$Transition01or10, "inFS", fs_state[1], "Terr", terr_state[1])
        if(key %in% names(coop_counts)) {
          coop_counts[key] <- coop_counts[key] + 1
        }
      }
    }
    
    ## (B) Count FS transitions
    if(nrow(FS_trans) > 0) {
      for (i in 1:nrow(FS_trans)) {
        row_i <- FS_trans[i, ]
        if (is.na(row_i$Transition01or10)) next
        # Look up Coop state at the moment of the FS change:
        coop_state <- Coop_seg$State[ Coop_seg$EdgeNumber == row_i$EdgeNumber & 
                                        Coop_seg$Start < row_i$DistanceFromRootwardNode & 
                                        Coop_seg$End   >= row_i$DistanceFromRootwardNode ]
        if(length(coop_state)==0) next
        # And Terr state:
        terr_state <- Terr_seg$State[ Terr_seg$EdgeNumber == row_i$EdgeNumber & 
                                        Terr_seg$Start < row_i$DistanceFromRootwardNode & 
                                        Terr_seg$End   >= row_i$DistanceFromRootwardNode ]
        if(length(terr_state)==0) next
        key <- paste0("FS", row_i$Transition01or10, "inCoop", coop_state[1], "Terr", terr_state[1])
        if(key %in% names(fs_counts)) {
          fs_counts[key] <- fs_counts[key] + 1
        }
      }
    }
    
    ## (C) Count Terr transitions
    if(nrow(Terr_trans) > 0) {
      for (i in 1:nrow(Terr_trans)) {
        row_i <- Terr_trans[i, ]
        if (is.na(row_i$Transition01or10)) next
        # Get the state of Coop and FS at the transition:
        coop_state <- Coop_seg$State[ Coop_seg$EdgeNumber == row_i$EdgeNumber & 
                                        Coop_seg$Start < row_i$DistanceFromRootwardNode & 
                                        Coop_seg$End   >= row_i$DistanceFromRootwardNode ]
        if(length(coop_state)==0) next
        fs_state <- FS_seg$State[ FS_seg$EdgeNumber == row_i$EdgeNumber & 
                                    FS_seg$Start < row_i$DistanceFromRootwardNode & 
                                    FS_seg$End   >= row_i$DistanceFromRootwardNode ]
        if(length(fs_state)==0) next
        key <- paste0("Terr", row_i$Transition01or10, "inCoop", coop_state[1], "FS", fs_state[1])
        if(key %in% names(terr_counts)) {
          terr_counts[key] <- terr_counts[key] + 1
        }
      }
    }
    
    ## Compute overall totals for each focal transition type (for expected counts)
    TotalCoop0to1 <- sum(coop_counts[grep("^Coop0to1", names(coop_counts))])
    TotalCoop1to0 <- sum(coop_counts[grep("^Coop1to0", names(coop_counts))])
    TotalFS0to1   <- sum(fs_counts[grep("^FS0to1", names(fs_counts))])
    TotalFS1to0   <- sum(fs_counts[grep("^FS1to0", names(fs_counts))])
    
    # For Terr transitions, get totals per transition type (e.g. "1to2" etc.)
    terr_totals <- list()
    for(tt in terr_transition_types) {
      key_pattern <- paste0("^Terr", tt)
      terr_totals[[tt]] <- sum(terr_counts[grep(key_pattern, names(terr_counts))])
    }
    
    ## Compute expected counts for each transition key 
    ## (expected = (total number of that focal transition) *
    ##              (proportion of time the other trait 1 is in state X) *
    ##              (proportion of time the other trait 2 is in state Y) ).
    
    # For Coop transitions:
    coop_expected <- coop_counts
    for (key in names(coop_counts)) {
      total_dir <- if (startsWith(key, "Coop0to1")) { TotalCoop0to1 } else
        if (startsWith(key, "Coop1to0")) { TotalCoop1to0 } else { NA }
      fs_state_val <- as.numeric(str_extract(key, "(?<=inFS)[0-9]"))
      terr_state_val <- as.numeric(str_extract(key, "(?<=Terr)[0-9]+"))
      if (!is.na(total_dir) && !is.na(fs_state_val) && !is.na(terr_state_val)) {
        overlap_key <- paste0("ObsPropFS", fs_state_val, "Terr", terr_state_val) #
        coop_expected[key] <- total_dir * ObsPropsFSTerr[overlap_key] #
      } else {
        coop_expected[key] <- NA
      }
    }
    
    # For FS transitions:
    fs_expected <- fs_counts
    for (key in names(fs_counts)) {
      total_dir <- if (startsWith(key, "FS0to1")) { TotalFS0to1 } else
        if (startsWith(key, "FS1to0")) { TotalFS1to0 } else { NA }
      coop_state_val <- as.numeric(str_extract(key, "(?<=inCoop)[0-9]"))
      terr_state_val <- as.numeric(str_extract(key, "(?<=Terr)[0-9]+"))
      if (!is.na(total_dir) && !is.na(coop_state_val) && !is.na(terr_state_val)) {
        overlap_key <- paste0("ObsPropCoop", coop_state_val, "Terr", terr_state_val) #
        fs_expected[key] <- total_dir * ObsPropsCoopTerr[overlap_key] #
      } else {
        fs_expected[key] <- NA
      }
    }
    
    # For Terr transitions:
    terr_expected <- terr_counts
    for (key in names(terr_counts)) {
      # Extract the Terr transition type (e.g., "Terr1to2")
      terr_trans_type <- str_extract(key, "Terr[0-9]+to[0-9]+")
      total_dir <- terr_totals[[str_remove(terr_trans_type, "Terr")]]
      coop_state_val <- as.numeric(str_extract(key, "(?<=inCoop)[0-9]"))
      fs_state_val <- as.numeric(str_extract(key, "(?<=FS)[0-9]"))
      if (!is.na(total_dir) && !is.na(coop_state_val) && !is.na(fs_state_val)) {
        overlap_key <- paste0("ObsPropCoop", coop_state_val, "FS", fs_state_val) #
        terr_expected[key] <- total_dir * ObsPropsCoopFS[overlap_key] #
      } else {
        terr_expected[key] <- NA
      }
    }
    
    ## Combine all results for this simulation into one named vector.
    # (We include the observed counts and, with "Expected" appended, the computed expectations.)
    outvec <- c(TreeNum = t,
                # Observed counts:
                coop_counts,
                fs_counts,
                terr_counts,
                # Expected counts:
                setNames(coop_expected, paste0(names(coop_expected), "Expected")),
                setNames(fs_expected,   paste0(names(fs_expected),   "Expected")),
                setNames(terr_expected, paste0(names(terr_expected), "Expected")),
                # Also include the overall proportions (for later checking):
                ObsPropCoop0 = ObsPropCoop0, ObsPropCoop1 = ObsPropCoop1,
                ObsPropFS0 = ObsPropFS0, ObsPropFS1 = ObsPropFS1,
                ObsPropTerr, 
                ObsPropsCoopTerr,
                ObsPropsFSTerr,
                ObsPropsCoopFS,
                ObsPropsAll
    )
    
    results_list[[t]] <- outvec
  }
  
  ## Combine all simulation results into one data.frame.
  results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  
  ## Convert any numeric columns (they may be character because of the rbind)
  results_df[] <- lapply(results_df, function(col) as.numeric(as.character(col)))
  
  # Add flag to indicate if this is 2-state or 3-state territory
  attr(results_df, "is_terr_2state") <- is_terr_2state
  attr(results_df, "terr_states") <- terr_states
  
  return(results_df)
}

processTransitionStats3 <- function(dfout) {
  trait1 = "HighConfidence_Coop"
  trait2 = "FemaleSong_Agg01"
  trait1StateLabels = getLabels(trait1)
  trait2StateLabels = getLabels(trait2)
  
  # Get info about territory states
  is_terr_2state <- attr(dfout, "is_terr_2state")
  if(is.null(is_terr_2state)) {
    # Try to detect if it's 2-state territory by checking column names
    is_terr_2state <- any(grepl("ObsPropTerr0", names(dfout)))
  }
  
  
  # Identify transition count columns:
  # (We assume that the observed transitions have names starting with
  #  "Coop", "FS", or "Terr" but do not include "Prop" or "TreeNum".)
  
  transition_cols <- grep("^(Coop|FS|Terr)(?!.*Prop)", names(dfout),
                          value = TRUE, perl = TRUE)
  transition_cols <- setdiff(transition_cols, "TreeNum")
  transition_cols <- grep("^Terr", transition_cols, invert = T, value = T, perl = T)
  
  
  
  # Gather observed and expected transition counts into long format.
  # (Note: In getTransitionStateCounts3() we appended expected counts
  #  by adding "Expected" to the observed key names.)
  
  dfMeltCounts <- dfout %>%
    gather(key = "Transition", value = "Count", all_of(transition_cols))
  
  # Create an indicator for whether a given row is from an observed column or an expected column.
  dfMeltCounts$ObservedVsExpected <- ifelse(grepl("Expected", dfMeltCounts$Transition),
                                            "Expected", "Observed")
  
  # Convert Count to numeric and replace any 0 with a small nonzero value (0.001)
  dfMeltCounts$Count <- as.numeric(dfMeltCounts$Count)
  dfMeltCounts$Count[dfMeltCounts$Count == 0] <- 0.001
  
  # Save the original column name as a label.
  dfMeltCounts$Label <- dfMeltCounts$Transition
  
  # Remove the string "Expected" from the Transition name so that the same
  # transition is identified regardless of whether it is observed or expected.
  dfMeltCounts$Transition <- gsub("Expected", "", dfMeltCounts$Transition)
  
  # Compute the natural log of the count (for LM analysis)
  dfMeltCounts$logCount <- log(dfMeltCounts$Count)
  
  
  # Linear model on log-transformed counts
  lmLogMult <- lm(logCount ~ ObservedVsExpected * Transition, data = dfMeltCounts)
  logLmANOVA <- anova(lmLogMult)
  pvalInteractionLog <- logLmANOVA$`Pr(>F)`[3]
  print(pvalInteractionLog)
  if (pvalInteractionLog < 0.0001) {
    pvalInteractionLogLabel <- "p < 0.0001"
  } else {
    pvalInteractionLogLabel <- paste("p =", round(pvalInteractionLog, digits = 5))
  }
  
  # Pairwise post-hoc tests using emmeans for the log-transformed model
  emmLog <- emmeans(lmLogMult, pairwise ~ ObservedVsExpected | Transition)
  contrastLog <- emmLog$contrasts
  logPairwisePostHoc <- summary(contrastLog, adjust = "tukey")
  
  
  # Linear model on raw (non-transformed) counts
  lmMult <- lm(Count ~ ObservedVsExpected * Transition, data = dfMeltCounts)
  LmANOVA <- anova(lmMult)
  pvalInteraction <- LmANOVA$`Pr(>F)`[3]
  print(pvalInteraction)
  if (pvalInteraction < 0.0001) {
    pvalInteractionLabel <- "p < 0.0001"
  } else {
    pvalInteractionLabel <- paste("p =", round(pvalInteraction, digits = 5))
  }
  
  emm <- emmeans(lmMult, pairwise ~ ObservedVsExpected | Transition)
  contrast <- emm$contrasts
  PairwisePostHoc <- summary(contrast, adjust = "tukey")
  
  
  # Count how many simulations (rows of dfout) had observed counts > expected.
  # We use the original dfout to compare for each transition.
  
  # Get the "base" keys (the ones without "Expected")
  observed_keys <- grep("Expected", names(dfout), invert = TRUE, value = TRUE)
  observed_keys <- observed_keys[grepl("^(Coop|FS|Terr)", observed_keys)]
  expected_keys <- paste0(observed_keys, "Expected")
  
  # In case some keys are missing, restrict to those present.
  observed_keys <- intersect(observed_keys, names(dfout))
  expected_keys <- intersect(expected_keys, names(dfout))
  
  observed_keys = observed_keys[which(observed_keys %in% gsub("Expected", "", expected_keys))]
  
  # Compare observed and expected counts (row‐by‐row)
  timesActualGreaterThanExpected <- dfout[, observed_keys] > dfout[, expected_keys]
  timesActualGreaterThanExpected <- as.data.frame(timesActualGreaterThanExpected)
  NtimesActualGreaterThanExpected <- colSums(timesActualGreaterThanExpected, na.rm = TRUE)
  
  
  # Bundle LM and posthoc results into a TransitionStats list
  TransitionStats <- list(
    logLmANOVA = logLmANOVA,
    logPairwisePostHoc = logPairwisePostHoc,
    LmANOVA = LmANOVA,
    PairwisePostHoc = PairwisePostHoc,
    NtimesActualGreaterThanExpected = NtimesActualGreaterThanExpected
  )
  
  # Set the factor levels for Transition in the melted data frame.
  dfMeltCounts$Transition <- factor(dfMeltCounts$Transition,
                                    levels = unique(dfMeltCounts$Transition))
  
  nsims_real <- nrow(dfout)
  
  
  # Create group labels for plotting: each label gives the transition
  # name and the number of sims where observed > expected.
  groupLabels <- sapply(levels(dfMeltCounts$Transition), function(x) {
    count_val <- NtimesActualGreaterThanExpected[x]
    if (is.na(count_val)) count_val <- 0
    paste0(x, "\nNsims Actual greater than\nExpected: ", count_val, "/", nsims_real)
  })
  
  
  # Additional processing for transition plots: make a data frame of p-values.
  # Here we use the log-transformed pairwise post-hoc results.
  # (This part may need tweaking depending on the exact structure of your contrast output.)
  pvaldf <- as.data.frame(logPairwisePostHoc)
  
  # Set significance labels based on p.value thresholds.
  pvaldf$SignificanceLabel <- "n.s."
  pvaldf$SignificanceLabel[pvaldf$p.value < 0.05] <- "p < 0.05"
  pvaldf$SignificanceLabel[pvaldf$p.value < 0.01] <- "p < 0.01"
  pvaldf$SignificanceLabel[pvaldf$p.value < 0.001] <- "p < 0.001"
  pvaldf$SignificanceLabel[pvaldf$p.value < 0.0001] <- "p < 0.0001"
  
  # Merge in the counts of how many simulations had observed > expected,
  # so that each transition (if present) carries that number.
  Ntimes_df <- data.frame(
    Transition = names(NtimesActualGreaterThanExpected),
    CountNActualGreaterThanExpected = as.integer(NtimesActualGreaterThanExpected),
    stringsAsFactors = FALSE
  )
  pvaldf <- merge(pvaldf, Ntimes_df, by = "Transition", all.x = TRUE)
  pvaldf$FractionActualGreaterThanExpected <- pvaldf$CountNActualGreaterThanExpected / nsims_real
  
  TransCountBoxplots <- ggplot(dfMeltCounts, aes(x = Transition, y = Count, fill = ObservedVsExpected)) +
    geom_boxplot(outlier.shape = NA) + # Exclude outliers
    theme_minimal() +
    labs(y = "Transition Counts", x = "", fill = "Observed/Expected") +
    scale_fill_manual(values = c("Observed" = "#762a83", "Expected" = "#1b7837")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), title = element_text(size = 8)) +
    ggtitle(paste(trait1, trait2, "\nnSimsObserved =", nsims_real, "    nSimsExpected =", nsims_real, "    Obs/Exp:TransCounts", pvalInteractionLogLabel)) +
    stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
    scale_x_discrete(labels = groupLabels)
  
  # Plot log
  logTransCountBoxplots <- ggplot(dfMeltCounts, aes(x = Transition, y = logCount, fill = ObservedVsExpected)) +
    geom_boxplot(outlier.shape = NA) + # Exclude outliers
    theme_minimal() +
    labs(y = "log(Transition Counts)", x = "", fill = "Observed/Expected") +
    scale_fill_manual(values = c("Observed" = "#762a83", "Expected" = "#1b7837")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), title = element_text(size = 8)) +
    ggtitle(paste(trait1, trait2, "\nnSimsObserved =", nsims_real, "    nSimsExpected =", nsims_real, "    Obs/Exp:TransCounts", pvalInteractionLogLabel)) +
    stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
    scale_x_discrete(labels = groupLabels)
  
  
  # Run transition plots for available territory states
  plot_list <- list()
  
  # Determine which territory suffixes to use based on the data
  if(is_terr_2state) {
    terr_suffixes <- c("Terr0", "Terr1")
  } else {
    terr_suffixes <- c("Terr1", "Terr2", "Terr3")
  }
  
  # Create plots for each territory state
  for(terr_suffix in terr_suffixes) {
    # Check if we have data for this territory state
    if(any(grepl(paste0(terr_suffix, "$"), names(dfout)))) {
      plot_title <- paste("Transition counts within", terr_suffix)
      plot_list[[terr_suffix]] <- make_transition_plot_for_terr(
        dfout = dfout, 
        pvaldf = pvaldf, 
        terrSuffix = terr_suffix, 
        trait1StateLabels = trait1StateLabels, 
        trait2StateLabels = trait2StateLabels, 
        plottitle = plot_title
      )
    }
  }
  
  
  # Return a list of outputs for further inspection and plotting.
  
  result_list <- list(
    dfMeltCounts = dfMeltCounts,
    TransitionStats = TransitionStats,
    groupLabels = groupLabels,
    pvaldf = pvaldf,
    TransCountBoxplots = TransCountBoxplots,
    logTransCountBoxplots = logTransCountBoxplots,
    is_terr_2state = is_terr_2state
  )
  
  # Add the territory transition plots
  for(terr_suffix in names(plot_list)) {
    result_list[[paste0(terr_suffix, "Transitions")]] <- plot_list[[terr_suffix]]
  }
  
  return(result_list)
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

make_transition_plot_for_terr <- function(
    dfout,               # your main data frame with transition counts
    pvaldf,              # your p-value data frame
    terrSuffix,          # e.g. "Terr1", "Terr2", or "Terr3"
    trait1StateLabels = c("0","1"),   # e.g. c("0","1") for Coop
    trait2StateLabels = c("0","1"),   # e.g. c("0","1") for FS
    plottitle            # title for the resulting plot
) {
  source("transition_plot.R")
  
  RateRef = as.data.frame(rbind(c("FS0to1inCoop0", "q12"),c("Coop0to1inFS0", "q13"),c("FS1to0inCoop0", "q21"), c("Coop0to1inFS1", "q24"),c("Coop1to0inFS0", "q31"), c("FS0to1inCoop1", "q34"), c("Coop1to0inFS1", "q42"),c("FS1to0inCoop1", "q43")))
  colnames(RateRef) <- c("FullTransitions", "Transitions")
  RateRef$qRate = RateRef$Transitions
  
  # 1) Create a copy of RateRef and append the Terr suffix to match dfout
  RateRefTerr <- RateRef
  RateRefTerr[["FullTransitions"]] <- paste0(RateRefTerr[["FullTransitions"]], terrSuffix)
  
  ratePvalsTerr = merge(RateRefTerr, pvaldf, by.x = "FullTransitions", by.y = "Transition")
  
  # remove these lines if want to use ANOVA pval labels instead
  ratePvalsTerr$PercentTrendingLabel = "<80%"
  ratePvalsTerr$PercentTrendingLabel[which(ratePvalsTerr$FractionActualGreaterThanExpected > .8 | ratePvalsTerr$FractionActualGreaterThanExpected < .2)] <- ">80%"
  ratePvalsTerr$PercentTrendingLabel[which(ratePvalsTerr$FractionActualGreaterThanExpected > .9 | ratePvalsTerr$FractionActualGreaterThanExpected < .1)] <- ">90%"
  ratePvalsTerr$PercentTrendingLabel[which(ratePvalsTerr$FractionActualGreaterThanExpected > .95 | ratePvalsTerr$FractionActualGreaterThanExpected < .05)] <- ">95%"
  ratePvalsTerr$PercentTrendingLabel[which(ratePvalsTerr$FractionActualGreaterThanExpected > .99 | ratePvalsTerr$FractionActualGreaterThanExpected < .01)] <- ">99%"
  
  
  # Find state percentages
  cols <- colnames(dfout)
  
  cols_to_keep <- cols[
    !startsWith(cols, "ObsProp") |                              # Keep if doesn't start with "ObsProp"
      (startsWith(cols, "ObsProp") &                              # OR if starts with "ObsProp" AND
         grepl("FS", cols) &                                        # contains "FS" AND
         grepl("Coop", cols) &                                      # contains "Coop" AND
         endsWith(cols, terrSuffix))                                   # ends with terrSuffix
  ]
  
  # Create new dataframe with only the kept columns
  dfTerr <- dfout[, cols_to_keep]
  
  # Define the desired order of the ObsProp columns
  desired_order <- paste0(c(
    "ObsPropFS0Coop0",
    "ObsPropFS1Coop0",
    "ObsPropFS0Coop1",
    "ObsPropFS1Coop1"), terrSuffix
  )
  
  # Get non-ObsProp columns
  other_cols <- names(dfTerr)[!grepl("^ObsProp", names(dfTerr))]
  
  # Reorder the dataframe
  dfTerr <- dfTerr[, c(other_cols, desired_order)]
  
  
  dfTerr[,paste0(ratePvalsTerr$FullTransitions,"DifferenceFromExpected")] = dfTerr[,ratePvalsTerr$FullTransitions] - dfTerr[,paste0(ratePvalsTerr$FullTransitions,"Expected")]
  
  old_names <- paste0(ratePvalsTerr$FullTransitions,"DifferenceFromExpected")
  new_names <- ratePvalsTerr$Transitions
  rename_map <- setNames(new_names, old_names)
  rename_map
  for (i in seq_along(colnames(dfTerr))) {
    col_i <- colnames(dfTerr)[i]
    if (col_i %in% names(rename_map)) {
      colnames(dfTerr)[i] <- rename_map[col_i]
    }
  }
  
  colnames(dfTerr)
  
  # 5) Call transition_plot function
  #    In particular, make sure to pass in 'ratePvals = ratePvalsTerr'
  plot_out <- transition_plot(
    df = dfTerr,
    trait1StateLabels = trait1StateLabels,
    trait2StateLabels = trait2StateLabels,
    scale_area_by = 1,
    offset = 0.2,
    lengthen = 0.2,
    ratePvals = ratePvalsTerr,
    plottitle = plottitle,
    center = "mean"
  )
  
  return(plot_out)
}

detect_unique_states <- function(simmap_collection) {
  # Check if we have a single simmap or a collection
  if (inherits(simmap_collection, "simmap")) {
    # Handle a single simmap
    simmap_list <- list(simmap_collection)
  } else {
    # Handle a collection (list, multiSimmap, multiPhylo)
    simmap_list <- simmap_collection
  }
  
  # Container for all unique states
  all_states <- character(0)
  
  # Process each simmap in the collection
  for (i in 1:length(simmap_list)) {
    simmap <- simmap_list[[i]]
    
    # Extract states from each edge in the maps
    for (edge_map in simmap$maps) {
      edge_states <- names(edge_map)
      all_states <- union(all_states, edge_states)
    }
    
    # Optionally check mapped.edge column names if they exist
    if (!is.null(simmap$mapped.edge)) {
      if (!is.null(colnames(simmap$mapped.edge))) {
        mapped_states <- colnames(simmap$mapped.edge)
        all_states <- union(all_states, mapped_states)
      }
    }
  }
  
  # Sort the states (preserving character/numeric distinctions)
  if (all(grepl("^[0-9]+$", all_states))) {
    # If all states are numeric strings, sort numerically
    all_states <- all_states[order(as.numeric(all_states))]
  } else {
    # Otherwise, sort alphabetically
    all_states <- sort(all_states)
  }
  
  # Print a verification message
  cat("Detected states:", paste(all_states, collapse=", "), "\n")
  
  return(all_states)
}


create_multidimensional_matrix <- function(states_per_trait) {
  # Check input
  if (!is.list(states_per_trait) || length(states_per_trait) == 0) {
    stop("states_per_trait must be a non-empty list of character vectors")
  }
  
  # Get dimensions for the array
  dims <- sapply(states_per_trait, length)
  
  # Create an empty array with the right dimensions
  result <- array(0, dim = dims)
  
  # Set proper dimnames for each dimension
  dimnames_list <- states_per_trait
  names(dimnames_list) <- paste0("trait", 1:length(states_per_trait))
  dimnames(result) <- dimnames_list
  
  # Print information about the created matrix
  dim_str <- paste(dims, collapse=" × ")
  cat("Created a", dim_str, "matrix for", length(dims), "traits\n")
  
  # Show state configuration
  for (i in 1:length(states_per_trait)) {
    cat("Trait", i, "states:", paste(states_per_trait[[i]], collapse=", "), "\n")
  }
  
  return(result)
}

