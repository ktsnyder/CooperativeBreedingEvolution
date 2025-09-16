## Get transitions from two sets of overlapping simmaps 
## Created 11/22/2023
## Kate Snyder
## Edited 11/27/2023 - calculate expected counts of transitions based on total transitions and time in each state

#transStateCounts = getTransitionStateCounts(Coopsimtrees = Coopsimtrees, FSsimtrees = FSsimtrees)

getTransitionStateCounts <- function(Coopsimtrees, FSsimtrees) {
  
  nTrees = min(c(length(FSsimtrees), length(Coopsimtrees)))
  
  countsdf = set.seed(10)
  
  for (t in 1:nTrees) {
    FSsimtree1 = FSsimtrees[[t]]
    Coopsimtree1 = Coopsimtrees[[t]]
    
    FStimePerState = describe.simmap(FSsimtree1)
    FSPropTime0= FStimePerState$times["prop","0"]
    FSPropTime1= FStimePerState$times["prop","1"]
    
    CooptimePerState = describe.simmap(Coopsimtree1)
    CoopPropTime0= CooptimePerState$times["prop","0"]
    CoopPropTime1= CooptimePerState$times["prop","1"]
    
    require(dplyr)
    
    Coopsegments = getSimmapSegments(Coopsimtree1)
    FSsegments = getSimmapSegments(FSsimtree1)
    
    CoopTransitions = Coopsegments$transitions_df
    CoopSegments = Coopsegments$segments_df
    FSTransitions = FSsegments$transitions_df
    FSSegments = FSsegments$segments_df
    
    # Initialize counts for both directions
    Coop0to1inFS0 <- Coop1to0inFS0 <- Coop0to1inFS1 <- Coop1to0inFS1 <- 0
    FS0to1inCoop0 <- FS1to0inCoop0 <- FS0to1inCoop1 <- FS1to0inCoop1 <- 0
    
    # Loop through each transition in CoopTransitions
    for (i in 1:nrow(CoopTransitions)) {
      transition <- CoopTransitions[i, ]
      fs_state_at_transition <- find_other_trait_state(transition$EdgeNumber, transition$DistanceFromRootwardNode, FSSegments)
      if (is.na(transition$Transition01or10)) {
        next
      } else {
        fs_state_at_transition = FSSegments[which(FSSegments$EdgeNumber == transition$EdgeNumber & FSSegments$Start < transition$DistanceFromRootwardNode & FSSegments$End >= transition$DistanceFromRootwardNode), "State"]
      }
      
      if (!is.na(fs_state_at_transition)) {
        if (transition$Transition01or10 == "0to1" && fs_state_at_transition == 0) {
          Coop0to1inFS0 <- Coop0to1inFS0 + 1
        } else if (transition$Transition01or10 == "1to0" && fs_state_at_transition == 0) {
          Coop1to0inFS0 <- Coop1to0inFS0 + 1
        } else if (transition$Transition01or10 == "0to1" && fs_state_at_transition == 1) {
          Coop0to1inFS1 <- Coop0to1inFS1 + 1
        } else if (transition$Transition01or10 == "1to0" && fs_state_at_transition == 1) {
          Coop1to0inFS1 <- Coop1to0inFS1 + 1
        }
        #print(c(Coop0to1inFS0, Coop1to0inFS0, Coop0to1inFS1, Coop1to0inFS1))
      }
    }
    # Process FS transitions
    for (i in 1:nrow(FSTransitions)) {
      transition <- FSTransitions[i, ]
      coop_state_at_transition <- find_other_trait_state(transition$EdgeNumber, transition$DistanceFromRootwardNode, CoopSegments)
      
      if (!is.na(coop_state_at_transition)) {
        if (transition$Transition01or10 == "0to1" && coop_state_at_transition == 0) {
          FS0to1inCoop0 <- FS0to1inCoop0 + 1
        } else if (transition$Transition01or10 == "1to0" && coop_state_at_transition == 0) {
          FS1to0inCoop0 <- FS1to0inCoop0 + 1
        } else if (transition$Transition01or10 == "0to1" && coop_state_at_transition == 1) {
          FS0to1inCoop1 <- FS0to1inCoop1 + 1
        } else if (transition$Transition01or10 == "1to0" && coop_state_at_transition == 1) {
          FS1to0inCoop1 <- FS1to0inCoop1 + 1
        }
       # print(c(FS0to1inCoop0, FS1to0inCoop0, FS0to1inCoop1, FS1to0inCoop1))
      }
    }
    outvec = c(t, Coop0to1inFS0, Coop1to0inFS0, Coop0to1inFS1, Coop1to0inFS1, FS0to1inCoop0, FS1to0inCoop0, FS0to1inCoop1, FS1to0inCoop1, CoopPropTime0, CoopPropTime1, FSPropTime0, FSPropTime1)
    names(outvec) <- c("TreeNum", "Coop0to1inFS0", "Coop1to0inFS0", "Coop0to1inFS1", "Coop1to0inFS1", "FS0to1inCoop0", "FS1to0inCoop0", "FS0to1inCoop1", "FS1to0inCoop1", "CoopPropTime0", "CoopPropTime1", "FSPropTime0", "FSPropTime1")
    
    countsdf = rbind(countsdf, outvec)
    
  } # end for t in 1:nTrees
  
  countsdf = as.data.frame(countsdf)
  
  countsdf$TotalFS0to1 = countsdf$FS0to1inCoop0 + countsdf$FS0to1inCoop1
  countsdf$TotalFS1to0 = countsdf$FS1to0inCoop0 + countsdf$FS1to0inCoop1
  countsdf$TotalCoop0to1 = countsdf$Coop0to1inFS0 + countsdf$Coop0to1inFS1
  countsdf$TotalCoop1to0 = countsdf$Coop1to0inFS0 + countsdf$Coop1to0inFS1
  
  countsdf$Coop0to1inFS0Expected = countsdf$TotalCoop0to1 * countsdf$FSPropTime0
  countsdf$Coop0to1inFS1Expected = countsdf$TotalCoop0to1 * countsdf$FSPropTime1
  countsdf$Coop1to0inFS0Expected = countsdf$TotalCoop1to0 * countsdf$FSPropTime0
  countsdf$Coop1to0inFS1Expected = countsdf$TotalCoop1to0 * countsdf$FSPropTime1
  countsdf$FS0to1inCoop0Expected = countsdf$TotalFS0to1 * countsdf$CoopPropTime0
  countsdf$FS0to1inCoop1Expected = countsdf$TotalFS0to1 * countsdf$CoopPropTime1
  countsdf$FS1to0inCoop0Expected = countsdf$TotalFS1to0 * countsdf$CoopPropTime0
  countsdf$FS1to0inCoop1Expected = countsdf$TotalFS1to0 * countsdf$CoopPropTime1
  
  return(countsdf)
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

# Function to find the state of the other trait during a transition
find_other_trait_state <- function(edge, transition_distance, other_segments) {
  matching_segment <- other_segments[other_segments$Edge == edge & 
                                       other_segments$Start <= transition_distance & 
                                       other_segments$End >= transition_distance, ]
  if (nrow(matching_segment) > 0) {
    return(matching_segment$State[1])
  }
  return(NA)
}


