# from chatGPT4 11/29/2023
# Kate Snyder
# Last edited: 12/04/2023 - fixed arrow identifiers/dataframe combination per "test simmap overlap weirdness.R", fixed right side arrow adjustments
# May still need to troubleshoot scale_area_by

library(ggplot2)

# Read the data
#df <- read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/BayesTraitsDiscreteML_MeanCoopTie2Coop-FemaleSong_Agg01 Jackknife Outputs/BayesTraitsDiscreteML_MeanCoopTie2Coop-FemaleSong_Agg01_removeAlaudidae.csv")



transition_plot <- function(df, trait1StateLabels = c("0","1"), trait2StateLabels = c("0","1"), scale_area_by = 0.5, offset = 0.15, lengthen = 0.4, ratePvals = NULL, plottitle = NULL) {
# Calculate the mean of each q__ column
transition_means <- colMeans(df[, grepl("^q[0-9]{2}$", names(df))])

# set beginning and end points for each arrow
q12 = c(2,4,3,4)
q13 = c(1,3,1,2)
q21 = c(3,4,2,4)
q24 = c(4,3,4,2)
q31 = c(1,2,1,3)
q34 = c(2,1,3,1)
q42 = c(4,2,4,3)
q43 = c(3,1,2,1)
transition_df = as.data.frame(rbind(q12, q13, q21, q24, q31, q34, q42, q43))
colnames(transition_df) = c("from_x", "from_y", "to_x", "to_y")

## these lines were causing rates to be assigned to incorrect arrows
#from_state = c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))
#to_state = c(2, 3, 1, 4, 1, 4, 2, 3)
#qColumns = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
#transition_df = cbind(qColumns,from_state, to_state, transition_means,transition_df)

# these 9 lines replace above 4 to correct arrow labeling
transition_means = cbind(names(transition_means), transition_means)
transition_means = as.data.frame(transition_means)
colnames(transition_means) = c("qRates", "transition_means")
transition_means
from_state = c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))
to_state = c(2, 3, 1, 4, 1, 4, 2, 3)
qColumns = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
transition_df_noMeans = cbind(qColumns,from_state, to_state, transition_df)
transition_df_noMeans
transition_df = merge(transition_df_noMeans, transition_means, by.x = "qColumns", by.y = "qRates")


# set state labels and positions of labels
label00 = paste(trait1StateLabels[1], trait2StateLabels[1], sep = "\n")
label01 = paste(trait1StateLabels[1], trait2StateLabels[2], sep = "\n")
label10 = paste(trait1StateLabels[2], trait2StateLabels[1], sep = "\n")
label11 = paste(trait1StateLabels[2], trait2StateLabels[2], sep = "\n")

point00 = c(1,4)
point01 = c(4,4)
point10 = c(1,1)
point11 = c(4,1)

statePoints = as.data.frame(rbind(point00, point01, point10, point11))
colnames(statePoints) <- c("x", "y")
stateTextDF = cbind(stateNum = c(1,2,3,4), stateLabel = c(label00, label01, label10, label11), statePoints)

transitions = merge(transition_df, stateTextDF, by.x = "from_state", by.y = "stateNum")

# Adjust the coordinates for horizontal and vertical orientation with offsets for parallel arrows
horizontalRates = c("q12", "q21", "q43", "q34")
verticalRates = c("q13", "q31", "q24", "q42")

# Apply the adjustment
transitions2list <- adjust_coordinates(transitions, horizontalRates, verticalRates, scale_area_by = scale_area_by, offset = offset, lengthen = lengthen)
transitions2 = transitions2list$dfTransition
transitions2$transition_means = as.numeric(transitions2$transition_means)

if ("p.value" %in% colnames(ratePvals)) {
  transitions2 = merge(ratePvals, transitions2, by.x = "qRate", by.y = "qColumns")
  transitions2$SignificanceLabel = paste(transitions2$SignificanceLabel) #, transitions2$qRate, sep="_")
}


# Define the color gradient
my_color_gradient <- scale_color_gradient2(low = "#2166ac", mid = "#f0f0f0", high = "#b2182b", midpoint = 0) #mid = "white", high = "#b2182b")

# Plot with corrected arrow orientations and offsets
transitionplot <- ggplot(data = transitions2) +
  geom_segment(aes(x = from_x, y = from_y, xend = to_x, yend = to_y, color = transition_means), 
               arrow = arrow(type = "closed", length = unit(0.01, "inches")), size = 7, linejoin = "mitre") +
  my_color_gradient +
  labs(y = "Transition Counts", x = "", color = "Mean Difference \nFrom Expected \nNumber of \nTransitions") +
  geom_text(aes(label = stateLabel, x = x, y = y), size = 4) +
  geom_text(aes(label = SignificanceLabel, x = to_x, y = to_y), size = 3) +
  theme_minimal() +
  xlim(c(0,5)) + 
  ylim(c(0.5,4.5)) + 
  ggtitle(label = plottitle) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )
transitionplot

output = list(transition_df = transitions2, transition_plot = transitionplot)

return(output)

}


# Function to adjust coordinates for parallel arrows
adjust_coordinates <- function(dfTransition, horizontalRates, verticalRates, offset = 0.1, lengthen = 0.5, scale_area_by = 1) {
  
  # adjust lengthen to make the arrow tips not so over-hangy
  lengthenTo = lengthen*0.5
  lengthenFrom = lengthen*1.3
  
  # Offset for parallel arrows
  
  # Offset the horizontal rates that need to be parallel
  for (rate in horizontalRates) {
    print(rate)
    if (rate %in% c("q12", "q34")) {
      dfTransition$from_y[dfTransition$qColumns == rate] <- dfTransition$from_y[dfTransition$qColumns == rate] - offset
      dfTransition$to_y[dfTransition$qColumns == rate] <- dfTransition$to_y[dfTransition$qColumns == rate] - offset
      dfTransition$from_x[dfTransition$qColumns == rate] <- dfTransition$from_x[dfTransition$qColumns == rate] - lengthenFrom
      dfTransition$to_x[dfTransition$qColumns == rate] <- dfTransition$to_x[dfTransition$qColumns == rate] + lengthenTo
    } else if (rate %in% c("q43", "q21")) {
      dfTransition$from_y[dfTransition$qColumns == rate] <- dfTransition$from_y[dfTransition$qColumns == rate] + offset
      dfTransition$to_y[dfTransition$qColumns == rate] <- dfTransition$to_y[dfTransition$qColumns == rate] + offset
      dfTransition$from_x[dfTransition$qColumns == rate] <- dfTransition$from_x[dfTransition$qColumns == rate] + lengthenFrom
      dfTransition$to_x[dfTransition$qColumns == rate] <- dfTransition$to_x[dfTransition$qColumns == rate] - lengthenTo
    }
  }
  # Offset the vertical rates that need to be parallel
  for (rate in verticalRates) {
    print(rate)
    if (rate %in% c("q13", "q24")) {
      dfTransition$from_x[dfTransition$qColumns == rate] <- dfTransition$from_x[dfTransition$qColumns == rate] - offset
      dfTransition$to_x[dfTransition$qColumns == rate] <- dfTransition$to_x[dfTransition$qColumns == rate] - offset
      dfTransition$from_y[dfTransition$qColumns == rate] <- dfTransition$from_y[dfTransition$qColumns == rate] + lengthenFrom 
      dfTransition$to_y[dfTransition$qColumns == rate] <- dfTransition$to_y[dfTransition$qColumns == rate] - lengthenTo
    } else if (rate %in% c("q42", "q31")) {
      dfTransition$from_x[dfTransition$qColumns == rate] <- dfTransition$from_x[dfTransition$qColumns == rate] + offset
      dfTransition$to_x[dfTransition$qColumns == rate] <- dfTransition$to_x[dfTransition$qColumns == rate] + offset
      dfTransition$from_y[dfTransition$qColumns == rate] <- dfTransition$from_y[dfTransition$qColumns == rate] - lengthenFrom 
      dfTransition$to_y[dfTransition$qColumns == rate] <- dfTransition$to_y[dfTransition$qColumns == rate] + lengthenTo
    }
    # if (rate %in% c("q13", "q31")) { # move vertical arrows slightly inward
    #   dfTransition$from_x[dfTransition$qColumns == rate] <- dfTransition$from_x[dfTransition$qColumns == rate] + offset
    #   dfTransition$to_x[dfTransition$qColumns == rate] <- dfTransition$to_x[dfTransition$qColumns == rate] + offset
    # } else if (rate %in% c("q24", "q42")) {
    #   dfTransition$from_x[dfTransition$qColumns == rate] <- dfTransition$from_x[dfTransition$qColumns == rate] - offset
    #   dfTransition$to_x[dfTransition$qColumns == rate] <- dfTransition$to_x[dfTransition$qColumns == rate] - offset
    # }
  } # end for rate in VerticalRates
  
  # Move state labels inwards a little
  #dfTransition$x[which(dfTransition$x == 1)] <- 1.2
  #dfTransition$x[which(dfTransition$x == 4)] <- 3.8
  
  dfTransitionUnscaled = dfTransition
  dfTransition[,c("from_x", "from_y", "to_x", "to_y", "x", "y")] <- dfTransition[,c("from_x", "from_y", "to_x", "to_y", "x", "y")]*scale_area_by
  
  outlist = list(dfTransition, dfTransitionUnscaled)
  names(outlist) = c("dfTransition", "dfTransitionUnscaled")
  
  return(outlist)
}


