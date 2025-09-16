# test dumb boxplots

RealDF = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ Griesser2023.Colonial01 FemaleSong_Agg01 REAL simmap overlap output nsim 500 HackettOscine .csv")
DummyDF = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ Griesser2023.Colonial01 FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap output nsim 2000 HackettOscine .csv")


require(reshape2)
dfRealMelt <- melt(RealDF, measure.vars = c("ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present"), variable.name = "ObservedState", value.name = "ObservedState.prop")
dfRealMelt$Which <- "Real"
dfDummyMelt <- melt(DummyDF, measure.vars = c("ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present"), variable.name = "ObservedState", value.name = "ObservedState.prop")
dfDummyMelt$Which <- "Dummy"
dfDummyMelt = rbind(dfRealMelt, dfDummyMelt)

# Create the basic boxplot without x-axis labels (xaxt = "n") and without outlines
boxplot(ObservedState.prop ~ Which*ObservedState, data = dfDummyMelt, xlab = "", ylab = "Observed State Proportion", xaxt = "n", outline = FALSE, boxwex = 0.5, col = "lightgray")

# Add points to the plot
#points(jitter(as.numeric(dfDummyMelt$ObservedState), amount = 0.05), dfDummyMelt$ObservedState.prop, pch = 16, col = "darkred", cex = 0.6)

# Add rotated x-axis labels
axis(1, at = 1:length(unique(dfDummyMelt$ObservedState)),
     labels = FALSE)
text(1:length(unique(dfDummyMelt$ObservedState)),
     par("usr")[3] - 0.001, srt = 45, adj = 1.2,
     labels = as.character(unique(dfDummyMelt$ObservedState)), xpd = TRUE, cex=0.8)




### ggplot version

library(tidyverse)

trait1 = RealDF$column1[1]
trait2 = RealDF$column2[1]

dfRealMelt <- RealDF %>%
  gather("ObservedState", "ObservedState.prop", ObsProp0Absent:ObsProp1Present) %>%
  mutate(Which = "Real")

dfDummyMelt <- DummyDF %>%
  gather("ObservedState", "ObservedState.prop", ObsProp0Absent:ObsProp1Present) %>%
  mutate(Which = "Dummy")

dfCombined <- rbind(dfRealMelt, dfDummyMelt)

# Create a combined label for x-axis
dfCombined <- dfCombined %>%
  mutate(Label = paste(ObservedState, Which, sep = "\n"))

# Create the boxplot
gg_plot <-  ggplot(dfCombined, aes(x = Label, y = ObservedState.prop, fill = Which)) +
  geom_boxplot(outlier.shape = NA) + # Exclude outliers
  theme_minimal() +
  labs(y = "Observed State Proportion", x = "") +
  scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste(trait1, trait2))


require(cowplot)
Huelplots <- calcHuel(dfout = RealDF, dfDummy = DummyDF, newplot = FALSE, otherlabel = "")

plot_grid(Huelplots[[1]], Huelplots[[2]], gg_plot, ncol = 1)
