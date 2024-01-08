## test simmap overlap scripts

# make a small phylogeny with dummy data
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
tree = read.nexus(treefile)

length(tree$tip.label)
droptips = sample(1:4685, 4665, replace = FALSE)
subtree = drop.tip(tree, tip = droptips)
plot(subtree)
species = subtree$tip.label
dummytrait1vec = c(rep(0, 10), rep(1, 10))
names(dummytrait1vec) = species
dummy1simmap = make.simmap(subtree, dummytrait1vec, model = "SYM", nsim = 6)
pdf("test subset simmap dummy1.pdf")
plotSimmap(dummy1simmap)
dev.off()
write.simmap(dummy1simmap, file = "test subset simmap dummy1.nex", format = "nexus")
dummytrait2vec = c(rep(0, 12), rep(1, 8))
names(dummytrait2vec) = species
dummy2simmap = make.simmap(subtree, dummytrait2vec, model = "SYM")
write.simmap(dummy2simmap, file = "test subset simmap dummy2.nex", format = "nexus")
pdf("test subset simmap dummy2.pdf")
plotSimmap(dummy2simmap)
dev.off()

# Manual Count:
# Dummy1 = Coop
# Dummy2 = FS
# N Dummy1 black to red in Dummy2 black: 0 
# N Dummy1 black to red in Dummy2 red: 1
# N Dummy1 red to black in Dummy2 black: 2
# N Dummy1 red to black in Dummy2 red: 0
# N Dummy2 black to red in Dummy1 black: 3
# N Dummy2 black to red in Dummy1 red: 1
# N Dummy2 red to black in Dummy1 black: 2
# N Dummy2 red to black in Dummy1 red: 1

Coopsimtree1 = dummy1simmap
FSsimtree1 = dummy2simmap

# getTransitionStateCounts() seems ok
# CharacterSimmaps, calcHuel seem correct


# Need to try doing full thing 
treefile = subtree
newdata = cbind(species, dummytrait1vec, dummytrait2vec)
newdata = as.data.frame(newdata)
colnames(newdata) = c("species","dummy1","dummy2")
dfout4 <- CharacterSimmaps(columns = c("dummy1", "dummy2"), df = newdata, tree =  treefile, dummy = FALSE, nsims = 3, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = TRUE)
dfDummy <- CharacterSimmaps(columns = c("dummy1", "dummy2"), df = newdata, tree =  treefile, dummy = TRUE, nsims = 3, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = TRUE, dummyMethod = "makeSimmap")
calcHuelout = calcHuel(dfout4, dfDummy4)
nPlots = length(calcHuelout)
plotname = paste0("Simmap Overlap Outputs/","TEST dummy1 dummy2", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts.pdf")
m3 <- marrangeGrob(calcHuelout[1:6], ncol = 1, nrow = 6)
ggsave(plotname, m3, width = 7.5, height = 3.5*nPlots, units = "in")


# testing end part of "find transition counts by state for two Discrete traits.R" getTransitionStateCounts()
countsdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/synthetic countsdf for test overlap_counts output nsim 3 HackettOscine .csv")
# [ ran through end part as described above ]
#write.csv(countsdf, "/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/synthetic countsdf for test overlap_counts output from getTransitionStateCounts nsim 3 HackettOscine .csv", row.names = FALSE)
countsdf = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/synthetic countsdf for test overlap_counts output from getTransitionStateCounts nsim 3 HackettOscine .csv")
transStateCounts = countsdf
# [ using this as intermediate output in test_trait_overlap_simmaps.R CharacterSimmaps(), along with modified dfout]
dfout = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ dummy1 dummy2 REAL simmap overlap_counts output nsim 3 HackettOscine .csv")
colnames(transStateCounts) %in% colnames(dfout)
dfout = dfout[,which(!colnames(dfout) %in% colnames(transStateCounts))]
# then start with line 223 in test_trait_overlap_simmaps.R: dfout = merge(dfout, transStateCounts, by.x = "treenum", by.y = "TreeNum")
write.csv(dfout, "/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/synthetic dfout for test overlap_counts output from CharacterSimmaps nsim 3 HackettOscine .csv")

# use dfDummy from above (though it only matters just so calcHuel runs and plots 3 and 4, but not concerned about those)
dfDummy = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ dummy1 dummy2 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 3 HackettOscine .csv")

# continue with Generate_Pub_Figs.R
#   tempdfDep = read.csv(DepFile)
#  tempdfInd = read.csv(IndFile)
tempdfDep = dfout
tempdfInd = dfDummy
calcHuelout = calcHuel(tempdfDep, tempdfInd)
Transitions = calcHuelout$TransitionStats$logPairwisePostHoc[2]
pvals = calcHuelout$TransitionStats$logPairwisePostHoc[7]
calcHuelout$p5 # looking good

# OUTPUT OF transition_plot() DOES NOT MATCH BOXPLOTS. Both transition_df and the plot
# testing more granularly... start slightly before transition_plot() in Generate_Pub_Figs.R
tempdfDep[,paste0(colnames(tempdfDep)[18:25],"DifferenceFromExpected")] = tempdfDep[,colnames(tempdfDep)[18:25]] - tempdfDep[,paste0(colnames(tempdfDep)[18:25],"Expected")]
#write.csv(tempdfDep, "/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/synthetic tempdfDep for test overlap_counts with DifferenceFromExpected columns nsim 3 HackettOscine .csv")


# Are these the incorrect lines?? From GeneratePubFIgs
RateRef = as.data.frame(rbind(c("FS0to1inCoop0", "q12"),c("Coop0to1inFS0", "q13"),c("FS1to0inCoop0", "q21"), c("Coop0to1inFS1", "q24"),c("Coop1to0inFS0", "q31"), c("FS0to1inCoop1", "q34"), c("Coop1to0inFS1", "q42"),c("FS1to0inCoop1", "q43")))
colnames(RateRef) <- c("Transitions", "qRate")
ratePvals = merge(RateRef, pvaldf, by.x = "Transitions", by.y = "Transition")

tempdfDep[,paste0(colnames(tempdfDep)[18:25],"DifferenceFromExpected")] = tempdfDep[,colnames(tempdfDep)[18:25]] - tempdfDep[,paste0(colnames(tempdfDep)[18:25],"Expected")] # this seems ok, or at least consistent with boxplots
colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop0DifferenceFromExpected")] <- "q12"
colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS0DifferenceFromExpected")] <- "q13"
colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop0DifferenceFromExpected")] <- "q21"
colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop0to1inFS1DifferenceFromExpected")] <- "q24"
colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS0DifferenceFromExpected")] <- "q31"
colnames(tempdfDep)[which(colnames(tempdfDep) == "FS0to1inCoop1DifferenceFromExpected")] <- "q34"
colnames(tempdfDep)[which(colnames(tempdfDep) == "Coop1to0inFS1DifferenceFromExpected")] <- "q42"
colnames(tempdfDep)[which(colnames(tempdfDep) == "FS1to0inCoop1DifferenceFromExpected")] <- "q43"


tempdfDep[,42:49] = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8)
transition_plot(df = tempdfDep, trait1StateLabels = trait1StateLabels, trait2StateLabels = trait2StateLabels, scale_area_by = 0.2, offset = 0.15, lengthen = 0.2, ratePvals = ratePvals)
# q12 color should be where q13 is (color level 1)
# q13 color should be where q31 is (color level 2)
# q21 color should be where q24 is (color level 3)
# q24 color should be where q42 is  (color level 4)
# q31 color should be where q12 is (color level 5)
# q34 color should be where q21 is (color level 6)
# q42 color should be where q34 is (color level 7)
# q43 color is where it should be


# THE ISSUE IS the rates are in a different order in df (tempdfDep), and when I just do cbind in transition plot below, it messes it up
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

#from_state = c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))
#to_state = c(2, 3, 1, 4, 1, 4, 2, 3)
#qColumns = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
#transition_df = cbind(qColumns,from_state, to_state, transition_means,transition_df)

## editing above commented out part
transition_means = cbind(names(transition_means), transition_means)
transition_means = as.data.frame(transition_means)
colnames(transition_means) = c("qRates", "transition_means")
transition_means
from_state = c(rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))
to_state = c(2, 3, 1, 4, 1, 4, 2, 3)
qColumns = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")
transition_df = cbind(qColumns,from_state, to_state, transition_df)
transition_df
transition_df_corrected = merge(transition_df, transition_means, by.x = "qColumns", by.y = "qRates")

# need to correct transition_plot.R with this^
