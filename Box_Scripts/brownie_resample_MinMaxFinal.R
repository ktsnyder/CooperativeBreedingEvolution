###### browniefunction wrapper for randomizing min/max/median values

boxpath = "/Users/kate/Library/CloudStorage/Box-Box"

setwd(file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/'))

source("browniefunction.R")
source("findQrates.R")
source("plotbrownie.R")

newdata = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/Data_R_Passerine_withTobias.csv')
treefile = file.path(boxpath, 'Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/ConsensusPasserineTreeHackett4_1000_OscineSubset.nex')

tree = read.nexus(treefile)
dfIn = read.csv(newdata)
rownames(dfIn) = dfIn$species

discreteCatLabels = c("Non-cooperative", "Cooperative")
nsim = 2 # due to the way each simmap is called in browniefunction() (as one element in a list of simmaps), we can't do just 1 sim
currentlabel = "ResampleMinMedMax"
nSeeds = 600

# get Q rates out here so we can use the pared-down dataframe in the loop
Qoutput <- findQrates(columns = "HighConfidence_Coop", plot=F, newtree = tree, newdata = dfIn)
qrates <- Qoutput$qrates


#### Run song features ####

for (song_base in c("Song.rep.", "Syllable.rep.")) {
#song_base = "Song.rep."

songcols <- grep(paste0("^",song_base), colnames(dfIn), value = TRUE)
songcols

# subset dataframe
complete_vars <- c("HighConfidence_Coop", songcols)
df <- dfIn[complete.cases(dfIn[,complete_vars]),]
df = df[,c("species", complete_vars)]


dfloop <- df
outdf = set.seed(10)
for (i in 1:nSeeds) {
  set.seed(i)
  sampleColName = paste0(song_base,"Sample",i)
  dfloop[[sampleColName]] <- apply(dfloop[, songcols], 1, sample, size = 1)
  brownieout = browniefunction(columns = c("HighConfidence_Coop", sampleColName), newdata = dfloop, newtree = treefile, nsim = nsim, islog = sampleColName, plotsimmaps = FALSE, setQrates = qrates)
  outdf = rbind(outdf, brownieout)
}

write.csv(outdf, paste0(Sys.Date(),"HighConfidence_Coop_",song_base, currentlabel, "_brownie",nsim,"sim", nSeeds, "resamples.csv"), row.names = F)
write.csv(dfloop, paste0(Sys.Date()," HighConfidence_Coop_",song_base, currentlabel, " ", nSeeds, "resampled columns.csv"), row.names = F)

plotbrownie(data = outdf, columns = c("HighConfidence_Coop",song_base), discreteCategoryLabels = discreteCatLabels, otherlabel = currentlabel, newpdf = TRUE, nsim = nSeeds, islog = TRUE)

}


source("brownie relative rates.R")
# compile results into one table
BrownieRelativeRates(BrownieOutputFolder = ".", otherlabel = currentlabel)
BrownieRelativeRates(BrownieOutputFolder = "OutputFiles", otherlabel = currentlabel)



