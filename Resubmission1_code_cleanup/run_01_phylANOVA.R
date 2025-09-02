# Run for phylANOVA results

source("subsettreedata.R")

newdata = "Data_R.csv"
df = read.csv(newdata)
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex" # to perform test using the alternative consensus tree for any given analysis, replace this with "ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex"
# treefile can also be set to any individual tree extracted from a BirdTree multiphylo object in order to test across many individual trees
songtraits = c("Song.rep.final","Syllable.rep.final", "Syll.song.final", "Duration.final", "Interval.final")
socialityMetrics = c("HighConfidence_Coop", "Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "grouping_Griesser2023", "social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "social_bonds_Griesser2023", "Territory_12vs3", "TerritorialityWeakVsStrong")

# Supplemental Table 3 - phylANOVA ----
phynovaDF = set.seed(10)
for (j in 1:length(songtraits)) {
  for (k in 1:length(socialityMetrics)) {
    songtrait = songtraits[j]
    tempgrouptrait = socialityMetrics[k]
    subsets = subsettreedata(columns = c(tempgrouptrait, songtrait), newdata = newdata, newtree = treefile)
    subsetdf = subsets$subsetdf
    discvec = subsetdf[,tempgrouptrait]
    names(discvec) = subsetdf$species
    contvec = subsetdf[,songtrait]
    #contvec = log(contvec)
    names(contvec) = subsetdf$species
    subsettree = subsets$subsettree
    Nspecies = length(subsettree$tip.label)
    Ngroups = length(unique(discvec))
    
    tryCatch({
      phylANOVAout = phylANOVA(subsettree, x = discvec, y = contvec, nsim = 50000, posthoc = TRUE)
      phylANOVAp = phylANOVAout$Pf
      temprow = c(tempgrouptrait, Ngroups, songtrait, Nspecies, phylANOVAp)
    }, error = function(e) {
      temprow = c(tempgrouptrait, Ngroups, songtrait, Nspecies, NA)
      message("Error in phylANOVA computation: ", e$message)
    })
    
    phynovaDF = rbind(phynovaDF, temprow)
    phynovaDF = as.data.frame(phynovaDF)
    colnames(phynovaDF) <- c("DiscreteTrait", "DiscreteNumGroups", "ContinuousTrait", "n_Species", "PhylANOVApval")
    print(paste(tempgrouptrait, songtrait))
    print(phylANOVAout)
  }
}
write.csv(phynovaDF, file = paste(Sys.Date(), "phylANOVA outputs Songs.csv"), row.names = F)

