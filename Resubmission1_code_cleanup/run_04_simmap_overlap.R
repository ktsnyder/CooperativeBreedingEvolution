### Overlapping stochastic character maps to assess co-occurrence of discrete trait states in evolutionary history -

newdata = "Data_R.csv"

if (!exists(treefile)) {
  treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"  
  print("Defaulting to treefile = ConsensusPasserineTreeHackett4_1000_OscineSubset.nex. Define treefile = ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex to use the alternative consensus tree.")
}

if (!exists("count_transitions")) {
  count_transitions = FALSE 
  print("Defaulting to count_transitions = FALSE. Define count_transitions = TRUE to count transitions between binary states.")
}

if (!exists("run_all_sociality_traits")) {
  run_all_sociality_traits = FALSE
  print("Defaulting to run_all_sociality_traits = FALSE. Define run_all_sociality_traits = TRUE to perform analyses for all sociality traits and alternative cooperative breeding classification methods.")
}

if (!exists("jackknife_families_above")) {
  jackknife_families_above = 65
  print("Defaulting to jackknifing only families with at least 65 species present in our dataset. Define jackknife_families_above = 3 to perform analyses jackknifing all families we tested removal of in the publication.")
}

# Figure 3A; Supplemental Table 9 - Co-occurrance of cooperative breeding, sociality traits, and female song  ----
source(file.path("Simmap_Overlap_functions", "CharacterSimmaps_modified.R"))
source(file.path("Simmap_Overlap_functions", "simmap_overlap_runner_helpers.R"))

if (!exists("nsims_real")) {
  nsims_real = 20 
  print("Defaulting to nsims_real = 20; to perform the analysis as in the publication, define nsims_real = 500.")
}
if (!exists("nsims_dummy")) {
  nsims_dummy = 20 
  print("Defaulting to nsims_dummy = 20; to perform the analysis as in the publication, define nsims_dummy = 500.")
}
if (!exists("nsims_real_jackknife")){
  nsims_real_jackknife = 10
  print("Defaulting to nsims_real_jackknife = 10; to perform the analysis as in the publication, define nsims_real_jackknife = 50.")
}
if (!exists("nsims_dummy_jackknife")){
  nsims_dummy_jackknife = 20
  print("Defaulting to nsims_dummy_jackknife = 20; to perform the analysis as in the publication, define nsims_real_jackknife = 200.")
}

if (run_all_sociality_traits) {
  targetMetrics = c("HighConfidence_Coop", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "BiagoliniCoop", "DowningCoop", "JetzCoopInclCockburn", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop", "HighConf_Coop_DefaultToCockburnInferred", "CockburnInferred")
} else {
  targetMetrics = c("HighConfidence_Coop") 
}

temptrait2 = "FemaleSong_Agg01" 
for (i in 1: length(targetMetrics)) {
  tempMetric = targetMetrics[i]
  print(i)
  print(tempMetric)
  runSimmapOverlapAnalysis(trait1 = tempMetric, trait2 = temptrait2, nsims_real = nsims_real, nsims_dummy = nsims_dummy, tree_file = treefile, data_file = newdata, save_outputs = T, calculate_transitions = count_transitions, plot_transitions = count_transitions)
}

# Supplemental Table 7 - Co-occurrance jackknife  ----
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
templabel = ""
columns = c(trait1, trait2)

source("findQrates.R")
Qout = findQrates(columns = "HighConfidence_Coop", newdata = newdata, newtree = treefile)
qrates= Qout$qrates

subset1 <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf

familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
# For sake of example, default to only doing families with at least 65 species

familyvec = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n >= jackknife_families_above)]
familyvec = c("None", familyvec)

for (j in 1:length(familyvec)) {
  
  familyToRemove = familyvec[j]
  templabel = paste0("removed",familyToRemove)
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2_AVONET != familyToRemove),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  print(paste("Simmap Overlap", trait1, trait2, familyToRemove, j, "out of", length(familyvec), "jacks. ", numSpecies, "species in this jackknifed tree."))
  
  HuelOut <- runSimmapOverlapAnalysis(trait1 = trait1, trait2 = trait2, nsims_real = nsims_real_jackknife, nsims_dummy = nsims_dummy_jackknife, tree_file = treefile, data_file = subsetdf, save_outputs = F, calculate_transitions = FALSE, plot_transitions = FALSE, other_label = templabel, setQratesTree = treefile, setQratesData = newdata)
  
} # end cycle through families for jackknife


# Extended Data Figure 5A, Extended Data Table 1 and Supplemental Table 10 - Co-occurrance of female song (or cooperative breeding) with sociality traits; if count_transitions = TRUE, also performs analyses for Supplemental Figure 2, Extended Data Figure 5A,B, ----
if (run_all_sociality_traits) {
  socialityMetrics = c("Griesser2017FamilialLiving", "Griesser2023.MoreThanTwoCaretakers", "Griesser2023.Asocial0VsSocial1", "Territory_12vs3", "TerritorialityWeakVsStrong","Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2023.Colonial01","Final.polygyny")
} else {
  socialityMetrics = c("Griesser2017FamilialLiving", "Griesser2023.MoreThanTwoCaretakers")
}

temptrait2 = "FemaleSong_Agg01"
#temptrait2 = "HighConfidence_Coop"  # uncomment to run analyses found in Supplemental Table 10

for (i in 1: length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(tempMetric)
  runSimmapOverlapAnalysis(trait1 = tempMetric, trait2 = temptrait2, nsims_real = nsims_real, nsims_dummy = nsims_dummy, tree_file = treefile, data_file = newdata, save_outputs = T, calculate_transitions = count_transitions, plot_transitions = count_transitions)
}

# Make summary table of all analyses examining two binary traits ----
source(file.path("Simmap_Overlap_functions", "extract_simmap_overlap_results.R"))

# Extended Data Figure 5B; Supplemental Table 11 - Co-occurrence between female song and multistate traits ----
multistateTraits = c("social_system_incl_nk_coop_Griesser2017", "grouping_Griesser2023")
othertraits = c("FemaleSong_Agg01", "FemaleSong_Agg01")
nsims = nsims_real
for (i in 1:length(multistateTraits)) {
  trait1 = multistateTrait = multistateTraits[i]
  trait2 = othertrait = othertraits[i]
  columns = c(multistateTrait, othertrait)
  subsetout = subsettreedata(columns = multistateTrait, newdata = newdata, newtree = treefile)
  subsetDisctree = subsetout$subsettree
  subsetDiscdf = subsetout$subsetdf
  
  discretetraitvecDisc = subsetDiscdf[,multistateTrait]
  names(discretetraitvecDisc) = subsetDiscdf$species
  
  ARDmodel <- ace(discretetraitvecDisc,subsetDisctree, type="discrete",model = "ARD")
  
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
  
  source("findQrates.R")
  FSrates <- findQrates(columns = othertrait, newdata = newdata, newtree = treefile)
  FSQ <- FSrates$qrates
  FSQAbsPres <- FSQ[3]
  FSQPresAbs <- FSQ[2]
  
  # Make simmaps from data subsetted to those with both
  subsetSong = subsettreedata(columns = c(multistateTrait, othertrait), newdata = newdata, newtree = treefile)
  subsetdf = subsetSong$subsetdf
  subsettree = subsetSong$subsettree
  
  discretetraitvec = subsetdf[,multistateTrait]
  names(discretetraitvec) = subsetdf$species
  othertraitvec = subsetdf[,othertrait]
  names(othertraitvec) = subsetdf$species
  
  simmapMultistate = make.simmap(subsettree, discretetraitvec, nsim = nsims, Q= rate_matrix, type = "discrete") 
  realDiscreteTraitVecList = list(discretetraitvec,discretetraitvec,discretetraitvec,discretetraitvec,discretetraitvec)
  
  simmapTrait2 = make.simmap(subsettree, othertraitvec, nsim = nsims, Q= FSQ, type = "discrete")
  realTrait2vecList = list(othertraitvec,othertraitvec,othertraitvec,othertraitvec,othertraitvec)
  
  ## DUMMY
  
  CoopsimtreesRand <- list()
  RanddiscretetraitvecList <- list()
  for (j in 1:nsims) { 
    Coopvec <- subsetdf[,columns[1]]
    CoopvecRandom <- sample(Coopvec)
    names(CoopvecRandom) <- subsetdf$species
    Coopsimtree <- make.simmap(tree = subsettree, x = CoopvecRandom, model = "ARD", nsim = 1, Q = rate_matrix)
    CoopsimtreesRand[[j]] <- Coopsimtree
    RanddiscretetraitvecList[[j]] <- CoopvecRandom
    if (j == 1) {
      CoopsimtreesMulti = Coopsimtree
    } else {
      CoopsimtreesMulti = c(CoopsimtreesMulti, Coopsimtree)
    }
    print(paste(j, Sys.time()))
  } # end for j in 1:nsims (Coop)
  Coopsimtrees <- CoopsimtreesRand
  
  # Make randomized versions of FemaleSong simmaps / DUMMY data
  FSsimtreesRand <- list()
  RandTrait2vecList <- list()
  print(paste("starting Dummy FemSong simmaps", Sys.time()))
  for (j in 1:nsims) {
    FSvec <- subsetdf[,columns[2]]
    FSvecRandom <- sample(FSvec)
    names(FSvecRandom) <- subsetdf$species
    FSsimtree <- make.simmap(tree = subsettree, x = FSvecRandom, model = "ARD", nsim = 1, Q = FSQ)
    FSsimtreesRand[[j]] <- FSsimtree
    RandTrait2vecList[[j]] <- FSvecRandom
    if (j == 1) {
      FSsimtreesMulti = FSsimtree
    } else {
      FSsimtreesMulti = c(FSsimtreesMulti, FSsimtree)
    }
    print(j)
    
  } # end for j in 1:nsims (FS)
  FSsimtrees<- FSsimtreesRand
  
  overlapdf = set.seed(10)
  for (k in 1:nsims) {
    # calculate overlap - real
    realOverlap = Map.Overlap(simmapMultistate[[k]], simmapTrait2[[k]])
    overlapVec = as.vector(realOverlap)
    new_names <- outer(rownames(realOverlap), colnames(realOverlap), paste, sep = "_FS")
    new_names <- as.vector(new_names)
    names(overlapVec) = paste0(new_names, "_REAL")
    overlapVec
    
    # calculate overlap - dummy
    dummyOverlap = Map.Overlap(Coopsimtrees[[k]], FSsimtrees[[k]])
    overlapVecDummy = as.vector(dummyOverlap)
    new_names2 <- outer(rownames(dummyOverlap), colnames(dummyOverlap), paste, sep = "_FS")
    new_names2 <- as.vector(new_names2)
    names(overlapVecDummy) = paste0(new_names2, "_DUMMY")
    overlapVecDummy
    
    temprow = c(k, multistateTrait, othertrait, overlapVec, overlapVecDummy)
    
    overlapdf = rbind(overlapdf, temprow)
    overlapdf = as.data.frame(overlapdf)
    colnames(overlapdf) = c("tree", "trait1", "trait2", names(overlapVec), names(overlapVecDummy))
  }
  write.csv(overlapdf, file = file.path("Simmap_Overlap_Outputs", paste("simmap overlap", multistateTrait, othertrait, nsims, "sims.csv")), row.names = F)
  
  calcHuelout = calcHuelflex(overlapdf)
  pdfname = file.path("Simmap_Overlap_Outputs", paste("simmap overlap", multistateTrait, othertrait, nsims, "sims.pdf"))
  require(gridExtra)
  single_page_plotBox <- grid.arrange(grobs = calcHuelout[1:3], ncol = 1)
  ggsave(pdfname, single_page_plotBox, width = 8, height = 12, units = "in", limitsize = FALSE)
  
  write.csv(calcHuelout$fraction_dummy_less_than_median_real, file = file.path("Simmap_Overlap_Outputs", paste("simmap overlap FractionDummyLessThanMedianReal", multistateTrait, othertrait, nsims,"sims.csv")))
}


# Figures 3C, 3D; Extended Data Figures 4A, 4B - Transition rates split by territoriality ----
source(file.path("Simmap_Overlap_functions", "TransitionCounts_3trait_flexTerr_fxns.R"))

if (count_transitions) {
  
  print(paste("Beginning transition-counts analysis for Cooperative breeding, Female song, and strength of territoriality for", nsims_real, "simulations. If nsims_real is large, this will take up to a few hours. Change nsims_real at the top of run_04_simmap_overlap.R to run for more or fewer simulations."))
  plot_transition_counts_3trait(Qdata = newdata, Qtree = treefile, columns = c("HighConfidence_Coop","FemaleSong_Agg01", "TerritorialityWeakVsStrong"), nsims = nsims_real)
  
  print(paste("Beginning transition-counts analysis for Cooperative breeding, Female song, and year-round territoriality for", nsims_real, "simulations. If nsims_real is large, this will take up to a few hours. Change nsims_real at the top of run_04_simmap_overlap.R to run for more or fewer simulations."))
  plot_transition_counts_3trait(Qdata = newdata, Qtree = treefile, columns = c("HighConfidence_Coop","FemaleSong_Agg01", "Territory_12vs3"), nsims = nsims_real)
  
}

