### Overlapping stochastic character maps to assess co-occurrence of discrete trait states in evolutionary history ----

newdata = "Data_R.csv"
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

# For the alternative tree, uncomment this line:
# treefile = "ConsensusPasserineTreeHackett4_1000_mean-edge_ignore-absent.nex

# Figure 3A; Supplemental Table 9 - Co-occurrance of cooperative breeding and female song  ----
source("test_trait_overlap_simmaps.R")
nsims_real = 10
nsims_dummy = 10

targetMetrics = c("HighConfidence_Coop", "Griesser2017FamilialLiving", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "BiagoliniCoop", "DowningCoop", "JetzCoopInclCockburn", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop", "HighConf_Coop_DefaultToCockburnInferred", "CockburnInferred")
temptrait2 = "FemaleSong_Agg01" 
for (i in 1: length(targetMetrics)) {
  tempMetric = targetMetrics[i]
  print(i)
  print(tempMetric)
  dfout <- CharacterSimmaps(columns = c(tempMetric,temptrait2), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = TRUE)
  dfDummy <- CharacterSimmaps(columns = c(tempMetric,temptrait2), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
  
  calcHuelout = calcHuel(dfout, dfDummy)
  require(gridExtra)
  plotname = file.path("Simmap Overlap Outputs", paste(tempMetric, temptrait2, nsims_real, nsims_dummy, treelabel, "withTransCounts.pdf"))
  
  calcHuelout2 = calcHuelout[c("p1", "p2", "p3", "p4","p5","p6")]
  nPlots = 6
  m3 <- marrangeGrob(calcHuelout2, ncol = 1, nrow = nPlots)
  ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
}

# Supplemental Table 7 - Co-occurrance jackknife  ----
trait1 = "HighConfidence_Coop"
trait2 = "FemaleSong_Agg01"
templabel = ""
columns = c(trait1, trait2)

nsims_real = 5
nsims_dummy = 10

source("findQrates.R")
Qout = findQrates(columns = "HighConfidence_Coop", newdata = newdata, newtree = treefile)
qrates= Qout$qrates

subset1 <- subsettreedata(columns = columns, newdata = newdata, newtree = treefile, skinnydata = FALSE)
subsettree1 <- subset1$subsettree
subsetdf1 <- subset1$subsetdf

familycounts = subsetdf1 %>% group_by(Family3_BirdtreeMatchSpecies2_AVONET) %>% count
familyvec = familycounts$Family3_BirdtreeMatchSpecies2_AVONET[which(familycounts$n > 65)]
familyvec = c("None", familyvec)

plotlist = list()
for (j in 1:length(familyvec)) {
  
  familyToRemove = familyvec[j]
  templabel = paste0(columns[1], " ", columns[2], " ", "removed",familyToRemove)
  tempdfIn = subsetdf1[which(subsetdf1$Family3_BirdtreeMatchSpecies2_AVONET != familyToRemove),]
  
  subsetbtw <- subsettreedata(columns = columns, newdata = tempdfIn, newtree = subsettree1, skinnydata = TRUE)
  subsettree <- subsetbtw$subsettree
  subsetdf <- subsetbtw$subsetdf
  subsetdf[,columns[1]] <- as.character(subsetdf[,columns[1]])
  subsetdf[,columns[2]] <- as.character(subsetdf[,columns[2]])
  numSpecies = length(subsetdf$species)
  
  print(paste("Simmap Overlap", trait1, trait2, familyToRemove, j, "out of", length(familyvec), "jacks. ", numSpecies, "species in this jackknifed tree."))
  
  dfout4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap", columnForGlobalQ = 1, columnGlobalQrates = qrates)
  dfDummy4 <- CharacterSimmaps(columns = columns, df = subsetdf, tree =  subsettree, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = templabel, dummyMethod = "makeSimmap", columnForGlobalQ = 1, columnGlobalQrates = qrates)
  
  HuelOut = calcHuel(dfout = dfout4, dfDummy = dfDummy4, newplot = FALSE, otherlabel = templabel)
  plotlist[[j]] <- HuelOut
} # end cycle through families for jackknife

# Plot jackknife output
require(gridExtra)
index = 0
plotlist2 = list()
for (i in 1:length(plotlist)) {
  tempplots = plotlist[[i]]
  for (k in 1:3) {
    index = index+1
    plotlist2[[index]] = tempplots[[k]]
  }
}
m1 <- marrangeGrob(plotlist2, ncol = 1, nrow = 3)
ggsave(paste(Sys.Date(), "jackknifed Simmap Overlaps Coop FemaleSong_Agg01.pdf"), m1, width = 8, height = 9, units = "in")


# Extended Data Figure 5A, Extended Data Table 1 and Supplemental Table 10 - Co-occurrance of female song (or cooperative breeding) with sociality traits ----
source("test_trait_overlap_simmaps.R")
nsims_real = 20
nsims_dummy = 20

socialityMetrics = c("Griesser2023.Asocial0VsSocial1", "Griesser2023.Colonial01" , "Griesser2023.MoreThanTwoCaretakers" ,"Griesser2023.LongSocialBonds","Griesser2023.GroupsLargerThanPair", "Griesser2023.TwoOrMoreCaretakers", "Griesser2023.LargestGroupSizes", "Griesser2023.SeasonOrLongerSocialBonds", "Griesser2017FamilialLiving","Final.polygyny" ,"MeanCoopTie2Noncoop", "MeanCoopTie2Coop", "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "Territory_12vs3", "TerritorialityWeakVsStrong")

temptrait2 = "FemaleSong_Agg01"
#temptrait2 = "HighConfidence_Coop"  # uncomment to run analyses found in Supplemental Table 10

for (i in 1: length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(i)
  print(tempMetric)
  dfout <- CharacterSimmaps(columns = c(tempMetric,temptrait2), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = NULL, plotSampleSimmaps = TRUE)
  dfDummy <- CharacterSimmaps(columns = c(tempMetric,temptrait2), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = NULL, dummyMethod = "makeSimmap")
  
  calcHuelout = calcHuel(dfout, dfDummy)
  require(gridExtra)
  plotname = file.path("Simmap Overlap Outputs", paste(tempMetric, temptrait2, nsims_real, nsims_dummy, treelabel, "withTransCounts.pdf"))
  
  calcHuelout2 = calcHuelout[c("p1", "p2", "p3", "p4","p5","p6")]
  nPlots = 6
  m3 <- marrangeGrob(calcHuelout2, ncol = 1, nrow = nPlots)
  ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
}


# Extended Data Figure 5B; Supplemental Table 11 - Co-occurrence between female song and multistate traits ----
multistateTraits = c("social_system_incl_nk_coop_Griesser2017", "grouping_Griesser2023")
othertraits = c("FemaleSong_Agg01", "FemaleSong_Agg01")
nsims = nsim
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
  write.csv(overlapdf, file = paste("simmap overlap", multistateTrait, othertrait, nsims, "sims.csv"), row.names = F)
  
  calcHuelout = calcHuelflex(overlapdf)
  pdfname = paste("simmap overlap", multistateTrait, othertrait, nsims, "sims.pdf")
  require(gridExtra)
  single_page_plotBox <- grid.arrange(grobs = calcHuelout[1:3], ncol = 1)
  ggsave(pdfname, single_page_plotBox, width = 8, height = 12, units = "in", limitsize = FALSE)
  
  write.csv(calcHuelout$fraction_dummy_less_than_median_real, file = paste("simmap overlap FractionDummyLessThanMedianReal", multistateTrait, othertrait, nsims,"sims.csv"))
}