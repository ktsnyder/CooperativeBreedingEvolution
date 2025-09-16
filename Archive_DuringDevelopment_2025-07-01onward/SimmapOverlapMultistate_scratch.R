### SimmapOverlap multistate trait Excerpted from Run_Analyses.R
newdata = "Data_R_2025-06-09.csv"
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

source("Archive_DuringDevelopment_2025-07-01onward/test_trait_overlap_simmaps.R")

df = read.csv(newdata)
df$TerritorialityWeakVsStrong[which(df$TerritorialityWeakVsStrong == 0)] <- "Weak"
df$TerritorialityWeakVsStrong[which(df$TerritorialityWeakVsStrong == 1)] <- "Strong"
df$HighConfidence_Coop[which(df$HighConfidence_Coop == 0)] <- "Noncooperative"
df$HighConfidence_Coop[which(df$HighConfidence_Coop == 1)] <- "Cooperative"
df$TerrWeakStrongXHighConfCoop <- paste(df$HighConfidence_Coop, df$TerritorialityWeakVsStrong, sep = "_")
require(stringr)
df$TerrWeakStrongXHighConfCoop[which(str_detect(df$TerrWeakStrongXHighConfCoop, "NA"))] <- NA
newdata = df

multistateTraits = c("social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "grouping_Griesser2023", "Territory")
multistateTraits = c("Territory", "Social.bond", "Social.bond")
othertraits = c("FemaleSong_Agg01", "HighConfidence_Coop", "FemaleSong_Agg01")
othertraits = "FemaleSong_Agg01"
multistateTraits = "TerrWeakStrongXHighConfCoop"
nsims = 500
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
  write.csv(overlapdf, file = paste("simmap overlap MultistateCustomQ", multistateTrait, othertrait, nsims, "sims.csv"), row.names = F)
  
  calcHuelout = calcHuelflex(overlapdf)
  pdfname = paste("simmap overlap MultistateCustomQ", multistateTrait, othertrait, nsims, "sims.pdf")
  require(gridExtra)
  single_page_plotBox <- grid.arrange(grobs = calcHuelout[1:3], ncol = 1)
  ggsave(pdfname, single_page_plotBox, width = 8, height = 18, units = "in", limitsize = FALSE)
  
  write.csv(calcHuelout$fraction_dummy_less_than_median_real, file = paste("simmap overlap MultistateCustomQ FractionDummyLessThanMedianReal", multistateTrait, othertrait, nsims,"sims.csv"))
}
