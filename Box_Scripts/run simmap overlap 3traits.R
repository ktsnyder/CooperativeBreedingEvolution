# CharacterSimmaps flex development for using multiple multi-state OR more than two states?
# 
# Last edited: 3/4/2025 - RandTrait3vecList[[j]] from RandTrait2vecList[[j]] where noted
# 
setwd('/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code') 
newdata = "Data_R_Passerine_withTobias.csv"
newdata = "scratchData_TerrWeakStrongYearRound_wJiayingDuet2025-03-04.csv"
newdata="Data_R_Passerine_withTobias_AVONET_JiayingDuet2025-03-04_WeakStrong2025-04-29.csv"
newdata = "Data_R_2025-06-09.csv"
df = read.csv(newdata)
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

source("subsettreedata.R")
source("/Users/kate/Library/CloudStorage/Box-Box/Kate_Nicole/CooperativeBreedingEvolutionOutputs/Nature Eco Evo - resubmission Code/test_trait_overlap_simmaps.R")
source("find transition counts by state for 2 Discrete traits.R")
source("transition_plot.R")
source("MapOverlapThree.R")

multistateTraits = c("social_system_incl_nk_coop_Griesser2017", "social_system_Griesser2017", "grouping_Griesser2023")
multistateTraits = c("Territory", "Social.bond") #, "Social.bond")
multistateTraits = "WeakStrongTerr.Rough"
multistateTraits = c("WeakStrongTerrYearRound.Rough", "Terr4Level.Rough")
multistateTraits = c("TerritorialityWeakVsStrongHighConf", "TerritorialityWeakVsStrong")
multistateTraits = c("Territory_12vs3")
othertraits = c("FemaleSong_Agg01", "FemaleSong_Agg01")
nsims = 500
for (i in 1:length(multistateTraits)) {
  trait1 = multistateTrait = multistateTraits[i]
  trait2 = othertrait = othertraits[i]
  trait3 = "HighConfidence_Coop"
  columns = c(multistateTrait, othertrait, trait3)
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
  
  source("findQrates.R")
  FSrates <- findQrates(columns = othertrait, newdata = newdata, newtree = treefile)
  FSQ <- FSrates$qrates
  FSQAbsPres <- FSQ[3]
  FSQPresAbs <- FSQ[2]
  
  CBrates <- findQrates(columns = trait3, newdata = newdata, newtree = treefile)
  CBQ <- CBrates$qrates
  CBQAbsPres <- CBQ[3]
  CBQPresAbs <- CBQ[2]
  
  # Make simmaps from data subsetted to those with all 3
  subsetSong = subsettreedata(columns = c(multistateTrait, othertrait, trait3), newdata = newdata, newtree = treefile)
  subsetdf = subsetSong$subsetdf
  subsettree = subsetSong$subsettree
  
  discretetraitvec = subsetdf[,multistateTrait]
  names(discretetraitvec) = subsetdf$species
  othertraitvec = subsetdf[,othertrait]
  names(othertraitvec) = subsetdf$species
  thirdtraitvec = subsetdf[,trait3]
  names(thirdtraitvec) = subsetdf$species
  
  simmapMultistate = make.simmap(subsettree, discretetraitvec, nsim = nsims, Q= rate_matrix, type = "discrete") 
  realDiscreteTraitVecList = list(discretetraitvec,discretetraitvec,discretetraitvec,discretetraitvec,discretetraitvec)
  
  simmapTrait2 = make.simmap(subsettree, othertraitvec, nsim = nsims, Q= FSQ, type = "discrete")
  realTrait2vecList = list(othertraitvec,othertraitvec,othertraitvec,othertraitvec,othertraitvec)
  
  simmapTrait3 = make.simmap(subsettree, thirdtraitvec, nsim = nsims, Q= CBQ, type = "discrete")
  realTrait3vecList = list(thirdtraitvec,thirdtraitvec,thirdtraitvec,thirdtraitvec,thirdtraitvec)
  
  pdf(file = file.path("Simmap Overlap Outputs", paste(trait1, trait2, trait3, "egSimmaps.pdf")),height=28,width=6)
  layout(matrix(1:9,nrow = 3,ncol=3, byrow = T))
  for (plottreeloop in 1:3) {
    multisimmap1 <- simmapMultistate[[plottreeloop]]
    simmapQ = multisimmap1$Q
    SimmapStates = colnames(multisimmap1$Q)
    py = c("orange", "#009E73", "blue", "#CC79A7")
    py4named <- py[1:length(SimmapStates)]
    names(py4named) <- SimmapStates
    tipcols = rep(NA, length(multisimmap1$tip.label))
    for (stateNum in 1:length(SimmapStates)) { # get color vector of tips
      tempstate = SimmapStates[stateNum]
      tipcols[which(multisimmap1$tip.label %in% names(which(discretetraitvec==tempstate)))] <- py4named[tempstate]
    }
    # Plot multistate (trait 1) simmap 
    plotSimmap(multisimmap1, fsize=0.1, lwd = 0.8, colors = py4named)
    tiplabels(pch=21,bg=tipcols, col = tipcols, cex=0.1)
    legend("bottomleft",legend = names(py4named), lwd=1,col=py4named, lty = c(rep(1,length(SimmapStates)),2), cex=1, title = multistateTrait) 
    
    # trait 2 simmap
    simmap1 <- simmapTrait2[[plottreeloop]]
    py = c("black","red")
    pynamed <- py
    names(pynamed) <- c(0,1)
    
    plotSimmap(simmap1,fsize=0.1, lwd = 0.8, colors = pynamed)
    #add tips
    treetiplabels <- simmap1$tip.label %in% names(othertraitvec[othertraitvec==1]) 
    tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
    numrates1 <- lapply(FSQ,round,digits=6)
    title(paste(" ", "\nmake.simmap","Qrates:", numrates1[2],numrates1[3], trait2),cex.main = 0.5)
    
    # trait 3 simmap
    simmapQset <- simmapTrait3[[plottreeloop]]
    plotSimmap(simmapQset,fsize=0.1,lwd=0.8, colors = pynamed)
    #add tips
    treetiplabels2 <- simmapQset$tip.label %in% names(thirdtraitvec[thirdtraitvec==1]) 
    tiplabels(pch=21,bg=py[as.numeric(treetiplabels2)+1], col = py[as.numeric(treetiplabels2)+1], cex=0.1)
    numrates2 <- lapply(CBQ,round,digits=6)
    title(main=paste(" ","\nARDmodel","Qrates:",numrates2[2],numrates2[3], trait3),cex.main = 0.5)
  }
  dev.off()
  
  ## DUMMY
  
  # Make randomized version of multistate trait simmaps / DUMMY data
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
  
  # Make randomized versions of Coop Breed simmaps / DUMMY data
  CBsimtreesRand <- list()
  RandTrait3vecList <- list()
  print(paste("starting Dummy CoopBreed simmaps", Sys.time()))
  for (j in 1:nsims) {
    CBvec <- subsetdf[,columns[3]]
    CBvecRandom <- sample(CBvec)
    names(CBvecRandom) <- subsetdf$species
    CBsimtree <- make.simmap(tree = subsettree, x = CBvecRandom, model = "ARD", nsim = 1, Q = CBQ)
    CBsimtreesRand[[j]] <- CBsimtree
    RandTrait3vecList[[j]] <- CBvecRandom # KTS corrected from RandTrait2vecList 3/4/2025
    if (j == 1) {
      CBsimtreesMulti = CBsimtree
    } else {
      CBsimtreesMulti = c(CBsimtreesMulti, CBsimtree)
    }
    print(j)
    
  } # end for j in 1:nsims (CB)
  CBsimtrees<- CBsimtreesRand
  
  overlapdf = set.seed(10)
  for (k in 1:nsims) {
    # calculate overlap - real
    realOverlap = Map.Overlap.Three(simmapMultistate[[k]], simmapTrait2[[k]], simmapTrait3[[k]])
    overlapVec = convert_overlap_3d_to_vector(realOverlap, "_REAL") # function in MapOverlapThree.R
    overlapVec
    
    # calculate overlap - dummy
    dummyOverlap = Map.Overlap.Three(Coopsimtrees[[k]], FSsimtrees[[k]], CBsimtrees[[k]])
    overlapVecDummy = convert_overlap_3d_to_vector(dummyOverlap, "_DUMMY") # function in MapOverlapThree.R
    overlapVecDummy
    
    temprow = c(k, multistateTrait, othertrait, trait3, overlapVec, overlapVecDummy)
    
    overlapdf = rbind(overlapdf, temprow)
    overlapdf = as.data.frame(overlapdf)
    colnames(overlapdf) = c("tree", "trait1", "trait2", "trait3", names(overlapVec), names(overlapVecDummy))
  }
  write.csv(overlapdf, file = paste("simmap overlap", multistateTrait, othertrait, trait3, nsims, "sims.csv"), row.names = F)
  
  calcHuelout = calcHuelflex_three(overlapdf)
  pdfname = paste("simmap overlap", multistateTrait, othertrait, trait3, nsims, "sims.pdf")
  require(gridExtra)
  single_page_plotBox <- grid.arrange(grobs = calcHuelout[1:4], ncol = 1)
  ggsave(pdfname, single_page_plotBox, width = 8, height = 16, units = "in", limitsize = FALSE)
  
  write.csv(calcHuelout$fraction_dummy_less_than_median_real, file = paste("simmap overlap FractionDummyLessThanMedianReal", multistateTrait, othertrait, trait3, nsims,"sims.csv"))
}

#Terr12vs3_FS_CB_simmaps_REALDUMMY <- list(Territory_12vs3_500simmaps_REAL = simmapMultistate, Territory_12vs3_500simmaps_DUMMY = Coopsimtrees, FemaleSong_Agg01_500simmaps_REAL = simmapTrait2, FemaleSong_Agg01_500simmaps_DUMMY = FSsimtrees, HighConfidence_Coop_500simmaps_REAL = simmapTrait3, HighConfidence_Coop_500simmaps_DUMMY = CBsimtrees)
#saveRDS(object = Terr12vs3_FS_CB_simmaps_REALDUMMY, file = "Simmap Overlap Outputs/Territory_12vs3 FemaleSong_Agg01 HighConfidence_Coop 500simmaps_REAL_DUMMY.rds")
