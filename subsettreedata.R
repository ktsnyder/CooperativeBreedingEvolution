########
#Coded by Kate T. Snyder
#Last Modified 9-14-2023
#Built using R Version 4.0.2
#
#ape_5.3  phytools_0.5-38   maps_3.1.0  btw_V1.0
########
########

################## Subset bird data and generate associated trees##################
 
#### Inputs: ----
## columns: vector of character values matching names of columns in dataframe. Will subset dataframe to only the subset of species which do not have NA values in any of the specified columns
## cladesubsetcolumn: a character object (aka a character vector of length 1) that is the column name of whichever column lists classification (e.g. Family) for a species
## cladesubsetvalue: a character object (aka a character vector of length 1) that is the clade/family (e.g. "Corvidae") that you want to analyze fdata for
## newdata: can be either a dataframe object (= the name of the dataframe in R environment, not in quotations) or a character object (= the name of the data file ending in .csv, in quotations)
## newtree: can be either an object of class "phylo" or a character object (= the name of the file ending in ".nex", in quotations)
## islog: defaults to FALSE, meaning that none of the data will be log-transformed. If we do want to log-transform some data, this should be a character object/vector (= the name of the column(s) of the data that you want to take the log() of, in quotations)
## suppresswarning: TRUE (default) or FALSE. Determines whether certain warnings that occur while running the code will be printed in the console or not
## skinnydata: TRUE or FALSE (default). Return dataframe subsetted just to columns mentioned plus species (including columns and cladesubsetcolumn)
 
#### Outputs (list of 2 items): ----
## subsettree - subsetted tree (.nex format)
## subsetdf - subsetted dataframe
## NonMatchedSpecies - vector of any species in data that are not in the tree

## 
#### Examples: ----
#outputall <- subsetbirddata(columns = c("PolygynyOrMonog", "EPP10threshold"), newdata = "20200123_SongDatabaseUpdate2019_R.csv", newtree = "2019-10-22matezilla2treeHack_nondicho.nex")
#outputnonpass <- subsetbirddata(columns = c("PolygynyOrMonog", "EPP10threshold"), cladesubsetcolumn = "oscine", cladesubsetvalue = "nonpasserine", newdata = "20200123_SongDatabaseUpdate2019_R.csv", newtree = "2019-10-22matezilla2treeHack_nondicho.nex")
#subsettreedata(columns = c("Final.polygyny", "Syllable.rep.final"), cladesubsetcolumn = "Family", cladesubsetvalue = c("Acrocephalidae", "Icteridae"), islog = "Syllable.rep.final", suppresswarning = FALSE, skinnydata = TRUE)

## 1/9/2020 update: make an input "columns" vector of all the columns to subset by - just for removing all rows with NAs
## 1/24/2020 update (v2.1): add arguments "cladesubsetcolumn" and "cladesubsetvalue" - can jackknife for oscines with cladesubsetcolumn = "oscine", cladesubsetvalue = "Oscine". Can jackknife by family similarly. Removed argument passeriformesonly. v2.1 completed 1/28/2020
## 2/11/2020 - changed data file to "2020-02-11_Song Database Update 2019.csv"
## 5/5/2020 - changed data file to "SnyderCreanza_NatComms2019_SupplementalData_R.csv", tree file to "birdzillatreeMaybeConsensus.nex", Removed "minmax = NA" from inputs, renamed file from "subsetbirddata2.1.R", Replace "matezillatree" with "fulltree", replace "$BirdtreeFormat" with [,1] (will make the first column the species column), add islog functionality
## 5/6/2020 - changed function from subsetbirddata to subsettreedata
## 5/12/2020 - now makes subfolder "OutputFiles" for output files
## 10-9-2020 - made ability to put in new data and tree directly rather than just the file name to be read; added code to remove any species from data that don't match a species in Phylogeny and tell the user that that happened - also ability to suppress that warning; nonmatchedspecies now outputted
## 6/4/2021 - can input vector of values as cladesubsetvalue to subset (rather than only one value); commented out dir.create; defaults columns = NULL - if not specified, will allow subsetting by just clade; change name of first column to "species"; added skinnydata option to return dataframe subsetted just to columns mentioned (including columns and cladesubsetcolumn); cladesubsetcolumn defaults to NULL instead of FALSE; commented out requiring packages "ape" and "maps")
## 9/13/2023 - changed order of newdata import to check whether it's dataframe first, then character, then FALSE because newdata == FALSE was throwing an error instead of a warning
## 
## 
## new objectives: 
##    ability to select clade based on a formula? like != "Passeriformes" or something
##    option to save file/tree? if not then probably get rid of subDir creation
##    return counts of each data type?
##    remove some # of species at random?
##    may have to change skinnydata to fit BayesTraits


#### Function ----
subsettreedata <- function(columns = NULL, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = FALSE, newtree = FALSE, islog = FALSE, suppresswarning = TRUE, skinnydata = FALSE) {
  require(phytools)
  #require(ape)
  #require(maps)
  
  # mainDir <- getwd()
  # subDir <- "OutputFiles"
  # if (!dir.exists(file.path(mainDir,subDir))) {
  #   dir.create(file.path(mainDir, subDir))
  # }
  
  output <- list()
  
  suppressWarnings(  # suppressing the warning that appears after last "if" if newdata is length > 1
    if (is.data.frame(newdata)) {
      alldatadfos <- newdata
    } else if (is.character(newdata)) {
      alldatadfos <- as.data.frame(read.csv(newdata))
    } else if (newdata == FALSE) {
      alldatadfos <- as.data.frame(read.csv("SnyderCreanza_NatComms2019_SupplementalData_R.csv"))
    }
  ) #end suppressWarnings
  
  #print(head(alldatadfos))
  
  
  if (is.character(newtree)) {
    fulltree <- read.nexus(newtree)
  } else if (is.list(newtree)) {
    fulltree <- newtree
  } else if (newtree == FALSE) {
    fulltree <- read.nexus("birdzillatreeMaybeConsensus.nex")
  } else {
    print("The tree is not a tree.")
  }
  
  alldatadfos[,1] <- as.character(alldatadfos[,1])
  colnames(alldatadfos)[1] <- "species"
  
  if (class(fulltree) == "multiPhylo") {
    subsettedMultiTree <- fullMultiTree <- fulltree
    fulltree <- fulltree[[1]]
  } else {
    subsettedMultiTree = NULL
  }
  
  
  ## check that all species in dataframe match a species in phylogeny
  matchtree <- alldatadfos[,1] %in% fulltree$tip.label
  alldatadfos <- alldatadfos[which(matchtree),]
  nonmatchedspecies <- alldatadfos[,1][!matchtree]
  numnonmatching <- length(nonmatchedspecies)
  
  if (suppresswarning == FALSE) {
    print(paste("There were", numnonmatching, "species in the data table that did not match a species in the phylogeny. These species were removed from the data:"))
    print(paste(nonmatchedspecies))
  }
  
  if (is.null(cladesubsetcolumn)) {
    alldatadf <- alldatadfos
  } else {
    includespecies <- alldatadfos[,1][which(alldatadfos[,cladesubsetcolumn] %in% cladesubsetvalue)]
    includetips <- match(includespecies, fulltree$tip.label)
    includetree <- keep.tip(fulltree, tip=includetips)
    fulltree <- includetree
    includerows <- alldatadfos[,1] %in% includespecies
    alldatadf <- alldatadfos[which(includerows),]
  }
  
  ####
  subsetdf <- alldatadf
  subsettree <- fulltree  # if multitree input, this is one tree pared down by clade subset, won't actually use this
  
  #print(subsettree)
  
  if (!is.null(columns)) {
    for (i in 1:length(columns)) {
      whichcol <- columns[i]
      subsetdf <- subsetdf[!is.na(subsetdf[,whichcol]),]
      havedatavec <- subsettree$tip.label %in% as.character(subsetdf[,1])
      matingtips <- which(havedatavec == TRUE)
      dropfortree <- which(havedatavec == FALSE)
      subsettree <- drop.tip(subsettree, tip = dropfortree)
      if (!is.null(subsettedMultiTree)) {
        keepspecies <- subsetdf[,1]
        OneTreeFromMultiTreeInput <- subsettedMultiTree[[1]]
        AllSpecies <- OneTreeFromMultiTreeInput$tip.label
        DropSpeciesMultiTree <- AllSpecies[which(!AllSpecies %in% keepspecies)]
        subsettedMultiTree <- lapply(subsettedMultiTree, drop.tip, tip = DropSpeciesMultiTree)
        class(subsettedMultiTree) <- "multiPhylo"
      }
      
      #   print(subsettree)
    }  # end for (i in length(columns))
  }  # end if !is.null(columns)
  
  if (islog != FALSE) {
    for (j in 1:length(islog)) {  #islog should be a character vector with column names of columns to be log-transformed
      whichcol <- islog[j]
      subsetdf[,whichcol] <- log(subsetdf[,whichcol])
    }
  } # end if (islog != FALSE)
  
  if (skinnydata == TRUE) {
    subsetdf <- subsetdf[, c("species", cladesubsetcolumn, columns)]
  }
  
  ####
  
  ditree <- multi2di(subsettree)
  ditree$edge.length[ditree$edge.length == 0] <- 0.0000000000000000001
  
  output$subsettree <- ditree
  output$subsetdf <- subsetdf
  output$NonMatchedSpecies <- nonmatchedspecies
  
  if (class(subsettedMultiTree) == "multiPhylo") {
    output$subsettree <- subsettedMultiTree
  }
  
  return(output)
}


