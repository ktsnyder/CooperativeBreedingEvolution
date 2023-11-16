
########
#Coded by Kate T. Snyder
#Last Modified 11-15-2023
#Built using RStudio Version 1.1.453
#R Version 3.4.2
#
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#BayesTraitsV2
########
# Modified 3/11/2022 - attempt to make it work with all the new code?
# Modified 3/16/2022 - reinstate ability to jackknife by all unique values in cladesubsetcolumn
# Modified 11/15/2023 - global Q rates now from just column 1 subset; changed column names in output to be DiscreteTrait and ContinuousTrait; add number of species to output; new output folder; more informative message on progress in loop

#Jackknifing brownie
# Do need to set cladesubsetcolumn to the column in the dataframe. cladeJackvalues can either be NULL, which will cause it to jackknife by each unique value in that column; otherwise, can specify which values to jack
# Must either set output to something (e.g. jackout <- jackbrowniefunction(...) or have allcsvs = TRUE)
# islog should be true or false, will apply to the continuous variable
# otherlabel should probably be the tree label

# e.g. jackbrowniefunction(columns = c("MeanCoopTie2Noncoop", "Song.rep.final"), islog = TRUE, matemodel = "ARD", matensim = 100, allcsvs = TRUE, plotsimmaps = FALSE, newtree = treefile, newdata = newdata, cladesubsetcolumn = "Family3_BirdtreeMatchSpecies2", cladeJackvalues = "Mimidae", otherlabel = "HackettOscine")

jackbrowniefunction <- function(columns, islog = FALSE, matemodel = "ARD", matensim = 10, allcsvs = FALSE, plotsimmaps = FALSE, newtree, newdata, cladesubsetcolumn = NULL, cladeJackvalues = NULL, otherlabel = NULL) {
  require(ape)
  require(phytools)
  require(base)
  require(mnormt)
  require(R.utils)
  source(file = "subsettreedata.R")
  source(file = "findQrates.R")
  output <- list()
  output$start <- Sys.time()
  print(Sys.time())
  
  MateParam = columns[1]
  SongParam = columns[2]
  
  subsetout <- subsettreedata(columns=columns, newtree = newtree, newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, islog = FALSE)  # these should be cladesubset__ = NULL since jackknifing happens later
  bothtree <- tree <- subsetout$subsettree
  songdf <- df <- subsetout$subsetdf
  #discretetraitvec <- df[,columns[1]]
  #names(discretetraitvec) <- df[,1]
  #continuoustraitvec <- df[,columns[2]]
  #names(continuoustraitvec) <- df[,1]  #species names 
  matecol <- MateParam
  songcol <- SongParam
  # subset <- subsetbirddata(MateParam = MateParam,SongParam = SongParam)
  # songdf <- subset$df
  # bothtree <- subset$ditree
  #  songdatavec <- subset$songcontvec
  #  matingdatavec <- subset$matevec
  # songcol <- subset$songcol
  # matecol <- subset$matecol
  
  ##Jackknife test removing each family
  if (is.null(cladeJackvalues)) {
    cladeJackvalues <- unique(songdf[,cladesubsetcolumn])
  }
  familyvec <- c("None",cladeJackvalues)
  brownielist <- list()
  
  Qoutput <- findQrates(columns = columns[1], plot=FALSE, newtree = newtree, newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = otherlabel)
  qrates <- Qoutput$qrates
  #print(qrates)
  # Qoutput <- findQrates(MateParam = MateParam, SongParam = "none")
  # qrates <- Qoutput$qrates
  
  for (k in 1:length(familyvec)) {
    jackedsongdf <- songdf  #have to do this so songdf doesn't get whittled down every time the for loop loops
    jackedsongdf <- jackedsongdf[which(jackedsongdf[, cladesubsetcolumn] != familyvec[k]),]
    
    numSpecies = length(jackedsongdf[,1])
    
    print(familyvec[k])
    
    notjackedvec <- bothtree$tip.label %in% as.character(jackedsongdf$species)
    dropforjack <- which(notjackedvec == FALSE)
    
    jacktree <- drop.tip(bothtree,tip = dropforjack)
    
    matevec <- as.character(jackedsongdf[,matecol])
    names(matevec) <- jackedsongdf$species
    songcontvec <- jackedsongdf[,songcol]
    names(songcontvec) <- jackedsongdf$species
    if (islog) {
      songcontvec = log(songcontvec)
      loglabel = "log"
    } else {
      loglabel = NULL
    }
    
    set.seed(10)
    print(paste("Beginning simmap for brownie",MateParam,SongParam)) 
    starttimebrownie <- Sys.time()
    
    if (plotsimmaps == TRUE) {
      print("Detour to plot simmaps...")
      findQrates(columns = columns, plot=plotsimmaps, newtree = newtree, newdata = jackedsongdf, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = paste0("Jacked",familyvec[k]), GlobalQrates = qrates)
    }
    
    simmappy <- make.simmap(jacktree,matevec,nsim=matensim, Q=qrates, message = FALSE) 
    print(class(simmappy))
    simmapsdone <- Sys.time()
    simmaptime <- simmapsdone - starttimebrownie
    print(paste("Simmaps generated. That step took this much time: ", simmaptime, ".  Starting for loop with ", matensim, " loops.", MateParam, SongParam, familyvec[k], k, "out of", length(familyvec), "jacks. ", numSpecies, "species in this jackknifed tree."))
    
    browniedata <- data.frame(DiscreteTrait=character(matensim),ContinuousTrait=character(matensim),Pval=numeric(matensim),ERRate=numeric(matensim),ERloglik=numeric(matensim),ERace=numeric(matensim),ARDRate0=numeric(matensim),ARDRate1=numeric(matensim),ARDloglik=numeric(matensim),ARDace=numeric(matensim),k2=numeric(matensim),convergence=character(matensim),simmapnumber=integer(matensim),jackedfam=character(matensim), stringsAsFactors = FALSE)
    
    
    for (i in 1:matensim) {
      simmapfor <- simmappy[[i]]
      brownieliteresults <- set.seed(10)
      
      tryCatch(
        expr = {
          withTimeout(expr={
            
            brownieliteresults <- brownie.lite(simmapfor,songcontvec,maxit=75000)
            browniedata[i,3] <- brownieliteresults$P.chisq
            browniedata[i,4] <- brownieliteresults$sig2.single
            browniedata[i,5] <- brownieliteresults$logL1
            browniedata[i,6] <- brownieliteresults$a.single
            browniedata[i,7] <- brownieliteresults$sig2.multiple[1]
            browniedata[i,8] <- brownieliteresults$sig2.multiple[2]
            browniedata[i,9] <- brownieliteresults$logL.multiple
            browniedata[i,10] <- brownieliteresults$a.multiple
            browniedata[i,11] <- brownieliteresults$k2
            browniedata[i,12] <- as.character(brownieliteresults$convergence)
            browniedata[i,13] <- i}, timeout = 16, cpu=Inf, onTimeout = "error")
        },
        TimeoutException = function(ex) {browniedata[i,3:13]<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,"timeout",i);
        print(paste("timeout",i));
        })
      
      browniedata[i,1] <- MateParam
      browniedata[i,2] <- SongParam
      browniedata[i,14] <- paste(familyvec[k])
      browniedata$numSpecies = numSpecies
      
      if (i %in% seq(0,2000,by=45)) {
        print(paste("End brownie loop iteration",i,Sys.time()))
      }
    }  #end for loop 1:matensim
    if (allcsvs == TRUE) {
      if (!dir.exists("BrownieJackknifeOutputs")) {
        dir.create("BrownieJackknifeOutputs")
      }
      outfile = paste0("BrownieJackknifeOutputs/", Sys.Date(), otherlabel, " Brownie", matensim, "sims ", MateParam, " ", loglabel, SongParam, " jacked", familyvec[k], ".csv")
      write.csv(file = outfile, x = browniedata, row.names = FALSE)
    }
    brownielist[[k]] <- browniedata
    
  } #end for loop going thru each family
  
  outputty <- list()
  outputty$brownielist <- brownielist
  outputty$familyvec <- familyvec
  return(outputty)
  
} #end jackfunction 
