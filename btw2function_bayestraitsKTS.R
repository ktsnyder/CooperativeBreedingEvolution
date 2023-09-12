########
#Coded by Kate T. Snyder
#Modified 10-31-2020 from Supplement2_btwfunction_newtree.R
#Modified 8/20/2021 - add newdata arg
# Modified 3/15/2022 - output nocorrD stuff too
# Modified 3/17/2022 - columns input instead of MateParam/SongParam, use BayesTraitsV4 and KTS version of Discrete, add treelabel
# Modified: 4/25/2022 - from btwfunction.R - started to replace Discrete() usage
#Last modified: 7/21/2023 - Split from btw2function_DiscreteKTS.R to make use bayestraitsKTS instead of DiscreteKTS - started 8/8/2023, mostly only tested for use with MCMC & Ace priors
#Built using RStudio Version 1.0.136
#R Version 3.4.1?
#
#ape_4.1  phytools_0.5-38   maps_3.1.0  btw_V1.0
#BayesTraitsV2 
########
########

# Split from btwfunction2.R to make it use btwDiscreteKTS.R functions

#Must set .BayesTraitsPath to location of program BayesTraitsV2, which must be located in your working directory. E.g.:
#.BayesTraitsPath <- "~/Documents/BayesTraits/BayesTraitsV2"
#.BayesTraitsPath <- "~/Documents/BayesTraitsV2"
.BayesTraitsPath <- "~/Documents/BayesTraitsV4.0.0-OSX/BayesTraitsV4"
#source("btwDiscreteKTS.R")
source("btwV2bayestraitsKTS.R")
# library(devtools)
# install_github("rgriff23/btw", ref="v1")
 library(btw)

btwfunction(columns = c("MeanCoopTie2Coop","Song.rep.final"), plot = FALSE, csvsout = TRUE, nsim = 10, newtreefile = "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex", newdata = "2023-06-20_CoopBreed-FemaleSong01HighConf-Song_Data_R.csv", treelabel = "Hackett", priors = "Ace", MCMCorML = "MCMC")
btwfunction(columns = c("MeanCoopTie2Noncoop","Song.rep.final"), plot = FALSE, csvsout = TRUE, nsim = 20, newtreefile = "/Users/kate/Desktop/CooperativeBreedingEvolution/2022-03-16ConsensusPasserineTreeHackett4_1000.nex", newdata = "2023-06-20_CoopBreed-FemaleSong01HighConf-Song_Data_R.csv", treelabel = "Hackett", priors = "Ace", MCMCorML = "MCMC")

btwfunction <- function(columns, plot=TRUE, jackknife = FALSE, csvsout = FALSE, nsim = 10, newtreefile = FALSE, newdata = FALSE, treelabel = NULL, priors = NULL, MCMCorML = "MCMC") {
  output <- list()
  output$start <- Sys.time()
  print(Sys.time())
  require(btw)
  require(phytools)
  require(base)
  #require(mnormt)
  require(geiger)
  require(R.utils)
  
  source(file = "subsettreedata.R")  #removed 10/31/2020 to test
  # source(file = "Supplement_BayesPlots.R")
  MateParam <- columns[1]
  SongParam <- columns[2]
  subset <- subsettreedata(columns = columns, newtree = newtreefile, newdata = newdata)
  subsettree <- bothtree <- subset$subsettree
  matecol <- MateParam
  songdf <- subset$subsetdf
  songcol <- SongParam
  OutputFolderPath = paste0("No-Priors_", columns[1], "-", columns[2])
  
  if (priors == "Ace") {
    source(file = "findQrates.R")
    column1qrates= tryCatch({
      findQrates(columns = columns[1], plot=FALSE, newtree = newtreefile, newdata = newdata, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL)
    }, 
    error = function(e) {
      message("An error occurred in finding q rates in column 1", e$message)
    })
    alpha1q = column1qrates$qrates[1,2]
    beta1q = column1qrates$qrates[2,1]
    OutputFolderPath = paste0("Priors-exp-AceQrates_", columns[1], "-", columns[2])
  }
  
  output$lengthsongdf <- length(songdf$species)
  
  if (jackknife ==TRUE) {
  ##Jackknife test removing each family
    familyvecNoNone <- unique(songdf$Family)
    familyvec <- c("None",as.character(familyvecNoNone))
  }
  
  
  subsetskinny <- subsettreedata(columns = columns, newtree = newtreefile, newdata = newdata, skinnydata = TRUE)
  subsetdf = subsetskinny$subsetdf
  subsetdf[which(subsetdf[,columns[2]] <= thresh),columns[2]] <- 0
  subsetdf[which(subsetdf[,columns[2]] > thresh),columns[2]] <- 1
  
  
  if (plot == TRUE) {
    pdf(paste(Sys.Date(),"Bayes",MateParam,SongParam,"Pval",nsim,"repsJacks.pdf", sep=""), width = 7, height = 10)
    par(mfrow=c(round(length(familyvec)/4)+1,4))
    #    layout(mat = matrix(nrow=20, ncol = 2, 1:40, byrow = TRUE))
    par(oma = c(1,3,3,1))
  }
  
  if (jackknife == TRUE) {
    
    transitions10reps <- data.frame()
    for (k in 1:length(familyvec)) {
      # for (k in 1:4) {
      jacksongdf <- songdf  #have to do this so songdf doesn't get whittled down every time the for loop loops
      jackedsongdf <- jacksongdf[which(jacksongdf$Family != familyvec[k]),]
      removedfamilydf <- jacksongdf[which(jacksongdf$Family == familyvec[k]),]
      
      notjackedvec <- bothtree$tip.label %in% as.character(jackedsongdf$species)
      dropforjack <- which(notjackedvec == FALSE)
      
      jacktree <- drop.tip(bothtree,tip = dropforjack)
      removedfamily <- familyvec[k]
      matevec <- as.character(jackedsongdf[,matecol])
      names(matevec) <- jackedsongdf$species
      songcontvec <- jackedsongdf[,songcol]
      songcontvec <- unique(songcontvec)
      print(paste("Without",familyvec[k],length(songcontvec)))
      LRstatall <- set.seed(10)
      LRpvalall <- set.seed(10)
      songcontvecall <- set.seed(10)
      transitions10reps <- data.frame(matrix(nrow=nsim*length(songcontvec), ncol=0))
      for (j in 1:nsim) {
        if (j == 1) {
          print(paste("Rep:",j, "Number of corrD points: 0", Sys.time()))
        } else if (j %in% c(10,20,40,60,80,100)) {
          print(paste("Rep:",j)) #,"Number of corrD points",length(transitions10reps[,3]), Sys.time()))
        }
        LRstat <- set.seed(10)
        LRpval <- set.seed(10)
        transitions <- data.frame(matrix(nrow=length(songcontvec), ncol=0)) 
        for (i in 1:length(songcontvec)) {
          transandp <- set.seed(10)
          songdiscvec <- jackedsongdf[,songcol]
          thresh = songcontvec[i]
          songdiscvec[songdiscvec <= thresh] <- 0
          songdiscvec[songdiscvec > thresh] <- 1
          names(songdiscvec) <- jackedsongdf$species
          btwdfjack <- as.data.frame(cbind(as.character(jackedsongdf$species),matevec,as.character(songdiscvec)))
          #nocorrD <- Discrete(jacktree, btwdfjack)
          commandIndML <- c("2","1", "Seed 10") # replace nocorrD
          IndMLout <- bayestraits(df,tree,commandIndML)
          
          #corrD <- Discrete(jacktree, btwdfjack, dependent=TRUE)
          commandDepML <- c("3","1", "Seed 10")#, "res q31 q42 1.5")  # replace corrD
          DepMLout <- bayestraits(df,tree,commandDepML, silent = FALSE)
          
          
          lrtestresults <- lrtest(corrD, nocorrD)
          LRstat[i] <- lrtestresults$LRstat
          LRpval[i] <- lrtestresults$pval
          transandp <- cbind(corrD,LRstat[i],LRpval[i],songcontvec[i], nocorrD)
          transitions <- rbind(transitions, transandp)
        } #end for i in 1:length(songcontvec)
        LRstatall <- c(LRstatall,LRstat)
        LRpvalall <- c(LRpvalall,LRpval)
        songcontvecall <- c(songcontvecall,songcontvec)
        transitions10reps <- rbind(transitions10reps,transitions)
      } #end for j in 1:nsim
      
      print(length(transitions10reps[,2]))
      
      if (csvsout == TRUE) {
        write.csv(transitions10reps,file=paste(Sys.Date(),"Bayes",MateParam,SongParam,nsim,"repsJack",familyvec[k],".csv",sep=""))
        csvfile <- paste(Sys.Date(),"Bayes",MateParam,SongParam,nsim,"repsJack",familyvec[k],".csv",sep="")
      }
      jackdf <- transitions10reps
      colnames(jackdf[,c(15:17)]) <- c("LRstat","LRpval","songcontvec")
      
      
      d = data.frame(songcontvecall,LRstatall, LRpvalall)
      
      if (plot == TRUE) {
        # par(mar = c(3,2,2,1))
        #   with(d, plot(songcontvecall, LRpvalall, pch=16, col="red3", 
        #                xlab=paste(SongParam, "threshold for binary (low/high) categorization"),ylab="p-value", main = paste(MateParam,"&",SongParam, "sign. of dependent correlation,","\nN =",length(songcontvec),"- Removed", paste(familyvec[k])), cex.main = 0.5, log = "x"
        #   ))
        plotBTjacks(MateParam,SongParam, d = jackdf, familysplit = removedfamily, nsim = nsim)
        
      } #end if plot ==TRUE
    } #end for loop covering family names
    #dev.off()
  } else { # if jackknife is anything but TRUE
    
    matevec <- as.character(songdf[,matecol])
    names(matevec) <- songdf$species
    songcontvec <- songdf[,songcol]
    songcontvec <- sort(unique(songcontvec))
    output$songcontvec <- songcontvec
    print(length(songcontvec))
    LRstatall <- set.seed(100)
    LRpvalall <- set.seed(100)
    
    songcontvecall <- set.seed(100)
    transitions100reps <- data.frame()
    outDepdfAll = set.seed(10)
    outInddfAll = set.seed(10)
    
    for (i in 1:length(songcontvec)) {
      LRstat <- set.seed(10)
      LRpval <- set.seed(10)
      transitions <- data.frame(matrix(nrow=length(songcontvec), ncol=0))
      transandp <- set.seed(10)
      songdiscvec <- songdf[,songcol]
      thresh = songcontvec[i]
      songdiscvec[songdiscvec <= thresh] <- 0
      songdiscvec[songdiscvec > thresh] <- 1
      names(songdiscvec) <- songdf$species
      btwdf <- as.data.frame(cbind(as.character(songdf$species),matevec,as.character(songdiscvec)))
      colnames(btwdf) = c("species", columns)
      
      if (priors == "Ace") {
        binarysongdf = read.csv(newdata)
        binarysongdf[which(binarysongdf[,columns[2]] <= thresh),columns[2]] <- 0
        binarysongdf[which(binarysongdf[,columns[2]] > thresh),columns[2]] <- 1
        column2qrates= tryCatch({
          findQrates(columns = columns[2], plot=FALSE, newtree = newtreefile, newdata = binarysongdf, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL)
        }, 
        error = function(e) {
          print("An error occurred in finding q rates in column 2, using prev threshold's q rates")
          return(useCol2Qrates)
        })
        useCol2Qrates <- column2qrates # if a threshold doesn't work with findQrates because of all being 0s or something, will use the previous threshold's Ace Q rates
        alpha2q = column2qrates$qrates[1,2]
        beta2q = column2qrates$qrates[2,1]
        AdditionalCommandsDep = c(paste("Prior q12 exp", alpha2q), paste("Prior q13 exp", alpha1q), paste("Prior q21 exp", beta2q), paste("Prior q24 exp", alpha1q), paste("Prior q31 exp", beta1q), paste("Prior q34 exp", alpha2q), paste("Prior q42 exp", beta1q), paste("Prior q43 exp", beta2q), "burnin 220000", "Stones 100 1000")
        AdditionalCommandsInd = c(paste("Prior alpha1 exp",alpha1q), paste("Prior beta1 exp", beta1q), paste("Prior alpha2 exp", alpha2q), paste("Prior beta2 exp", beta2q), "burnin 220000", "Stones 100 1000")
      } else {
        AdditionalCommandsDep = NULL
        AdditionalCommandsInd = NULL
        alpha2q = NULL
        beta2q = NULL
        alpha1q = NULL
        beta1q = NULL
      }
      
      outDepdf = set.seed(10)
      outInddf = set.seed(10)
      for (j in 1:nsim) {
        print(paste("Rep:",j, Sys.time()))
        # if (j == 1) {
        #   print(paste("Rep:",j, "Number of corrD points: 0", Sys.time()))
        # } else if (j %in% c(10,20,30,40,50,60,70,80,90)) {
        #   print(paste("Rep:",j, Sys.time())) #,"Number of corrD points",length(transitions100reps[,3]),Sys.time()))
        # }
        
        
        if (MCMCorML == "MCMC") {
          MLorMCMC = "2"
        } else if (MCMCorML == "ML") {
          MLorMCMC = "1"
        }
        
        nocorrCommands = c("2", MLorMCMC, AdditionalCommandsInd)
        corrCommands = c("3", MLorMCMC, AdditionalCommandsDep)
        
        outPriorAllInd = bayestraitsKTS(data = subsetdf, tree = subsettree, commands = nocorrCommands, remove_files = F, BTdirpath = "~/Documents", silent = T, OutputFolderPath = OutputFolderPath, TestPrior = TestPrior)
        outPriorAllDep =  bayestraitsKTS(data = subsetdf, tree = subsettree, commands = corrCommands, remove_files = F, BTdirpath = "~/Documents", silent = T, OutputFolderPath = OutputFolderPath, TestPrior = TestPrior)
        
        if (MCMCorML == "ML") { 
          lrtestresults <- lrtestV1(corrD, nocorrD)
          # LRstat[i] <- lrtestresults$LRstat
          # LRpval[i] <- lrtestresults$pval
          # transandp <- cbind(corrD,LRstat[i],LRpval[i],songcontvec[i], nocorrD)
          LRstat <- lrtestresults$LRstat
          LRpval <- lrtestresults$pval
          transandp <- cbind(corrD,LRstat,LRpval,thresh, nocorrD)
          transitions <- rbind(transitions, transandp)
        } else if (MCMCorML == "MCMC") {
          require(stringr)
          
          # Independent
          outPriorAllIndOptions <- outPriorAllInd$Log$options
          outPriorAllIndResults <- outPriorAllInd$Log$results
          outPriorAllIndStonesLh <- outPriorAllInd$Stones$logMarLH
          Model = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Model")]), "Model: ")
          #Iterations = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Iterations")]), "Iterations: ")
          #BurnIn = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Burn in")]), "Burn in: ")
          if (sum(str_detect(outPriorAllIndOptions, "Iterations")) == 1) {
            Iterations = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Iterations")]), "Iterations: ")
          } else {Iterations = NA}
          if (sum(str_detect(outPriorAllIndOptions, "Burn in")) == 1) {
            BurnIn = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Burn in")]), "Burn in: ")
          } else {BurnIn = NA}
          Seed = str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Seed")]), "Seed: ")
          ScheduleFile = NA #str_remove(str_squish(outPriorAllIndOptions[str_detect(outPriorAllIndOptions, "Schedule File:")]), "Schedule File: ")
          #StonesLh = outPriorAllIndStonesLh
          if (!is.null(outPriorAllIndStonesLh)) {
            StonesLh = outPriorAllIndStonesLh
          } else {StonesLh = NA}
          meanResultsInd = apply(outPriorAllIndResults[,c("Lh", "alpha1", "beta1", "alpha2", "beta2")], 2, mean)
          names(meanResultsInd)[which(names(meanResultsInd) == "Lh")] <- "SamplingMeanLh"
          outIndrow = c(i, Model, Iterations, BurnIn, Seed, ScheduleFile, StonesLh, meanResultsInd, thresh)
          names(outIndrow) = c("Sim", "Model", "Iterations", "BurnIn", "Seed", "ScheduleFile", "StonesLh", "SamplingMeanLh", "alpha1", "beta1", "alpha2", "beta2", "Threshold")
          outInddf = rbind(outInddf, outIndrow)
          outInddf = as.data.frame(outInddf)
          SamplingMeanLhInd = meanResultsInd["SamplingMeanLh"] # for calculating bayesfactor w/ dep results
          StonesLhInd = StonesLh # for calculating bayesfactor w/ dep results
          
          # Dependent
          outPriorAllDepOptions <- outPriorAllDep$Log$options
          outPriorAllDepResults <- outPriorAllDep$Log$results
          outPriorAllDepStonesLh <- outPriorAllDep$Stones$logMarLH
          Model = str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Model")]), "Model: ")
          if (sum(str_detect(outPriorAllDepOptions, "Iterations")) == 1) {
            Iterations = str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Iterations")]), "Iterations: ")
          } else {Iterations = NA}
          if (sum(str_detect(outPriorAllDepOptions, "Burn in")) == 1) {
            BurnIn = str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Burn in")]), "Burn in: ")
          } else {BurnIn = NA}
          Seed = str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Seed")]), "Seed: ")
          ScheduleFile = NA #str_remove(str_squish(outPriorAllDepOptions[str_detect(outPriorAllDepOptions, "Schedule File:")]), "Schedule File: ")
          if (!is.null(outPriorAllDepStonesLh)) {
            StonesLh = outPriorAllDepStonesLh
          } else {StonesLh = NA}
          meanResults = apply(outPriorAllDepResults[,c("Lh", "q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43")], 2, mean)
          names(meanResults)[which(names(meanResults) == "Lh")] <- "SamplingMeanLh"
          
          BayesFactor_SamplingMeanLh = 2*(meanResults["SamplingMeanLh"] - SamplingMeanLhInd)
          BayesFactor_StonesLh = 2*(StonesLh - StonesLhInd)
          
          outDeprow = c(i, Model, Iterations, BurnIn, Seed, ScheduleFile, StonesLh, meanResults, thresh, alpha1q, beta1q, alpha2q, beta2q, StonesLhInd, SamplingMeanLhInd,  BayesFactor_SamplingMeanLh, BayesFactor_StonesLh)
          names(outDeprow) = c("Sim", "Model", "Iterations", "BurnIn", "Seed", "ScheduleFile", "StonesLh", "SamplingMeanLh", "q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43", "Threshold", "ACEpriorAlpha1", "ACEpriorBeta1", "ACEpriorAlpha2", "ACEpriorBeta2", "StonesLhInd", "SamplingMeanLhInd", "BayesFactor_SamplingMeanLh", "BayesFactor_StonesLh")
          outDepdf = rbind(outDepdf, outDeprow)
          outDepdf = as.data.frame(outDepdf)
          
        } # end   else if (MCMCorML == "MCMC")
      } # end j in 1:nsims
      
      # if (i %in% seq(1,100, by = 9)) {
      #   filenameDep = paste0(OutputFolderPath, "/",columns[1], "-", columns[2], "_", "Dependent_", Sys.Date(),".csv")
      #   write.csv(outDepdf, file = filenameDep, row.names = FALSE)
      #   filenameInd = paste0(OutputFolderPath, "/", columns[1], "-", columns[2], "_", "Independent_", Sys.Date(),".csv")
      #   write.csv(outInddf, file = filenameInd, row.names = FALSE)
      #   print(paste("Finished threshold loop #", i, "with columns", columns[1], columns[2], "at", Sys.time()))
      # } # end if i in seq(1,100, by = 9))
      
      if (MCMCorML == "ML") {
        LRstatall <- c(LRstatall,LRstat)
        LRpvalall <- c(LRpvalall,LRpval)
        songcontvecall <- c(songcontvecall,songcontvec)
        transitions100reps <- rbind(transitions100reps,transitions)
        print(length(transitions100reps[,2]))
        colnames(transitions100reps)[which(colnames(transitions100reps) == "thresh")] <- "songcontvec"
        df <-  transitions100reps
        d = data.frame(songcontvecall,LRstatall, LRpvalall)
        if (csvsout == TRUE) {
          write.csv(d,file=paste(Sys.Date(),"Bayes",MateParam,SongParam,treelabel,nsim,"reps.csv",sep=""))
          csvfile <- paste(Sys.Date(),"Bayes",MateParam,SongParam,treelabel,nsim,"reps.csv",sep="")
          output$nojackdf<- transitions100reps
        }
      } else if (MCMCorML == "MCMC") {
        outInddfAll = rbind(outInddfAll, outInddf)
        outDepdfAll = rbind(outDepdfAll, outDepdf)
        csvfileInd = paste(MateParam, SongParam, "BayesIndependent", MCMCorML, priors, treelabel, nsim, "reps", Sys.Date(), ".csv")
        csvfileDep = paste(MateParam, SongParam, "BayesDependent", MCMCorML, priors, treelabel, nsim, "reps", Sys.Date(), ".csv")
        write.csv(outInddfAll, file = csvfileInd)
        write.csv(outDepdfAll, file = csvfileDep)
      } # end else if (MCMCorML == "MCMC")
      print(paste("Finished threshold loop #", i, "with columns", columns[1], columns[2], "at", Sys.time()))
    } # end i in 1:length(songcontvec)

    #colnames(df[,c(15:17)]) <- c("LRstat","LRpval","songcontvec")
    
    if (plot == TRUE) {
      transitionBinplots(MateParam,SongParam,df = df,newpdf = TRUE, nsim = nsim)
    } #end if plot == True 
    
  } #end for loop (jackknife == FALSE)
  # dev.off()
  output$finish <- Sys.time()
  output$runtime <- output$finish - output$start
  return(output)
} #end btwfunction


cyclematesongs <- function(jackknife = FALSE, csvsout = FALSE, plotbt = TRUE, nsim = 10, treevec = trees) {
  for (n in 1:length(MateParams)) {
    for (m in 1:length(SongParams)) {
      MateParam = MateParams[n]
      SongParam = SongParams[m]
      for (i in treevec) {
        temptree <- paste0("samplematezillaHack",i,".nex")
      btwfunction(MateParam,SongParam,plot = plotbt,jackknife = jackknife,csvsout = csvsout, newtreefile = temptree, nsim = nsim)
      } #end for i in 10:29 going through trees
      #    transitionplots(MateParams[n],SongParams[n],csvfiles[n]) 
    } #end for m in 1:songparams
  } #end for n in 1:MateParams
} #end function BTfunc

