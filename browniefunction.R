#Coded by Kate T. Snyder
#Last Modified 8-18-2020
#Built using RStudio Version 1.1.453
#R Version 3.5.2?
#
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#BayesTraitsV2
########
# 
# Jan2020 update: make compatible with findQratesNewTree2.0.R and subsetbirddata2.0.R
# 5/12/2020: matensim--> nsim, remove plot arg, added args
# 8/18/2021: add arg plotsimmaps (TRUE/FALSE) to go into findQrates; add to tryCatch ability to accommadate compute error; changed default arg "cladesubsetcolumn = FALSE" to "= NULL"
# 8/26/2021: added otherlabel arg - for csv output and to feed into findQrates
# 4/26/2024: changed object returned from browniedf (no longer used) to browniedata
# 4/27/2024: added findQrates outputs to browniedata columns
# 
# e.g.
# brownieout <- browniefunction(c("Final.polygyny", "Syllable.rep.final"), islog = "Syllable.rep.final", nsim = 500)

browniefunction <- function(columns, newtree = FALSE, newdata = FALSE, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, nsim = 500, islog = FALSE, phylanovaP = "not evaluated", plotsimmaps = FALSE, otherlabel = NULL) {
  require(R.utils)
  require(phytools)
 # require(ape)
#  require(base)
  
  source(file = "findQrates.R")
 
  print("Beginning simmaps for brownie") 
  starttimebrownie <- Sys.time()
  
  
  subsetout <- subsettreedata(columns=columns, newtree = newtree, newdata = newdata, cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue, islog = islog)
  tree <- subsetout$subsettree
  df <- subsetout$subsetdf
  discretetraitvec <- df[,columns[1]]
  names(discretetraitvec) <- df[,1]
  continuoustraitvec <- df[,columns[2]]
  names(continuoustraitvec) <- df[,1]  #species names 
  
  #phylanova <- phylANOVA(tree,discretetraitvec,continuoustraitvec, nsim=nsim)
  #phylanovaP = phylanova[[2]] #pval 
  
  Qoutput <- findQrates(columns, plot=plotsimmaps, newtree = newtree, newdata = newdata, cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue, otherlabel = otherlabel)
  qrates <- Qoutput$qrates
  print(qrates)
  
  # make OutputFiles folder if not present
  mainDir <- getwd()
  subDir <- "OutputFiles"
  if (!dir.exists(file.path(mainDir,subDir))) {
    dir.create(file.path(mainDir, subDir))
  }
  
  simmappy <- make.simmap(tree,discretetraitvec,nsim=nsim,Q=qrates) 
  write.simmap(simmappy, file=paste(getwd(),"/OutputFiles/",Sys.Date(),columns[1], columns[2], otherlabel, nsim,"simmaps",".txt",sep=""))
  simmapsdone <- Sys.time()
  simmaptime <- simmapsdone - starttimebrownie
  
  print(paste("Simmaps generated. That step took this much time: ", simmaptime, ".  Starting for loop with ", nsim, " loops."))
  
  browniedata <- data.frame(DiscreteTrait=character(nsim),ContinuousTrait=character(nsim),Pval=numeric(nsim),ERRate=numeric(nsim),ERloglik=numeric(nsim),ERace=numeric(nsim),ARDRate0=numeric(nsim),ARDRate1=numeric(nsim),ARDloglik=numeric(nsim),ARDace=numeric(nsim),k2=numeric(nsim),convergence=character(nsim),simmapnumber=integer(nsim),phylanovaP=numeric(nsim),stringsAsFactors = FALSE)
  brownied <- set.seed(10)
  browniedf <- set.seed(10)
  
  for (i in 1:nsim) {
    simmapfor <- simmappy[[i]]
    brownieliteresults <- set.seed(10)
    
    tryCatch(
      expr = {
        withTimeout(expr={
          
          brownieliteresults <- brownie.lite(simmapfor,continuoustraitvec,maxit=75000)
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
      TimeoutException = function(ex) {browniedata[i,2:13]<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"timeout",i);
      print(paste("timeout",i));
      },
      error = function(e) {browniedata[i,2:13]<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"error",i);
      print(paste("compute error",i));
      })
    
    browniedata[i,1] <- columns[1] #column name of Discrete trait
    browniedata[i,2] <- columns[2] #column name of Continuous trait
    browniedata[i,14] <- phylanovaP
    
    if (i %in% c(50,100,160,200,400,600,800,1000,1200,1400)) {
      write.csv(browniedata, file = paste(getwd(),"/OutputFiles/",Sys.Date(),columns[1], columns[2],otherlabel, "_brownie",nsim,"sim", cladesubsetvalue,".csv",sep=""), row.names = F) #cumulative brownie data results, saved during long process
      print(paste("Saved data - Loop", i))
    }
    if (i %in% seq(0,2000,by=50)) {
      print(paste("End brownie loop iteration",i,Sys.time()))
    }
    
  }  #end for loop 1:nsim
  print("End brownie loop")
  endtimebrownie <- Sys.time()
  looptime = endtimebrownie-starttimebrownie
  
  browniedata$ERsimmapQ = gsub("ERrates ", "", Qoutput$ERrates)
  browniedata$ARDsimmapQ0to1 = gsub("ARDrates ", "", Qoutput$ARDrates[2])
  browniedata$ARDsimmapQ1to0 = gsub("ARDrates ", "", Qoutput$ARDrates[1])
  browniedata$ERsimmapQ.LogLik = Qoutput$anovaERARD$`Log lik.`[1]
  browniedata$ARDsimmapQ.LogLik = Qoutput$anovaERARD$`Log lik.`[2]
  browniedata$ARDvERsimmapQ.LRtestPval = Qoutput$anovaERARD$`Pr(>|Chi|)`[2]
  
  write.csv(browniedata, file = paste(getwd(),"/OutputFiles/",Sys.Date(),columns[1], columns[2],otherlabel,"_brownie",nsim,"sim",cladesubsetvalue,".csv",sep=""), row.names = F)
  
  return(browniedata)
}
