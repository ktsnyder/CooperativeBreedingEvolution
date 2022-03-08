########
#Coded by Kate T. Snyder
#Last Modified 6-8-2021
#Built using RStudio Version 1.1.453
#R Version 3.5.1
#
#phytools_0.6-44
#maps_3.3.0        ape_5.1 
########
#
#Jan2020 update: make compatible with new subsetbirddata function taking column vector rather than SongParam/MateParam
#V2.0 --> v2.1 - change "passeriformesonly" input to be able to do oscine, suboscine, or passeriformes  
#1/28/2020 (v2.2): #ACTUALLY change to cladesubsetcolumn, cladesubsetvalue; also remove "matemodel" arg
#5/12/2020 - remove matensim arg, change subset source from subsetbirddata2.2.R, added subdirectory to pdf output
#6/4/2021 - if plotting, checks for or makes OutputFiles subdirectory; made all text in titles of subplots be on 2nd line because margin weirdness
#8/27/2021 - add otherlabel arg to file name if plotting simmaps; added tip labels (points); added named colors for plotting simmap
#8/27/2021 - findQrates seems to calculate Q for the data subsetted by both columns, rather than just the discrete column... it should be computing Q based on whole set of discrete data! Granted, this is the case if columns input is only the discrete column... but we still want to plot simmaps of double-subsetted trees probably. Solution: add another subset within the plotting statement, make original subset only subset based on columns[1]. Done.
# 3/8/2022 - this version does not contain the setmodel parameter added in ~/Desktop/Phylobiology/findQrates.R


#findQrates - to be used within matingfunction to set the rates of transition between states for the building of simmaps for brownie
#Can also be used to generate simmaps with the computed rates with plot = TRUE

#Example:
#qoutputpass <- findQrates(columns = c("MonogamyOrNot"), plot = TRUE, cladesubsetcolumn = "oscine", cladesubsetvalue = "nonpasserine", newdata = FALSE, newtree = "2019-10-22matezilla2treeHack_nondicho.nex")
#qoutputNonpass <- findQrates(columns = c("MonogamyOrNot"), plot = TRUE, cladesubsetcolumn = "oscine", cladesubsetvalue = "Nonpasserine", newtree = "2019-10-22matezilla2treeHack_nondicho.nex")

findQrates <- function(columns, plot=FALSE, newtree = FALSE, newdata = FALSE, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL) {
require(phytools)
source(file = "subsettreedata.R")
  
    output <- list()
  

subsetoutput <- subsettreedata(columns[1], cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue, newdata = newdata, newtree=newtree)  

  tree <- subsetoutput$subsettree  
  df <- subsetoutput$subsetdf
  print(length(tree$tip.label))
  discretetraitvec <- df[,columns[1]]
  names(discretetraitvec) <- df[,1] 
  
  ##########
    ERmodel <- ace(discretetraitvec,tree, type="discrete",model = "ER")
    ARDmodel <- ace(discretetraitvec,tree, type="discrete",model = "ARD")
    anovaERARD <- anova(ERmodel,ARDmodel)
    testrates=c(ERmodel$rates,ERmodel$lik.anc[1,], ARDmodel$rates,ARDmodel$lik.anc[1,])
    output$anovaERARD <- anovaERARD
    output$ERrates <- paste("ERrates",ERmodel$rates)
    output$ERlikanc <- paste("ERlik.anc",ERmodel$lik.anc[1,])
    output$ARDrates <- paste("ARDrates",ARDmodel$rates)
    output$ARDlikanc <- paste("ARDlik.anc",ARDmodel$lik.anc[1,])
  
    
    matsize <- length(ARDmodel$rates)
    if(matsize == 1){matsize<- matsize+1}
    #create matrix of correct size; make diag negative; multiply by rates
    blankqrates <- matrix(rep(1,matsize^2),matsize,matsize)
    diag(blankqrates) <- -1
    qrates <- c(ARDmodel$rates[2],ARDmodel$rates[1])*blankqrates  #kts reversed 2/21/2020
    #pulls the state names and sets as col and row names
    rownames(qrates) <- dimnames(ARDmodel$lik.anc)[[2]]
    colnames(qrates) <- dimnames(ARDmodel$lik.anc)[[2]]
    
    if (plot == TRUE) {
      subsetoutput <- subsettreedata(columns, cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue, newdata = newdata, newtree=newtree) 
      df <- subsetoutput$subsetdf
      tree <- subsetoutput$subsettree
      # make OutputFiles folder if not present
      mainDir <- getwd()
      subDir <- "OutputFiles"
      if (!dir.exists(file.path(mainDir,subDir))) {
        dir.create(file.path(mainDir, subDir))
      }
      
      discretetraitsimmap <- make.simmap(tree,discretetraitvec,model = "ARD", nsim = 3) #makes 3 simmaps for viewing purposes
      discretetraitsimmapQset <- make.simmap(tree,discretetraitvec,model = "ARD", nsim = 3, Q = qrates)
      pdf(file = paste(mainDir,"/OutputFiles/",Sys.Date(),columns[1], columns[2],cladesubsetvalue, otherlabel, "egSimmaps.pdf",sep=""),height=12,width=6)
      layout(matrix(1:6,nrow = 2,ncol=3))
      for (i in 1:3) {
        simmap <- discretetraitsimmap[[i]]
        py = c("black","red")
        pynamed <- py
        names(pynamed) <- c(0,1)
        
        # Plot simmap 
        plotSimmap(simmap,fsize=0.2, lwd = 0.8, colors = pynamed)
        #add tips
        treetiplabels <- simmap$tip.label %in% names(discretetraitvec[discretetraitvec==1]) 
        tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
        numrates1 <- lapply(simmap[[7]],round,digits=9)
        title(paste(" ", "\nmake.simmap","Qrates:", numrates1[2],numrates1[3]),cex.main = 0.5)
        
        simmapQset <- discretetraitsimmapQset[[i]]
        plotSimmap(simmapQset,fsize=0.2,lwd=0.8, colors = pynamed)
        #add tips
        treetiplabels <- simmapQset$tip.label %in% names(discretetraitvec[discretetraitvec==1]) 
        tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
        numrates <- lapply(qrates,round,digits=9)
        title(main=paste(" ","\nARDmodel","Qrates (output for Brownie):",numrates[2],numrates[3]),cex.main = 0.5)
      }
      dev.off()
    }
    
  output$qrates <- qrates
  return(output)

}
