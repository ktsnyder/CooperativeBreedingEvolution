########
#Coded by Kate T. Snyder
#Last Modified 11/15/2023
#Built using RStudio Version 1.1.453
#R Version 4.0
#
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#BayesTraitsV2
########

# 8/18/2021 - change MateParam to DiscreteTrait, $species_in_birdtree to $species
# 9-12-2023 - outline
# 11/15/2023 - trycatch around wilcox test

#Plots merged boxplot/scatterplot for each song characteristic for each mating classification

 require(phytools)
# 
# datafile <- "2020-10-30_CoopSong_AnyCoop_NatCommsSubset.csv"
# ourdf <- read.csv(file = datafile)
# birdtree <- read.nexus("birdzillatreeMaybeConsensus.nex")
# columnnames <- "CoopBreed"
# passertree <- read.nexus("2020-10-11ConsensusPasserineTreeHack100.nex")
# treefile <- "birdzillatreeMaybeConsensus.nex"
# 
#### Example usage: ----
# scatterboxes(DiscreteTrait = "CoopBreed", newdata = datafile, newtree = passertree)


scatterboxes <- function(DiscreteTrait = "CoopBreed", newdata = FALSE, newtree = FALSE, discreteCategoryLabels = c("Non-cooperative","Cooperative"), otherlabel = NULL) {
  require(phytools)
  source(file = "subsettreedata.R")
  

  pdf(file=paste(Sys.Date(),DiscreteTrait,otherlabel,"Scatterbox.pdf",sep=""),width = 8, height = 4)
  
  
  layout(mat=matrix(1:4,nrow=1,ncol=4, byrow=TRUE))
  
  SongParams <- c("Syllable.rep.final","Syll.song.final","Song.rep.final","Duration.final","Interval.final","Song.rate","Continuity")
  SongNames <- c("Syllable repertoire","Syllables per song","Song repertoire","Song duration","Intersong interval","Song rate","Continuity")
  
  
  for (m in 1:length(SongParams)) {
    SongParam <- SongParams[m]
    subset <- subsettreedata(columns = c(DiscreteTrait,SongParam), newdata = newdata, newtree = newtree, islog = SongParam)  # this line calls a separate 
    songcol <- SongParam
    matecol <- DiscreteTrait
    subsetdf <- subset$subsetdf
    songdatavec <- subsetdf[,SongParam]
    matingdatavec <- subsetdf[,DiscreteTrait]
    names(songdatavec) <- subsetdf$species
    names(matingdatavec) <- subsetdf$species
    tree <- subset$subsettree
    print(songdatavec)
    
    set.seed(10)
    phylanovaresults <- phylANOVA(tree,matingdatavec,songdatavec,nsim=50000)
    df <- dfwilcox <- subsetdf[,c(matecol,songcol)]
    datalist <- list()
    
 #   if (DiscreteTrait == "EPP") { 
      datalist <- list() 
      datalistnames <- set.seed(10)
      df[,DiscreteTrait][which(df[,DiscreteTrait] == 0)] <- discreteCategoryLabels[1]
      df[,DiscreteTrait][which(df[,DiscreteTrait] == 1)] <- discreteCategoryLabels[2]
      factors <- discreteCategoryLabels 
      for (i in 1:2) {
        onefactordf <- df[which(df[,DiscreteTrait] == factors[i]),]
        tempfactordata <- onefactordf[,songcol]
        datalist[[i]] <- tempfactordata
        datalistnames[i] <- as.character(factors[i]) 
      } #end for i in 1:length(factors) 
      names(datalist) <- datalistnames
 #   } #end else if DiscreteTrait == "EPP"
    
      wilcoxresults <- tryCatch({wilcox.test(datalist[[1]],datalist[[2]])}, error=function(e) {
        message('An Error Occurred')
        print(e)
        return(NULL)
      })
      
      #wilcoxresults <- wilcox.test(datalist[[1]],datalist[[2]])
      
      if (is.null(wilcoxresults)) {
        wilcoxresults = NULL
        wilcoxresults$p.value = NA
      }

    labels <- c(paste(datalistnames[1]," \nN =",length(datalist[[1]])), paste(datalistnames[2]," \nN =",length(datalist[[2]])))
    main <- paste("Wilcoxon rank-sum p =",round(wilcoxresults$p.value,6) ," \nphylANOVA p =" , phylanovaresults$Pf)
    
    
    
    ScatterBox(data = datalist, main = main, ylab = SongNames[m], labels = labels)
  } # end for m in 1:length(SongParams)
  dev.off()
} # end function scatterboxes



ScatterBox <- function(data,main="",sub="",xlab="",ylab="",labels=FALSE,col1=NA,col2=rgb(0,0,0,.7),log=""){
  #make data a list, create empty plot, get quantiles and get offset
  if(is.list(data)==FALSE){
    data <- list(data)
  }
  offset <- 1:length(data)-1
  plot(1,type='n',xlim=c(0,2),ylim=c(0,max(unlist(data))),xaxt='n',
       xlab=xlab, ylab=ylab, main=main, sub=sub, cex.main = 0.8) #,log="y",cex.main = 1)
  axis(1, at=.5+offset, labels = labels, las = 1, cex.axis = 0.6)
  quan <- sapply(data, quantile, c(.25,.5,.75),type=8)
  #box
  rect(.25+offset,quan[1,],.5+offset,quan[3,],col = col1)
  #mean
  segments(.25+offset,quan[2,],.5+offset,quan[2,], lwd=2)
  
  
  #Scatter colors
  if(length(col2)!=length(data)){
    col2 <- rep(col2,length(data))
  }
  
  #get data for whiskers, scatter, and outliers
  up <- vector("numeric",length(data))
  down <- vector("numeric",length(data))
  index <- vector("list",length(data))
  jity <- vector("list",length(data))
  for(i in seq_along(data)){
    up[i] <- min(max(data[[i]]), quan[3,i] + 1.5*(quan[3,i]-quan[1,i]))
    down[i] <- max(min(data[[i]]), quan[1,i] - 1.5*(quan[3,i]-quan[1,i]))
    index <- which(data[[i]] > up[i] | data[[i]] < down[i])
    jity <- jitter(rep(.625,length(data[[i]])),8)+offset[i]
    type <- rep(20,length(data[[i]]))
    points(jity, data[[i]], pch=type, cex = 0.2, col="black")
  }
  
  #whiskerlength and block
  segments(.5+offset,c(quan[3,],quan[1,]),.5+offset,c(up,down),lty = 2)
  segments(.5+offset,c(up,down),.375+offset,c(up,down),lwd = 1)
}
