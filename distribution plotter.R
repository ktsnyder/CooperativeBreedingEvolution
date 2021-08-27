# Adapted from /Creanza Lab/Mating Systems/OC 2019 Misc/distribution plotter_kate (2020_03_04 22_03_48 UTC).R

# DistributionPlots("OC","Syllsong")

#alldatadf <- as.data.frame(read.csv("~/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/2021-08-18CoopSong_MeanCoop_NatCommsSubset.csv"))
#nodups <-alldatadf[which(duplicated(alldatadf$Syllable.rep.final)==FALSE),]
#nodups <- alldatadf #use when want all values
#alldatadf <- nodups

DistributionPlots <- function(MateParam=c("Polygyny","EPP", "CoopBreed"),SongParam=c("Song","Syllsong","Syllrep","Interval","Duration","Rate","Continuity")){
  #pull data, combine, kill dups, sort, separate out closed-ended learners
  
  # if (MateParam == "Polygyny") {
  #   nodups <- alldatadf[!is.na(alldatadf$Final.polygyny),]
  # } else if (MateParam == "EPP") {
  #   nodups <- alldatadf[!is.na(alldatadf$Final.EPP),] 
  # } else if (MateParam=="none") {
  #   nodups <- alldatadf
  # } else {
  #   return(print("Not a valid value for MateParam"))}
  # output$lengthnodups <- length(nodups$BirdtreeFormat)
  nodups <- alldatadf
  if (SongParam == "Syllrep") {
    songdf <- nodups[!is.na(nodups$Syllable.rep.final),]
    nodups <- songdf[order(songdf$Syllable.rep.final),]   #sort the reps smallest to largest
    maxsong <- max(songdf$Syllable.rep.final)
    dfcol <- "Syllable.rep.final"
  } else if (SongParam == "Syllsong") {
    songdf <- nodups[!is.na(nodups$Syll.song.final),]
    nodups <- songdf[order(songdf$Syll.song.final),]
    maxsong <- max(songdf$Syll.song.final)
    dfcol <- "Syll.song.final"
  } else if (SongParam == "Song") {
    songdf <- nodups[!is.na(nodups$Song.rep.final),]
    nodups <- songdf[order(songdf$Song.rep.final),]
    maxsong <- max(songdf$Song.rep.final)
    dfcol <- "Song.rep.final"
  } else if (SongParam == "Interval") {
    songdf <- nodups[!is.na(nodups$Interval.final),]
    nodups <- songdf[order(songdf$Interval.final),]
    maxsong <- max(songdf$Interval.final)
    dfcol <- "Interval.final"
  } else if (SongParam == "Duration") {
    songdf <- nodups[!is.na(nodups$Duration.final),]
    nodups <- songdf[order(songdf$Duration.final),]
    maxsong <- max(songdf$Duration.final)
    dfcol <- "Duration.final"
  } else if (SongParam == "Rate") {
    songdf <- nodups[!is.na(nodups$Song.rate),]
    nodups <- songdf[order(songdf$Song.rate),]
    maxsong <- max(songdf$Song.rate)
    dfcol <- "Song.rate"
  } else if (SongParam == "Continuity") {
    songdf <- nodups[!is.na(nodups$Continuity),]
    nodups <- songdf[order(songdf$Continuity),]
    maxsong <- max(songdf$Continuity)
    dfcol <- "Continuity"
  } else if (SongParam == "none") {
    songdf <- nodups
  } else {return(print("Not a valid value for SongParam"))}

  #logic = plot everything in one color and then plot over that a subset in a different color
  #this part pulls out that subset and then gets their index for positional things
if (MateParam == "CoopBreed") {
  opens <- which(nodups$CoopBreed==1)
  closeds <- which(nodups$CoopBreed==0)
  nahs <- which(is.na(nodups$CoopBreed))
  polygyny <-nodups[opens,]
  monogamy <- nodups[closeds,]
}

  #stuff for graphing
  ticks <- length(nodups$species)
#k  loc <- which(nodups$Syllable.rep.final == 38)-.5 #this was for my dotted line.  I don't think you need it
  
  #I made this to got into a set of four.  I assume you are not doing that, so just remove the
  #new=T and thd fig=blah stuff and you are good to go for a side-by-side of lin and log
#  par(mfrow=c(1,2),mar = c(5,3.5,.5,1))
  pdf(paste(Sys.Date(),MateParam,SongParam,"_dist2.pdf", sep=""),width=12,height=4) 
  layout(cbind(1,2))
  par(mar=c(4,4,1,1))
  
  #legend position info
  legtop <- maxsong - .025*maxsong
  
  #plot for linear graph
  #plot everything NOTHING; make the window (there is some extra garbage in this line of code)
  plot(0, type="n", xlim = c(1,ticks), ylim = c(0, maxsong),
       las = 2, xaxt = "n", xlab = "", ylab = dfcol, pch=16, col = "red",  cex = .5,
       cex.axis = .75
       #the below line is again for my dotted line with you prolly don't need
#k       panel.first = {segments(loc, 0, loc, max(nodups$Syllable.rep.final), lty = 2, lwd = 1.5, col = "grey60")}
)
  #connects your lines in BLACK
  segments(2:ticks-1, nodups[,dfcol][2:ticks-1], 2:ticks, nodups[,dfcol][2:ticks],
           lty = 1, lwd = 1.5, col = "black")
  #plot all points in RED
  points(1:ticks, nodups[,dfcol], pch=16, col = "black", cex=0.7)
  #plot your subste in BLUE
  points(closeds, monogamy[,dfcol], pch=16, col = "blue", cex=0.7)
  points(opens, polygyny[,dfcol], pch=16, col = "red", cex=0.7)

  
  #self-explanatory <3
  axis(1, 1:ticks, labels = nodups$species, cex.axis = .3, las=2, font = 3)
  legend(1, legtop, legend=c("Non-cooperative", "Cooperative"),col=c("blue", "red"), pch = 16, cex=.8)
  nodups[,dfcol] <- log(nodups[,dfcol])
  polygyny[,dfcol] <- log(polygyny[,dfcol])
  monogamy[,dfcol] <- log(monogamy[,dfcol])
  plot(0, type="n", xlim = c(1,ticks), ylim = c(0, log(maxsong)),
       las = 2, xaxt = "n", xlab = "", ylab =paste("ln",dfcol), pch=16, col = "red",  cex = .5,cex.axis = .75)

  points(1:ticks, nodups[,dfcol], pch=16, col = "black", cex=0.7)
  segments(2:ticks-1, nodups[,dfcol][2:ticks-1], 2:ticks, nodups[,dfcol][2:ticks], lty = 1, lwd = 1.5, col = "black")
  #plot your subste in BLUE
  points(closeds, monogamy[,dfcol], pch=16, col = "blue", cex=0.7)
  points(opens, polygyny[,dfcol], pch=16, col = "red", cex=0.7)
  axis(1, 1:ticks, labels = nodups$species, cex.axis = .3, las=2, font  = 3)
  legtop <- maxsong - .025*maxsong
  legend(1, legtop, legend=c("Non-cooperative", "Cooperative"),
         col=c("blue", "red"), pch = 16, cex=.8)
 # grid(10,NA)

  dev.off()
}
