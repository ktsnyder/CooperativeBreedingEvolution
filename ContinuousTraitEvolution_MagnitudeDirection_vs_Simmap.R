## Get direction and magnitude of continuous character change
## 


nsims = 100

newdata = "2022-03-10CoopSong_All.csv"
df <- read.csv(newdata)
columns <- c("MeanCoopTie2Noncoop", "Syll.song.final")

source("findQrates.R")
cooprates <- findQrates(columns = columns[1], newdata = df, newtree = tree)
coopQ <- cooprates$qrates

subsetout <- subsettreedata(columns = columns, newdata = df, newtree = HackConsensus, islog = columns[2])
subsetdf <- subsetout$subsetdf
subsettree <- subsetout$subsettree
contvec <- subsetdf[,columns[2]]
names(contvec) <- subsetdf$species
ancOut<- fastAnc(subsettree, contvec)
ancOutdf <- as.data.frame(cbind(names(ancOut), ancOut))
colnames(ancOutdf) <- c("label", "value")


tipnum <- 1:length(subsettree$tip.label)
tipsdf <- as.data.frame(cbind(tipnum, subsettree$tip.label))
subsetdftips <- merge(subsetdf, tipsdf, by.x = "species", by.y = "V2")
subsetdftips <- subsetdftips[,c("tipnum",columns[2])]
colnames(subsetdftips) <- c("label", "value")

contdf <- rbind(ancOutdf, subsetdftips)

xDisc <- subsetdf[,columns[1]]
names(xDisc) <- subsetdf$species

simmap <- make.simmap(tree = subsettree, x = xDisc, model = "ARD", nsim = nsims)
plotSimmap(simmap, fsize = 0.8)
  # tip label numbers are in order of simmap$tip.labels
  # tip.label is the same (and in same order) for both simmap and subsettree
  # in simmap$mapped.edge, the second number corresponds to tip.label number 
  # ancOut provides values at nodes only, would need to add in tip values
head(simmap$mapped.edge)
class(rownames(simmap$mapped.edge))

## Next: make below into a loop through each simmap in simmap
require(stringr)
require(tidyr)
nodes <- str_split(rownames(simmap$mapped.edge), ",")
simmap_edges <- as.data.frame(simmap$mapped.edge)
simmap_edges <- cbind(rownames(simmap_edges),simmap_edges)
colnames(simmap_edges) <- c("nodes","timeIn0","timeIn1")
simmap_edgesdf <- separate(simmap_edges, "nodes", into = c("node1", "node2"), sep = ",")

simmap_edgesdf <- merge(simmap_edgesdf, contdf, by.x = "node1", by.y = "label")
simmap_edgesdf <- merge(simmap_edgesdf, contdf, by.x = "node2", by.y = "label", suffixes = c("_node1","_node2"), sort = FALSE)

delta_contval <- as.numeric(simmap_edgesdf$value_node1) - as.numeric(simmap_edgesdf$value_node2)
edgesdf <- cbind(simmap_edgesdf, delta_contval)

changeDuring0 <- edgesdf$delta_contval*edgesdf$timeIn0
changeDuring1 <- edgesdf$delta_contval*edgesdf$timeIn1
edgesdf <- cbind(edgesdf, changeDuring0, changeDuring1)

netChangeDuring0 <- sum(changeDuring0)
totalTimeIn0 <- sum(edgesdf$timeIn0)
MeanChangeDuring0 <- netChangeDuring0/totalTimeIn0
netChangeDuring1 <- sum(changeDuring1)
totalTimeIn1 <- sum(edgesdf$timeIn1)
MeanChangeDuring1 <- netChangeDuring1/totalTimeIn1

totalPositiveChangeDuring0 <- sum(changeDuring0[which(changeDuring0 > 0)])
totalNegativeChangeDuring0 <- sum(changeDuring0[which(changeDuring0 < 0)])
totalPositiveChangeDuring1 <- sum(changeDuring1[which(changeDuring1 > 0)])
totalNegativeChangeDuring1 <- sum(changeDuring1[which(changeDuring1 < 0)])

row <- c(totalPositiveChangeDuring0, totalNegativeChangeDuring0, netChangeDuring0, totalTimeIn0, MeanChangeDuring0, totalPositiveChangeDuring1, totalNegativeChangeDuring1, netChangeDuring1, totalTimeIn1, MeanChangeDuring1)
