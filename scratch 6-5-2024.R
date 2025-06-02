library(clootl)

ex1 <- extractTree(species=c("amerob", "canwar", "reevir1", "yerwar", "gockin"),
                   output.type="code", taxonomy.year="current", version="current")


trees = dataStore$trees

trees10.2021 = trees$Aves_1.2$tree.samples$year2021
trees10.2021[[1]]
trees10.2022 = trees$Aves_1.2$tree.samples$year2022
trees10.2022[[1]]
trees10.2023 = trees$Aves_1.2$tree.samples$year2023
trees10.2023[[1]]


setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")
data = read.csv("2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv")
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
OscineBirdTree = read.nexus(treefile)

dataSub = data[which(data$species %in% OscineBirdTree$tip.label),]
dataSub %>% group_by(HighConfidence_Coop) %>% count

sum(OscineBirdTree$tip.label %in% trees10.2023[[1]]$tip.label)
length(OscineBirdTree$tip.label)

dataSub = data[which(!is.na(data$HighConfidence_Coop) & !is.na(data$FemaleSong_Agg01)),]
sum(dataSub$species %in% trees10.2023[[1]]$tip.label)


#newdata = "2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
newdata = "2024-06-06_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

socialityMetrics = c("BiagoliniCoop", "DowningCoop", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop") # "JetzCoopInclCockburn"

treelabel = "HackettOscine_GlobalQrates"
datalabel = NULL
nsims_real = 500
nsims_dummy = 500

findQout = findQrates(columns = "HighConfidence_Coop", newtree = treefile, newdata = newdata)
GlobalQrates = findQout$qrates

source("test_trait_overlap_simmaps.R")
socialityPlots = list()
socialityStats = list()
mediansdfGlobalQ = set.seed(10)
for (i in 1: length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(i)
  print(tempMetric)
  print("Starting Real")
  dfout4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = treelabel, datalabel = datalabel, plotSampleSimmaps = FALSE, columnForGlobalQ = 1, columnGlobalQrates = GlobalQrates)
  print("Starting Dummy")
  dfDummy4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = treelabel, datalabel = datalabel, dummyMethod = "makeSimmap", columnForGlobalQ = 1, columnGlobalQrates = GlobalQrates)
  
  print(paste("Starting calcHuel", tempMetric))
  calcHuelout = calcHuel(dfout4, dfDummy4)
  require(gridExtra)
  plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts.pdf")

  calcHuelout2 = calcHuelout[c("p1", "p2", "p3", "p4","p5","p6")]
  nPlots = 6
  m3 <- marrangeGrob(calcHuelout2, ncol = 1, nrow = nPlots)
  ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
  #socialityPlots[[i]] <- calcHuelout
  temprow = calcHuelout$mediansRow
  mediansdfGlobalQ = rbind(mediansdfGlobalQ, temprow)
}
write.csv(mediansdfGlobalQ, "SingleSourceCB FemaleSong_Agg01 simmap overlap medians_GlobalQ.csv")



treelabel = "HackettOscine_LocalQrates"
datalabel = NULL
nsims_real = 500
nsims_dummy = 500
source("test_trait_overlap_simmaps.R")
socialityPlots = list()
socialityStats = list()
mediansdfLocalQ = set.seed(10)
for (i in 1: length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(i)
  print(tempMetric)
  print("Starting Real")
  dfout4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = treelabel, datalabel = datalabel, plotSampleSimmaps = FALSE)
  print("Starting Dummy")
  dfDummy4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = treelabel, datalabel = datalabel, dummyMethod = "makeSimmap")
  
  print(paste("Starting calcHuel", tempMetric))
  calcHuelout = calcHuel(dfout4, dfDummy4)
  require(gridExtra)
  plotname = paste0("Simmap Overlap Outputs/",tempMetric, " FemaleSong_Agg01 ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts.pdf")
  #nPlots = length(calcHuelout)-5
  #m3 <- marrangeGrob(calcHuelout, ncol = 1, nrow = nPlots)
  #ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
  calcHuelout2 = calcHuelout[c("p1", "p2", "p3", "p4","p5","p6")]
  nPlots = 6
  m3 <- marrangeGrob(calcHuelout2, ncol = 1, nrow = nPlots)
  ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
  
  calcHuelout$TransitionStats
  
  temprow = calcHuelout$mediansRow
  mediansdfLocalQ = rbind(mediansdfLocalQ, temprow)
  #socialityPlots[[i]] <- calcHuelout
}
write.csv(mediansdfLocalQ, "SingleSourceCB FemaleSong_Agg01 simmap overlap medians_LocalQ.csv")

JetzLocalDummy = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ JetzCoopInclCockburn FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 500 HackettOscine_LocalQrates .csv")
JetzLocalReal = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ JetzCoopInclCockburn FemaleSong_Agg01 REAL simmap overlap_counts output nsim 500 HackettOscine_LocalQrates .csv")
calcHuelout = calcHuel(JetzLocalReal, JetzLocalDummy)

JetzGlobalDummy = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ JetzCoopInclCockburn FemaleSong_Agg01 DUMMYResampledMkSimmap-CoopFS simmap overlap_counts output nsim 500 HackettOscine_GlobalQrates .csv")
JetzGlobalReal = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Simmap Overlap Outputs/ JetzCoopInclCockburn FemaleSong_Agg01 REAL simmap overlap_counts output nsim 500 HackettOscine_GlobalQrates .csv")
calcHuelout = calcHuel(JetzGlobalReal, JetzGlobalDummy)

socialityMetrics = c("BiagoliniCoop", "DowningCoop", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop", "JetzCoopInclCockburn")
filelist = list.files("Simmap Overlap Outputs", pattern = "Qrates .csv")

mediansdfout = set.seed(10)
for (i in 1:length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  metricfiles = filelist[which(str_detect(filelist, tempMetric))]
  LocalFiles = metricfiles[which(str_detect(metricfiles, "Local"))]
  dummyFile = LocalFiles[which(str_detect(LocalFiles, "DUMMY"))]
  realFile = LocalFiles[which(str_detect(LocalFiles, "REAL"))]
  dummyLocal = read.csv(paste0("Simmap Overlap Outputs/", dummyFile))
  realLocal = read.csv(paste0("Simmap Overlap Outputs/", realFile))
  calcHueloutLocal = calcHuel(realLocal, dummyLocal, otherlabel = "LocalQ")
  Localrow = calcHueloutLocal$mediansRow
  Qrates = "Local"
  Localrow = cbind(Localrow, Qrates, realFile, dummyFile)
  print(Localrow)
  mediansdfout = rbind(mediansdfout, Localrow)
  
  GlobalFiles = metricfiles[which(str_detect(metricfiles, "Global"))]
  dummyFile = GlobalFiles[which(str_detect(GlobalFiles, "DUMMY"))]
  realFile = GlobalFiles[which(str_detect(GlobalFiles, "REAL"))]
  dummyGlobal = read.csv(paste0("Simmap Overlap Outputs/", dummyFile))
  realGlobal = read.csv(paste0("Simmap Overlap Outputs/", realFile))
  calcHueloutGlobal = calcHuel(realGlobal, dummyGlobal, otherlabel = "GlobalQ")
  Globalrow = calcHueloutGlobal$mediansRow
  Qrates = "Global"
  Globalrow = cbind(Globalrow, Qrates, realFile, dummyFile)
  print(Globalrow)
  mediansdfout = rbind(mediansdfout, Globalrow)
  
}
write.csv(mediansdfout, "single source simmap overlap medians summary LocalVsGlobalQ.csv")



#### High Conf Coop and High Conf FS
tempMetric = "HighConfidence_Coop"
treelabel = "HackettOscine_LocalQrates"
nsims_real = 500
nsims_dummy = 500
dfout4 <- CharacterSimmaps(columns = c(tempMetric,"HighConfidence_FemaleSong"), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = treelabel, datalabel = NULL, plotSampleSimmaps = FALSE)
print("Starting Dummy")
dfDummy4 <- CharacterSimmaps(columns = c(tempMetric,"HighConfidence_FemaleSong"), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = treelabel, datalabel = NULL, dummyMethod = "makeSimmap")

print(paste("Starting calcHuel", tempMetric))
calcHuelout = calcHuel(dfout4, dfDummy4)
require(gridExtra)
plotname = paste0("Simmap Overlap Outputs/",tempMetric, " HighConfidence_FemaleSong ", nsims_real, " ", nsims_dummy, " ", treelabel, " withTransCounts.pdf")

calcHuelout2 = calcHuelout[c("p1", "p2", "p3", "p4","p5","p6")]
nPlots = 6
m3 <- marrangeGrob(calcHuelout2, ncol = 1, nrow = nPlots)
ggsave(plotname, m3, width = 7.5, height = 3.8*nPlots, units = "in")
#socialityPlots[[i]] <- calcHuelout
temprow = calcHuelout$mediansRow
mediansdfGlobalQ = rbind(mediansdfGlobalQ, temprow)



data0607 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/2024-06-07_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv")
data0606 = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/2024-06-06_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv")
all.equal(data0606, data0607)

data = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Sandbox/Data_R.csv")
all.equal(data0607, data) #TRUE
sum(is.na(data0606$AnyNoncoopEqualsNoncoop)) # 0 = bad
