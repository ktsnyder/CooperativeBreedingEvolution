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


newdata = "2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

socialityMetrics = c("BiagoliniCoop", "DowningCoop", "JetzCoop", "CockburnCoop" , "Griesser2017Coop" , "DaleCoop", "CornwallisCoop")

treelabel = "HackettOscine_GlobalQrates"
nsims_real = 500
nsims_dummy = 500

findQout = findQrates(columns = "HighConfidence_Coop", newtree = treefile, newdata = newdata)
GlobalQrates = findQout$qrates

source("test_trait_overlap_simmaps.R")
socialityPlots = list()
socialityStats = list()
for (i in 2: length(socialityMetrics)) {
  tempMetric = socialityMetrics[i]
  print(i)
  print(tempMetric)
  dfout4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = FALSE, nsims = nsims_real, treelabel = "HackettOscine", datalabel = datalabel, plotSampleSimmaps = TRUE, columnForGlobalQ = 1, columnGlobalQrates = GlobalQrates)
  dfDummy4 <- CharacterSimmaps(columns = c(tempMetric,"FemaleSong_Agg01"), df = newdata, tree =  treefile, dummy = TRUE, nsims = nsims_dummy, treelabel = "HackettOscine", datalabel = datalabel, dummyMethod = "makeSimmap", columnForGlobalQ = 1, columnGlobalQrates = GlobalQrates)
  
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
  #socialityPlots[[i]] <- calcHuelout
}
