# PhyloPath scratch

library(phylopath)

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")

source("subsettreedata.R")

newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"

tree = read.nexus(treefile)
dfIn = read.csv(newdata)
rownames(dfIn) = dfIn$species

dfIn[which(dfIn$species == "Vermivora_ruficapilla"),]

phylo_path()


df = read.csv("/Users/kate/Downloads/M-FsongOrdinal&AllNatHist_ForFinalAnalysis_v4_ForPublication.csv")
unique(df$FemSongFinal)
sapply(df, unique)

df %>% group_by(FemSongFinal, prs_abs) %>% count

df %>% group_by(DimorphSongElb_ord, prs_abs) %>% count
1309-223

df %>% group_by(DimorphSongOcc_ord, prs_abs) %>% count
1309-273

df %>% group_by(DimorphSongLen_ord, prs_abs) %>% count
1309-324


