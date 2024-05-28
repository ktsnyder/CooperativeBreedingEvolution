# scratch find discrepancies in song data
# Vermivora_ruficapilla and Oporornis_philadelphia
# 5/9/2024

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB")

require(readxl)
require(stringr)
filelist = list.files()

outdf = set.seed(10)
for (i in 2:length(filelist)) {
  tempfile = filelist[i]
  readable = FALSE
  if (str_detect(tempfile, ".xlsx")) {
    df = read_xlsx(tempfile)
    readable = TRUE
  } else if (str_detect(tempfile, ".xls")) {
    df = read_xls(tempfile)
    readable = TRUE
  } else if (str_detect(tempfile, ".csv")) {
    df = read.csv(tempfile)
    readable = TRUE
  }
  
  if (readable) {
    if ("Syllable.rep.final" %in% colnames(df)) {
      if (!"Oporornis_philadelphia" %in% df[[1]]) {
        OphiladelphiaSyllrep = OphiladelphiaSyllsong = OphiladelphiaSongrep = "species not present"
      
      } else {
        OphiladelphiaSyllrep = df[which(df[,1] == "Oporornis_philadelphia"),"Syllable.rep.final"][[1]]
        OphiladelphiaSyllsong = df[which(df[,1] == "Oporornis_philadelphia"),"Syll.song.final"][[1]]
        OphiladelphiaSongrep = df[which(df[,1] == "Oporornis_philadelphia"),"Song.rep.final"][[1]]
      }
      
      if (!"Vermivora_ruficapilla" %in% df[[1]]) {
        VruficapillaSyllrep = VruficapillaSyllsong = VruficapillaSongrep = "species not present"
      } else {
        VruficapillaSyllrep = df[which(df[,1] == "Vermivora_ruficapilla"),"Syllable.rep.final"][[1]]
        VruficapillaSyllsong = df[which(df[,1] == "Vermivora_ruficapilla"),"Syll.song.final"][[1]]
        VruficapillaSongrep = df[which(df[,1] == "Vermivora_ruficapilla"),"Song.rep.final"][[1]]
      }
      temprow = c(tempfile, OphiladelphiaSyllrep, OphiladelphiaSyllsong, OphiladelphiaSongrep, VruficapillaSyllrep,VruficapillaSyllsong, VruficapillaSongrep )
  
    } else {
      temprow = c(tempfile, "column not present", "column not present", "column not present","column not present","column not present","column not present")
     # temprow = as.data.frame(temprow)
     # 
    }
    names(temprow) = c("tempfile", "OphiladelphiaSyllrep", "OphiladelphiaSyllsong", "OphiladelphiaSongrep", "VruficapillaSyllrep", "VruficapillaSyllsong", "VruficapillaSongrep" )
    temprow = c(temprow, sum(is.na(temprow)))
    outdf = rbind(outdf, temprow)
  }
}
outdf = as.data.frame(outdf)
write.csv(outdf,"Check song data for V_ruficapilla and O_philadelphia4.csv")


"Oporornis_philadelphia"
"Vermivora_ruficapilla"
Rdata = read.csv("SongData_R_Update.csv")
Refsdata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/SupplementDataRefs_Update.csv")
Refsdata = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/SupplementDataRefs_Update_2024-05-13.csv")
which(Refsdata$BirdtreeSpecies == "Vermivora_ruficapilla")
Refsdata$BirdtreeSpecies = Refsdata$BirdtreeFormat
Refsdata = Refsdata[which(Refsdata$BirdtreeFormat != ""),]
Refsdata[870,]
Rdatasubset = Rdata[which(Rdata$BirdtreeFormat %in% Refsdata$BirdtreeSpecies),]
Refsdata$Syllable.rep.final == Rdatasubset$Syllable.rep.final
Rdatasubset$Syllable.rep.final == Refsdata$Syllable.rep.final
sum(is.na(Rdatasubset$Syllable.rep.final))
sum(is.na(Refsdata$Syllable.rep.final))

are_different <- !((Rdatasubset$Syllable.rep.final == Refsdata$Syllable.rep.final) | (is.na(Rdatasubset$Syllable.rep.final) & is.na(Refsdata$Syllable.rep.final)))
sum(are_different, na.rm = T)

colSums(!is.na(Refsdata[,c("Syllable.rep.final","Syll.song.final", "Song.rep.final", "Interval.final", "Duration.final")]))
colSums(!is.na(Rdata[,c("Syllable.rep.final","Syll.song.final", "Song.rep.final", "Interval.final", "Duration.final")]))


Rdatasubset$BirdtreeFormat == Refsdata$BirdtreeSpecies

which(is.na(Rdatasubset$Syllable.rep.final) & !is.na(Refsdata$Syllable.rep.final))
Rdatasubset[c(522,870),]
Refsdata[c(522,870),]
which(is.na(Rdatasubset$Syll.song.final) & !is.na(Refsdata$Syll.song.final))
Rdatasubset[c(233),]
Refsdata[c(233),]
which(is.na(Rdatasubset$Song.rep.final) & !is.na(Refsdata$Song.rep.final)) #none
which(is.na(Rdatasubset$Duration.final) & !is.na(Refsdata$Duration.final)) #none
which(is.na(Rdatasubset$Interval.final) & !is.na(Refsdata$Interval.final)) #none


sum(Rdatasubset$Syll.song.final == Refsdata$Syll.song.final, na.rm = T)
sum(Rdatasubset$Syllable.rep.final == Refsdata$Syllable.rep.final, na.rm = T)
sum(Rdatasubset$Song.rep.final == Refsdata$Song.rep.final, na.rm = T)
sum(Rdatasubset$Interval.final == Refsdata$Interval.final, na.rm = T)
sum(Rdatasubset$Duration.final == Refsdata$Duration.final, na.rm = T)

which(Rdatasubset$Syllable.rep.final != Refsdata$Syllable.rep.final)
which(Rdatasubset$Syll.song.final != Refsdata$Syll.song.final)
which(Rdatasubset$Song.rep.final != Refsdata$Song.rep.final)
which(Rdatasubset$Interval.final != Refsdata$Interval.final)
Rdatasubset[c(16,18,320),]
Refsdata[c(16,18,320),]
which(Rdatasubset$Duration.final != Refsdata$Duration.final)
Rdatasubset[c(237,266),]
Refsdata[c(237,266),]


# after all data re-compiled
AggNew = read.csv("2024-05-13_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv")
AggOld = read.csv("2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv")

AggNew = read.csv("2024-05-13_Aggregate_CBSource_Data_AllColumns.csv")
AggOld = read.csv("2024-02-24_Aggregate_CBSource_Data_AllColumns.csv")


sum(AggNew$species == AggOld$species)
sum(AggNew$HighConfidence_Coop == AggOld$HighConfidence_Coop, na.rm = T)
sum(!is.na(AggNew$HighConfidence_Coop))
sum(!is.na(AggOld$HighConfidence_Coop))
sum(AggNew$HighConfidence_FemaleSong == AggOld$HighConfidence_FemaleSong, na.rm = T)
sum(!is.na(AggNew$HighConfidence_FemaleSong))
sum(!is.na(AggOld$HighConfidence_FemaleSong))
sum(AggNew$FemaleSong_Agg01 == AggOld$FemaleSong_Agg01, na.rm = T)
sum(!is.na(AggNew$FemaleSong_Agg01))
sum(!is.na(AggOld$FemaleSong_Agg01))
sum(AggNew$MeanCoopTie2Noncoop == AggOld$MeanCoopTie2Noncoop, na.rm = T)
sum(!is.na(AggNew$MeanCoopTie2Noncoop))
sum(!is.na(AggOld$MeanCoopTie2Noncoop))
sum(AggNew$Syll.song.final == AggOld$Syll.song.final, na.rm = T)
sum(!is.na(AggNew$Syll.song.final))
sum(!is.na(AggOld$Syll.song.final))
sum(AggNew$Syllable.rep.final == AggOld$Syllable.rep.final, na.rm = T)
sum(!is.na(AggNew$Syllable.rep.final))
sum(!is.na(AggOld$Syllable.rep.final))
sum(AggNew$Song.rep.final == AggOld$Song.rep.final, na.rm = T)
sum(!is.na(AggNew$Song.rep.final))
sum(!is.na(AggOld$Song.rep.final))
sum(AggNew$Duration.final == AggOld$Duration.final, na.rm = T)
sum(!is.na(AggNew$Duration.final))
sum(!is.na(AggOld$Duration.final))
sum(AggNew$Interval.final == AggOld$Interval.final, na.rm = T)
sum(!is.na(AggNew$Interval.final))
sum(!is.na(AggOld$Interval.final))

### redo brownie ----

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution/")
source("subsettreedata.R")
source("browniefunction.R")
source("plotbrownie.R")

songfeatures <- c("Song.rep.final", "Syllable.rep.final", "Syll.song.final", "Song.rep.max")
CBcolumn <- "HighConfidence_Coop"
discreteCatLabels = c("Non-cooperative", "Cooperative")
nsim = 500
currentlabel <- "AddVruficapillaAndOphiladelphia"

newdata = "2024-02-24_CoopBreed-FemaleSong-Song-Sociality01_PasseriformesData_HighConfCoopCol_R.csv"
treefile = "2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
df = read.csv(newdata)
# 5/13/2024
df$Syllable.rep.final[which(df$species == "Oporornis_philadelphia")] = df$Syll.song.final[which(df$species == "Oporornis_philadelphia")]
df$Syllable.rep.final[which(df$species == "Vermivora_ruficapilla")] = 2.7
df$Syll.song.final[which(df$species == "Vermivora_ruficapilla")] = 2.7
df$Song.rep.final[which(df$species == "Oporornis_philadelphia")] = 1.5
df$Song.rep.max[which(df$species == "Oporornis_philadelphia")] = 2

for (k in 1:length(songfeatures)) {
  print(Sys.time())
  feature <- songfeatures[k]
  print(feature)
  browniefunction(columns = c(CBcolumn, feature), newdata = df, newtree = treefile, nsim = nsim, islog = feature, plotsimmaps = TRUE, otherlabel = currentlabel)
  plotbrownie(data = paste0(Sys.Date(),CBcolumn,feature, currentlabel, "_brownie",nsim,"sim.csv"), columns = c(CBcolumn,feature), discreteCategoryLabels = discreteCatLabels, otherlabel = currentlabel, newpdf = TRUE, nsim = nsim, islog = TRUE)
} # end for k



df = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/2024-05-13HighConfidence_CoopSong.rep.finalUpdatedSongData_brownie500sim.csv")
dfold = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/Brownie/2024-02-25HighConfidence_CoopSong.rep.finalHackettOscine_brownie500sim.csv")
sum(df$Pval < 0.05)/500
df$ARDRate1 == dfold$ARDRate1
df$ARDRate0 == dfold$ARDRate0
df$ERRate == dfold$ERRate

df = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/2024-05-13HighConfidence_CoopSyllable.rep.finalUpdatedSongData_brownie500sim.csv")
dfold = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/OutputFiles/Brownie/2024-02-25HighConfidence_CoopSyllable.rep.finalHackettOscine_brownie500sim.csv")
sum(df$Pval < 0.05)/500
sum(dfold$Pval < 0.05)/500
df$ARDRate1 == dfold$ARDRate1
df$ARDRate0 == dfold$ARDRate0
df$ERRate == dfold$ERRate

head(df)
head(dfold)
