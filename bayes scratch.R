# Scratch - Bayes Discrete


simplebtwOut = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/BayesTraitsDiscrete_Hackett-Tie2Coop-FSHighConf mlt100_noRes.csv")
simplebtwOut = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/BayesTraitsDiscrete_Hackett-TieNoncoop-FSAgg mlt100_noRes.csv")
simplebtwOut = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/BayesTraitsDiscrete_Hackett-OmitTies-FSHighConf mlt100_noRes.csv")
simplebtwOut = read.csv("/Users/kate/Desktop/CooperativeBreedingEvolution/BayesTraitsDiscrete_Hackett-Tie2Noncoop-FSHighConf mlt100_noRes.csv")
columns = c("HighConfidence_FemaleSong", "MeanCoopTie2Noncoop")

nsim = length(simplebtwOut$Seed)
plotDiscreteBayes(columns = columns, simplebtwOut = simplebtwOut, nocorrDdf = NULL, newpdf = TRUE, nsim = nsim, treelabel = "Hackett")
