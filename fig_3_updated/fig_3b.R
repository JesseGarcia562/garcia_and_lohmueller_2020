library(tidyverse)
library(data.table)
library(glue)

prepMatchedDataForPermutation <- function(MatchedLDData){
  MatchedLDData$Variation <- as.character(MatchedLDData$Variation)
  MatchedLDData$Variation[MatchedLDData$Variation == "=nonsynonymous_SNV"] <- "Nonsynonymous SNV"
  MatchedLDData$Variation[MatchedLDData$Variation == "=synonymous_SNV"] <- "Synonymous SNV"
  return(MatchedLDData)
}

permuteLDStat <- function(preppedMatchedData, SummaryStat){
  nonsynonymousStats <- preppedMatchedData[[SummaryStat]][preppedMatchedData$Variation == "Nonsynonymous SNV" ]
  synonymousStats <- preppedMatchedData[[SummaryStat]][preppedMatchedData$Variation == "Synonymous SNV" ]
  permute(nonsynonymousStats=nonsynonymousStats, synonymousStats=synonymousStats, statistic = SummaryStat )
  
}

plottingDPermutations <- function(MatchedLDData, LDStat="D"){
  preppedMatchedData <- prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + theme(legend.position="none")
}


oneTest <- function(nonsynonymousStats, synonymousStats){
  differences <- nonsynonymousStats - synonymousStats
  #For Randomization test on matched samples you need to randomly change the sign of the differences. 
  #Sample 1 and -1, multiply by differences
  signs <- sample(size=length(differences), c(1,-1), replace = TRUE)
  newsample <- signs*differences
  randomizedMean <- mean(newsample)
  return(randomizedMean)
}

permute<-function(nonsynonymousStats,synonymousStats, replications = 10000, statistic="Some LD Statistic" ){
  actualDiff<-mean(nonsynonymousStats - synonymousStats)
  permutations<-replicate(replications, oneTest(nonsynonymousStats, synonymousStats))
  if (statistic == "rSquare"){
    title <- bquote( ~ R^2)
  } else if (statistic == "DPrime") {
    title <- bquote( ~ "D' ")
  } else if (statistic == "D"){
    title <- bquote(~ "D ")
  } else if (statistic == "FreqAB"){
    title <- bquote(~ "Count of AB Haplotypes")
  }
  permDF<- as.data.frame(permutations)
  names(permDF)<-"MeanDif"
  ggplot(permDF, aes(x=MeanDif)) + 
    geom_histogram() + 
    labs(title=title, x="Mean Difference \n Between Matched Pairs", y="Frequency") + 
    geom_vline(aes(xintercept = actualDiff, color = "Observed Mean Difference"), linetype=2, size=2) + 
    theme_bw()
}


prepMatchedDataForPermutation<-function(MatchedLDData){
  MatchedLDData$Variation <- as.character(MatchedLDData$Variation)
  MatchedLDData$Variation[MatchedLDData$Variation == "=nonsynonymous_SNV"] <- "Nonsynonymous SNV"
  MatchedLDData$Variation[MatchedLDData$Variation == "=synonymous_SNV"] <- "Synonymous SNV"
  return(MatchedLDData)
}

permuteLDStat<-function(preppedMatchedData, SummaryStat){
  nonsynonymousStats <- preppedMatchedData[[SummaryStat]][preppedMatchedData$Variation == "Nonsynonymous SNV" ]
  synonymousStats <- preppedMatchedData[[SummaryStat]][preppedMatchedData$Variation == "Synonymous SNV" ]
  permute(nonsynonymousStats=nonsynonymousStats, synonymousStats=synonymousStats, statistic = SummaryStat )
  
}

plottingDPermutations <- function(MatchedLDData, LDStat="D"){
  preppedMatchedData <- prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + theme(legend.position="none")
}

plottingDPrimePermutations <- function(MatchedLDData, LDStat="DPrime"){
  preppedMatchedData <- prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + labs(color="") + theme(legend.position = c(0.89,.8)) + theme(legend.position="none")
  
}

plottingRSquarePermutations <- function(MatchedLDData, LDStat="rSquare"){
  preppedMatchedData <- prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + theme(legend.position="none")
}

plottingFreqABPermutations <- function(MatchedLDData, LDStat="FreqAB"){
  preppedMatchedData <- prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + theme(legend.position="none")
}



plotPermutations3<-function(MatchedLDData, population){
  DPermGraph <- plottingDPermutations(MatchedLDData = MatchedLDData) +theme(text = element_text(size=14)) +   scale_x_continuous(labels = scales::scientific, breaks = scales::pretty_breaks(n=3)) + geom_histogram( fill="slategray3") + xlab(NULL)
  DPrimePermGraph <- plottingDPrimePermutations(MatchedLDData = MatchedLDData) +theme(text = element_text(size=14))+   scale_x_continuous(labels = scales::scientific, breaks = scales::pretty_breaks(n=3))  + geom_histogram( fill="slategray3") + xlab(NULL)
  RSquarePermGraph <- plottingRSquarePermutations(MatchedLDData = MatchedLDData)+theme(text = element_text(size=14))+   scale_x_continuous(labels = scales::scientific, breaks = scales::pretty_breaks(n=3)) + geom_histogram( fill="slategray3") + xlab("Mean Difference")
  p<-plot_grid(DPermGraph  +
                 theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")), 
               DPrimePermGraph +
                 theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")), 
               RSquarePermGraph +
                 theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")), ncol=1)
  p
}




theme_set((theme_bw(base_size=22)))



set.seed(1)
yri_df_b<-as_tibble(read_rds("../data2/YRI_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="YRI") 




plotPermutations3(yri_df_b, population="YRI")
