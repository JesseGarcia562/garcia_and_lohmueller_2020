library(tidyverse)

oneTest<-function(nonsynonymousStats, synonymousStats){
differences<-nonsynonymousStats - synonymousStats
#For Randomization test on matched samples you need to randomly change the sign of the differences. 
#Sample 1 and -1, multiply by differences
signs <- sample(size=length(differences), c(1,-1), replace = TRUE)
newsample <- signs*differences
randomizedMean<-mean(newsample)
return(randomizedMean)
}

permute<-function(nonsynonymousStats,synonymousStats, replications = 10000, statistic="Some LD Statistic" ){

  actualDiff<-mean(nonsynonymousStats - synonymousStats)
  permutations<-replicate(replications, oneTest(nonsynonymousStats, synonymousStats))
  
  
  
     
   if (statistic == "rSquare"){
     title <- bquote("Randomization Test on Matched Samples" ~ R^2)
   } else if (statistic == "DPrime") {
     title <- bquote( ~ "Randomization Test on Matched Samples D' ")
   } else if (statistic == "D"){
     title <- bquote(~ "Randomization Test on Matched Samples D ")
   } else if (statistic == "FreqAB"){
     title <- bquote(~ "Randomization Test on Matched Samples Count of AB Haplotypes ")
   }
  
  
  
  
  #hist(permutations, main = title, xlim=c(actualDiff ,max(permutations)   ), abline(v=actualDiff, col =  "blue", lty = 3))
       
  permDF<- as.data.frame(permutations)
  names(permDF)<-"MeanDif"
       ggplot(permDF, aes(x=MeanDif)) + geom_histogram() + labs(title=title, x="Mean Difference Between Matched Pairs", y="Frequency") + geom_vline(aes(xintercept = actualDiff, color = "Observed Mean Difference"), linetype=2)+ theme_bw()



}


prepMatchedDataForPermutation<-function(MatchedLDData){
  MatchedLDData$Variation<-as.character(MatchedLDData$Variation)
MatchedLDData$Variation[MatchedLDData$Variation == "=nonsynonymous_SNV"] <- "Nonsynonymous SNV"
MatchedLDData$Variation[MatchedLDData$Variation == "=synonymous_SNV"] <- "Synonymous SNV"
return(MatchedLDData)
}

permuteLDStat<-function(preppedMatchedData, SummaryStat){
  nonsynonymousStats<-preppedMatchedData[[SummaryStat]][preppedMatchedData$Variation == "Nonsynonymous SNV" ]
  
  
   synonymousStats<-preppedMatchedData[[SummaryStat]][preppedMatchedData$Variation == "Synonymous SNV" ]
   

   
   
permute(nonsynonymousStats=nonsynonymousStats, synonymousStats=synonymousStats, statistic = SummaryStat )

}

plottingDPermutations<-function(MatchedLDData, LDStat="D"){
  preppedMatchedData<-prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + theme(legend.position="none")
}

plottingDPrimePermutations<-function(MatchedLDData, LDStat="DPrime"){
  preppedMatchedData<-prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + labs(color="") + theme(legend.position = c(0.89,.8)) + theme(legend.position="none")

}

plottingRSquarePermutations<-function(MatchedLDData, LDStat="rSquare"){
  preppedMatchedData<-prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + theme(legend.position="none")
}

plottingFreqABPermutations<-function(MatchedLDData, LDStat="FreqAB"){
  preppedMatchedData<-prepMatchedDataForPermutation(MatchedLDData)
  permuteLDStat(preppedMatchedData=preppedMatchedData,SummaryStat= LDStat) + theme(legend.position="none")
}



plotPermutations<-function(MatchedLDData){
DPermGraph<-plottingDPermutations(MatchedLDData = MatchedLDData) 
DPrimePermGraph<-plottingDPrimePermutations(MatchedLDData = MatchedLDData) 
RSquarePermGraph<-plottingRSquarePermutations(MatchedLDData = MatchedLDData)
FreqABPermGraph<-plottingFreqABPermutations(MatchedLDData = MatchedLDData)
plot_grid(DPermGraph, DPrimePermGraph, RSquarePermGraph, FreqABPermGraph, labels = "AUTO", ncol=2)
}





#No Selection New DFE Recom 1e-8 (Scaled 1e-7)
matchedLDTableNoSelectionAvgRecombDataNewDFE<-readRDS("/Users/jessegarcia/Documents/LD Simulations Backup/burninVCF/h_0.5/matchedLDTableNoSelectionAvgRecombNewDFE.rds")
si_fig_14 <- plotPermutations(matchedLDTableNoSelectionAvgRecombDataNewDFE)
si_fig_14
ggsave(filename="figures/si_figure_14_no_selection_permutation.tiff", plot=si_fig_14, width=20, height=12)
