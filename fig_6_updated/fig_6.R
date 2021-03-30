library(data.table)
library(tidyverse)
library(glue)
library(binr)
library(mltools)
library(janitor)

## Colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



midpoints <- function(x, dp=2){
lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
return(round(lower+(upper-lower)/2, dp))
}



getCutBreaks <- function(x, dp=30){
lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
upper <- tail(as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x))),n=1)
return(c(lower, upper))
}


getBins <- function(x, dp=30){
lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))




bins<-tibble(
  lower=lower,
  upper=upper
)

return(bins)
}


findTheNumberOfVariantsPerAC <- function(empiricalDf) {
  
  counts <- empiricalDf %>% 
    group_by(ac, geneticDistanceBreaks, variation) %>%
    count()
  
  return(counts)
  
}


binSimulatedIntoEmpiricalBreaks <- function(binnedEmpiricalBreaks, simulatedGeneticDistance) {
  bins <- getBins(levels(binnedEmpiricalBreaks))
  simulatedGeneticDistanceBreaks <- bin_data(simulatedGeneticDistance, boundaryType = "lcro]", bins=bins) 

  return(simulatedGeneticDistanceBreaks)
  
  
}


sampleFromCounts <- function(binnedSimulation, empiricalCounts){
  samples <- list()
  for (i in 1:nrow(empiricalCounts) ){
    
    AC <- empiricalCounts$ac[i]
    geneticDistanceBreaks <- empiricalCounts$geneticDistanceBreaks[i]
    variation <- empiricalCounts$variation[i]
    sampleSize <- empiricalCounts$n[i]
    
    
    test <- binnedSimulation %>%
      filter(ac == AC, geneticDistanceBreaks == !!geneticDistanceBreaks, variation == !!variation) %>%
      sample_n(sampleSize)
    
    samples[[i]] <- test
  }
  
  samples <- bind_rows(samples)
  return(samples)
}


differenceInD <- function(df){
  nsD <- df %>% filter(variation == "Nonsynonymous") %>% pull(mean_d)
  sD <- df %>% filter(variation == "Synonymous" ) %>% pull(mean_d)
  
  difference <- nsD-sD
  ratio <- nsD/sD
  normalized_difference <- ( mean(nsD) - mean(sD) ) / mean(sD)
  
  ldDf <- tibble(
    difference=difference,
    ratio=ratio,
    normalized_difference=normalized_difference,
    nsD=nsD, 
    sD=sD  )
  
  return(ldDf)
  
}


binEmpiricalIntoSpecificEmpiricalBreaks <- function(binnedEmpiricalBreaks, empiricalGeneticDistance) {

   
  bins <- getBins(levels(binnedEmpiricalBreaks))
  empiricalGeneticDistanceBreaks <- bin_data(empiricalGeneticDistance, boundaryType = "lcro]", bins=bins) 

  return(empiricalGeneticDistanceBreaks)
  
  
}


mean_normalized_difference_in_d <- function(df, population){
  df %>% 
  ungroup() %>% 
  filter(genetic_distance != 0) %>%
  mutate(geneticDistanceBreaks=bin_data(genetic_distance, bins=5,binType = "quantile" )) %>%
  mutate(variation=case_when(variation=="synonymous_SNV" ~ "Synonymous", variation=="nonsynonymous_SNV" ~ "Nonsynonymous"    )) %>%
  group_by(variation, geneticDistanceBreaks) %>% add_tally() %>%
  summarise(mean_d=mean(d),count=unique(n) )  %>%  
  ungroup() %>% 
  group_by(geneticDistanceBreaks) %>%
  group_map( ~ differenceInD(.x) ) %>%
  ggplot(aes(x=geneticDistanceBreaks, y=normalized_difference, group=1)) +
  geom_line() +
  labs(x="Genetic distance quantiles (cM)", y="Mean normalized difference in D (NS - S) / S ", colour="", subtitle=population) +
  theme_bw() +
  geom_hline(yintercept = 0, colour="grey", size=2, linetype="dashed")
  
  
}

bin_one_population_into_another_compute_mean_difference <- function(ld_df,population_of_ld_df, empiricalDf) {
  
  ld_df$geneticDistanceBreaks <- binEmpiricalIntoSpecificEmpiricalBreaks( binnedEmpiricalBreaks = empiricalDf$geneticDistanceBreaks, 
                                           empiricalGeneticDistance = ld_df$genetic_distance)
  
  ld_df %>% 
    ungroup() %>% 
    filter(genetic_distance != 0) %>%
    mutate(variation=case_when(variation=="synonymous_SNV" ~ "Synonymous", variation=="nonsynonymous_SNV" ~ "Nonsynonymous"    )) %>%
    group_by(variation, geneticDistanceBreaks) %>% 
    add_tally() %>%
    summarise(mean_d=mean(d), count_of_pairs=unique(n) )  %>%  
    ungroup() %>% 
    group_by(geneticDistanceBreaks) %>%
    group_map( ~ differenceInD(.x) ) %>%
    mutate(population=population_of_ld_df, type=parse_character(glue("Empirical ({population_of_ld_df})"))) %>% 
    ungroup()

}


differenceInD_2 <- function(df){
  nsD <- df %>% filter(variation == "Nonsynonymous") %>% pull(mean_d)
  sD <- df %>% filter(variation == "Synonymous" ) %>% pull(mean_d)
  
  difference <- nsD-sD
  ratio <- nsD/sD
  normalized_difference <- ( mean(nsD) - mean(sD) ) / mean(sD)
  count_of_ns_pairs <- df %>% filter(variation == "Nonsynonymous") %>% pull(count_of_pairs)
  count_of_s_pairs <- df %>% filter(variation == "Synonymous") %>% pull(count_of_pairs)

  ldDf <- tibble(
    difference=difference,
    ratio=ratio,
    normalized_difference=normalized_difference,
    nsD=nsD, 
    sD=sD,
    count_of_ns_pairs=count_of_ns_pairs,
    count_of_s_pairs=count_of_s_pairs
  )
  
  return(ldDf)
  
}

ACFilteredLD <- readRDS("ACFilteredLDWithGeneticDistanceJan18.RDS")
matchedPhysicalBreaks50Data <- readRDS("matchedPhysicalBreaks50Data.rds")

mergedEmpirical <- as_tibble(merge(ACFilteredLD, matchedPhysicalBreaks50Data)) %>% 
  filter(AC <= 5, GeneticDistance != 0) %>% 
  clean_names() %>% 
  mutate(variation=case_when(
  variation == "=nonsynonymous_SNV" ~ "Nonsynonymous", 
  variation == "=synonymous_SNV" ~ "Synonymous" 
))

  

empiricalDf <- mergedEmpirical %>%
  mutate(geneticDistanceBreaks=bin_data(genetic_distance, bins=5, binType="quantile"))


empiricalCounts <- findTheNumberOfVariantsPerAC(empiricalDf = empiricalDf)
set.seed(125)

samples <- read_rds("resamples_small_recomb_rates_ac_1_5_june_29_2019.rds")





test <- samples %>%
  group_by(resample, variation, geneticDistanceBreaks) %>%
  summarise(mean_d=mean(d)) %>%
  group_by(resample, geneticDistanceBreaks) %>%
  group_modify( ~ differenceInD(.x) )

test$midpoints <- midpoints(test$geneticDistanceBreaks, dp=30)


summaryOfDifference <- test %>%
  ungroup() %>%
  group_by(geneticDistanceBreaks) %>%
  summarise(median_difference=median(difference), 
            median_ratio=median(ratio),
            mean_difference=mean(difference),
            mean_ratio=mean(ratio),
            mean_normalized_difference=mean(normalized_difference)) %>%
  ungroup() 

empiricalDifferences <- empiricalDf  %>%
  group_by(variation, geneticDistanceBreaks) %>%
  summarise(mean_d=mean(d))  %>%  
  ungroup() %>% 
  group_by(geneticDistanceBreaks) %>%
  group_modify( ~ differenceInD(.x) )

neutralModel <- tibble(
  difference=0, 
  type="Neutral", 
  resample=4444, 
  alpha=1,
  geneticDistanceBreaks=summaryOfDifference$geneticDistanceBreaks,
  size=1, 
  colour="#999999",
  ratio=1,
  normalized_difference=0
)



graphDf <- bind_rows(test %>% mutate(type="Resample mean", alpha=0.4, size=1, colour="#E69F00") , 
                   summaryOfDifference %>% mutate(difference=mean_difference,ratio=mean_ratio,normalized_difference=mean_normalized_difference, type="Simulated", resample=9999, alpha=1, size=1.1, colour="#E69F00") ,
                   empiricalDifferences %>% mutate(type="Empirical (YRI)" , resample=10000, alpha=1, size=1.1,colour="#009E73" ),
                   neutralModel
                   ) 




cols <- c("Simulated" = "#E69F00", "Empirical (YRI)" = "#009E73", "Resample mean" ="#E69F00", "Neutral" ="#999999")
lines <- c("Simulated" = "solid", "Empirical (YRI)" = "solid", "Resample mean" ="solid", "Neutral" ="dashed")
alphas <- c("Simulated" = 1, "Empirical (YRI)" = 1, "Resample mean" =0.19, "Neutral" =1)
sizes<-c("Simulated" = 3, "Empirical (YRI)" = 3, "Resample mean" =1, "Neutral" =2)
xbreaks<-c("[1.62011535221e-11, 6.35875037262e-05)" = "[1.62e-11, 6.36e-05)",
           "[6.35875037262e-05, 0.000597266333334)" = "[6.36e-05, 5.97e-04)",
           "[0.000597266333334, 0.00252735883636)" ="[5.97e-04, 2.53e-03)",
            "[0.00252735883636, 0.00928393376262)" = "[2.53e-03, 9.28e-03)",
            "[0.00928393376262, 0.732906302888]" =  "[9.28e-03, 0.733]"
           )



graphDf %>%
  ggplot(aes(x=geneticDistanceBreaks, y=normalized_difference, group=resample, color=type, alpha=type, size=type, linetype=type)) +
  geom_line() +
  labs(x="Genetic distance quantiles (cM)", y="Mean normalized difference in D (NS - S) / S ", colour="") +
  scale_color_manual(values=cols, breaks=c("Empirical (YRI)", "Neutral", "Simulated")) +
  scale_linetype_manual(values=lines, guide=F) +
  scale_alpha_manual(values=alphas, guide=F) +
  scale_size_manual(values=sizes, guide=F) + 
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  scale_x_discrete(labels=xbreaks) + 
  theme(legend.position = "bottom")



