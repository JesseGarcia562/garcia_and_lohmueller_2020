library(tidyverse)
library(glue)
library(mltools)

midpoints <- function(x, dp=2){
lower <- as.numeric(gsub(",.*","",gsub("\\(|\\[|\\)|\\]","", x)))
upper <- as.numeric(gsub(".*,","",gsub("\\(|\\[|\\)|\\]","", x)))
return(round(lower+(upper-lower)/2, dp))
}

# Analyzing Lower recombination rate doubletons

lowRecombSLD<-read_rds("gravelMigrationAndNoMigrationSimulations_attemptFour_smallerRecombinations_June10_SLD.rds")



NSLD<-lowRecombSLD %>% select(recombinationRate, seed, migration, NSLD) %>% unnest(NSLD) %>% unnest(unbiasedD)
SLD<-lowRecombSLD %>% select(recombinationRate, seed, migration, SLD) %>% unnest(SLD)
SLD$type<-SLD$unbiasedD %>% map(~ typeof(.x)) %>% unlist()
SLD<-SLD %>% 
  filter(type != "NULL") %>%
  unnest(unbiasedD)

lowRecombLD<-bind_rows(SLD %>% mutate(variation="Synonymous"), NSLD %>% mutate(variation = "Nonsynonymous")) %>% 
  mutate(dprime=case_when(
  d > 0 ~ d/(.02*.98),
  d < 0 ~ d/(.02*.02),
  d == 0 ~ 0
  
  
  
)) %>%
  mutate(recombinationRateNumber=recombinationRate) %>%
  mutate(recombinationRate = glue("r={recombinationRate}")) %>%
  mutate(recombinationRate=fct_rev(recombinationRate)) %>%
  mutate(geneticDistance=distance*recombinationRateNumber*(1/.01)) 


si_fig_6 <- lowRecombLD %>% 
  filter(migration=="YesMigration") %>% 
  group_by(recombinationRateNumber ) %>%
  mutate(geneticDistanceBreaks=as.character(bin_data(geneticDistance, bins=5, binType="quantile"))) %>% 
  ungroup() %>%
  group_by(recombinationRate, variation, geneticDistanceBreaks) %>%
  summarise(completeRepulsion=mean( dprime == -1)) %>%
  mutate(geneticDistanceMidpoints=midpoints(geneticDistanceBreaks, dp = 30)) %>%
  ggplot(aes(x=geneticDistanceMidpoints, y=completeRepulsion, colour=variation)) +
  geom_point() +
  facet_grid(~ recombinationRate, scales = "free") +
  geom_line()+ 
  labs(x="Genetic Distance Midpoints", y="Fraction of Doubletons with complete repulsion (D' = -1)", colour="Variation" )  + 
  theme_bw() +
  theme(text=element_text(size=18))



si_fig_6
ggsave(filename="figures/si_figure_6_fraction_repulsion.tiff", plot=si_fig_6, width=20, height=12)
