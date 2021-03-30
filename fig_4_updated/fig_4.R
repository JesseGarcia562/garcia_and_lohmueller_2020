library(GenomicRanges)
library(tidyverse)
library(janitor)
library(data.table)
library(testthat)
library(glue)
library(vroom)
library(cowplot)
library(viridis)
library(binr)
library(mltools)
library(FSA)
library(knitr)
library(kableExtra)
library(ggbeeswarm)
library(janitor)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



ceu_df_b <- as_tibble(read_rds(path = "ceu_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="CEU")
chb_df_b <- as_tibble(read_rds(path ="chb_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="CHB")
mxl_df_b <- as_tibble(read_rds(path = "MXL_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="MXL")
jpt_df_b <- as_tibble(read_rds(path ="JPT_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="JPT")
yri_df_b <- as_tibble(read_rds(path ="YRI_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="YRI")


#yri_df_b
empirical_matched <- bind_rows(ceu_df_b, chb_df_b,mxl_df_b , jpt_df_b)


plot_ld_across_genetic_distance <- function(ACFilteredLD,AlleleCountMax, bins, removeZeroCM){
  
  if(removeZeroCM){
    doubletonsLD<- na.omit(ACFilteredLD %>% filter( AC <= AlleleCountMax & GeneticDistance != 0))
  } else {
  doubletonsLD<- na.omit(ACFilteredLD %>% filter( AC <= AlleleCountMax))
}

doubletonsLD$GeneticDistanceBreaks <- bin_data(doubletonsLD$GeneticDistance, bins=bins, binType="quantile")


percentileNames <- c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%")


doubletonsLD <- doubletonsLD %>% 
  mutate(Variation=case_when(
    Variation=="=nonsynonymous_SNV" ~ "Nonsynonymous",
    Variation=="=synonymous_SNV" ~ "Synonymous"
  ))

return(doubletonsLD)

}






  
empirical_matched$GeneticDistanceBreaks <- bin_data(empirical_matched$genetic_distance, bins=10, binType="quantile")

empirical_matched <- empirical_matched %>%  
  mutate(Variation=case_when(
  Variation == "=nonsynonymous_SNV" ~ "Nonsynonymous", 
  Variation == "=synonymous_SNV" ~ "Synonymous" 
))

empirical_matched %>% 
  group_by(GeneticDistanceBreaks, Variation, population) %>%
  summarise(mean_d=mean(D)) %>%
  ggplot(aes(x=GeneticDistanceBreaks, y=mean_d, colour=Variation, group=Variation)) +
  geom_line(size=1.5) + 
  theme_bw() +
  facet_wrap(~population, ncol=1) + 
  labs(y="Mean D of Matched Pairs", x="Genetic Distance Bins") + 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2")) +
  theme(text=element_text(size=16), axis.text.x = element_text(angle=45, hjust=1, size=10)) 
  



empirical_matched %>% 
  group_by(GeneticDistanceBreaks, Variation, population) %>%
  summarise(mean_d=mean(D)) %>%
  ggplot(aes(x=GeneticDistanceBreaks, y=mean_d, colour=Variation, group=Variation)) +
  geom_point(size=1.5) + 
  geom_line(size=1) +
  theme_bw() +
  facet_wrap(~population, ncol = 1) + 
  labs(y="Mean D of Matched Pairs", x="Genetic Distance Bins (cM)") + 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2")) +
  scale_x_discrete(labels=c("3.47e-05", "3.32e-04", "1.56e-03", "6.46e-03", "6.28e-01")) +
  theme(text=element_text(size=16)) 


empirical_matched %>% 
  group_by(GeneticDistanceBreaks, Variation, population) %>%
  summarise(mean_d=mean(DPrime)) %>%
  ggplot(aes(x=GeneticDistanceBreaks, y=mean_d, colour=Variation, group=Variation)) +
  geom_point(size=1.5) + 
  geom_line(size=1) +
  theme_bw() +
  facet_wrap(~population, ncol=1) + 
  labs(y="Mean D' of Matched Pairs", x="Genetic Distance Bins (cM)") + 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2")) +
  scale_x_discrete(labels=c("3.47e-05", "3.32e-04", "1.56e-03", "6.46e-03", "6.28e-01")) +
  theme(text=element_text(size=18))
