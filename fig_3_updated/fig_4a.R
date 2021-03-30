library(tidyverse)
library(data.table)
library(glue)

se <- function(x) sqrt(var(x)/length(x))


yri_df_b<-as_tibble(read_rds(file ="../data2/YRI_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="YRI") %>%  mutate(Variation=case_when(
  Variation == "=nonsynonymous_SNV" ~ "Nonsynonymous", 
  Variation == "=synonymous_SNV" ~ "Synonymous" 
)) 

yri_df_b$GeneticDistanceBreaks<-bin_data(yri_df_b$genetic_distance, bins=10, binType="quantile")

yri_plot<-yri_df_b %>%
    group_by(GeneticDistanceBreaks, Variation) %>%
  summarise(mean_d=mean(DPrime), 
            sd_d=se(DPrime)) %>%
  ggplot(aes(x=GeneticDistanceBreaks, y=mean_d, colour=Variation, group=Variation)) +
  geom_errorbar(aes(ymin=mean_d-sd_d, ymax=mean_d+sd_d), width=.2) +
  geom_point(size=1.5) + 
  geom_line(size=1) +
  theme_bw() +
  labs(y="Mean D' of Matched Pairs", x="Genetic Distance Bins (cM)") + 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2")) +
  scale_x_discrete(labels=c("1st quantile", "2nd quantile", "3rd quantile", "4th quantile", "5th quantile", "6th quantile"))+
  theme(text=element_text(size=22))





yri_plot + theme(legend.position = "bottom")
