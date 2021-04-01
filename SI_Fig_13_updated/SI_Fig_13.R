library(tidyverse)


yri_df_b<-as_tibble(read_rds(file ="YRI_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="YRI") %>%  mutate(Variation=case_when(
  Variation == "=nonsynonymous_SNV" ~ "Nonsynonymous SNV", 
  Variation == "=synonymous_SNV" ~ "Synonymous SNV" 
)) 

yri_df_b$GeneticDistanceBreaks<-bin_data(yri_df_b$genetic_distance, bins=10, binType="quantile")


pos_1<-yri_df_b %>% filter(variation=="synonymous_SNV") %>% dplyr::rename(POS1_syn = POS1, POS2_syn=POS2)


pos_2<-yri_df_b %>% filter(variation == "nonsynonymous_SNV") %>% dplyr::rename(POS1_nonsyn = POS1, POS2_nonsyn=POS2)


matches<-bind_cols(pos_1, pos_2,)



## Bind columns changes the names of columns. physical_distance...4 is the name for the distance between the synonymous pair and 
## physical_distance...27 is the name of the distance between nonsynonymous pairs.
plot_1<-matches %>%
  ggplot(aes(x=physical_distance...4, y=physical_distance...27)) +
  geom_point() +
  labs(x="Distance (bp) between pair of NS SNPs", y="Distance (bp) between pair of S SNPs")


plot_2<-yri_df_b %>%
  ggplot(aes(x=allele_count, fill=Variation)) +
  geom_bar(position = "dodge_2") +
  labs(x="Allele count", y="Count", fill="Variation")


plot_3<-yri_df_b %>%
  ggplot(aes(x=sqrt(genetic_distance), color=Variation)) +
  geom_density() +
  labs(x="Square Root of Genetic Distance (cM)", color="Variation") 



plot_4<-yri_df_b %>%
  ggplot(aes(x=(BValuePOS1+BValuePOS2)/2, color=Variation)) +
  geom_density() +
  labs(x="Mean B Value", colour="Variation")

plot_5<-yri_df_b %>%
  ggplot(aes(x=chromosome, fill=Variation)) +
  geom_bar(position = "dodge_2") +
  labs(x="Chromosome", y="Count", fill="Variation")

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  plot_2 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(plot_1,plot_2 + theme(legend.position = "none"), plot_3 + theme(legend.position = "none"), plot_4 + theme(legend.position = "none"), plot_5+theme(legend.position = "none"),legend,labels = c("A","B", "C", "D", "E", "") )

