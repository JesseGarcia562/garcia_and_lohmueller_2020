library(tidyverse)

yri_df_b<-as_tibble(read_rds(path ="../data2/YRI_b_value_breaks_5_physical_breaks_200.rds")) %>%
  dplyr::rename(D=X6,variation=variation_type, r_square=genetic_distance, genetic_distance=X3, DPrime=X7) %>%
  mutate(Variation=glue("={variation}"), rSquare=r_square) %>% mutate(population="YRI") %>%  mutate(Variation=case_when(
  Variation == "=nonsynonymous_SNV" ~ "Nonsynonymous", 
  Variation == "=synonymous_SNV" ~ "Synonymous" 
)) 


yri_df_b %>%
  ggplot(aes(x=genetic_distance, y=DPrime, colour=Variation)) +
  geom_smooth(method='gam', formula = y ~ s(x, bs = "cs")) +
  theme_bw() +
  labs(y="Mean D' of Matched Pairs", x="Genetic Distance") + 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2")) + 
  theme(legend.position = "bottom") +
  theme(text=element_text(size=22))

ggsave("figure_4_a_august_24_2020.png", width=15, height =8)
