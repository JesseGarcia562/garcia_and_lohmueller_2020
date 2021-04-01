library(data.table)
library(tidyverse)
library(glue)
library(mltools)
# Path on computer /Users/jessegarcia/Documents/SLiM_ParallelRProject/data/doubletonsLD_gravel_migration_April19.csv
ld<-fread("doubletonsLD_gravel_migration_April19.csv")
ld<-as_tibble(ld) %>% mutate(migration=case_when(
  migration=="NoMigration" ~ "No Migration", 
  migration == "YesMigration" ~ "Migration"
)) 

ld<-ld %>% mutate(dprime=case_when(
  d > 0 ~ d/(.02*.98),
  d < 0 ~ d/(.02*.02)
  
  
  
  
)) %>%
  mutate(recombinationRate = glue("r={recombinationRate}")) %>%
  mutate(recombinationRate=fct_rev(recombinationRate)) %>%
  as_tibble()



si_fig_8 <- ld %>% 
  ggplot(aes(x=distance, y=d, colour=variation, linetype=migration)) +
  geom_smooth() +
  facet_wrap(~ recombinationRate) + 
  xlim(0,100000) +
  labs(x="Distance (bp)", y="D", linetype="Migration", colour="Variation") + 
  scale_color_manual(values=c("Synonymous"="dodgerblue1", "Nonsynonymous"="darkorchid2")) +
   theme(strip.background =element_rect(fill="white"), text = element_text(size=18)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3), labels = scales::scientific) +
  theme_bw()
si_fig_8
ggsave(filename="figures/si_figure_8_d_migration_no_migration.tiff", plot=si_fig_8, width=20, height=12)
