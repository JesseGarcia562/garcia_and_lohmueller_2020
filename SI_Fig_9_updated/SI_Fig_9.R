library(tidyverse)
library(data.table)
library(glue)


label_translation<-c("NoMigration" = "Model 3 \n(No migration)", "YesMigration" = "Model 2 \n(Migration)")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbbPalette_labels <- c("African"="#56B4E9", "African (no migration)"="#009E73","East Asian"="#F0E442", "European" = "#0072B2")
# Path on Computer
infoOfVariants<-read_rds("gravel_migration_no_migration_variant_metadata.rds") %>% ungroup() 


si_fig_9 <- infoOfVariants %>% 
  mutate(recombinationRate = case_when(
    recombinationRate == 1e-09 ~ "r=1e-09", 
    recombinationRate == 1e-08 ~ "r=1e-08"
  )) %>%
  mutate(recombinationRate=factor(recombinationRate, levels=c("r=1e-09", "r=1e-08"))) %>%
  filter(ac==2) %>%
  group_by(migration, seed,extra_po, recombinationRate ) %>%
  summarise(mean_sel=mean(selection_coefficient)) %>% 
  ungroup() %>%
  mutate(extra_po = case_when(
    extra_po == "No Migration" ~ "African (no migration)", 
    TRUE ~ extra_po
  )) %>%
  mutate(pretty_migration = factor(label_translation[migration]) )%>%
  mutate(pretty_migration = fct_rev(pretty_migration) )%>%
  ggplot(aes(x=pretty_migration, y=mean_sel, fill=extra_po)) +
           geom_boxplot(position = position_dodge(preserve = "single")) + 
  labs(x="", y="Mean selection coefficient", fill="Population origin") +
  facet_wrap(~recombinationRate) +
  scale_fill_manual(values=cbbPalette_labels) +
    theme_bw() +
  theme(text=element_text(size=22))
si_fig_9

ggsave(filename="figures/si_figure_9_selection_migration_no_migration.tiff", plot=si_fig_9, width=20, height=12)

