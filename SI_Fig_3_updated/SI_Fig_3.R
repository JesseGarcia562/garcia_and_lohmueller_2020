library(tidyverse)

ldTable<-as_tibble(fread("decomposeGravel/LDTableOfAllACExceptIntergenic.csv"))


ldTable$Scenario<-case_when(

ldTable$Selection == "../data//SelectionYes" ~ "Full Gravel Model",

ldTable$Selection == "../data/decomposeGravel//DFESelectionGravelOnlyYRI" ~ "Population Size of 7,310 with Expansion to 14,474",

ldTable$Selection == "../data/decomposeGravel//DFESelectionGravelNoMigration" ~ "Gravel Model Without Migration",

ldTable$Selection == "../data/decomposeGravel//DFESelectionConstantPopGravel" ~ "Constant Population Size of 14,474"

)


si_fig_3 <- ldTable %>% 
  filter(AC < 100,Scenario == "Constant Population Size of 14,474") %>% 
  group_by(Scenario, AC, Variation) %>% summarise(meanDist=mean(Distance),  meanSel=mean( ((SelectionCoefficient.x+SelectionCoefficient.y)/2)) ) %>% 
  ggplot(aes(x=AC, y=meanSel, colour=Variation)) + geom_line() + geom_point() + theme_bw(base_size=16) + 
  labs(y="Mean Selection Coefficient", x="Allele Count")   

si_fig_3

ggsave(filename="../figures/si_figure_3_mean_selection.tiff", plot=si_fig_3, width=20, height=12)
