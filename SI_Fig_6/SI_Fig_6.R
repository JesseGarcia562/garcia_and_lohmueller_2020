library(data.table)
library(tidyverse)
library(glue)

ldTable<-as_tibble(fread("../data/decomposeGravel/LDTableOfAllACExceptIntergenic.csv"))


ldTable$Scenario<-case_when(

ldTable$Selection == "../data//SelectionYes" ~ "Full Gravel Model",

ldTable$Selection == "../data/decomposeGravel//DFESelectionGravelOnlyYRI" ~ "Population Size of 7,310 with Expansion to 14,474",

ldTable$Selection == "../data/decomposeGravel//DFESelectionGravelNoMigration" ~ "Gravel Model Without Migration",

ldTable$Selection == "../data/decomposeGravel//DFESelectionConstantPopGravel" ~ "Constant Population Size of 14,474"

) 
ldTable %>% filter(AC < 100,Scenario == "Constant Population Size of 14,474" ) %>% group_by(Scenario, AC, Variation) %>% summarise(meanDist=mean(Distance),  meanSel=mean( ((SelectionCoefficient.x+SelectionCoefficient.y)/2)) ) %>% 
  ggplot(aes(x=AC, y=meanSel, colour=Variation)) + geom_line() + geom_point() + theme_bw(base_size=16) + 
  labs(y="Mean Selection Coefficient Across Pair")   +
  facet_grid(~Scenario )
