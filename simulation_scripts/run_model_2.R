library(tidyverse)
library(glue)


df<-read_rds("../data/gravelMigrationAndNoMigrationSimulations_attemptThree_round2.rds")


SGETaskID <- parse_integer(Sys.getenv("SGE_TASK_ID")) 

glue("recomb: {df$recombinationRate[SGETaskID]} seed: {df$seed[SGETaskID]} migration: {df$migration[SGETaskID]}")

system( df$Script[SGETaskID] )
