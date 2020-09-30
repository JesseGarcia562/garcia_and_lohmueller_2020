library(data.table)
library(tidyverse)
library(glue)
library(furrr)
library(broom)
plan(multiprocess)
source("../code/functions/annotateVCF_Variation_AC.R")
functionToSource<-list.files(path="../code/functions/",pattern=".*",  full.names=T)
functionToSource %>% map(~source(.x))
possibly_computeLD<-possibly(computeLD, otherwise=NA)
possibly_annotateVCF_Variation_AC<-possibly(annotateVCF_Variation_AC, otherwise=NA)
possibly_fread<-possibly(fread, otherwise=NA)

analyze<-function(outputPath, distanceThreshold = 100000, AC=2, Variation=1){
	commandToFilter<-glue("grep -E '#|MT={Variation};AC={AC};DP=1000' '{outputPath}'")	
	
	AC<-enquo(AC)
	Variation<-enquo(Variation)
	vcf<-as_tibble(fread(cmd=commandToFilter, header=T)) 
	annotatedVCF<-annotateVCF_Variation_AC(vcf) %>%
	filter(AC == !!AC, Variation == !!Variation) %>%
	mutate(selection_coefficient=str_extract(INFO, pattern="(?<=S=).*;DOM")) %>%
	mutate(selection_coefficient=parse_number(selection_coefficient))
	LD<-computeLD(annotatedVCF,distanceThreshold = distanceThreshold) %>% select(item1, item2, distance, unbiasedD)
	LD<-left_join(LD,annotatedVCF %>% select(POS, selection_coefficient), by=c("item1"="POS")) 
	LD<-left_join(LD, annotatedVCF %>% select(POS,selection_coefficient), by=c("item2"="POS"))
}
possibly_analyze<-possibly(analyze, otherwise=as_tibble(NA))


simulation_df<-read_rds("../data/higher_frequency_ld_model_1.rds")


SGETaskID<-parse_integer(Sys.getenv("SGE_TASK_ID"))


SGETaskID=2

simulation_df<-simulation_df %>%
	mutate(fileExists=file.exists(vcfOutpath)) %>%
	mutate(migration="no") %>%
	filter(fileExists==TRUE)


AC<-1:100


simulation_df<-simulation_df %>%
	nest(recombinationRate, seed, genomeLength, Script, fileExists)


simulation_df<-crossing(simulation_df, AC)


set.seed(1)
subSample<-simulation_df %>%
	slice(SGETaskID)





 
extractFrequencyFilteredLDStats<-function(simulationDf,part){ 

AC=simulationDf$AC

simulationDf$NSLD<-simulationDf$vcfOutpath%>%
        map(~possibly_analyze(.x, AC = AC, Variation =1))
        
        simulationDf$SLD<-simulationDf$vcfOutpath %>%
        map(~possibly_analyze(.x, AC =AC, Variation = 2))






NSLD<-simulationDf %>%
	unnest(data) %>%
        select( recombinationRate, seed, migration, NSLD) %>%
	unnest(NSLD) %>% 
	mutate(Variation = "Nonsynonymous") %>%
	mutate(AC = AC)


SLD<-simulationDf %>%
	unnest(data) %>%
        select(recombinationRate, seed, migration, SLD) %>%
        unnest(SLD) %>%
	mutate(Variation = "Synonymous") %>%
	mutate(AC = AC)

LD<-bind_rows(NSLD, SLD)


write_rds(LD, glue("/u/home/j/jessegar/project-klohmuel/higher_frequency_ld/higher_frequency_ld_AverageRecombinationRate_September29_2020_AC_{AC}_part_{part}.rds"))

}



extractFrequencyFilteredLDStats(simulationDf=subSample, part=SGETaskID)
