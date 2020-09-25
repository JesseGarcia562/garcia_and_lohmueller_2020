library(dplyr)
library(GenomicRanges)
library(tidyverse)
library(glue)
nameBValueColumns<-function(bValues){
  names(bValues) <- c("Chr", "Start", "End", "BValue")
  return(bValues)
}

turnLDTableIntoBed<-function(ACFilteredLD){
LDBedStartPos1<-ACFilteredLD$POS1
LDBedStartPos2<-ACFilteredLD$POS2
chr<-ACFilteredLD$Chromosome


bedStarts<-c(LDBedStartPos1, LDBedStartPos2)
bedEnds<-c(LDBedStartPos1, LDBedStartPos2)
bedChr<-chr

ldBed<-data.frame(Chr=chr, Start=bedStarts, End=bedEnds)

#Removes duplicated rows
uniqueLDBed<-ldBed %>% distinct
return(uniqueLDBed)
}


ldBedIntoGRange<-function(uniqueLDBed){
ranges = IRanges(start=uniqueLDBed$Start, end =uniqueLDBed$Start )

ldGRange<-GRanges(seqnames = uniqueLDBed$Chr,ranges)
return(ldGRange)
}


bValueBedIntoGRange<-function(bValues){
bValueGRange<-GRanges(seqnames=bValues$Chr,IRanges(start=bValues$Start, end=bValues$End), BValue=bValues$BValue)
return(bValueGRange)
}


intersectLDAndBValue<-function(ldGRange, bValueGRange,uniqueLDBed, bValues){
hits<-findOverlaps(ldGRange, bValueGRange)


indexOfBValueIntersections<-subjectHits(hits)
indexOfLDIntersections<-queryHits(hits)

intersectedLD<-uniqueLDBed[indexOfLDIntersections,]

intersectedB<-bValues[indexOfBValueIntersections,]

#Removing rownames for debugging. 
rownames(intersectedB)<-c()

intersectedLD$BValues<-intersectedB$BValue

return(intersectedLD)
}



append_bvalues<-function(
			 empirical_ld_table,
			 hg19_bvalue_bed
			 ){
ACFilteredLD<-read_rds(empirical_ld_table) %>% dplyr::rename(POS1 = X1, POS2 = X2) %>% mutate(Chromosome=glue("chr{chromosome}"))
bValuesFile<-hg19_bvalue_bed
bValues<-read.delim(bValuesFile, header = F, stringsAsFactors=F)
bValues<- as_tibble(nameBValueColumns(bValues))


#Changing LD Table into BED Format

uniqueLDBed<-turnLDTableIntoBed(ACFilteredLD = ACFilteredLD)
ldGRange<-ldBedIntoGRange(uniqueLDBed = uniqueLDBed )

bValueGRange<-bValueBedIntoGRange(bValues=bValues)
intersectedLD<-intersectLDAndBValue(ldGRange = ldGRange, bValueGRange = bValueGRange, uniqueLDBed=uniqueLDBed, bValues=bValues)
joinedPOS1<-inner_join(x=ACFilteredLD, y=intersectedLD, by= c("POS1"="Start", "Chromosome"="Chr"))
ACFilteredLDPOS1<-joinedPOS1 %>% rename("BValues" = "BValuePOS1")
joinedPOS2<-inner_join(x=ACFilteredLDPOS1, y=intersectedLD, by=c("POS2"= "Start", "Chromosome" = "Chr"))
ACFilteredLDBValues<-joinedPOS2 %>% rename("BValues" = "BValuePOS2")

ACFilteredLDBValues<-ACFilteredLDBValues %>% select(-c(End.x,End.y))


}
