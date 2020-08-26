library(data.table)
library(stringr)
library(tidyr)
library(reshape2)
library(plyr)
library(glue)
library(dplyr)
library(future)
library(furrr)
library(tidyverse)


frequencyFilteredLD<-function(VCF,distanceLimit=10000){
  SNPS<-VCF$POS
  Variation<-VCF$Variation[1]
  #IF only ONE Snp in this AC Category, return NULL
  if(length(SNPS) == 1){
    return(as_tibble(NULL))
  }
  
  #Calculate the comparisons that could be made.
  Comparisons<-combn(SNPS,2, simplify = TRUE)
  
  #Applying Distance Filter
  distances<-diff(Comparisons)
  Comparisons<- Comparisons[,abs(distances) < distanceLimit]
  
  
  AC=VCF$AC[1]
  #Removing AC Column
  VCF$AC<-NULL
  
  calculateLD<-function(Comparisons){

    pos1 = Comparisons[1]
    pos2 = Comparisons[2]
    
    
    
    gen1 = t(VCF[VCF$POS == pos1,][-c(1:7)])
    gen2 = t(VCF[VCF$POS == pos2,][-c(1:7)])
    
    # Get total number of alleles. Since Humans are diploids, multiply genotypes by 2
    total_alleles = nrow(gen1) * 2
    #Getting the count of the reference Alleles. 
    p_count = sum(str_count(gen1, "0"))
    q_count = sum(str_count(gen2, "0"))
    
    if (sum(str_count(gen1, fixed("..")) > 0) | sum(str_count(gen2, fixed("..")) > 0)){
      stop("Remove all uncallable variants and rerun analysis.")
    }
    
    if (p_count != q_count) {
      cat(pos1, "Has an allele frequency of ", p_count, pos2, " Has an allele frequency of ", q_count)
      stop("Your variants are not filtered by Allele count")
    }
    
    #Should not be happening in my data. Could happen in unfiltered Data
    if (total_alleles == 0){
      Dr = NA
      r2 = NA
      D_matrix[as.character(pos1), as.character(pos2)] = NA
      r_matrix[as.character(pos1), as.character(pos2)] = NA
    }
    
    #Getting the allele frequencies. p1 and q1 are Reference Alleles. 
    p1 = p_count / total_alleles
    p2 = 1 - p1
    q1 = q_count / total_alleles
    q2 = 1 - q1
    
    
    # Get observed genotypes. First split genotypes into Haplotypes
    gen1Split<-str_split_fixed(gen1, "\\|",2)
    colnames(gen1Split)<-c("Haplo1", "Haplo2")
    gen2Split<-str_split_fixed(gen2, "\\|",2)
    colnames(gen2Split)<-c("Haplo1", "Haplo2")
    #Recombine Haplotypes. Each haplotype1 element refers to the binary decision if Variant1 and Variant 2 Exist at Position1 and Position2 within 1 Individual on the left of the "|"
    haplotype1<-paste(gen1Split[,"Haplo1"],gen2Split[,"Haplo1"], sep ="")
    haplotype2<-paste(gen1Split[,"Haplo2"],gen2Split[,"Haplo2"], sep ="")
    
    gen_counts<-table(c(haplotype1, haplotype2))
    
    # Get expected genotypes counts
    #0 Refers to Reference, 1 Refers to Alternate
    exp_gen=c("00"=NA, "01"=NA, "10"=NA, "11"=NA)
    exp_gen["00"] = p1 * q1 * total_alleles
    exp_gen["01"] = p1 * q2 * total_alleles
    exp_gen["10"] = p2 * q1 * total_alleles
    exp_gen["11"] = p2 * q2 * total_alleles
    
    # Quantification of LD
    # Frequencies of genotypes
    freq_00 = gen_counts["00"] / total_alleles
    freq_01 = gen_counts["01"] / total_alleles
    freq_10 = gen_counts["10"] / total_alleles
    freq_11 = gen_counts["11"] / total_alleles
    
    AB=gen_counts["11"]
    Ab=gen_counts["10"]
    aB=gen_counts["01"]
    ab=gen_counts["00"]
    #Often times a specific genotype is never observed. This means its frequency should be 0.
    #However, because of the table function above, it causes the frequency to be NA. 
    
    fixZeroObserved <- function(observedFrequency){
      if (is.na(observedFrequency)){
        observedFrequency = 0
      } else{
        observedFrequency = observedFrequency
      }
      return(observedFrequency)
    }
    
    freq_00=fixZeroObserved(freq_00)
    freq_01=fixZeroObserved(freq_01)
    freq_10=fixZeroObserved(freq_10)
    freq_11=fixZeroObserved(freq_11)
    
    AB=fixZeroObserved(AB)
    Ab=fixZeroObserved(Ab)
    aB=fixZeroObserved(aB)
    ab=fixZeroObserved(ab)
    
    
    D = (freq_00*freq_11) - (freq_01*freq_10)
    
    
    # Compute r2
    
    r2 = D^2 / (q1 * q2 * p1 * p2)
    #Catching Divide by 0 errors
    if(is.infinite(r2)){
      r2=NA
    }
    
    
    
    # Compute Dr (D')
    if (D < 0) {
      Dr = D / min(p1*q1, p2*q2)
    } else {
      Dr = D / min(p1*q2, p2*q1)
    }
    if(is.infinite(Dr)){
      Dr=NA
    }
    
    Distance=pos2-pos1
    
    
    line<-c(POS1=pos1, POS2=pos2, DPrime=Dr, D=D, rSquare=r2, AB=AB, Ab=Ab, aB=aB, ab=ab, Distance=Distance, Variation=Variation, AC=AC)
    
  }
  
  #Checks if theres only two SNPS to calculate LD For, if so, then apply will NOT work. Need to just call the function calculateLD
  if(is.null(ncol(Comparisons))) {
    frequencyFilteredLD<-t(calculateLD(Comparisons)) 
  } else if ( ncol(Comparisons) == 0 ){
    #Sometimes No Comparisons can be Made. 
    return(as_tibble(NULL) )
  } else {
    frequencyFilteredLD<-t(apply(Comparisons, MARGIN = 2, FUN = calculateLD))
  }
  
  
  
  colnames(frequencyFilteredLD) <-c("POS1", "POS2", "DPrime", "D", "rSquare", "FreqAB","FreqAb", "FreqaB", "Freqab", "Distance", "Variation", "AC")
  
  return(as_tibble(frequencyFilteredLD))
}




vcfDirectory<-"../data/"
fileExtension <- "PolarizedAndACFix.vcf$"




computeFrequencyFilteredLD<-function(vcfDirectory="../data/", fileExtension="PolarizedAndACFix.vcf$"){
fileLocation <- list.files(path=vcfDirectory, pattern = fileExtension, full.names = TRUE)
fileNames<- list.files(path=vcfDirectory, pattern = fileExtension, full.names = FALSE)


df<-tibble(
  fileLocation=fileLocation,
  fileNames=fileNames
)


df<-df %>% separate(fileNames, c("Step", "Chromosome", "Population"),remove=F, extra="drop") %>% select( -Step, -fileNames) 

df<-df %>% expand(nesting(fileLocation,Chromosome, Population), Variation=c("=nonsynonymous_SNV", "=synonymous_SNV"))

df <-df %>% mutate(command=glue("grep -v  '^##' {fileLocation} | sed s/^#// | grep {Variation}")) 

df$VCF <- df$command %>% future_map(~as_tibble(fread(.x,stringsAsFactors=F,header=F,sep="\t", data.table = FALSE), .progress=TRUE))

# 

df$VCF<-df$VCF %>% future_map(~.x %>% separate(V8, c("AncestralAllele", "AAChange", "AC", "AF", "AFRAF"), sep=";" , remove=T, extra="drop") %>% mutate(AC=parse_number(AC), AF=parse_number(AF),AFRAF=parse_number(AFRAF)) %>% select(-AncestralAllele,-AAChange, -V3, -V4,-V5,-V6,-V7,-V9), .progress=TRUE)
df$VCF<-df$VCF %>% map(~.x %>% rename(Chr="V1", POS="V2"))

oldNames<-names(df$VCF[[1]])[str_detect(names(df$VCF[[1]]), "V")]
individualNumber<-seq_along(1:length(oldNames))
newNames<-glue("Individual{individualNumber}")

df$VCF <- df$VCF %>% map(~.x %>% rename_at( vars(oldNames), ~ newNames  )) 
               
df<-df%>% unnest(VCF) %>% filter(AC != 0, AC != 100) %>% select(-command, -fileLocation)

ld<-df %>% group_by(Chr, Variation, AC) %>% do(frequencyFilteredLD(.,distanceLimit=10000)) %>% ungroup() %>% mutate_at(vars(-Variation), parse_number)



return(ld)

}
