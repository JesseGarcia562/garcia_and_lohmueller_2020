library(glue)
library(tidyverse)
library(furrr)

library(pryr)
plan(multiprocess)

chromosome<-c(1:22)




acquireRawData<-function(individuals="NA18907,NA19235,NA18867,NA18933,NA18853,NA19108,NA18917,NA18923,NA18908,NA19247,NA18489,NA18510,NA18915,NA18881,NA19149,NA19131,NA19092,NA19225,NA18504,NA19239,NA18502,NA18916,NA19095,NA18912,NA19222,NA19213,NA18
871,NA19146,NA19198,NA18520,NA19118,NA18877,NA19209,NA19153,NA18517,NA18488,NA19184,NA18523,NA19214,NA19119,NA18873,NA19138,NA18879,NA19117,NA18861,NA18924,NA19207,NA19238,NA18909,NA19130",
                         bcftools="/u/local/apps/bcftools/1.2/gcc-4.4.7/bin/bcftools",
                         population="YRI",
                         chrNum){
  
  outputVCF<-glue("../data/01_chr{chrNum}_{population}_Individuals.vcf")
  command<-glue("{bcftools} view -s {individuals} /u/project/klohmuel/DataRepository/Human/Variants/1000G_Phase3_WGS_zippedVCFs/ALL.chr{chrNum}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > {outputVCF}")
  system(command)
}

polarizeVariants<-function(epoFastaDir="/u/project/klohmuel/DataRepository/Human/ReferenceGenomeFiles/EPOAncestral_release75_fasta/homo_sapiens_ancestor_GRCh37_e71_EPO",
                   population="YRI",
                   polarizer="AssignEPOAncestral_Allele_rmLowQual.py",
                   chrNum){
  
  
  inVCF<-glue("../data/01_chr{chrNum}_{population}_Individuals.vcf")
  ancestralAllele<-glue("{epoFastaDir}/homo_sapiens_ancestor_{chrNum}.fa")
  intermediateOut=glue("../data/Intermediate_{population}_EPORecode_Chr{chrNum}.tmp")
  
  
  
  outVCF=glue("../data/02_{population}_Chr{chrNum}_EPORefAlleleRecode.vcf")
  
  
  
  command<-glue("python {polarizer} {ancestralAllele} {inVCF} {intermediateOut} {outVCF} ")
  
  results<-system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  
}

annotateVariants<-function(annovar="../../annovar/table_annovar.pl", 
                           database='../../annovar/',
                           chrNum,
			   population){
  

  outVCF<-glue("../data/03_chr{chrNum}_{population}_Annotated.vcf")
  inVCF<-glue("../data/02_{population}_Chr{chrNum}_EPORefAlleleRecode.vcf")
  
  command<-glue("{annovar} {inVCF} {database}/humandb/ -buildver hg19 -out {outVCF} -remove -protocol refGene -operation g -nastring . -vcfinput")
  
  
  
  system(command)
  
  
}

fillAC<-function(fillAnACScript="/u/project/klohmuel/tanya_data/softwares/vcftools_perl/src/perl/fill-an-ac",
		 chrNum,
		 population){

inVCF<-glue("../data/03_chr{chrNum}_{population}_Annotated.vcf.hg19_multianno.vcf")
outVCF<-glue("../data/04_chr{chrNum}_{population}_PolarizedAndACFix.vcf")
command<-glue("{fillAnACScript} {inVCF} > {outVCF}")
system(command)
}


create_data_set<-function(individuals,chromosome, population){
chromosome %>% future_map(~acquireRawData(individuals=individuals 
					  ,chrNum = .x, population=population), .progress = T)

chromosome %>% future_map(~polarizeVariants(population=population, chrNum=.x), .progress=T)


chromosome %>% future_map(~annotateVariants(population=population, chrNum=.x), .progress=T)

chromosome %>% future_map(~fillAC(chrNum=.x, population=population), .progress=T)
}


create_data_set(individuals = "NA19785,NA19681,NA19658,NA19670,NA19722,NA19676,NA19664,NA19717,NA19777,NA19661,NA19756,NA19731,NA19788,NA19774,NA19752,NA19678,NA19749,NA19762,NA19735,NA19729,NA19725,NA19740,NA19761,NA19776,NA19764,NA19684,NA19792,NA19795,NA19682,NA19786,NA19789,NA19657,NA19755,NA19
770,NA19794,NA19782,NA19719,NA19758,NA19720,NA19771,NA19654,NA19651,NA19723,NA19663,NA19734,NA19746,NA19759,NA19747,NA19779,NA19783",
		chromosome=1:22,
		population="MXL"
		)


create_data_set(individuals = "NA18907,NA19235,NA18867,NA18933,NA18853,NA19108,NA18917,NA18923,NA18908,NA19247,NA18489,NA18510,NA18915,NA18881,NA19149,NA19131,NA19092,NA19225,NA18504,NA19239,NA18502,NA18916,NA19095,NA18912,NA19222,NA19213,NA18871,NA19146,NA19198,NA18520,NA19118,NA18877,NA19209,NA19
153,NA18517,NA18488,NA19184,NA18523,NA19214,NA19119,NA18873,NA19138,NA18879,NA19117,NA18861,NA18924,NA19207,NA19238,NA18909,NA19130",
		chromosome=1:22,
		population="YRI"
		)



create_data_set(individuals = "NA18975,NA18945,NA18977,NA18980,NA19079,NA19076,NA19084,NA18986,NA19009,NA18979,NA19058,NA19063,NA18956,NA18988,NA19005,NA19066,NA18962,NA19011,NA18994,NA19001,NA18947,NA18967,NA18966,NA18985,NA18954,NA19087,NA18989,NA19086,NA19006,NA19055,NA19007,NA18942,NA18997,NA19
064,NA18961,NA18959,NA19078,NA18999,NA19090,NA18974,NA19082,NA18952,NA18943,NA19081,NA19068,NA19083,NA18963,NA18973,NA19056,NA19062",
		chromosome=1:22,
		population="JPT"
		)
   
   
create_data_set(individuals= "NA12827,NA11843,NA07037,NA11829,NA11992,NA11831,NA10851,NA11920,NA12814,NA07056,NA12413,NA12154,NA12004,NA12815,NA12414,NA11892,NA12400,NA12776,NA12283,NA12156,NA12144,NA12347,NA12340,NA12777,NA12751,NA11994,NA07051,NA12872,NA12750,NA11918,NA1
1894,NA12749,NA11930,NA11933,NA11995,NA12044,NA12286,NA12046,NA12763,NA07347,NA06994,NA12341,NA06989,NA11931,NA12348,NA12716,NA11932,NA12874,NA12890,NA12830" 
					  ,chromosome=1:22, population="CEU")
            )
            
create_data_set(individuals= "NA18567,NA18533,NA18570,NA18574,NA18642,NA18639,NA18647,NA18591,NA18618,NA18573,NA18624,NA18630,NA18544,NA18593,NA18614,NA18632,NA18550,NA18619,NA18599,NA18611,NA18535,NA18559,NA18555,NA18582,NA18745,NA18545,NA18595,NA18558,NA18615,NA18621,NA
18616,NA18530,NA18603,NA18643,NA18549,NA18547,NA18592,NA18608,NA18749,NA18566,NA18645,NA18542,NA18531,NA18644,NA18634,NA18617,NA18638,NA18564,NA18622,NA18641" 
					  ,chromosome = 1:22, population="CHB")
