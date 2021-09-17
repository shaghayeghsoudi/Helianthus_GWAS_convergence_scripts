library(stringr)
library(plyr)
library(dplyr)

### useful sources
# https://www.biostars.org/p/196393/
# https://stackoverflow.com/questions/30079720/union-and-intersection-of-intervals
# genomicRanges tutorial: 
# genomicRanges tutorial: 
#https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
#https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

 


setwd("~/Documents/input_data/tests_analysis_for_paper/null_w/clusters/sep15")

#total_genome<-#(153889008+180944244+168477679+178891894+219039782+103879295+103863077+153157343+209807207+246277002+168412038+166464245+197256366+174421972+171235706+188567792+214708607)
#2999293257


### stackoverflow link
### genomerange: https://stackoverflow.com/questions/30079720/union-and-intersection-of-intervals
############
data_all<-read.table(file ="convergent_clustering_0.9_1_cM_summary", header = TRUE)

### add unique direction id 
ann_arg <- data_all[which(data_all$direction=="Annuus_Argophyllus" | data_all$direction=="Argophyllus_Annuus"),]
ann_arg$direction_ID<-rep("Annuus_Argophyllus")

ann_pf <- data_all[which(data_all$direction=="Annuus_petfal" |data_all$direction=="petfal_Annuus" ),]
ann_pf$direction_ID<-rep("Annuus_petfal")

ann_pp <- data_all[which(data_all$direction=="Annuus_petpet" | data_all$direction=="petpet_Annuus"),]
ann_pp$direction_ID<-rep("Annuus_petpet")

arg_pf <- data_all[which(data_all$direction=="Argophyllus_petfal" | data_all$direction=="petfal_Argophyllus"),]
arg_pf$direction_ID<-rep("Argophyllus_petfal")

arg_pp <- data_all[which(data_all$direction=="Argophyllus_petpet" | data_all$direction=="petpet_Argophyllus"),]
arg_pp$direction_ID<-rep("Argophyllus_petpet")

pf_pp <- data_all[which(data_all$direction=="petfal_petpet" | data_all$direction=="petpet_petfal"),]
pf_pp$direction_ID<-rep("petfal_petpet")

data_all_withID<-do.call("rbind", list(ann_arg,ann_pf,ann_pp,arg_pf,arg_pp,pf_pp))
data_all_withID$start<-as.numeric(gsub(":.*", "", data_all_withID$range))
data_all_withID$end<-as.numeric(gsub(".*:", "", data_all_withID$range))

gr <- GRanges(seqnames = rep(1,nrow(data_all_withID)),IRanges(start = data_all_withID$start,end =      data_all_withID$end))


known<-read.table(file = "known_genes_FT.txt", header = TRUE)
colnames(known)<-c("chrom","start","end","ID")

gr2 <- GRanges(seqnames = rep(1,nrow(known)),IRanges(start = known$start,end =      known$end))











