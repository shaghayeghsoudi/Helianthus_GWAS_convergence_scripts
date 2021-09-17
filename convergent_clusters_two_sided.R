
# rm(list = ls())
### analyzing convergent clusters (numbers and length) ####
require(RColorBrewer)

## load libraries
library(Biostrings)   ## for tilegenome
#### to download biostrings:  https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(GenomicRanges)   ## for tilegenome, to download :https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
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


### break the main corercted file into each comparison type
chroms<-c("Ha412HOChr01","Ha412HOChr02","Ha412HOChr03","Ha412HOChr04","Ha412HOChr05","Ha412HOChr06","Ha412HOChr07","Ha412HOChr08","Ha412HOChr09","Ha412HOChr10","Ha412HOChr11","Ha412HOChr12","Ha412HOChr13","Ha412HOChr14","Ha412HOChr15","Ha412HOChr16","Ha412HOChr17")
compars<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","Argophyllus_petpet","petfal_petpet")

data_baypass<-data_all_withID[data_all_withID$analysis=="baypass",]
data_spearman<-data_all_withID[data_all_withID$analysis=="spearman",]

data_gwas_corr<-data_all_withID[data_all_withID$analysis=="GWAS_corrected",]
data_gwas_Uncorr<-data_all_withID[data_all_withID$analysis=="GWAS_uncorrected",]



### calulate size and number of convergent clusters per chromosome
out_res<-NULL
for (i in 1:length(compars)){
  
  for (j in 1:length(chroms)){
    data_focal<-data_baypass[data_baypass$direction_ID==compars[i] & data_baypass$chromosome==chroms[j] , ]
    data_focal$start<-as.numeric(gsub(":.*", "", data_focal$range))
    data_focal$end<-as.numeric(gsub(".*:", "", data_focal$range))
    
    data_sizenumber<-data_focal[,c(7,8)]
    data_sizenumberSum<-colSums(data_sizenumber)
    
    sub_good<-data.frame(t(colSums(data_sizenumber)))
    
    gr <- GRanges(seqnames = rep(1,nrow(data_focal)),IRanges(start = data_focal$start,end =      data_focal$end))
    gr_cluster_unique<-gr %>% unique
    union<-as.data.frame(reduce(gr_cluster_unique))[,2:3]
    union$size<-(union$end-union$start)
    aa<-data.frame(t(colSums(union)))
    
    sub_good<-cbind(aa,sub_good)
    colnames(sub_good)<-c("start_union","end_union","size_union","size", "N_convergent_window")
    
    sub_good$direction <- compars[i]
    sub_good$chromosome <- chroms[j]
    sub_good$N_convergent_cluster<-nrow(data.frame(gr_cluster_unique))
    out_res <- rbind (out_res,sub_good)
    
    
    #outliers_count <- tapply (data_baypass_focal$N_convergent, list(as.character (data_baypass_focal$chromosome=="chroms[j]")),sum)
    
  }
}

write.table(out_res,file = "out_res_convergent_CLUSTERS_2binning_Summed_Size_number_per_chromosome_Climate_BayPass_union_percomparison.table", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


###############################################################
##### plot the size and number convergent cluster results #####
###############################################################

setwd("~/Documents/input_data/tests_analysis_for_paper/null_w/clusters/sep15")
null_phenUnocorr<-read.table(file = "out_res_convergent_CLUSTERS_2binning_Summed_Size_number_per_chromosome_Climate_Soil_Spearman_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]


all_comparisonSize   <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     #### overlapped fixed (union)
#all_comparisonSize   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])    ### with overlaps
colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)



all_comparisonCount   <- cbind(ann_arg[,8], ann_pf[,8],ann_pp[,8],arg_pf[,8],arg_pp[,8],pf_pp[,8])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)



pdf (file = "barplots_convergent_clusters_perChrom_recomadjusted_BayPass_size_Union_number.pdf",width = 8, height= 12)
par(mar = c(5,4.5,4.1,2.1))
par(mfrow=c(2,1))

barplot(all_comparisonSize, 
        
        #beside=TRUE, ylim=c(0,2000000000),ylab = "size of convergent clusters", xlab = "chromosome", cex.lab=1.5,    ## with overlaps
        #cex.axis=1.2,width=c(0.08,0.08,0.08,0.08,0.08,0.08),cex = 0.9,las=0,col=brewer.pal(n = 6, name = "RdBu"))
        
        beside=TRUE, ylim=c(0,200000000),ylab = "size of convergent clusters", xlab = "chromosome", cex.lab=1.5,     cex.axis=1.2,width=c(0.08,0.08,0.08,0.08,0.08,0.08),cex = 0.9,las=0,col=brewer.pal(n = 6, name = "RdBu"))

legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

barplot(all_comparisonCount, 
        
        beside=TRUE, ylim=c(0,1000),ylab = "number of convergent clusters", xlab = "chromosome", cex.lab=1.5, cex.axis=1.2,width=c(0.08,0.08,0.08,0.08,0.08,0.08),cex = 0.9,las=0,col=brewer.pal(n = 6, name = "RdBu"))

legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

dev.off()

############################################################
############################################################
####### two sided plots (total size and number)#######
############################################################
############################################################

rm(list = ls())
## phenotype uncorrected
## two sided bar plot (size and number of convergent clusters)

null_phenUnocorr<-read.table(file = "out_res_convergent_CLUSTERS_2binning_Summed_Size_number_per_chromosome_Phenotype_GWAS_Uncorrected_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]


all_comparisonSize   <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     #### overlapped fixed (union)
#all_comparisonSize   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])    ### with overlaps


colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)


all_comparisonCount   <- cbind(ann_arg[,8], ann_pf[,8],ann_pp[,8],arg_pf[,8],arg_pp[,8],pf_pp[,8])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)


pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_Phenotype_Uncorrceted_Union_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")  #ylim=c(-2000000000,2000000000) ## withoverlap


barplot(-all_comparisonCount*500000,  #-all_comparisonCount*2000000 with overlap
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-200000000,200000000,50000000),
     labels = c(rev(seq(0,400,100)), seq(50000000,200000000,50000000)),las = 0,cex.axis=0.9)


# now add y-axis labels
mtext("size of convergent clusters", 2, line = 3, at = 100000000,cex.lab=3)
mtext("number of convergent clusters", 2, line = 3, at = -100000000,cex.lab=3)


legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

dev.off()

############################################################
#### two-sided phenotype corrected GWAS (union)
rm(list = ls())

null_phenUnocorr<-read.table(file = "out_res_convergent_CLUSTERS_2binning_Summed_Size_number_per_chromosome_Phenotype_GWAS_Corrected_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]


all_comparisonSize   <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     #### overlapped fixed (union)
#all_comparisonSize   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])    ### with overlaps


colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)



all_comparisonCount   <- cbind(ann_arg[,8], ann_pf[,8],ann_pp[,8],arg_pf[,8],arg_pp[,8],pf_pp[,8])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)


pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_Phenotype_Corrceted_Union_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")


barplot(-all_comparisonCount*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-200000000,200000000,50000000),
     labels = c(rev(seq(0,400,100)), seq(50000000,200000000,50000000)),las = 0,cex.axis=0.9)


# now add y-axis labels
mtext("size of convergent clusters", 2, line = 3, at = 100000000,cex.lab=3)
mtext("number of convergent clusters", 2, line = 3, at = -100000000,cex.lab=3)


legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.falax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

dev.off()

###################################

#### two-sided phenotype corrected (withoverlap)
rm(list = ls())

null_phenUnocorr<-read.table(file = "out_res_convergent_CLUSTERS_2binning_Summed_Size_number_per_chromosome_Phenotype_GWAS_corrected_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]


#all_comparisonSize   <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     #### overlapped fixed (union)
all_comparisonSize   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])    ### with overlaps


colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)



all_comparisonCount   <- cbind(ann_arg[,8], ann_pf[,8],ann_pp[,8],arg_pf[,8],arg_pp[,8],pf_pp[,8])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)


pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_Phenotype_Corrceted_with_overlap_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-2000000000,2000000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")


barplot(-all_comparisonCount*2000000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-2000000000,2000000000,500000000),
     labels = c(rev(seq(0,1200,300)), seq(500000000,2000000000,500000000)),las = 0,cex.axis=0.9)


# now add y-axis labels
mtext("size of convergent clusters", 2, line = 3, at = 900000000,cex.lab=3)
mtext("number of convergent clusters", 2, line = 3, at = -800000000,cex.lab=3)


legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.falax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

dev.off()


######################################
######################################
###### Climate ######
#### two-sided climate corrected (with overlap)

rm(list = ls())
null_phenUnocorr<-read.table(file = "out_res_convergent_CLUSTERS_2binning_Summed_Size_number_per_chromosome_Climate_BayPass_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]


all_comparisonSize   <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     #### overlapped fixed (union)
#all_comparisonSize   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])    ### with overlaps


colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)



all_comparisonCount   <- cbind(ann_arg[,8], ann_pf[,8],ann_pp[,8],arg_pf[,8],arg_pp[,8],pf_pp[,8])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)




pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_ClimateSoil_BayPass_withOverlaps_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-400000000,400000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")


barplot(-all_comparisonCount*2000000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-400000000,400000000,100000000),
     labels = c(rev(seq(0,400,100)), seq(100000000,400000000,100000000)),las = 0,cex.axis=0.9)


# now add y-axis labels
mtext("size of convergent clusters", 2, line = 3, at = 200000000,cex.lab=3)
mtext("number of convergent clusters", 2, line = 3, at = -200000000,cex.lab=3)


legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.falax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

dev.off()


###### with union
all_comparisonSize   <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     #### overlapped fixed (union)
   ### with overlaps


colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)



all_comparisonCount   <- cbind(ann_arg[,8], ann_pf[,8],ann_pp[,8],arg_pf[,8],arg_pp[,8],pf_pp[,8])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)




pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_ClimateSoil_BayPass_AdjustedAxis_Union_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-100000000,100000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")


barplot(-all_comparisonCount*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-100000000,100000000,25000000),
     labels = c(rev(seq(0,200,50)), seq(25000000,100000000,25000000)),las = 0,cex.axis=0.9)

# now add y-axis labels
mtext("size of convergent clusters", 2, line = 3, at = 60000000,cex.lab=3)
mtext("number of convergent clusters", 2, line = 3, at = -50000000,cex.lab=3)


legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.falax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

dev.off()





#########################################
#### two-sided climate Uncorrected (with overlpa)

null_phenUnocorr<-read.table(file = "out_res_convergent_CLUSTERS_2binning_Summed_Size_number_per_chromosome_Climate_Spearman_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]


all_comparisonSize   <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     #### overlapped fixed (union)
#all_comparisonSize   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])    ### with overlaps


colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)



all_comparisonCount   <- cbind(ann_arg[,8], ann_pf[,8],ann_pp[,8],arg_pf[,8],arg_pp[,8],pf_pp[,8])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)



pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_ClimateSoil_Spearman_withOverlap_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")


barplot(-all_comparisonCount*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-200000000,200000000,50000000),
     labels = c(rev(seq(0,400,100)), seq(50000000,200000000,50000000)),las = 0,cex.axis=0.9)

# now add y-axis labels
mtext("size of convergent clusters", 2, line = 3, at = 100000000,cex.lab=3)
mtext("number of convergent clusters", 2, line = 3, at = -100000000,cex.lab=3)


legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

dev.off()




####### with Union 


all_comparisonSize   <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])    ### with UNion


colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)



all_comparisonCount   <- cbind(ann_arg[,8], ann_pf[,8],ann_pp[,8],arg_pf[,8],arg_pp[,8],pf_pp[,8])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)



pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_ClimateSoil_Spearman_Union_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")


barplot(-all_comparisonCount*200000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-200000000,200000000,50000000),
     labels = c(rev(seq(0,1200,300)), seq(50000000,200000000,50000000)),las = 0,cex.axis=0.9)

# now add y-axis labels
mtext("size of convergent clusters", 2, line = 3, at = 110000000,cex.lab=3)
mtext("number of convergent clusters", 2, line = 3, at = -80000000,cex.lab=3)


legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.8, cex = 0.8)

dev.off()





##############################################################################################################
#######################################################
#### two-sided climate and phenotype (Uncorrected) ####

setwd("~/Documents/input_data/tests_analysis_for_paper/null_w/clusters")
null_climSP<-read.table(file = "out_res_convergent_CLUSTERS_2binning_soil_climate_Spearman_unique_percomparison_perchromosome.table", header = TRUE)
null_climSP$test<-rep("spearman")
ann_arg <- null_climSP[which(null_climSP$direction=="Annuus_Argophyllus"),]
ann_pf <- null_climSP[which(null_climSP$direction=="Annuus_petfal"),]
ann_pp <- null_climSP[which(null_climSP$direction=="Annuus_petpet"),]
arg_pf <- null_climSP[which(null_climSP$direction=="Argophyllus_petfal"),]
arg_pp <- null_climSP[which(null_climSP$direction=="Argophyllus_petpet"),]
pf_pp <- null_climSP[which(null_climSP$direction=="petfal_petpet"),]


all_comparisonSize_sp   <- cbind(ann_arg[,1], ann_pf[,1],ann_pp[,1],arg_pf[,1],arg_pp[,1],pf_pp[,1])
colnames(all_comparisonSize_sp) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize_sp) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize_sp<-t(all_comparisonSize_sp)


null_phenUnocorr<-read.table(file = "out_res_convergent_CLUSTERS_2binning_phenotype_Uncorrected_unique_percomparison_perchromosome.table", header = TRUE)
null_phenUnocorr$test<-rep("gwas_uncorrected")
ann_arg_gwas <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf_gwas <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp_gwas <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf_gwas <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp_gwas <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp_gwas <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]


all_comparisonSize_gwasun   <- cbind(ann_arg_gwas[,1], ann_pf_gwas[,1],ann_pp_gwas[,1],arg_pf_gwas[,1],arg_pp_gwas[,1],pf_pp_gwas[,1])
colnames(all_comparisonSize_gwasun) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize_gwasun) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize_gwasun<-t(all_comparisonSize_gwasun)


### plot results 
pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_Spearman_GWASUncorrected_size_number.pdf",width = 8, height= 8)
par(mar = c(5.7,4.5,4.1,2.1))
barplot(all_comparisonSize_gwasun, 
        
        beside=TRUE, ylim=c(-1200000000,1200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")


barplot(-all_comparisonSize_sp, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-1200000000,1200000000,600000000),
     labels = c(rev(seq(0,1200000000,600000000)), seq(600000000,1200000000,600000000)),las = 0,cex.axis=1.3)


# now add y-axis labels
mtext("size of cluster (phenotype uncorrected)", 2, line = 3, at = 700000000,cex.lab=3)
mtext("size of cluster (climate uncorrected)", 2, line = 3, at = -700000000,cex.lab=3)

legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.falax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.7, cex = 0.7)

dev.off()

###################################################################
#######################################################
#### two-sided climate (uncorercted) and phenotype (ncorrected) ####

setwd("~/Documents/input_data/tests_analysis_for_paper/null_w/clusters")
null_climSP<-read.table(file = "out_res_convergent_CLUSTERS_2binning_soil_climate_Spearman_unique_percomparison_perchromosome.table", header = TRUE)
null_climSP$test<-rep("spearman")
ann_arg <- null_climSP[which(null_climSP$direction=="Annuus_Argophyllus"),]
ann_pf <- null_climSP[which(null_climSP$direction=="Annuus_petfal"),]
ann_pp <- null_climSP[which(null_climSP$direction=="Annuus_petpet"),]
arg_pf <- null_climSP[which(null_climSP$direction=="Argophyllus_petfal"),]
arg_pp <- null_climSP[which(null_climSP$direction=="Argophyllus_petpet"),]
pf_pp <- null_climSP[which(null_climSP$direction=="petfal_petpet"),]


all_comparisonSize_sp   <- cbind(ann_arg[,1], ann_pf[,1],ann_pp[,1],arg_pf[,1],arg_pp[,1],pf_pp[,1])
colnames(all_comparisonSize_sp) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize_sp) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize_sp<-t(all_comparisonSize_sp)


null_phencorr<-read.table(file = "out_res_convergent_CLUSTERS_2binning_phenotype_Corrected_unique_percomparison_perchromosome.table", header = TRUE)
null_phencorr$test<-rep("gwas_uncorrected")
ann_arg_gwas <- null_phencorr[which(null_phencorr$direction=="Annuus_Argophyllus"),]
ann_pf_gwas <- null_phencorr[which(null_phencorr$direction=="Annuus_petfal"),]
ann_pp_gwas <- null_phencorr[which(null_phencorr$direction=="Annuus_petpet"),]
arg_pf_gwas <- null_phencorr[which(null_phencorr$direction=="Argophyllus_petfal"),]
arg_pp_gwas <- null_phencorr[which(null_phencorr$direction=="Argophyllus_petpet"),]
pf_pp_gwas <- null_phencorr[which(null_phencorr$direction=="petfal_petpet"),]


all_comparisonSize_gwascorr   <- cbind(ann_arg_gwas[,1], ann_pf_gwas[,1],ann_pp_gwas[,1],arg_pf_gwas[,1],arg_pp_gwas[,1],pf_pp_gwas[,1])
colnames(all_comparisonSize_gwascorr) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize_gwascorr) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize_gwascorr <-t(all_comparisonSize_gwascorr )


### plot results 
pdf (file = "Two_barplots_convergent_clusters_perChrom_recomadjusted_Spearman_GWASCorrected_size_number.pdf",width = 8, height= 8)
par(mar = c(5.7,4.5,4.1,2.1))
barplot(all_comparisonSize_gwascorr, 
        
        beside=TRUE, ylim=c(-1200000000,1200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")


barplot(-all_comparisonSize_sp, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

axis(2, lwd = 1,  at = seq(-1200000000,1200000000,600000000),
     labels = c(rev(seq(0,1200000000,600000000)), seq(600000000,1200000000,600000000)),las = 0,cex.axis=1.3)


# now add y-axis labels
mtext("size of cluster (phenotype corrected)", 2, line = 3, at = 700000000,cex.lab=3)
mtext("size of cluster (climate uncorrected)", 2, line = 3, at = -700000000,cex.lab=3)

legend("topleft", c("H.annuus-H.argophyllus","H.annuus-H.pet.falax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris"), pch=15, 
       col=brewer.pal(n = 6, name = "RdBu"), 
       bty="n",y.intersp= 0.7, cex = 0.7)

dev.off()