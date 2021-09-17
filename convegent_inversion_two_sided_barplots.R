
## This script calulates what proportion of convergent clusters overlap with inversions and plot them as two-dided bar plots ##
## load required libraries
library(GenomicRanges)
library(RColorBrewer)
library("magrittr")
library("dplyr")

### useful sources
# https://www.biostars.org/p/196393/
# https://stackoverflow.com/questions/30079720/union-and-intersection-of-intervals
# genomicRanges tutorial: 
#https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
#https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html


setwd("~/Documents/input_data/tests_analysis_for_paper/inversion_convergence")
con_inv<-read.table(file = "convergent_cluster_inversion_overlap_LD0.9_1cM_ranges_verbose", header = TRUE)

#########################################

chroms<-c("Ha412HOChr01","Ha412HOChr02","Ha412HOChr03","Ha412HOChr04","Ha412HOChr05","Ha412HOChr06","Ha412HOChr07","Ha412HOChr08","Ha412HOChr09","Ha412HOChr10","Ha412HOChr11","Ha412HOChr12","Ha412HOChr13","Ha412HOChr14","Ha412HOChr15","Ha412HOChr16","Ha412HOChr17")
compars<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","Argophyllus_petpet","petfal_petpet")

data_baypass<-con_inv[con_inv$analysis=="baypass",]
data_spearman<-con_inv[con_inv$analysis=="spearman",]

data_gwas_corr<-con_inv[con_inv$analysis=="GWAS_corrected",]
data_gwas_Uncorr<-con_inv[con_inv$analysis=="GWAS_uncorrected",]



out_res<-NULL
for (i in 1:length(compars)){
  
  for (j in 1:length(chroms)){
    data_focal<-data_baypass[data_baypass$comparison==compars[i] & data_baypass$chromosome==chroms[j] , ]
    
    #### get total cluster size ####
    gr_cluster <- GRanges(seqnames = rep(1,nrow(data_focal)),IRanges(start = data_focal$cluster_start,end =      data_focal$cluster_end))
    
    gr_cluster_unique<-gr_cluster %>% unique
    union<-as.data.frame(reduce(gr_cluster_unique))[,2:3]
    union$UniqueclusterSize<-(union$end -union$start)
    total_cluster_size<-data.frame(t(colSums(union)))
    colnames(total_cluster_size)<-c("cluster_start","cluster_end","total_cluster_length")
    total_cluster_number<-data.frame("total_N.cluster"=nrow(data.frame(gr_cluster_unique)))
    all_clusters<-cbind(total_cluster_size,total_cluster_number)
    
    ### find total range of overlaps between unique cluster and total inversions
    inversion<-na.omit(data_focal)
    gr_inversion <- GRanges(seqnames = rep(1,nrow(inversion)),IRanges(start = inversion$inversion_start,end =      inversion$inversion_end))
    CI_intersect<-data.frame(intersect(gr_cluster_unique,gr_inversion))[,2:4]
    total_overlap<-data.frame(t(colSums(CI_intersect)))
    total_overlap_size<-data.frame("total_overlap_size"=total_overlap[,3])
    
    ### unique number of clusters overlapping with inversion ###
    gr_cluster_uniqueNumber <- data.frame(GRanges(seqnames = rep(1,nrow(inversion)),IRanges(start = inversion$cluster_start,end =        inversion$cluster_end)))
    total_number_overlap<-data.frame("total_N.overlapped"=nrow(gr_cluster_uniqueNumber%>%unique))
    all_overlaps<-cbind(total_overlap_size,total_number_overlap)
    
    sub_good<-cbind(all_clusters,all_overlaps)
    
    sub_good$direction <- compars[i]
    sub_good$chromosome <- chroms[j]
    out_res <- rbind (out_res,sub_good)
    
    
    #outliers_count <- tapply (data_baypass_focal$N_convergent, list(as.character (data_baypass_focal$chromosome=="chroms[j]")),sum)
    
  }
}


write.table(out_res,file = "out_res_Inversion_Convergent_CLUSTERS_overlapping_Summed_Size_number_per_chromosome_Climate_BayPass_union_percomparison.table", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


######################################################################################################################
### PLOT THE RESULTS ####
###### two-sided phenotype Uncorrected GWAS (union) #######
#rm(list = ls())

null_phenUnocorr<-read.table(file = "out_res_Inversion_Convergent_CLUSTERS_overlapping_Summed_Size_number_per_chromosome_Phenotype_GWAS_Uncorrected_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]

## size 
all_comparisonSize  <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     
all_comparisonSize_overlap<-cbind(ann_arg[,5], ann_pf[,5],ann_pp[,5],arg_pf[,5],arg_pp[,5],pf_pp[,5])     

colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)

colnames(all_comparisonSize_overlap) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize_overlap) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize_overlap<-t(all_comparisonSize_overlap)


## count
all_comparisonCount   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)


all_comparisonCount_overlap   <- cbind(ann_arg[,6], ann_pf[,6],ann_pp[,6],arg_pf[,6],arg_pp[,6],pf_pp[,6])
colnames(all_comparisonCount_overlap) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount_overlap) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount_overlap<-t(all_comparisonCount_overlap)





pdf (file = "Two_barplots_CI_convergent_clusters_Inversions_perChrom_recomadjusted_Phenotype_UnCorrceted_Union_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")

par (new = TRUE)
barplot(all_comparisonSize_overlap, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),axes = FALSE, xlab= "Chromosome",xaxt = "n", yaxt = "n", bty = "n", angle = 45, density = 40, col = "black")



barplot(-all_comparisonCount*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

par (new = TRUE)
barplot(-all_comparisonCount_overlap*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08), add= TRUE,axes = FALSE,xaxt = "n", yaxt = "n", bty = "n", angle = 45, density = 40, col = "black")



axis(2, lwd = 1,  at = seq(-200000000,200000000,50000000),
     labels = c(rev(seq(0,400,100)), seq(50000000,200000000,50000000)),las = 0,cex.axis=0.9)


# now add y-axis labels
mtext("Total size convergent clusters (bp)", 2, line = 3, at = 100000000,cex.lab=3)
mtext("Total number convergent clusters", 2, line = 3, at = -100000000,cex.lab=3)


legend("topleft", 
legend= c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris","overlaps with inversion"),  
       #col=brewer.pal(n = 6, name = "RdBu"), 
       
       fill = c("#B2182B","#EF8A62","#FDDBC7","#D1E5F0","#67A9CF","#2166AC","black"),density=c(NA, NA,NA,NA,NA,NA,40),bty="n",
       y.intersp= 0.8, cex = 0.8)

dev.off()


###########################################################
###########################################################
#### Corrected GWAS ####

#rm(list = ls())

null_phenUnocorr<-read.table(file = "out_res_Inversion_Convergent_CLUSTERS_overlapping_Summed_Size_number_per_chromosome_Phenotype_GWAS_Corrected_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]

## size 
all_comparisonSize  <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     
all_comparisonSize_overlap<-cbind(ann_arg[,5], ann_pf[,5],ann_pp[,5],arg_pf[,5],arg_pp[,5],pf_pp[,5])     

colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)

colnames(all_comparisonSize_overlap) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize_overlap) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize_overlap<-t(all_comparisonSize_overlap)


## count
all_comparisonCount   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)


all_comparisonCount_overlap   <- cbind(ann_arg[,6], ann_pf[,6],ann_pp[,6],arg_pf[,6],arg_pp[,6],pf_pp[,6])
colnames(all_comparisonCount_overlap) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount_overlap) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount_overlap<-t(all_comparisonCount_overlap)





pdf (file = "Two_barplots_CI_convergent_clusters_Inversions_perChrom_recomadjusted_Phenotype_Corrceted_Union_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")

par (new = TRUE)
barplot(all_comparisonSize_overlap, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),axes = FALSE, xlab= "Chromosome",xaxt = "n", yaxt = "n", bty = "n", angle = 45, density = 40, col = "black")



barplot(-all_comparisonCount*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

par (new = TRUE)
barplot(-all_comparisonCount_overlap*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08), add= TRUE,axes = FALSE,xaxt = "n", yaxt = "n", bty = "n", angle = 45, density = 40, col = "black")



axis(2, lwd = 1,  at = seq(-200000000,200000000,50000000),
     labels = c(rev(seq(0,400,100)), seq(50000000,200000000,50000000)),las = 0,cex.axis=0.9)


# now add y-axis labels
mtext("Total size convergent clusters (bp)", 2, line = 3, at = 100000000,cex.lab=3)
mtext("Total number convergent clusters", 2, line = 3, at = -100000000,cex.lab=3)


legend("topleft", 
legend= c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris","overlaps with inversion"),  
       #col=brewer.pal(n = 6, name = "RdBu"), 
       
       fill = c("#B2182B","#EF8A62","#FDDBC7","#D1E5F0","#67A9CF","#2166AC","black"),density=c(NA, NA,NA,NA,NA,NA,40),bty="n",
       y.intersp= 0.8, cex = 0.8)

dev.off()



###########################################################
###########################################################
#### climate and soil (Spearman)

#rm(list = ls())

null_phenUnocorr<-read.table(file = "out_res_Inversion_Convergent_CLUSTERS_overlapping_Summed_Size_number_per_chromosome_Climate_Spearman_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]

## size 
all_comparisonSize  <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     
all_comparisonSize_overlap<-cbind(ann_arg[,5], ann_pf[,5],ann_pp[,5],arg_pf[,5],arg_pp[,5],pf_pp[,5])     

colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)

colnames(all_comparisonSize_overlap) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize_overlap) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize_overlap<-t(all_comparisonSize_overlap)


## count
all_comparisonCount   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)


all_comparisonCount_overlap   <- cbind(ann_arg[,6], ann_pf[,6],ann_pp[,6],arg_pf[,6],arg_pp[,6],pf_pp[,6])
colnames(all_comparisonCount_overlap) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount_overlap) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount_overlap<-t(all_comparisonCount_overlap)





pdf (file = "Two_barplots_CI_convergent_clusters_Inversions_perChrom_recomadjusted_Climate_Spearman_Union_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")

par (new = TRUE)
barplot(all_comparisonSize_overlap, 
        
        beside=TRUE, ylim=c(-200000000,200000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),axes = FALSE, xlab= "Chromosome",xaxt = "n", yaxt = "n", bty = "n", angle = 45, density = 40, col = "black")



barplot(-all_comparisonCount*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

par (new = TRUE)
barplot(-all_comparisonCount_overlap*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08), add= TRUE,axes = FALSE,xaxt = "n", yaxt = "n", bty = "n", angle = 45, density = 40, col = "black")



axis(2, lwd = 1,  at = seq(-200000000,200000000,50000000),
     labels = c(rev(seq(0,400,100)), seq(50000000,200000000,50000000)),las = 0,cex.axis=0.9)


# now add y-axis labels
mtext("Total size convergent clusters (bp)", 2, line = 3, at = 100000000,cex.lab=3)
mtext("Total number convergent clusters", 2, line = 3, at = -100000000,cex.lab=3)


legend("topleft", 
legend= c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris","overlaps with inversion"),  
       #col=brewer.pal(n = 6, name = "RdBu"), 
       
       fill = c("#B2182B","#EF8A62","#FDDBC7","#D1E5F0","#67A9CF","#2166AC","black"),density=c(NA, NA,NA,NA,NA,NA,40),bty="n",
       y.intersp= 0.8, cex = 0.8)

dev.off()




###########################################################
###########################################################
#### climate and soil (Corrected_BayPass)

#rm(list = ls())

null_phenUnocorr<-read.table(file = "out_res_Inversion_Convergent_CLUSTERS_overlapping_Summed_Size_number_per_chromosome_Climate_BayPass_union_percomparison.table", header = TRUE)

ann_arg <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_Argophyllus"),]
ann_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petfal"),]
ann_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Annuus_petpet"),]
arg_pf <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petfal"),]
arg_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="Argophyllus_petpet"),]
pf_pp <- null_phenUnocorr[which(null_phenUnocorr$direction=="petfal_petpet"),]

## size 
all_comparisonSize  <- cbind(ann_arg[,3], ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])     
all_comparisonSize_overlap<-cbind(ann_arg[,5], ann_pf[,5],ann_pp[,5],arg_pf[,5],arg_pp[,5],pf_pp[,5])     

colnames(all_comparisonSize) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize<-t(all_comparisonSize)

colnames(all_comparisonSize_overlap) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonSize_overlap) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonSize_overlap<-t(all_comparisonSize_overlap)


## count
all_comparisonCount   <- cbind(ann_arg[,4], ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])
colnames(all_comparisonCount) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount<-t(all_comparisonCount)


all_comparisonCount_overlap   <- cbind(ann_arg[,6], ann_pf[,6],ann_pp[,6],arg_pf[,6],arg_pp[,6],pf_pp[,6])
colnames(all_comparisonCount_overlap) <- c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")
rownames(all_comparisonCount_overlap) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")
all_comparisonCount_overlap<-t(all_comparisonCount_overlap)



pdf (file = "Two_barplots_CI_convergent_clusters_Inversions_perChrom_recomadjusted_Climate_BayPass_AdjustedAxis_Union_size_number.pdf",width = 8, height= 8)
par(mar = c(5.1,4.5,4.1,2.1))
barplot(all_comparisonSize, 
        
        beside=TRUE, ylim=c(-100000000,100000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), axes = FALSE, xlab= "Chromosome")

par (new = TRUE)
barplot(all_comparisonSize_overlap, 
        
        beside=TRUE, ylim=c(-100000000,100000000),width=c(0.08,0.08,0.08,0.08,0.08,0.08),axes = FALSE, xlab= "Chromosome",xaxt = "n", yaxt = "n", bty = "n", angle = 45, density = 40, col = "black")



barplot(-all_comparisonCount*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08),col=brewer.pal(n = 6, name = "RdBu"), add= TRUE,axes = FALSE)

par (new = TRUE)
barplot(-all_comparisonCount_overlap*500000, 
        
        beside=TRUE, width=c(0.08,0.08,0.08,0.08,0.08,0.08), add= TRUE,axes = FALSE,xaxt = "n", yaxt = "n", bty = "n", angle = 45, density = 40, col = "black")



axis(2, lwd = 1,  at = seq(-100000000,100000000,25000000),
     labels = c(rev(seq(0,200,50)), seq(25000000,100000000,25000000)),las = 0,cex.axis=0.9)


# now add y-axis labels
mtext("Total size convergent clusters (bp)", 2, line = 3, at = 60000000,cex.lab=3)
mtext("Total number convergent clusters", 2, line = 3, at = -50000000,cex.lab=3)


legend("topleft", 
legend= c("H.annuus-H.argophyllus","H.annuus-H.pet.fallax","H.annuus-H.pet.petiolaris","H.argophyllus-H.pet.fallax","H.argophyllus-H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris","overlaps with inversion"),  
       #col=brewer.pal(n = 6, name = "RdBu"), 
       
       fill = c("#B2182B","#EF8A62","#FDDBC7","#D1E5F0","#67A9CF","#2166AC","black"),density=c(NA, NA,NA,NA,NA,NA,40),bty="n",
       y.intersp= 0.8, cex = 0.8)

dev.off()

