#!/usr/bin/env Rscript

#### find overlaps between all 5K  windows and genic-intergenic regions. Script use parallel mode


#rm(list = ls())

## with chi-square test- overlap with genic intergenic 
library(stringr)
library(GenomicRanges)
library(foreach)
library(doParallel)

#setwd("~/Dropbox/input_data/coding_noncoding_overlap_picmin")

windows_5k<-read.table(file = "bins_good_5000.txt", header = FALSE)[,c(1:4)]

colnames(windows_5k)<-c("window_id","chrom","start","end")
windows_5k$window_id<-paste("window",windows_5k$window_id,sep = "_")
windows_5k$chrom<-gsub("chr","",windows_5k$chrom)

windows_5k$chrom<-sprintf("%02d", as.numeric(as.character(windows_5k$chrom)))
windows_5k$chrom<-paste("Ha412HOChr",windows_5k$chrom, sep = "")



gff<-read.table(file = "HAN412_Eugene_curated_v1_1.gff3", header= FALSE)[,c(1:5,9)]
gff_good<-gff[!grepl("Ha412HOChr00",gff$V1),]

gff_good$Vn<-str_split_fixed(gff_good$V9, ";", 3)[,2]
gff_good$Vq<-str_split_fixed(gff_good$Vn, ":", 2)[,2]

gff_goodfin<-gff_good[,c("V1","V3","V4","V5","Vq")]
colnames(gff_goodfin)<-c("chrom","position_id","start","end","ensemble_gene_id")


## subset genic regions
gff_goodfin<-gff_goodfin[gff_goodfin$position_id=="gene",]
gff_genic_gr<-GRanges(IRanges(start = gff_goodfin$start,end =  gff_goodfin$end),seqnames = gff_goodfin$chrom, ensemble_gene_id=gff_goodfin$ensemble_gene_id)
gff_intergenic_gr<-gaps(gff_genic_gr)
windows<-unique(windows_5k$window_id)

parallel::detectCores()
n.cores <- parallel::detectCores() - 30

### **create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(my.cluster)

## **register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()


chroms<-unique(windows_5k$chrom)

foreach(j = 1:length(chroms))%dopar%{   

    library(GenomicRanges) 
    focal_chroms_window<-windows_5k[windows_5k$chrom==chroms[j],]
    windows<-unique(focal_chroms_window$window_id)
    

    out_res_b<-NULL
    for (i in 1:length(windows)){
   
    #for (i in 1:200){
    focal_window<-focal_chroms_window[focal_chroms_window$window_id==windows[i],]
    focal_window_gr<-GRanges(IRanges(start=focal_window$start , end = focal_window$end), seqnames = focal_window$chrom)
    focal_intersect_genic<-data.frame(intersect(focal_window_gr,gff_genic_gr))
    focal_intersect_intergenic<-data.frame(intersect(focal_window_gr,gff_intergenic_gr))

    
    if (nrow(focal_intersect_genic) >=1  && nrow(focal_intersect_intergenic) == 0){
        focal_window$genic_nongenis_status<-"genic"
    } 
    
    if (nrow(focal_intersect_genic) == 0 && nrow(focal_intersect_intergenic) >= 1) {
        focal_window$genic_nongenis_status<-"intergenic"
    }

    if (nrow(focal_intersect_genic) >=1 && nrow(focal_intersect_intergenic) >= 1) {
        focal_window$genic_nongenis_status<-"genic_intergenic"
    }
   
   out_res_b<-rbind(focal_window,out_res_b)
  }

write.table(out_res_b, file = paste("outputs/overlapps_5Kwindows_with_genic_intergenic_regions_",chroms[j],".table",sep = "" ),col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)
#write.table(out_res_b, file = "overlapps_5Kwindows_with_genic_intergenic_regions_Updated.table", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

}


parallel::stopCluster(cl = my.cluster)
