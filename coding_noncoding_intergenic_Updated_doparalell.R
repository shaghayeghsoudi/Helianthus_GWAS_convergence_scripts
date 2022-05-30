#!/usr/bin/env Rscript


#### find overlaps between all 5K  windows and coding-noncoding-intergenic regions. Script use parallel mode using "foreach" and "doParalell" packages in R to accelerate the process

## Load required packages
library(stringr)
library(GenomicRanges)
library(foreach)
library(doParallel)


## Laod gff file and adjust required formats
gff<-read.table(file = "HAN412_Eugene_curated_v1_1.gff3", header= FALSE)[,c(1:5,9)]
gff_good<-gff[!grepl("Ha412HOChr00",gff$V1),]

gff_good$Vn<-str_split_fixed(gff_good$V9, ";", 3)[,2]
gff_good$Vq<-str_split_fixed(gff_good$Vn, ":", 2)[,2]

gff_goodfin<-gff_good[,c("V1","V3","V4","V5","Vq")]
colnames(gff_goodfin)<-c("chrom","position_id","start","end","ensemble_gene_id")

gff_exon<-gff_goodfin[gff_goodfin$position_id=="CDS",]
gff_non<-gff_goodfin[(gff_goodfin$position_id=="five_prime_UTR" |gff_goodfin$position_id=="CDS"| gff_goodfin$position_id=="three_prime_UTR"), ]
genes<-unique(gff_non$ensemble_gene_id)


curated_bh_Shaghayegh<-read.table(file = "HAN412_Eugene_curated_v1_1_introns_curated_bh_Shaghayegh.txt", header = TRUE)
out_res<-curated_bh_Shaghayegh


out_res_good_intron<-out_res[out_res$start!=1,]
out_res_good_intron<-out_res_good_intron[,c("seqnames","position_id","start","end","ensemble_gene_id")]
colnames(out_res_good_intron)<-c("chrom","position_id","start","end", "ensemble_gene_id")

gff_non_5utrs<-gff_non[gff_non$position_id=="five_prime_UTR", ]
gff_non_5utrs$start<-(gff_non_5utrs$start-500)

gff_non_3utrs<-gff_non[gff_non$position_id=="three_prime_UTR", ]
gff_non_3utrs$end<-(gff_non_3utrs$end+500)

gff_CDS<-gff_non[gff_non$position_id=="CDS", ]
allnoncoding<-rbind(out_res_good_intron,gff_non_5utrs,gff_non_3utrs)

## 5k nonoverlapping windows covering the whole genome
windows_5k<-read.table(file = "bins_good_5000.txt", header = FALSE)[,c(1:4)]
#windows_5k<-read.table(file = "bins_good_5000_random10K.txt", header = FALSE)[,c(1:4)]
colnames(windows_5k)<-c("window_id","chrom","start","end")
windows_5k$window_id<-paste("window",windows_5k$window_id,sep = "_")
windows_5k$chrom<-gsub("chr","",windows_5k$chrom)

windows_5k$chrom<-sprintf("%02d", as.numeric(as.character(windows_5k$chrom)))
windows_5k$chrom<-paste("Ha412HOChr",windows_5k$chrom, sep = "")
#windows<-unique(windows_5k$window_id)


parallel::detectCores()
n.cores <- parallel::detectCores() - 27

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
    out_res_c<-NULL
    
    for (i in 1:length(windows)){
 
    focal_window<-focal_chroms_window[focal_chroms_window$window_id==windows[i],]
    focal_window_gr<-GRanges(IRanges(start=focal_window$start , end = focal_window$end), seqnames = focal_window$chrom)
    
    gff_CDS_gr<-GRanges(IRanges(start=gff_CDS$start , end = gff_CDS$end), seqnames = gff_CDS$chrom)
    allnoncoding_gr<-GRanges(IRanges(start=allnoncoding$start , end = allnoncoding$end), seqnames = allnoncoding$chrom)
    
    
    focal_intersect_coding<-data.frame(intersect(focal_window_gr,gff_CDS_gr))
    focal_intersect_noncoding<-data.frame(intersect(focal_window_gr,allnoncoding_gr))

    
    if (nrow(focal_intersect_coding) == 0 &&  nrow(focal_intersect_noncoding) == 0){
        focal_window$coding_noncoding_status<-"intergenic"

    } 
    
    if (nrow(focal_intersect_coding) >= 1 &&  nrow(focal_intersect_noncoding) == 0) {
        focal_window$coding_noncoding_status<-"coding"
    } 
    
    if (nrow(focal_intersect_coding) == 0 &&  nrow(focal_intersect_noncoding) >= 1) {
        focal_window$coding_noncoding_status<-"noncoding"
    }
    
    if (nrow(focal_intersect_coding) >= 1 &&  nrow(focal_intersect_noncoding) >= 1) {
        focal_window$coding_noncoding_status<-"coding_noncoding"
    }
   out_res_c<-rbind(focal_window,out_res_c)
    }

   write.table(out_res_c, file = paste("overlapps_5Kwindows_with_coding_noncoding_intergenic_regions_",chroms[j],".table",sep = "" ),col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

}


parallel::stopCluster(cl = my.cluster)
