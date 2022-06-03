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


#### check picmin windows overlap with genic or intergenic regions


picmin_unclust<-read.delim(file = "picmin_fdr20_results_unclustered_4way_sh.txt", header = FALSE)
picmin_unclust$V9<-gsub("llevation" ,"Elevation",picmin_unclust$V9 )
picmin_unclust_selected_cols<-picmin_unclust[,c(5,8,9)]
picmin_unclust_selected_cols<-picmin_unclust_selected_cols[!duplicated(picmin_unclust_selected_cols$V8),]
colnames(picmin_unclust_selected_cols)<-c("chrom","window","variable")

### load 5k windows with overlap information

files<-list.files("overlapps_5Kwindows_with_coding_noncoding_intergenic_regions", pattern = "overlapps_5Kwindows_with_coding_noncoding_intergenic_regions_*" , full.names= TRUE)

windows<-lapply(files, function(x) {
    read.table(x, header = TRUE, sep = "\t")
})

aa <- do.call("rbind", windows) 

aa$coding_noncoding_status<-gsub("coding_noncoding", "genic",aa$coding_noncoding_status)
aa$coding_noncoding_status<-gsub("noncoding","genic",aa$coding_noncoding_status)
aa$coding_noncoding_status<-gsub("coding","genic",aa$coding_noncoding_status)

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


chroms<-unique(aa$chrom)


foreach(j = 1:length(chroms))%dopar%{   


    #library(GenomicRanges)
  
    focal_chroms_window<-aa[aa$chrom==chroms[j],]
    windows<-unique(focal_chroms_window$window_id)
    out_res_n<-NULL

    for (k in 1:length(windows)){

     #for (k in 1:100){
     focal_window_genic_intergenic<-focal_chroms_window[k,]

        if (focal_window_genic_intergenic$window_id %in% picmin_unclust_selected_cols$window){

        focal_window_genic_intergenic$picmin_status<-"picmin_window"
        } else {

        focal_window_genic_intergenic$picmin_status<-"Not_picmin_window"
        }

    out_res_n<-rbind(focal_window_genic_intergenic,out_res_n)
    }


write.table(out_res_n, file = paste("chi/5k_windows_with_genic_intergenic_status_checked_picmin_or_notpicmin_",chroms[j,".table", sep = ""),col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

### perform chisq.test 
#table(out_res_n$coding_noncoding_status, out_res_n$picmin_status)
#chitest<-data.frame(chisq.test(out_res_n$coding_noncoding_status, out_res_n$picmin_status,simulate.p.value = TRUE)[c("statistic","parameter","p.value")])

#write.table(chitest, "out_resk_chisq.test_5k_windows_with_genic_intergenic_status_association_with_picmin_or_notpicmin.table", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


wra_files<-list.files("chi", pattern = "5k_windows_with_genic_intergenic_status_checked_picmin_or_notpicmin_*", full.names = TRUE)


attackStats <- lapply(wra_files,function(x) {
    read.delim(x,  header=TRUE, sep = "\t")
    })

out_res_n <- do.call("rbind", attackStats) 
write.table(out_res_n ,file = "5k_windows_with_genic_intergenic_status_checked_picmin_or_notpicmin.table", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

### perform chisq.test 
table(out_res_n$coding_noncoding_status, out_res_n$picmin_status)
chitest<-data.frame(chisq.test(out_res_n$coding_noncoding_status, out_res_n$picmin_status,simulate.p.value = TRUE)[c("statistic","parameter","p.value")])

write.table(chitest, "out_resk_chisq.test_5k_windows_with_genic_intergenic_status_association_with_picmin_or_notpicmin.table", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



##
