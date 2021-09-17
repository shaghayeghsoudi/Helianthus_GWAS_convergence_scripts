

#!/usr/bin/env Rscript
###############################################
### This script finds what proportion of convergent windows and control background windows  overlaps with CDS and regulatory regions 
### script has six inputs ###
###############################################


library(GenomicRanges)
library(dplyr)
setwd("/data/users/soudi/paper_sunflowers/coding_noncoding/inputs")


### gff files
gff_gene<-read.table(file = "gene.gff.table", header = FALSE)
gff_gene$V6<-gsub("ID=","",gff_gene$V6)
gff_gene$down<-rep("down")
gff_gene$up<-rep("up")

gff_gene$downStart<-(gff_gene$V4-500)
gff_gene$downEnd<-(gff_gene$V4-1)

gff_gene$upStart<-(gff_gene$V5+1)
gff_gene$upEnd<-(gff_gene$V5+500)

gff_down<-gff_gene[,c(1,2,7,9,10,6)]
colnames(gff_down)<-c("chrom","emv","genomic_position","start","end","gene_id")
gff_up<-gff_gene[,c(1,2,8,11,12,6)]
colnames(gff_up)<-c("chrom","emv","genomic_position","start","end","gene_id")



gff_cds_utr<-read.table(file = "cds_utr", header = FALSE)
gff_cds_utr$V6<-gsub("ID=three_prime_utr","",gff_cds_utr$V6) 
gff_cds_utr$V6<-gsub("ID=five_prime_utr","",gff_cds_utr$V6)
gff_cds_utr$V6<-gsub("ID=","",gff_cds_utr$V6)
gff_cds_utr$V6<-gsub("-RA","",gff_cds_utr$V6)
gff_cds_utr$V6<-gsub("-RB","",gff_cds_utr$V6)
colnames(gff_cds_utr)<-c("chrom","emv","genomic_position","start","end","gene_id")

### final cds and utr files 
cds<-gff_cds_utr[(!grepl("three_prime_UTR",gff_cds_utr$genomic_position)) & (!grepl("five_prime_UTR",gff_cds_utr$genomic_position)) & (!grepl("UTR",gff_cds_utr$genomic_position)) &(!grepl("mRNA",gff_cds_utr$genomic_position)) & (!grepl("gene",gff_cds_utr$genomic_position)),]

gffutr<-gff_cds_utr[(grepl("three_prime_UTR",gff_cds_utr$genomic_position)) | (grepl("five_prime_UTR",gff_cds_utr$genomic_position)) | (grepl("UTR",gff_cds_utr$genomic_position)),]

regulatory <- do.call("rbind", list(gffutr, gff_down, gff_up)) 


###############################################
###### load recombination rate bin file #######
#recom<-read.table(file ="recombination_rate_bins/out_res_5k_with_5quantile.table", header = TRUE)
recom<-read.table(file ="out_res_5k_with_5quantile.table", header = TRUE)

recom<-recom[!(duplicated(recom$window)),]


### load winow5k coordinate
#window5k<-read.table(file = "recombination_rate_bins/bins_good_5000", header = FALSE)
window5k<-read.table(file = "bins_good_5000", header = FALSE)


colnames(window5k)<-c("window_id","chrom","start","end","length","non")
window5k$window_id<-sub("^","window_",window5k$window_id)
window5k$chrom<-sub("chr","Ha412HOChr",window5k$chrom)


#### load 5k windows with assigned overlapped genes (output of bed intersect) ####
window5k_genes_bed<-read.table(file = "intersection_5k_genes.table", header = FALSE)
window5k_genes_bed<-window5k_genes_bed[,c(4,8)]
window5k_genes_bed$chrom<-substr(window5k_genes_bed$V8,0,12)
colnames(window5k_genes_bed)<-c("window_id","gene_id","chrom")

#assign bin ID and window start and end to 5kwindows
window5k_genes_bed<-merge(window5k_genes_bed,window5k, by.x = "window_id", by.y="window_id" )[,c(1,2,3,5,6)]
colnames(window5k_genes_bed)<-c("window_id",  "gene_id", "chrom", "start" , "end")
recom_window5k_coor<-merge(window5k_genes_bed, recom, by.x ="window_id" , by.y="window" )[,c(1:5,9)]   ### control windows wih assigned quantile


#### all comparisons convergent windows file #####
con_table<-read.table(file = "CONV_WIN_GENE_overlap", header = TRUE)
comparisons<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","petfal_petpet","Argophyllus_petpet")
#aa<-con_table[con_table$comparison== "Annuus_Argophyllus",]
#aa<-aa[!(duplicated(aa$gene_id)),]
#aa$chrom<-gsub(":.*$","",aa$window5)
#gene_id_foc<-unique(aa$gene_id)
assigned_recom_con<-merge(con_table, recom, by.x="window_name", by.y = "window")[,-c(16,17,18)]  ### convergent genes with assigned quantiles
assigned_recom_con<-assigned_recom_con[assigned_recom_con$gene_id%in%gff_gene$V6,]

#############################

comparisons<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","petfal_petpet","Argophyllus_petpet")


out_res<-NULL

for (ii in 1:length(comparisons)){	
	aa<-assigned_recom_con[(assigned_recom_con$comparison==comparisons[ii]),]
	aa_foc<-aa[(aa$analysis=="spearman")  | (aa$analysis=="GWAS_corrected"),]
	aa_foc<-aa_foc[!(duplicated(aa_foc$gene_id)),]
    #aa_foc$chrom<-gsub(":.*$","",aa_foc$window5)
    #gene_id_foc<-unique(aa_foc$gene_id)
    
	
	for (jj in 1:5){      ## loop through each 5 recombination bin quantiles
		
		assigned_recom_focal<-aa_foc[aa_foc$quantile==jj,]          ## convergent windows
		window_id_foc<-unique(assigned_recom_focal$window_name)
				
				
        for (zz in 1:length(window_id_foc)){   ### loop through each convergent window

			
         foc<-assigned_recom_focal[assigned_recom_focal$window_name==window_id_foc[zz],]   
         cds_focal<-cds[cds$gene_id==foc$gene_id,]
         reg_focal<-regulatory[regulatory$gene_id==foc$gene_id,]
    
    
         gr_aa <- GRanges(IRanges(start = foc$window_start,end =      foc$window_end),seqnames = foc$Chr)
         gr_cds<-GRanges(IRanges(start = cds_focal$start,end =      cds_focal$end),seqnames = cds_focal$chrom,pos= cds_focal$genomic_position )
         gr_regulatory<-GRanges(IRanges(start = reg_focal$start,end =      reg_focal$end),seqnames = reg_focal$chrom,pos=     reg_focal$genomic_position )
         coding_overlap<-data.frame(intersect(gr_aa, gr_cds))[,1:4]
  

  
      proporlength<-sum(coding_overlap$width)/5000
      length<-sum(coding_overlap$width)
      window_id<-as.character(window_id_foc[zz])
      genomic_pos<-"CDS"
      recombin<-jj
      type<-comparisons[ii]     
      status<-"convergent" 
      coding_overlap_con<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type, status))
  
     
      reg_overlap<-data.frame(intersect(gr_aa,gr_regulatory))[,1:4]        ### overlaping of convergent windows with regulatory regions
      proporlength<-sum(reg_overlap$width)/5000
      length<-sum(reg_overlap$width)
      window_id<-as.character(window_id_foc[zz])
      genomic_pos<-"regulatory"
      recombin<-jj
      type<-comparisons[ii]   
      status<-"convergent"   
      reg_overlap_con<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type, status))
      


cds_reg<-do.call("rbind", list(coding_overlap_con,reg_overlap_con)) 
out_res <- rbind (out_res,cds_reg)

  } ### zz loop (convergent window)
  
  
   
     ###########################################################
     #### control random windows for the NULL distribution #####
     recom_window5k_coor_foc<-recom_window5k_coor[recom_window5k_coor$quantile==jj,]
     qq<-recom_window5k_coor_foc[!(recom_window5k_coor_foc$window_id %in% assigned_recom_focal$window_name),] 
   
     control_window<-unique(qq$window_id)
 
 
 for (rr in 1:length(control_window)){   ## loop through each random control window
	
 	
    control_gene_foc<-qq[qq$window_id==control_window[rr],]
    cds_focal_control<-cds[cds$gene_id==control_gene_foc$gene_id ,]
    reg_focal_control<-regulatory[regulatory$gene_id==control_gene_foc$gene_id ,]
    
    
    gr_ww <- GRanges(IRanges(start = control_gene_foc$start,end =      control_gene_foc$end),seqnames = control_gene_foc$chrom)
    gr_cds_ww<-GRanges(IRanges(start = cds_focal_control$start,end =      cds_focal_control$end),seqnames = cds_focal_control$chrom,pos= cds_focal_control$genomic_position )
    gr_regulatory_ww<-GRanges(IRanges(start = reg_focal_control$start,end =      reg_focal_control$end),seqnames = reg_focal_control$chrom,pos=     reg_focal_control$genomic_position )
    coding_overlap_ww<-data.frame(intersect(gr_ww, gr_cds_ww))[,1:4]
    #reg_overlap_ww<-data.frame(intersect(gr_ww, gr_regulatory_ww))[,1:4]

   
  proporlength<-sum(coding_overlap_ww$width)/5000
  length<-sum(coding_overlap_ww$width)
  window_id<-as.character(control_window[rr])
  genomic_pos<-"CDS"
  recombin<-jj
  type<-comparisons[ii] 
  status<-"random_control"  
  coding_overlap_null<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type,status))
  
   
  ### regulatory random windows
  reg_overlap_ww<-data.frame(intersect(gr_ww, gr_regulatory_ww))[,1:4]
  

  proporlength<-sum(reg_overlap_ww$width)/5000
  length<-sum(reg_overlap_ww$width)
  window_id<-as.character(control_window[rr])
  genomic_pos<-"regulatory"
  recombin<-jj
  type<-comparisons[ii] 
  status<-"random_control"  
  reg_overlap_null<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type,status))
  
cds_reg<-do.call("rbind", list(coding_overlap_null,reg_overlap_null)) 
out_res <- rbind (out_res,cds_reg)
  
         } ## rr loop (randomcontrol window)
 
     } ## jj loop
 
} ## ii loop

write.table(out_res, file = "test_cds_reg_script.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


#########
#rm(list=ls())
setwd("~/Documents/input_data/tests_analysis_for_paper/coding_noncoding")
cds_reg<-read.table(file = "test_cds_reg_script.txt", header = TRUE)
comparisons<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","petfal_petpet","Argophyllus_petpet")



out_res<-NULL
for (yy in 1:length(comparisons)){
	for (pp in 1:5){
		
		cds_reg_foc<- cds_reg[(cds_reg$type ==comparisons[yy]) & (cds_reg$recombin==pp),] 
		
		### coding
		sample_cds<-cds_reg_foc[(cds_reg_foc$genomic_pos == "CDS") & (cds_reg_foc$status == "convergent"),]
		mean_sample_cds<-mean(cds_reg_foc[(cds_reg_foc$genomic_pos == "CDS") & (cds_reg_foc$status == "convergent"),1])
		sd_sample_cds<-sd(cds_reg_foc[(cds_reg_foc$genomic_pos == "CDS") & (cds_reg_foc$status == "convergent"),1])
		
		null_mean_cds<-mean(cds_reg_foc[(cds_reg_foc$genomic_pos == "CDS") & (cds_reg_foc$status == "random_control"),1])
		zscore_cds<-(mean_sample_cds-null_mean_cds)/(sd_sample_cds/sqrt(nrow(sample_cds)))
		
		#pvalue_2sided_cds<-2*pnorm(-abs(zscore))  ### two sided
		pvalue_Onesided_cds<-1-(pnorm(abs(zscore_cds)))
		cds<-cbind(mean_sample_cds,sd_sample_cds,null_mean_cds,zscore_cds,pvalue_Onesided_cds)
		
		
		### regulatory
		sample_reg<-cds_reg_foc[(cds_reg_foc$genomic_pos == "regulatory") & (cds_reg_foc$status == "convergent"),]
		mean_sample_reg<-mean(cds_reg_foc[(cds_reg_foc$genomic_pos == "regulatory") & (cds_reg_foc$status == "convergent"),1])
		sd_sample_reg<-sd(cds_reg_foc[(cds_reg_foc$genomic_pos == "regulatory") & (cds_reg_foc$status == "convergent"),1])
		
		null_mean_reg<-mean(cds_reg_foc[(cds_reg_foc$genomic_pos == "regulatory") & (cds_reg_foc$status == "random_control"),1])
		zscore_reg<-(mean_sample_reg-null_mean_reg)/(sd_sample_reg/sqrt(nrow(sample_reg)))
		#pvalue_2sided_reg<-2*pnorm(-abs(zscore))
		pvalue_Onesided_reg<-1-(pnorm(abs(zscore_reg)))
		reg<-data.frame(cbind(mean_sample_reg,sd_sample_reg,null_mean_reg,zscore_reg,pvalue_Onesided_reg))
        
        
        cds_reg_table<-cbind(cds,reg)
        cds_reg_table$comparison<-comparisons[yy]
        cds_reg_table$recombin<-pp
        out_res<-rbind (out_res,cds_reg_table)
		
	}
}

write.table(out_res, file = "out_convergent_windows_coding_noncoding_overlaps_stats.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")




















