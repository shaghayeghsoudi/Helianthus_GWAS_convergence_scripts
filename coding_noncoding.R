
library(GenomicRanges)
library(dplyr)
setwd("~/Documents/input_data/tests_analysis_for_paper/coding_noncoding")


### gff files
gff_gene<-read.table(file = "gene", header = FALSE)
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
recom<-read.table(file ="recombination_rate_bins/out_res_5k_with_5quantile.table", header = TRUE)
recom<-recom[!(duplicated(recom$window)),]


### load winow5k coordinate
window5k<-read.table(file = "recombination_rate_bins/bins_good_5000", header = FALSE)
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


for (ii in 1:length(comparisons)){	
	aa<-assigned_recom_con[(assigned_recom_con$comparison==comparisons[ii]),]
	aa_foc<-aa[(aa$analysis=="spearman")  | (aa$analysis=="GWAS_corrected"),]
	aa_foc<-aa_foc[!(duplicated(aa_foc$gene_id)),]
    #aa_foc$chrom<-gsub(":.*$","",aa_foc$window5)
    #gene_id_foc<-unique(aa_foc$gene_id)
    
	
	for (jj in 1:5){ ## 5 recombination bin quantiles
		
		assigned_recom_focal<-aa_foc[aa_foc$quantile==jj,]          ## convergent windows
		window_id_foc<-unique(assigned_recom_focal$window_name)
				
#out_res<-NULL

		for (zz in 1:length(window_id_foc)){
			
         foc<-assigned_recom_focal[assigned_recom_focal$window_name==window_id_foc[zz],]   
         cds_focal<-cds[cds$gene_id==foc$gene_id,]
         reg_focal<-regulatory[regulatory$gene_id==foc$gene_id,]
    
    
         gr_aa <- GRanges(IRanges(start = foc$window_start,end =      foc$window_end),seqnames = foc$Chr)
         gr_cds<-GRanges(IRanges(start = cds_focal$start,end =      cds_focal$end),seqnames = cds_focal$chrom,pos= cds_focal$genomic_position )
         gr_regulatory<-GRanges(IRanges(start = reg_focal$start,end =      reg_focal$end),seqnames = reg_focal$chrom,pos=     reg_focal$genomic_position )
         coding_overlap<-data.frame(intersect(gr_aa, gr_cds))[,1:4]
  
       if (nrow(coding_overlap)>=1){
  
      proporlength<-sum(coding_overlap$width)/5000
      length<-sum(coding_overlap$width)
      window_id<-as.character(window_id_foc[zz])
      genomic_pos<-"CDS"
      recombin<-jj
      type<-comparisons[ii]     
      status<-"convergent" 
      coding_overlap_con<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type, status))
  
     
     
      }else{
      	
      	 coding_overlap_con <- setNames(data.frame(t(c(0, 0, 0, 0, 0, 0, 0))),
    c("proporlength","length","window_id","genomic_pos","recombin","type","status"))
          coding_overlap_con$window_id_con<-as.character(window_id_foc[zz])   	
  	      coding_overlap_con$genomic_pos_con<-"CDS"
    	  coding_overlap_con$recombin<-jj
          coding_overlap_con$type<-comparisons[ii]
          coding_overlap_con$status<-"convergent"
          coding_overlap_con<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type, status))
      	
      }
  	      	
  
  reg_overlap<-data.frame(intersect(gr_aa,gr_regulatory))[,1:4]        ### overlaping of convergent windows with regulatory regions
  
  if (nrow(reg_overlap)>=1){
      proporlength<-sum(reg_overlap$width)/5000
      length<-sum(reg_overlap$width)
      window_id<-as.character(window_id_foc[zz])
      genomic_pos<-"regulatory"
      recombin<-jj
      type<-comparisons[ii]   
      status<-"convergent"   
      reg_overlap_con<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type, status))
      
        
  }else{
  	      
   reg_overlap_con <- setNames(data.frame(t(c(0, 0, 0, 0,0,0,0))),
    c("proporlength","length","window_id","genomic_pos","recombin","type","status"))
           reg_overlap_con$window_id_con<-as.character(window_id_foc[zz])   	
  	       reg_overlap_con$genomic_pos_con<-"regulatory"
           recombin<-jj
           type<-comparisons[ii]      
           status<-"convergent"
           reg_overlap_con<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type, status))    	
     	} ## else loop

		
	# } ## zz loop (focal convergent genes)
	
  
cds_reg<-do.call("rbind", list(coding_overlap_con,reg_overlap_con)) 
out_res <- rbind (out_res,cds_reg)

   
   ###########################################################
   #### control random windows for the NULL distribution #####
   recom_window5k_coor_foc<-recom_window5k_coor[recom_window5k_coor$quantile==jj,]
   qq<-recom_window5k_coor_foc[!(recom_window5k_coor_foc$window_id %in% assigned_recom_focal$window_name),] 
   
   control_window<-unique(qq$window_id)
 
 
 for (rr in 1:length(control_window)){
	
    control_gene_foc<-qq[qq$window_id==control_window[rr],]
    cds_focal_control<-cds[cds$gene_id==control_gene_foc$gene_id ,]
    reg_focal_control<-regulatory[regulatory$gene_id==control_gene_foc$gene_id ,]
    
    
    gr_ww <- GRanges(IRanges(start = control_gene_foc$start,end =      control_gene_foc$end),seqnames = control_gene_foc$chrom)
    gr_cds_ww<-GRanges(IRanges(start = cds_focal_control$start,end =      cds_focal_control$end),seqnames = cds_focal_control$chrom,pos= cds_focal_control$genomic_position )
    gr_regulatory_ww<-GRanges(IRanges(start = reg_focal_control$start,end =      reg_focal_control$end),seqnames = reg_focal_control$chrom,pos=     reg_focal_control$genomic_position )
    coding_overlap_ww<-data.frame(intersect(gr_ww, gr_cds_ww))[,1:4]
    #reg_overlap_ww<-data.frame(intersect(gr_ww, gr_regulatory_ww))[,1:4]

 
 if (nrow(coding_overlap_ww)>=1){
  
  proporlength<-sum(coding_overlap_ww$width)/5000
  length<-sum(coding_overlap_ww$width)
  window_id<-as.character(control_window[rr])
  genomic_pos<-"CDS"
  recombin<-jj
  type<-comparisons[ii] 
  status<-"random_control"  
  coding_overlap_null<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type,status))
  
  }else{
  	
  	
  	coding_overlap_null <- setNames(data.frame(t(c(0, 0, 0, 0,0,0,0))),
    c("proporlength","length","window_id","genomic_pos","recombin","type","status"))
           coding_overlap_null$window_id<-as.character(control_window[rr])   	
  	       coding_overlap_null$genomic_pos<-"CDS"
  	       coding_overlap_null$recombin<-jj
  	       coding_overlap_null$type<-comparisons[ii] 
  	       coding_overlap_null$status<-"random_control"
  	       coding_overlap_null<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type,status))  	      
  	}
  
  
  ### regulatory random windows
  reg_overlap_ww<-data.frame(intersect(gr_ww, gr_regulatory_ww))[,1:4]
  
  if (nrow(reg_overlap_ww)>=1){
  proporlength<-sum(reg_overlap_ww$width)/5000
  length<-sum(reg_overlap_ww$width)
  window_id<-as.character(control_window[rr])
  genomic_pos<-"regulatory"
  recombin<-jj
  type<-comparisons[ii] 
  status<-"random_control"  
  reg_overlap_null<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type,status))
  
  }else{
  	      
   
   reg_overlap_null <- setNames(data.frame(t(c(0, 0, 0, 0,0,0,0))),
    c("proporlength","length","window_id","genomic_pos","recombin","type","status"))
           reg_overlap_null$window_id<-as.character(control_window[rr])    
           reg_overlap_null$genomic_pos<-"regulatory"   
           reg_overlap_null$recombin<-jj
           reg_overlap_null$type<-comparisons[ii] 
           reg_overlap_null$status<-"random_control"  
           reg_overlap_null<-cbind(data.frame(proporlength,length,window_id,genomic_pos,recombin,type,status))
              	
         	    
     	    } ## else loop
  
         #} ## rr loop
 



cds_reg<-do.call("rbind", list(coding_overlap_null,reg_overlap_null)) 
out_res <- rbind (out_res,cds_reg)


   } ### jj loop ### 1-5 quantile loop


} ## ii loop  ### comparison loop




















### convergent windows ###

out_res<-NULL
for (i in 1:length(gene_id_foc)){
	#for (i in 1:100){
  aa_foc<-aa[aa$gene_id==gene_id_foc[i],]
  aa_foc<-aa_foc[1,]
  cds_focal<-cds[cds$gene_id==gene_id_foc[i],]
  reg_focal<-regulatory[regulatory$gene_id==gene_id_foc[i],]

  
  gr_aa <- GRanges(IRanges(start = aa_foc$window_start,end =      aa_foc$window_end),seqnames = aa_foc$chrom)
  gr_cds<-GRanges(IRanges(start = cds_focal$start,end =      cds_focal$end),seqnames = cds_focal$chrom,pos= cds_focal$genomic_position )
  gr_regulatory<-GRanges(IRanges(start = reg_focal$start,end =      reg_focal$end),seqnames = reg_focal$chrom,pos=     reg_focal$genomic_position )
  coding_overlap<-data.frame(intersect(gr_aa, gr_cds))[,1:4]
  
  if (nrow(coding_overlap)>=1){
  
  coding_overlap$proporlength<-(coding_overlap$end-coding_overlap$start)/5000
  coding_overlap$length<-(coding_overlap$end-coding_overlap$start)
  coding_overlap$gene_id<-gene_id_foc[i]
  coding_overlap$genomic_pos<-"CDS"
  
  }else{
  	
  	
  	coding_overlap <- setNames(data.frame(t(c(0, 0, 0, 0, 0, 0, 0,"CDS"))),
    c("seqnames","start","end","width","proporlength","length","gene_id","genomic_pos"))
           coding_overlap$gene_id<-gene_id_foc[i]    	
  	
  	}
  
  reg_overlap<-data.frame(intersect(gr_aa,gr_regulatory))[,1:4]
  if (nrow(reg_overlap)>=1){
  reg_overlap$proporlength<-(reg_overlap$end-reg_overlap$start)/5000
  reg_overlap$length<-(reg_overlap$end-reg_overlap$start)
  reg_overlap$gene_id<-gene_id_foc[i]
  reg_overlap$genomic_pos<-"regulatory"
  
  }else{
  	      
   reg_overlap <- setNames(data.frame(t(c(0, 0, 0, 0, 0, 0, 0,"regulatory"))),
    c( "seqnames","start","end","width","proporlength","length","gene_id","genomic_pos"))
           reg_overlap$gene_id<-gene_id_foc[i]                     	
  	}
  	
  cds_reg<-do.call("rbind", list(coding_overlap,reg_overlap)) 
  out_res <- rbind (out_res,cds_reg)
  
}
 
 
 
  
  