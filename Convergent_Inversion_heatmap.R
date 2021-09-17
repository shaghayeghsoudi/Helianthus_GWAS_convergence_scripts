#### analyzing convergent wondows ####
rm(list = ls())
setwd("/data/users/soudi/paper_sunflowers/null_W/inversion_conversion_overlap")

### load climate and soil vars Corrected #####
## load corrected table
null_corr<-read.table(file= "convergent_inversion_overlap_merged_P_LD0.9_1cM", header = TRUE)

### assign vars
climate_vars<-c("latitude","longitude","llevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH")

soil_vars<-c("BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS")


env_soil<-c("latitude","longitude","llevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","RH","BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS")


comparison_type<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petpet","petfal_petpet","Argophyllus_petfal")


### null_corr focal test
#null_corr<-null_corr[grep("baypass",null_corr$analysis),]
null_corr<-null_corr[grep("spearman",null_corr$analysis),]


out_res<-NULL
for (i in 1:length(comparison_type)){
		
	for (j in 1:length(env_soil)){
	
	foc<-null_corr[null_corr$direction==comparison_type[i]& null_corr$variable==env_soil[j],]
	if(nrow (foc)>0){
	
	#foc_com<-foc[!duplicated(foc$window_name),]
	
	convergent_count<-foc[,c(5,9,11,12,14)]
	out_res_good<-data.frame(convergent_count)
    out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    #out_res_good$type<-as.character(unique(foc$type))   
    out_res_good$type<-"environment"
    out_res_good$analysis<-foc$analysis
    out_res <- rbind (out_res_good,out_res)
    
    
    }else{
     
   convergent_count<-data.frame(N_cluster = c(0),porportion_number_cluster_overlap_inversion = c(0), porportion_number_cluster_overlap_inversion_P = c(0),porportion_length_cluster_overlap_inversion=c(0),porportion_length_cluster_overlap_inversion_P=c(0))
    
	out_res_good<-data.frame(convergent_count)
	out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    out_res_good$type<-"environment"
    out_res_good$analysis<-"spearman"               #### CHANGE here according to the analysis ####
    out_res <- rbind (out_res_good,out_res)
            } ## else statement
        }
}

absent<-c("PAS")
comp<-c("Argophyllus_petfal","Argophyllus_petpet","Annuus_Argophyllus")
out_res$porportion_number_cluster_overlap_inversion[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"
out_res$porportion_number_cluster_overlap_inversion_P[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"
out_res$porportion_length_cluster_overlap_inversion[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"
out_res$porportion_length_cluster_overlap_inversion_P[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"
out_res$N_cluster[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"


write.table(out_res, file = "convergent_clusters_overlaps_with_inversions_pervariable.Spearman.table", col.names = TRUE, row.names = FALSE, quote= FALSE)




##################################
############ phenotypes ##########
rm(list = ls())
setwd("/data/users/soudi/paper_sunflowers/null_W/inversion_conversion_overlap")

### load climate and soil vars Corrected #####
## load corrected table
null_corr<-read.table(file= "convergent_inversion_overlap_merged_P_LD0.9_1cM", header = TRUE)


env_soil<-c("DTF","Leaf_height_mid-width","Leaf_circular","Leaf_area","Leaf_ellipsoid","Leaf_eccentricity","Leaf_curved_shape_index","Leaf_curvedHeight_maxWidth","Leaf_distal_eccentricity","Leaf_C_N_ratio","Leaf_eccentricity_area_index","Leaf_obovoid","Leaf_width_widest_pos","Leaf_shape_index_internal","Leaf_total_C","Leaf_total_N","Leaf_shape_index_external_I","Leaf_rectangular","Leaf_maximum_height","Leaf_perimeter","Leaf_proximal_eccentricity","Leaf_shape_index_external_II","Leaf_width_mid-height","Leaf_maximum_width","RGB_proportion_red","RGB_proportion_blue","Total_RGB","RGB_proportion_green","SLA","Guides_3_petals","Internode_length","Distance_of_first_branching_from_ground","Days_to_budding","Guides_individual_petal","Flower_FHDD_ratio","Inflorescence_diameter","Disk_diameter","Primary_branches","Ligule_width","Ligule_LW_ratio","Plant_height_at_flowering","Ligule_length","Ligules_number","LIR","Seed_shape_index_external_I","Seed_HW_ratio","Seed_eccentricity","Seed_maximum_width","Seed_perimeter","Seed_curved_height","Seed_ellipsoid","Seed_area","Seed_eccentricity_area_index","Seed_curved_shape_index","Seed_maximum_height","Seed_rectangular","Seed_circular","Seed_ovoid","Seed_distal_eccentricity","Seed_height_mid_width","Seed_width_mid_height","TLN","Seed_proximal_eccentricity","Stem_diameter_at_flowering","Seed_width_widest_pos","Stem_diameter_final_after_5th_node","Seed_shape_index_external_II","Stem_diameter_final_before_1st_node","Seed_shape_index_internal","Peduncle_length_of_first_flower")

comparison_type<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petpet","petfal_petpet","Argophyllus_petfal")


### null_corr focal test
null_corr<-null_corr[grep("GWAS_corrected",null_corr$analysis),]
#null_corr<-null_corr[grep("GWAS_uncorrected",null_corr$analysis),]





out_res<-NULL
for (i in 1:length(comparison_type)){
		
	for (j in 1:length(env_soil)){
	
	foc<-null_corr[null_corr$direction==comparison_type[i]& null_corr$variable==env_soil[j],]
	if(nrow (foc)>0){
	
	#foc_com<-foc[!duplicated(foc$window_name),]
	
	convergent_count<-foc[,c(5,9,11,12,14)]
	out_res_good<-data.frame(convergent_count)
    out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    #out_res_good$type<-as.character(unique(foc$type))   
    out_res_good$type<-"phenotype"
    out_res_good$analysis<-foc$analysis
    out_res <- rbind (out_res_good,out_res)

    
    
    }else{
     
    convergent_count<-data.frame(N_cluster = c(0),porportion_number_cluster_overlap_inversion = c(0), porportion_number_cluster_overlap_inversion_P = c(0),porportion_length_cluster_overlap_inversion=c(0),porportion_length_cluster_overlap_inversion_P=c(0))
    

	out_res_good<-data.frame(convergent_count)
	out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    out_res_good$type<-"phenotype"
    out_res_good$analysis<-"GWAS_corrected"               #### CHANGE here according to the analysis ####
    out_res <- rbind (out_res_good,out_res)

            } ## else statement
        }
}

## vars which are not in argophyllus
absent<-c("TLN","LIR","DTF","Stem_diameter_at_flowering","Plant_height_at_flowering" ,"Internode_length" ,"Stem_diameter_final_before_1st_node" ,"Stem_diameter_final_after_5th_node","Distance_of_first_branching_from_ground","Primary_branches","Disk_diameter","Ligule_width","Ligule_length","Ligule_LW_ratio","Flower_FHDD_ratio" ,"Ligules_number","Seed_area","Seed_maximum_width","Seed_width_mid_height","Seed_height_mid_width","Seed_maximum_height","Seed_curved_height","Seed_HW_ratio","Seed_shape_index_external_I","Seed_shape_index_external_II","Seed_curved_shape_index" ,"Seed_ellipsoid","Seed_circular","Seed_rectangular","Seed_ovoid","Seed_width_widest_pos","Seed_eccentricity","Seed_shape_index_internal","Seed_eccentricity_area_index","Inflorescence_diameter","Seed_proximal_eccentricity","Seed_distal_eccentricity","Peduncle_length_of_first_flower")

comp<-c("Argophyllus_petfal","Argophyllus_petpet","Annuus_Argophyllus")
out_res$porportion_number_cluster_overlap_inversion[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"
out_res$porportion_number_cluster_overlap_inversion_P[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"
out_res$porportion_length_cluster_overlap_inversion[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"
out_res$porportion_length_cluster_overlap_inversion_P[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"
out_res$N_cluster[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"


write.table(out_res, file = "convergent_clusters_overlaps_with_inversions_pervariable.GWAS_corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


##############################################
###### plot results in pannel heatmap ########
##############################################

rm(list = ls())
#### Speearma and GWAS corrected
library(pheatmap)
library(RColorBrewer)

setwd("~/Documents/input_data/tests_analysis_for_paper/inversion_convergence/convergent_inversion_overlaps_per_variable_heatmap")

### load phenotypes
pheno<-read.table(file = "convergent_clusters_overlaps_with_inversions_pervariable.GWAS_corrected.table", header = TRUE)
pheno$porportion_number_cluster_overlap_inversion[pheno$N_cluster==0]<-NA
pheno$porportion_number_cluster_overlap_inversion_P[pheno$N_cluster==0]<-NA
pheno$porportion_length_cluster_overlap_inversion[pheno$N_cluster==0]<-NA
pheno$porportion_length_cluster_overlap_inversion_P[pheno$N_cluster==0]<-NA
pheno$N_cluster[pheno$N_cluster==0]<-NA


## load environmental variables
env<-read.table(file = "convergent_clusters_overlaps_with_inversions_pervariable.Spearman.table", header = TRUE)
env$porportion_number_cluster_overlap_inversion[env$N_cluster==0]<-NA
env$porportion_number_cluster_overlap_inversion_P[env$N_cluster==0]<-NA
env$porportion_length_cluster_overlap_inversion[env$N_cluster==0]<-NA
env$porportion_length_cluster_overlap_inversion_P[env$N_cluster==0]<-NA
env$N_cluster[env$N_cluster==0]<-NA

##
data_all<-do.call("rbind", list(env,pheno))
data_all$porportion_number_cluster_overlap_inversion_P_TRUE<-(1-(data_all$porportion_number_cluster_overlap_inversion_P))
data_all$porportion_length_cluster_overlap_inversion_p_TRUE<-(1-(data_all$porportion_length_cluster_overlap_inversion_P))
data_all<-data_all[gsub("llevation","elevation",data_all)]


ann_arg<-data_all[data_all$comparison=="Annuus_Argophyllus",]
ann_pf<-data_all[data_all$comparison=="Annuus_petfal",]
ann_pp<-data_all[data_all$comparison=="Annuus_petpet",]
arg_pf<-data_all[data_all$comparison=="Argophyllus_petfal",]
arg_pp<-data_all[data_all$comparison=="Argophyllus_petpet",]
pf_pp<-data_all[data_all$comparison=="petfal_petpet",]


## proportion number overlap
all_comparisonnumber   <- cbind(ann_arg[,2],ann_pf[,2],ann_pp[,2],arg_pf[,2],arg_pp[,2],pf_pp[,2])
all_comparisonnumber<-t(all_comparisonnumber)
colnames(all_comparisonnumber)<-ann_arg$var

rownames(all_comparisonnumber)<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")


test_labels_number <- cbind(ann_arg[,10],ann_pf[,10],ann_pp[,10],arg_pf[,10],arg_pp[,10],pf_pp[,10]) 
test_labels_number<-t(test_labels_number)

colnames(test_labels_number)<-ann_arg$var

rownames(test_labels_number)<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")

test_labels_number[test_labels_number < 0.05 ] <- "*"
test_labels_number[test_labels_number > 0.05] <- ""
test_labels_number[is.na(test_labels_number)] <- ""


### add gaps ####
par(mar = c(1, 0.5, 1, 1))
pdf(file = "Convergence_Inversion_Overlaps_GWAS_Corrected_Spearman_proportion_number.pdf", height= 4, width = 15)

aa<-data.frame(category = c(rep("Soil", 15), rep("Climate", 25),rep("Phenotype", 70)))
row.names(aa)<-colnames(all_comparisonnumber)


pheatmap(all_comparisonnumber,cluster_cols = F, cluster_rows = F,cex= 1, col=brewer.pal(11,"RdBu"),cellwidth=8,cellheight = 8,fontsize = 7,gaps_col = c(15,40,110),annotation_colors = anno_colors,na_col = "gray70",display_numbers = test_labels_number,number_color = "black",ontsize_number=200)

dev.off()


###################################
#### proportion length overlap ####
all_comparisonlength   <- cbind(ann_arg[,4],ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])
all_comparisonlength <-t(all_comparisonlength )
colnames(all_comparisonlength )<-ann_arg$var

rownames(all_comparisonlength )<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")


test_labels_length <- cbind(ann_arg[,11],ann_pf[,11],ann_pp[,11],arg_pf[,11],arg_pp[,11],pf_pp[,11]) 
test_labels_length<-t(test_labels_length)

colnames(test_labels_length)<-ann_arg$var

rownames(test_labels_length)<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")

test_labels_length[test_labels_length < 0.05 ] <- "*"
test_labels_length[test_labels_length > 0.05] <- ""
test_labels_length[is.na(test_labels_length)] <- ""


### add gaps ####
par(mar = c(1, 0.5, 1, 1))
pdf(file = "Convergence_Inversion_Overlaps_GWAS_Corrected_Spearman_proportion_length.pdf", height= 4, width = 15)

aa<-data.frame(category = c(rep("Soil", 15), rep("Climate", 25),rep("Phenotype", 70)))
row.names(aa)<-colnames(test_labels_length)


pheatmap(all_comparisonlength,cluster_cols = F, cluster_rows = F,cex= 1, col=brewer.pal(11,"RdBu"),cellwidth=8,cellheight = 8,fontsize = 7,gaps_col = c(15,40,110),annotation_colors = anno_colors,na_col = "gray70",display_numbers = test_labels_length,number_color = "black",ontsize_number=200)

dev.off()


