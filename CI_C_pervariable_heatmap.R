
#### This Rscript plots heatmaps for the convergence-CI overlaps ####
#### Shaghayegh October 2020 ####
rm(list = ls())
setwd("/data/users/soudi/paper_sunflowers/null_W/C_CI_overlap")

null_corr<-read.table(file= "C_CI_overlap_merged_P_LD0.9_1cM", header = TRUE)

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
	if(nrow(foc)>0){
	
	#foc_com<-foc[!duplicated(foc$window_name),]
	
	convergent_count<-foc[,c(5,7,8,9,11,12,14)]
	out_res_good<-data.frame(convergent_count)
    out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    #out_res_good$type<-as.character(unique(foc$type))   
    out_res_good$type<-"environment"
    out_res_good$analysis<-foc$analysis
    out_res <- rbind (out_res_good,out_res)
    
    
    #}
    #if (foc$N_original_cluster_overlab_CI==0){
    	
    #convergent_count<-foc[,c(5,7,8,9,11,12,14)]
	#out_res_good<-data.frame(convergent_count)
	#out_res_good$var<-env_soil[j]
    #out_res_good$comparison<-comparison_type[i]
    #out_res_good$type<-"environment"
    #out_res_good$analysis<-"spearman"               #### CHANGE here according to the analysis ####
    #out_res <- rbind (out_res_good,out_res)
    	
    	
    }else{
    
    
        
       convergent_count<-data.frame(N_original_cluster  = NA,N_original_cluster_overlab_CI = NA, length_C_overlap_CI = NA,porportion_number_original_cluster_overlap_CI=NA,porportion_number_original_cluster_overlap_CI_P=NA,porportion_length_C_overlap_CI=NA,porportion_length_C_overlap_CI_P=NA)
    

    
	out_res_good<-data.frame(convergent_count)
	out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    out_res_good$type<-"environment"
    out_res_good$analysis<-"BayPass"               #### CHANGE here according to the analysis ####
    out_res <- rbind (out_res_good,out_res)
            } ## else statement
        }
}


write.table(out_res, file = "C_CI_per_variable.BayPass.table", col.names = TRUE, row.names = FALSE, quote= FALSE)

###############################################
##################################
############ phenotypes ##########


rm(list = ls())
setwd("/data/users/soudi/paper_sunflowers/null_W/C_CI_overlap")

null_corr<-read.table(file= "C_CI_overlap_merged_P_LD0.9_1cM", header = TRUE)


env_soil<-c("DTF","Leaf_height_mid-width","Leaf_circular","Leaf_area","Leaf_ellipsoid","Leaf_eccentricity","Leaf_curved_shape_index","Leaf_curvedHeight_maxWidth","Leaf_distal_eccentricity","Leaf_C_N_ratio","Leaf_eccentricity_area_index","Leaf_obovoid","Leaf_width_widest_pos","Leaf_shape_index_internal","Leaf_total_C","Leaf_total_N","Leaf_shape_index_external_I","Leaf_rectangular","Leaf_maximum_height","Leaf_perimeter","Leaf_proximal_eccentricity","Leaf_shape_index_external_II","Leaf_width_mid-height","Leaf_maximum_width","RGB_proportion_red","RGB_proportion_blue","Total_RGB","RGB_proportion_green","SLA","Guides_3_petals","Internode_length","Distance_of_first_branching_from_ground","Days_to_budding","Guides_individual_petal","Flower_FHDD_ratio","Inflorescence_diameter","Disk_diameter","Primary_branches","Ligule_width","Ligule_LW_ratio","Plant_height_at_flowering","Ligule_length","Ligules_number","LIR","Seed_shape_index_external_I","Seed_HW_ratio","Seed_eccentricity","Seed_maximum_width","Seed_perimeter","Seed_curved_height","Seed_ellipsoid","Seed_area","Seed_eccentricity_area_index","Seed_curved_shape_index","Seed_maximum_height","Seed_rectangular","Seed_circular","Seed_ovoid","Seed_distal_eccentricity","Seed_height_mid_width","Seed_width_mid_height","TLN","Seed_proximal_eccentricity","Stem_diameter_at_flowering","Seed_width_widest_pos","Stem_diameter_final_after_5th_node","Seed_shape_index_external_II","Stem_diameter_final_before_1st_node","Seed_shape_index_internal","Peduncle_length_of_first_flower")

comparison_type<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petpet","petfal_petpet","Argophyllus_petfal")


### null_corr focal test
null_corr<-null_corr[grep("GWAS_corrected",null_corr$analysis),]
#null_corr<-null_corr[grep("GWAS_uncorrected",null_corr$analysis),]


out_res<-NULL
for (i in 1:length(comparison_type)){
		
	for (j in 1:length(env_soil)){
	
	foc<-null_corr[null_corr$direction==comparison_type[i]& null_corr$variable==env_soil[j],]
	if(nrow(foc)>0){
	
	#foc_com<-foc[!duplicated(foc$window_name),]
	
	convergent_count<-foc[,c(5,7,8,9,11,12,14)]
	out_res_good<-data.frame(convergent_count)
    out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    #out_res_good$type<-as.character(unique(foc$type))   
    out_res_good$type<-"phenotype"
    out_res_good$analysis<-foc$analysis
    out_res <- rbind (out_res_good,out_res)
    
    
    #}
    #if (foc$N_original_cluster_overlab_CI==0){
    	
    #convergent_count<-foc[,c(5,7,8,9,11,12,14)]
	#out_res_good<-data.frame(convergent_count)
	#out_res_good$var<-env_soil[j]
    #out_res_good$comparison<-comparison_type[i]
    #out_res_good$type<-"environment"
    #out_res_good$analysis<-"spearman"               #### CHANGE here according to the analysis ####
    #out_res <- rbind (out_res_good,out_res)
    	
    	
    }else{
    
    
       convergent_count<-data.frame(N_original_cluster  = NA,N_original_cluster_overlab_CI = NA, length_C_overlap_CI =  NA,porportion_number_original_cluster_overlap_CI=NA,porportion_number_original_cluster_overlap_CI_P=NA,porportion_length_C_overlap_CI=NA,porportion_length_C_overlap_CI_P=NA)
    

    
	out_res_good<-data.frame(convergent_count)
	out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    out_res_good$type<-"phenotype"
    out_res_good$analysis<-"GWAS_corrected"               #### CHANGE here according to the analysis ####
    out_res <- rbind (out_res_good,out_res)
            } ## else statement
        }
}


write.table(out_res, file = "C_CI_per_variable.GWAS_corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


###################################
####### plot the results ##########
###################################

rm(list = ls())
#### Speearma and GWAS corrected
library(pheatmap)
library(RColorBrewer)

setwd("~/Documents/input_data/tests_analysis_for_paper/null_w/clusters/C_CI_per_variable")

pheno<-read.table(file = "C_CI_per_variable.GWAS_Uncorrected.table", header = TRUE)
env<-read.table(file = "C_CI_per_variable.BayPass.table", header = TRUE)

data_all<-do.call("rbind", list(env,pheno))
data_all$porportion_number_original_cluster_overlap_CI_P_TRUE<-(1-(data_all$porportion_number_original_cluster_overlap_CI_P))
data_all$porportion_length_C_overlap_CI_P_TRUE<-(1-(data_all$porportion_length_C_overlap_CI_P))

ann_arg<-data_all[data_all$comparison=="Annuus_Argophyllus",]
ann_pf<-data_all[data_all$comparison=="Annuus_petfal",]
ann_pp<-data_all[data_all$comparison=="Annuus_petpet",]
arg_pf<-data_all[data_all$comparison=="Argophyllus_petfal",]
arg_pp<-data_all[data_all$comparison=="Argophyllus_petpet",]
pf_pp<-data_all[data_all$comparison=="petfal_petpet",]



## proportion size and length overlap
all_comparisonSize   <- cbind(ann_arg[,4],ann_pf[,4],ann_pp[,4],arg_pf[,4],arg_pp[,4],pf_pp[,4])    #### number ###
#all_comparisonSize   <- cbind(ann_arg[,6],ann_pf[,6],ann_pp[,6],arg_pf[,6],arg_pp[,6],pf_pp[,6])    #### length ###

all_comparisonSize<-t(all_comparisonSize)
colnames(all_comparisonSize)<-ann_arg$var

rownames(all_comparisonSize)<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")

### proportion length ###
test_labels <- cbind(ann_arg[,12],ann_pf[,12],ann_pp[,12],arg_pf[,12],arg_pp[,12],pf_pp[,12])  #### proportion number ###
#test_labels <- cbind(ann_arg[,13],ann_pf[,13],ann_pp[,13],arg_pf[,13],arg_pp[,13],pf_pp[,13])  #### proportion length ###


test_labels<-t(test_labels)
colnames(test_labels)<-ann_arg$var

rownames(test_labels)<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")

test_labels[test_labels <= 0.05 ] <- "*"
test_labels[test_labels > 0.05] <- ""
test_labels[is.na(test_labels)] <- ""


### add gaps ####

pdf(file ="Convergence_CI_GWAS_UnCorrected_BayPass_proportion_number.pdf", height= 4, width= 15)
par(mar = c(1, 0.5, 1, 1))
aa<-data.frame(category = c(rep("Soil", 15), rep("Climate", 25),rep("Phenotype", 70)))
row.names(aa)<-colnames(all_comparisonSize)


pheatmap(all_comparisonSize,cluster_cols = F, cluster_rows = F,cex= 1, col=brewer.pal(11,"RdBu"),cellwidth=8,cellheight = 8,fontsize = 7,gaps_col = c(15,40,110),annotation_colors = anno_colors,na_col = "gray70",display_numbers = test_labels,number_color = "black",ontsize_number=200)

dev.off()


#####################################
##### size and length overlap #######

## proportion size and length overlap
all_comparisonlength   <- cbind(ann_arg[,3],ann_pf[,3],ann_pp[,3],arg_pf[,3],arg_pp[,3],pf_pp[,3])    #### length ###
all_comparisonnumber  <- cbind(ann_arg[,2],ann_pf[,2],ann_pp[,2],arg_pf[,2],arg_pp[,2],pf_pp[,2])    #### number ###

all_comparisonlength<-t(all_comparisonlength)
colnames(all_comparisonlength)<-ann_arg$var
rownames(all_comparisonlength)<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")


all_comparisonnumber<-t(all_comparisonnumber)
colnames(all_comparisonnumber)<-ann_arg$var
rownames(all_comparisonnumber)<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")


### add gaps ####

pdf(file ="Convergence_CI_GWAS_Corrected_Spearman_length.pdf", height= 4, width= 15)
par(mar = c(1, 0.1, 1, 1))
aa<-data.frame(category = c(rep("Soil", 15), rep("Climate", 25),rep("Phenotype", 70)))
row.names(aa)<-colnames(all_comparisonlength)
pheatmap(all_comparisonlength,cluster_cols = F, cluster_rows = F,cex= 1, col=brewer.pal(11,"RdBu"),cellwidth=7.5,cellheight = 8,fontsize = 7,gaps_col = c(15,40,110),annotation_colors = anno_colors,na_col = "gray70")

dev.off()



pdf(file ="Convergence_CI_GWAS_Corrected_Spearman_number.pdf", height= 4, width= 14)
par(mar = c(1, 0.5, 1, 1))
aa<-data.frame(category = c(rep("Soil", 15), rep("Climate", 25),rep("Phenotype", 70)))
row.names(aa)<-colnames(all_comparisonnumber)
pheatmap(all_comparisonnumber,cluster_cols = F, cluster_rows = F,cex= 1, col=brewer.pal(11,"RdBu"),cellwidth=7.5,cellheight = 8,fontsize = 7,gaps_col = c(15,40,110),annotation_colors = anno_colors,na_col = "gray70")

dev.off()
