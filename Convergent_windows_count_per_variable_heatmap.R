#### analyzing convergent wondows ####
rm(list = ls())
setwd("/data/users/soudi/paper_sunflowers/null_W")

### load climate and soil vars Corrected #####
## load corrected table
null_corr<-read.table(file= "convergent_2binning_soil_climate_BayPass.table", header = FALSE)
colnames(null_corr)<-c("data","type","analysis","Taxa1","Taxa2","comparison","test_name","window_name","window5","Emp_p")

### assign vars
climate_vars<-c("latitude","longitude","llevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH")

soil_vars<-c("BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS")


env_soil<-c("latitude","longitude","llevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","RH","BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS")


comparison_type<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petpet","petfal_petpet","Argophyllus_petfal")


#### count total number of convergent windows per variable for each pair

out_res<-NULL
for (i in 1:length(comparison_type)){
		
	for (j in 1:length(env_soil)){
	
	foc<-null_corr[null_corr$comparison==comparison_type[i]& null_corr$test_name==env_soil[j],]
	if(nrow (foc)>0){
	
	#foc_com<-foc[!duplicated(foc$window_name),]
	
	convergent_count<-nrow(foc)
	out_res_good<-data.frame(convergent_count)
    out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    #out_res_good$type<-as.character(unique(foc$type))   
    out_res_good$type<-"environment"
    out_res_good$analysis<-"baypass"
    out_res <- rbind (out_res_good,out_res)
    
    
    }else{
     
    convergent_count<-0
	out_res_good<-data.frame(convergent_count)
	out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    out_res_good$type<-"environment"
    out_res_good$analysis<-"baypass"
    out_res <- rbind (out_res_good,out_res)
            } ## else statement
        }
}

absent<-c("MAP")
comp<-c("Argophyllus_petfal","Argophyllus_petpet","Annuus_Argophyllus")
out_res$convergent_count[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"

write.table(out_res, file = "convergent_5Kwinows_count_per_variable_BayPass.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


##################################
############ phenotypes ##########

#### analyzing convergent wondows ####
rm(list = ls())
setwd("/data/users/soudi/paper_sunflowers/null_W")

### load climate and soil vars Corrected #####
## load corrected table
null_corr<-read.table(file= "convergent_2binning_phenotype_GWAS_corrected", header = FALSE)
colnames(null_corr)<-c("data","type","analysis","Taxa1","Taxa2","comparison","test_name","window_name","window5","Emp_p")


### assign vars
env_soil<-c("DTF","Leaf_height_mid-width","Leaf_circular","Leaf_area","Leaf_ellipsoid","Leaf_eccentricity","Leaf_curved_shape_index","Leaf_curvedHeight_maxWidth","Leaf_distal_eccentricity","Leaf_C_N_ratio","Leaf_eccentricity_area_index","Leaf_obovoid","Leaf_width_widest_pos","Leaf_shape_index_internal","Leaf_total_C","Leaf_total_N","Leaf_shape_index_external_I","Leaf_rectangular","Leaf_maximum_height","Leaf_perimeter","Leaf_proximal_eccentricity","Leaf_shape_index_external_II","Leaf_width_mid-height","Leaf_maximum_width","RGB_proportion_red","RGB_proportion_blue","Total_RGB","RGB_proportion_green","SLA","Guides_3_petals","Internode_length","Distance_of_first_branching_from_ground","Days_to_budding","Guides_individual_petal","Flower_FHDD_ratio","Inflorescence_diameter","Disk_diameter","Primary_branches","Ligule_width","Ligule_LW_ratio","Plant_height_at_flowering","Ligule_length","Ligules_number","LIR","Seed_shape_index_external_I","Seed_HW_ratio","Seed_eccentricity","Seed_maximum_width","Seed_perimeter","Seed_curved_height","Seed_ellipsoid","Seed_area","Seed_eccentricity_area_index","Seed_curved_shape_index","Seed_maximum_height","Seed_rectangular","Seed_circular","Seed_ovoid","Seed_distal_eccentricity","Seed_height_mid_width","Seed_width_mid_height","TLN","Seed_proximal_eccentricity","Stem_diameter_at_flowering","Seed_width_widest_pos","Stem_diameter_final_after_5th_node","Seed_shape_index_external_II","Stem_diameter_final_before_1st_node","Seed_shape_index_internal","Peduncle_length_of_first_flower")



comparison_type<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petpet","petfal_petpet","Argophyllus_petfal")
#### count total number of convergent windows per variable for each pair

out_res<-NULL
for (i in 1:length(comparison_type)){
		
	for (j in 1:length(env_soil)){
	
	foc<-null_corr[null_corr$comparison==comparison_type[i]& null_corr$test_name==env_soil[j],]
	if(nrow (foc)>0){
	
	#foc_com<-foc[!duplicated(foc$window_name),]
	
	convergent_count<-nrow(foc)
	out_res_good<-data.frame(convergent_count)
    out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    #out_res_good$type<-as.character(unique(foc$type))   
    out_res_good$type<-"phenotype"
    out_res_good$analysis<-"GWAS_corrected"
    out_res <- rbind (out_res_good,out_res)
    
    
    }else{
     
    convergent_count<-0
	out_res_good<-data.frame(convergent_count)
	out_res_good$var<-env_soil[j]
    out_res_good$comparison<-comparison_type[i]
    out_res_good$type<-"phenotype"
    out_res_good$analysis<-"GWAS_corrected"
    out_res <- rbind (out_res_good,out_res)
            } ## else statement
        }
}

## vars which are not in argophyllus
absent<-c("TLN","LIR","DTF","Stem_diameter_at_flowering","Plant_height_at_flowering" ,"Internode_length" ,"Stem_diameter_final_before_1st_node" ,"Stem_diameter_final_after_5th_node","Distance_of_first_branching_from_ground","Primary_branches","Disk_diameter","Ligule_width","Ligule_length","Ligule_LW_ratio","Flower_FHDD_ratio" ,"Ligules_number","Seed_area","Seed_maximum_width","Seed_width_mid_height","Seed_height_mid_width","Seed_maximum_height","Seed_curved_height","Seed_HW_ratio","Seed_shape_index_external_I","Seed_shape_index_external_II","Seed_curved_shape_index" ,"Seed_ellipsoid","Seed_circular","Seed_rectangular","Seed_ovoid","Seed_width_widest_pos","Seed_eccentricity","Seed_shape_index_internal","Seed_eccentricity_area_index","Inflorescence_diameter","Seed_proximal_eccentricity","Seed_distal_eccentricity","Peduncle_length_of_first_flower")

comp<-c("Argophyllus_petfal","Argophyllus_petpet","Annuus_Argophyllus")

out_res$convergent_count[out_res$var%in%absent & out_res$comparison%in%comp]<-"NA"

write.table(out_res, file = "convergent_5Kwinows_count_per_variable_Phenotype_GWAS_ccorrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)

##############################################
###### plot results in pannel heatmap ########
##############################################
rm(list = ls())
#### Speearma and GWAS corrected
library(pheatmap)
library(RColorBrewer)

setwd("~/Documents/input_data/tests_analysis_for_paper/null_w/count_convergent_windows_per_variable")

env<-read.table(file = "convergent_5Kwinows_count_per_variable_Spearman.table", header = TRUE)
pheno<-read.table(file = "convergent_5Kwinows_count_per_variable_Phenotype_GWAS_corrected.table", header = TRUE)

data_all<-do.call("rbind", list(env,pheno))

ann_arg<-data_all[data_all$comparison=="Annuus_Argophyllus",]
ann_pf<-data_all[data_all$comparison=="Annuus_petfal",]
ann_pp<-data_all[data_all$comparison=="Annuus_petpet",]
arg_pf<-data_all[data_all$comparison=="Argophyllus_petfal",]
arg_pp<-data_all[data_all$comparison=="Argophyllus_petpet",]
pf_pp<-data_all[data_all$comparison=="petfal_petpet",]


all_comparisonSize   <- cbind(ann_arg[,1],ann_pf[,1],ann_pp[,1],arg_pf[,1],arg_pp[,1],pf_pp[,1])
all_comparisonSize<-t(all_comparisonSize)
colnames(all_comparisonSize)<-ann_arg$var

rownames(all_comparisonSize)<-c("H.annuus-H.argophyllus","H.annuus-H.petiolaris.fallax","H.annuus-H.petiolaris.petiolaris","H.argophyllus-H.petiolaris.fallax","H.argophyllus-H.petiolaris.petiolaris","H.petiolaris.fallax-H.petiolaris.petiolaris")

### add gaps ####
par(mar = c(1, 0.5, 1, 1))
pdf(file = "heatmap_convergent_windows_count_Spearman_GWAS_Corrected.pdf", height= 4, width = 15)

aa<-data.frame(category = c(rep("Soil", 15), rep("Climate", 25),rep("Phenotype", 70)))
row.names(aa)<-colnames(all_comparisonSize)
pheatmap(all_comparisonSize,cluster_cols = F, cluster_rows = F,cex= 1, col=brewer.pal(11,"RdBu"),cellwidth=8,cellheight = 8,fontsize = 7,gaps_col = c(15,40,110),annotation_colors = anno_colors,na_col = "gray70")
dev.off()

##############################################
