setwd("/data/users/soudi/paper_sunflowers/")
#rm(list = ls())

## this script count top candidate windows per variable for corercted and uncorrected methods and plot then 
###################################################
## load corrected files (climate and soil)
###################################################

env_vars<-c("latitude","longitude","llevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH","BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS")



## annuus
ann_Uncorr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_spearman_H.annuus", header = FALSE)
colnames(ann_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ann_Uncorr<-ann_Uncorr[!(ann_Uncorr$outlier_count==0),]
ann_Uncorr_top<-ann_Uncorr[(ann_Uncorr$outlier_count>ann_Uncorr$p4),]
#ann_UncorrUnique<-ann_Uncorr_top[!duplicated(ann_Uncorr_top$windows), ]
#ann_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(env_vars)){
	ann_Uncorr_foc<-ann_Uncorr_top[ann_Uncorr_top$test_name==env_vars[i],]
	top_candidate_count<-nrow(ann_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-env_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_climate/total_top_candidate_winows_count_per_variable_Ordered_H.annuus_Spearman.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


#### H.argophyllus ####
arg_Uncorr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_spearman_H.argophyllus", header = FALSE)
colnames(arg_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
arg_Uncorr<-arg_Uncorr[!(arg_Uncorr$outlier_count==0),]
arg_Uncorr_top<-arg_Uncorr[(arg_Uncorr$outlier_count>arg_Uncorr$p4),]
#arg_Uncorr_topUnique<-arg_Uncorr_top[!duplicated(arg_Uncorr_top$windows), ]
#arg_Uncorr_topUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(env_vars)){
	arg_Uncorr_foc<-arg_Uncorr_top[arg_Uncorr_top$test_name==env_vars[i],]
	top_candidate_count<-nrow(arg_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-env_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_climate/total_top_candidate_winows_count_per_variable_Ordered_H.argophyllus_BayPass_Corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


#####
ppfal_Uncorr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_baypass_H.petiolaris.fallax", header = FALSE)
colnames(ppfal_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppfal_Uncorr<-ppfal_Uncorr[!(ppfal_Uncorr$outlier_count==0),]
ppfal_Uncorr_top<-ppfal_Uncorr[(ppfal_Uncorr$outlier_count>ppfal_Uncorr$p4),]
#ppfal_Uncorr_topUnique<-ppfal_Uncorr_top[!duplicated(ppfal_Uncorr_top$windows), ]
#ppfal_Uncorr_topUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(env_vars)){
	ppfal_Uncorr_foc<-ppfal_Uncorr_top[ppfal_Uncorr_top$test_name==env_vars[i],]
	top_candidate_count<-nrow(ppfal_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-env_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]

}
write.table(out_res, file = "top_candidates_Mojtaba_climate/total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.fallx_BayPass_Corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)




#####
ppet_Uncorr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_baypass_H.petiolaris.Pet", header = FALSE)
colnames(ppet_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppet_Uncorr<-ppet_Uncorr[!(ppet_Uncorr$outlier_count==0),]
ppet_Uncorr_top<-ppet_Uncorr[(ppet_Uncorr$outlier_count>ppet_Uncorr$p4),]
#ppet_Uncorr_topUnique<-ppet_Uncorr_top[!duplicated(ppet_Uncorr_top$windows), ]
#ppet_Uncorr_topUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(env_vars)){
	ppet_Uncorr_foc<-ppet_Uncorr_top[ppet_Uncorr_top$test_name==env_vars[i],]
	top_candidate_count<-nrow(ppet_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-env_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]

}
write.table(out_res, file = "top_candidates_Mojtaba_climate/total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.petiolaris_BayPass_Corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)

####################################
#### phenotypes ####

pheno_vars<-c("Days_to_budding","Disk_diameter","Distance_of_first_branching_from_ground","DTF","Flower_FHDD_ratio","Guides_3_petals","Guides_individual_petal","Inflorescence_diameter","Internode_length","Leaf_area","Leaf_circular","Leaf_C_N_ratio","Leaf_curved_height","Leaf_curvedHeight_maxWidth","Leaf_curved_shape_index","Leaf_distal_eccentricity","Leaf_eccentricity_area_index","Leaf_ellipsoid","Leaf_height_mid-width","Leaf_maximum_height","Leaf_maximum_width","Leaf_obovoid","Leaf_perimeter","Leaf_proximal_eccentricity","Leaf_rectangular","Leaf_shape_index_external_I","Leaf_shape_index_external_II","Leaf_shape_index_internal","Leaf_total_C","Leaf_total_N","Leaf_width_mid-height","Leaf_width_widest_pos","Ligule_length","Ligule_LW_ratio","Ligules_number","Ligule_width","LIR","Plant_height_at_flowering","Primary_branches","RGB_proportion_blue","RGB_proportion_green","RGB_proportion_red","Seed_area","Seed_circular","Seed_curved_height","Seed_curved_shape_index","Seed_distal_eccentricity","Seed_eccentricity","Seed_eccentricity_area_index","Seed_ellipsoid","Seed_height_mid_width","Seed_HW_ratio","Seed_maximum_height","Seed_maximum_width","Seed_ovoid","Seed_perimeter","Seed_proximal_eccentricity","Seed_rectangular","Seed_shape_index_external_I","Seed_shape_index_external_II","Seed_shape_index_internal","Seed_width_mid_height","Seed_width_widest_pos","SLA","Stem_diameter_at_flowering","Stem_diameter_final_after_5th_node","Stem_diameter_final_before_1st_node","TLN","Total_RGB")


##### annuus (GWAS corrected) ######
ann_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_corrected_H.annuus", header = FALSE)
colnames(ann_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ann_Uncorr<-ann_Uncorr[!(ann_Uncorr$outlier_count==0),]
ann_Uncorr_top<-ann_Uncorr[(ann_Uncorr$outlier_count>ann_Uncorr$p4),]
#ann_UncorrUnique<-ann_Uncorr_top[!duplicated(ann_Uncorr_top$windows), ]
#ann_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(pheno_vars)){
	ann_Uncorr_foc<-ann_Uncorr_top[ann_Uncorr_top$test_name==pheno_vars[i],]
	top_candidate_count<-nrow(ann_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-pheno_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_phenotype/total_top_candidate_winows_count_per_variable_Ordered_H.annuus_Phenotype_Corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


###########################################
###### argophyllus (GWAS corrected) #######

arg_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_corrected_H.argophyllus", header = FALSE)
colnames(arg_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
arg_Uncorr<-arg_Uncorr[!(arg_Uncorr$outlier_count==0),]
arg_Uncorr_top<-arg_Uncorr[(arg_Uncorr$outlier_count>arg_Uncorr$p4),]
#arg_UncorrUnique<-arg_Uncorr_top[!duplicated(arg_Uncorr_top$windows), ]
#arg_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(pheno_vars)){
	arg_Uncorr_foc<-arg_Uncorr_top[arg_Uncorr_top$test_name==pheno_vars[i],]
	top_candidate_count<-nrow(arg_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-pheno_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_phenotype/total_top_candidate_winows_count_per_variable_Ordered_H.argophyllus_Phenotype_Corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)

##################################################
####### petiolaris.fallax (GWAS corrected) #######

pf_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_corrected_H.petiolaris.fallax", header = FALSE)
colnames(pf_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
pf_Uncorr<-pf_Uncorr[!(pf_Uncorr$outlier_count==0),]
pf_Uncorr_top<-pf_Uncorr[(pf_Uncorr$outlier_count>pf_Uncorr$p4),]
#pf_UncorrUnique<-pf_Uncorr_top[!duplicated(pf_Uncorr_top$windows), ]
#pf_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(pheno_vars)){
	pf_Uncorr_foc<-pf_Uncorr_top[pf_Uncorr_top$test_name==pheno_vars[i],]
	top_candidate_count<-nrow(pf_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-pheno_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_phenotype/total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.fallax_Phenotype_Corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


#### petiolaris.petiolaris #####
pp_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_corrected_H.petiolaris.petiolaris", header = FALSE)
colnames(pp_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
pp_Uncorr<-pp_Uncorr[!(pp_Uncorr$outlier_count==0),]
pp_Uncorr_top<-pp_Uncorr[(pp_Uncorr$outlier_count>pp_Uncorr$p4),]
#pp_UncorrUnique<-pp_Uncorr_top[!duplicated(pp_Uncorr_top$windows), ]
#pp_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(pheno_vars)){
	pp_Uncorr_foc<-pp_Uncorr_top[pp_Uncorr_top$test_name==pheno_vars[i],]
	top_candidate_count<-nrow(pp_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-pheno_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_phenotype/total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.petiolaris_Phenotype_Corrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)




############### GWAS Uncorrected #############

##### H. annuus (GWAS Uncorrected) ######
ann_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_uncorrected_H.annus", header = FALSE)
colnames(ann_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ann_Uncorr<-ann_Uncorr[!(ann_Uncorr$outlier_count==0),]
ann_Uncorr_top<-ann_Uncorr[(ann_Uncorr$outlier_count>ann_Uncorr$p4),]
#ann_UncorrUnique<-ann_Uncorr_top[!duplicated(ann_Uncorr_top$windows), ]
#ann_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(pheno_vars)){
	ann_Uncorr_foc<-ann_Uncorr_top[ann_Uncorr_top$test_name==pheno_vars[i],]
	top_candidate_count<-nrow(ann_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-pheno_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_phenotype/total_top_candidate_winows_count_per_variable_Ordered_H.annuus_Phenotype_Uncorrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


###########################################
###### H. argophyllus (GWAS Uncorrected) ##

arg_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_uncorrected_H.argophyllus", header = FALSE)
colnames(arg_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
arg_Uncorr<-arg_Uncorr[!(arg_Uncorr$outlier_count==0),]
arg_Uncorr_top<-arg_Uncorr[(arg_Uncorr$outlier_count>arg_Uncorr$p4),]
#arg_UncorrUnique<-arg_Uncorr_top[!duplicated(arg_Uncorr_top$windows), ]
#arg_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(pheno_vars)){
	arg_Uncorr_foc<-arg_Uncorr_top[arg_Uncorr_top$test_name==pheno_vars[i],]
	top_candidate_count<-nrow(arg_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-pheno_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_phenotype/total_top_candidate_winows_count_per_variable_Ordered_H.argophyllus_Phenotype_Uncorrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)


##################################################
####### petiolaris.fallax (GWAS Uncorrected) #####

pf_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_uncorrected_H.petiolaris.fallax", header = FALSE)
colnames(pf_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
pf_Uncorr<-pf_Uncorr[!(pf_Uncorr$outlier_count==0),]
pf_Uncorr_top<-pf_Uncorr[(pf_Uncorr$outlier_count>pf_Uncorr$p4),]
#pf_UncorrUnique<-pf_Uncorr_top[!duplicated(pf_Uncorr_top$windows), ]
#pf_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(pheno_vars)){
	pf_Uncorr_foc<-pf_Uncorr_top[pf_Uncorr_top$test_name==pheno_vars[i],]
	top_candidate_count<-nrow(pf_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-pheno_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_phenotype/total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.fallax_Phenotype_Uncorrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)

##################################################
#### petiolaris.petiolaris (GWAS Uncorrected)#####
pp_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_uncorrected_H.petiolaris.petiolaris", header = FALSE)
colnames(pp_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
pp_Uncorr<-pp_Uncorr[!(pp_Uncorr$outlier_count==0),]
pp_Uncorr_top<-pp_Uncorr[(pp_Uncorr$outlier_count>pp_Uncorr$p4),]
#pp_UncorrUnique<-pp_Uncorr_top[!duplicated(pp_Uncorr_top$windows), ]
#pp_UncorrUnique$status<-rep("Uncorrected")

out_res<-NULL
for (i in 1:length(pheno_vars)){
	pp_Uncorr_foc<-pp_Uncorr_top[pp_Uncorr_top$test_name==pheno_vars[i],]
	top_candidate_count<-nrow(pp_Uncorr_foc)
	out_res_good<-data.frame(top_candidate_count)
    out_res_good$var<-pheno_vars[i]
    out_res <- rbind (out_res_good,out_res)
    out_res<-out_res[order(out_res$top_candidate_count),]
}
write.table(out_res, file = "top_candidates_Mojtaba_phenotype/total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.petiolaris_Phenotype_Uncorrected.table", col.names = TRUE, row.names = FALSE, quote= FALSE)



#-------------------------------------------------------------------
###################################################################
######## plot the results (corrected vs. uncorrected dotplots)#####
###################################################################
#-------------------------------------------------------------------


rm(list = ls())

library(ggplot2)
library(forcats)
setwd("~/Documents/input_data/tests_analysis_for_paper/top_candidate_windows/top_cadidates_count_per_variable")

## oreder of vars
env_vars<-c("latitude","longitude","elevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH","BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS","Days_to_budding","Disk_diameter","Distance_of_first_branching_from_ground","DTF","Flower_FHDD_ratio","Guides_3_petals","Guides_individual_petal","Inflorescence_diameter","Internode_length","Leaf_area","Leaf_circular","Leaf_C_N_ratio","Leaf_curved_height","Leaf_curvedHeight_maxWidth","Leaf_curved_shape_index","Leaf_distal_eccentricity","Leaf_eccentricity_area_index","Leaf_ellipsoid","Leaf_height_mid-width","Leaf_maximum_height","Leaf_maximum_width","Leaf_obovoid","Leaf_perimeter","Leaf_proximal_eccentricity","Leaf_rectangular","Leaf_shape_index_external_I","Leaf_shape_index_external_II","Leaf_shape_index_internal","Leaf_total_C","Leaf_total_N","Leaf_width_mid-height","Leaf_width_widest_pos","Ligule_length","Ligule_LW_ratio","Ligules_number","Ligu le_width","LIR","Plant_height_at_flowering","Primary_branches","RGB_proportion_blue","RGB_proportion_green","RGB_proportion_red","Seed_area","Seed_circular","Seed_curved_height","Seed_curved_shape_index","Seed_distal_eccentricity","Seed_eccentricity","Seed_eccentricity_area_index","Seed_ellipsoid","Seed_height_mid_width","Seed_HW_ratio","Seed_maximum_height","Seed_maximum_width","Seed_ovoid","Seed_perimeter","Seed_proximal_eccentricity","Seed_rectangular","Seed_shape_index_external_I","Seed_shape_index_external_II","Seed_shape_index_internal","Seed_width_mid_height","Seed_width_widest_pos","SLA","Stem_diameter_at_flowering","Stem_diameter_final_after_5th_node","Stem_diameter_final_before_1st_node","TLN","Total_RGB")


## annus
ann_env_correc<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.annuus_BayPass_Corrected.table", header = TRUE)
ann_env_correc$analysis<-rep("corrected")
ann_env_correc$species<-rep("H.annuus")
ann_env_correc$var<-gsub("llevation","elevation",ann_env_correc$var)

ann_env_uncorrec<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.annuus_Spearman.table", header = TRUE)
ann_env_uncorrec$analysis<-rep("uncorrected")
ann_env_uncorrec$species<-rep("H.annuus")
ann_env_uncorrec$var<-gsub("llevation","elevation",ann_env_uncorrec$var)

all_annuus_env<-do.call("rbind", list(ann_env_correc,ann_env_uncorrec)) 
#all_annuus_env_ordered<-all_annuus_env[order(all_annuus_env$top_candidate_count),]
#row.names(all_annuus_env)<- all_annuus_env$var
#all_annuus_env$var<-NULL


ann_pheno_correc<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.annuus_Phenotype_Corrected.table", header = TRUE)
ann_pheno_correc$analysis<-rep("corrected")
ann_pheno_correc$species<-rep("H.annuus")
#ann_pheno_correc<-ann_pheno_correc[!grepl("LIR", ann_pheno_correc$var),]


ann_pheno_uncorrec<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.annuus_Phenotype_Uncorrected.table", header = TRUE)
ann_pheno_uncorrec$analysis<-rep("uncorrected")
ann_pheno_uncorrec$species<-rep("H.annuus")
#ann_pheno_uncorrec<-ann_pheno_uncorrec[!grepl("LIR", ann_pheno_uncorrec$var),]


all_annuus <- do.call("rbind", list(ann_env_correc, ann_env_uncorrec,ann_pheno_correc,ann_pheno_uncorrec)) 


all_annuus_good<-all_annuus[order(unlist(sapply(all_annuus$var, function(x) which(env_vars == x)))),]
#all_annuus_good<-na.omit(all_annuus_good)

ggsave("dotplot_H.annuus_nnumeber_top_candidate_windows_all_vars_analyses.pdf", width = 18, height = 8)
p<-ggplot(data= all_annuus_good, mapping = aes(x = reorder(var, top_candidate_count), top_candidate_count)) +
geom_point(aes(color = analysis), size = 2.7, alpha = 0.4)+scale_x_discrete(guide = guide_axis(angle = 60))+theme(text = element_text(size=10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 scale_color_manual(values=c("#D95F02", "#1B9E77"))

p + aes(x = fct_inorder(var))
print(p)
dev.off()



ggsave("dotplot_H.annuus_asending_Order_nnumeber_top_candidate_windows_all_vars_analyses.pdf", width = 18, height = 8)
p<-ggplot(data= all_annuus_good, mapping = aes(x = reorder(var, top_candidate_count), top_candidate_count)) +
geom_point(aes(color = analysis), size = 2.7, alpha = 0.4)+scale_x_discrete(guide = guide_axis(angle = 60))+theme(text = element_text(size=10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 scale_color_manual(values=c("#D95F02", "#1B9E77"))

#p + aes(x = fct_inorder(var))
print(p)
dev.off()


############################                                                                                                                                 
###### H.argophyllus #######
############################

rm(list = ls())

library(ggplot2)
library(forcats)
setwd("~/Documents/input_data/tests_analysis_for_paper/top_candidate_windows/top_cadidates_count_per_variable")

## oreder of vars
env_vars<-c("latitude","longitude","elevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH","BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS","Days_to_budding","Disk_diameter","Distance_of_first_branching_from_ground","DTF","Flower_FHDD_ratio","Guides_3_petals","Guides_individual_petal","Inflorescence_diameter","Internode_length","Leaf_area","Leaf_circular","Leaf_C_N_ratio","Leaf_curved_height","Leaf_curvedHeight_maxWidth","Leaf_curved_shape_index","Leaf_distal_eccentricity","Leaf_eccentricity_area_index","Leaf_ellipsoid","Leaf_height_mid-width","Leaf_maximum_height","Leaf_maximum_width","Leaf_obovoid","Leaf_perimeter","Leaf_proximal_eccentricity","Leaf_rectangular","Leaf_shape_index_external_I","Leaf_shape_index_external_II","Leaf_shape_index_internal","Leaf_total_C","Leaf_total_N","Leaf_width_mid-height","Leaf_width_widest_pos","Ligule_length","Ligule_LW_ratio","Ligules_number","Ligu le_width","LIR","Plant_height_at_flowering","Primary_branches","RGB_proportion_blue","RGB_proportion_green","RGB_proportion_red","Seed_area","Seed_circular","Seed_curved_height","Seed_curved_shape_index","Seed_distal_eccentricity","Seed_eccentricity","Seed_eccentricity_area_index","Seed_ellipsoid","Seed_height_mid_width","Seed_HW_ratio","Seed_maximum_height","Seed_maximum_width","Seed_ovoid","Seed_perimeter","Seed_proximal_eccentricity","Seed_rectangular","Seed_shape_index_external_I","Seed_shape_index_external_II","Seed_shape_index_internal","Seed_width_mid_height","Seed_width_widest_pos","SLA","Stem_diameter_at_flowering","Stem_diameter_final_after_5th_node","Stem_diameter_final_before_1st_node","TLN","Total_RGB")




arg_env_correc<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.argophyllus_BayPass_Corrected.table", header = TRUE)
arg_env_correc$analysis<-rep("corrected")
#ann_env_correc$species<-rep("H.annuus")
arg_env_correc$var<-gsub("llevation","elevation",arg_env_correc$var)

arg_env_uncorrec<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.argophyllus_Spearman.table", header = TRUE)
arg_env_uncorrec$analysis<-rep("uncorrected")
#ann_env_uncorrec$species<-rep("H.annuus")
arg_env_uncorrec$var<-gsub("llevation","elevation",arg_env_uncorrec$var)

all_arg_env<-do.call("rbind", list(arg_env_correc,arg_env_uncorrec)) 
#all_annuus_env_ordered<-all_annuus_env[order(all_annuus_env$top_candidate_count),]

#row.names(all_arg_env_ordered)<- all_arg_env_ordered$var
#all_arg_env_ordered$var<-NULL


arg_pheno_correc<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.argophyllus_Phenotype_Corrected.table", header = TRUE)
arg_pheno_correc$analysis<-rep("corrected")
#arg_pheno_correc$species<-rep("H.annuus")
#arg_pheno_correc<-arg_pheno_correc[!grepl("LIR", arg_pheno_correc$var),]


arg_pheno_uncorrec<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.argophyllus_Phenotype_Uncorrected.table", header = TRUE)
arg_pheno_uncorrec$analysis<-rep("uncorrected")
#arg_pheno_uncorrec$species<-rep("H.annuus")
#arg_pheno_uncorrec<-arg_pheno_uncorrec[!grepl("LIR", arg_pheno_uncorrec$var),]


all_arg <- do.call("rbind", list(arg_env_correc, arg_env_uncorrec,arg_pheno_correc,arg_pheno_uncorrec)) 

all_arg_good<-all_arg[order(unlist(sapply(all_arg$var, function(x) which(env_vars == x)))),]
all_arg_good<-na.omit(all_arg_good)

ggsave("dotplot_H.argophyllus_nnumeber_top_candidate_windows_all_vars_analyses.pdf", width = 18, height = 8)

p<-ggplot(data= all_arg_good, mapping = aes(x = reorder(var, top_candidate_count), top_candidate_count)) +
geom_point(aes(color = analysis), size = 2.7, alpha= 0.4)+scale_x_discrete(guide = guide_axis(angle = 60))+theme(text = element_text(size=10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 scale_color_manual(values=c("#D95F02", "#1B9E77"))

p + aes(x = fct_inorder(var))
print(p)
dev.off()



ggsave("dotplot_H.argophyllus_asending_Order_nnumeber_top_candidate_windows_all_vars_analyses.pdf", width = 18, height = 8)

p<-ggplot(data= all_arg_good, mapping = aes(x = reorder(var, top_candidate_count), top_candidate_count)) +
geom_point(aes(color = analysis), size = 2.7, alpha= 0.4)+scale_x_discrete(guide = guide_axis(angle = 60))+theme(text = element_text(size=10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 scale_color_manual(values=c("#D95F02", "#1B9E77"))

#p + aes(x = fct_inorder(var))
print(p)
dev.off()


##################################                                                                                                                                 
###### H.petiolaris.fallax #######
##################################
rm(list = ls())

library(ggplot2)
library(forcats)
setwd("~/Documents/input_data/tests_analysis_for_paper/top_candidate_windows/top_cadidates_count_per_variable")

## oreder of vars
env_vars<-c("latitude","longitude","elevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH","BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS","Days_to_budding","Disk_diameter","Distance_of_first_branching_from_ground","DTF","Flower_FHDD_ratio","Guides_3_petals","Guides_individual_petal","Inflorescence_diameter","Internode_length","Leaf_area","Leaf_circular","Leaf_C_N_ratio","Leaf_curved_height","Leaf_curvedHeight_maxWidth","Leaf_curved_shape_index","Leaf_distal_eccentricity","Leaf_eccentricity_area_index","Leaf_ellipsoid","Leaf_height_mid-width","Leaf_maximum_height","Leaf_maximum_width","Leaf_obovoid","Leaf_perimeter","Leaf_proximal_eccentricity","Leaf_rectangular","Leaf_shape_index_external_I","Leaf_shape_index_external_II","Leaf_shape_index_internal","Leaf_total_C","Leaf_total_N","Leaf_width_mid-height","Leaf_width_widest_pos","Ligule_length","Ligule_LW_ratio","Ligules_number","Ligu le_width","LIR","Plant_height_at_flowering","Primary_branches","RGB_proportion_blue","RGB_proportion_green","RGB_proportion_red","Seed_area","Seed_circular","Seed_curved_height","Seed_curved_shape_index","Seed_distal_eccentricity","Seed_eccentricity","Seed_eccentricity_area_index","Seed_ellipsoid","Seed_height_mid_width","Seed_HW_ratio","Seed_maximum_height","Seed_maximum_width","Seed_ovoid","Seed_perimeter","Seed_proximal_eccentricity","Seed_rectangular","Seed_shape_index_external_I","Seed_shape_index_external_II","Seed_shape_index_internal","Seed_width_mid_height","Seed_width_widest_pos","SLA","Stem_diameter_at_flowering","Stem_diameter_final_after_5th_node","Stem_diameter_final_before_1st_node","TLN","Total_RGB")



## annus
pf_env_correc<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.fallx_BayPass_Corrected.table", header = TRUE)
pf_env_correc$analysis<-rep("corrected")
pf_env_correc$species<-rep("H.pet.fallax")
pf_env_correc$var<-gsub("llevation","elevation",pf_env_correc$var)


pf_env_uncorrec<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.fallx_Spearman.table", header = TRUE)
pf_env_uncorrec$analysis<-rep("uncorrected")
pf_env_uncorrec$species<-rep("H.pet.fallax")
pf_env_uncorrec$var<-gsub("llevation","elevation",pf_env_correc$var)

all_annuus_env<-do.call("rbind", list(pf_env_correc,pf_env_uncorrec)) 
#all_annuus_env_ordered<-all_annuus_env[order(all_annuus_env$top_candidate_count),]
#row.names(all_annuus_env)<- all_annuus_env$var
#all_annuus_env$var<-NULL


pf_pheno_correc<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.fallax_Phenotype_Corrected.table", header = TRUE)
pf_pheno_correc$analysis<-rep("corrected")
pf_pheno_correc$species<-rep("H.pet.fallax")
#pf_pheno_correc<-pf_pheno_correc[!grepl("LIR", pf_pheno_correc$var),]


pf_pheno_uncorrec<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.fallax_Phenotype_Uncorrected.table", header = TRUE)
pf_pheno_uncorrec$analysis<-rep("uncorrected")
pf_pheno_uncorrec$species<-rep("H.pet.fallax")
#pf_pheno_uncorrec<-pf_pheno_uncorrec[!grepl("LIR", pf_pheno_uncorrec$var),]


all_pf <- do.call("rbind", list(pf_env_correc, pf_env_uncorrec,pf_pheno_correc,pf_pheno_uncorrec)) 


all_pf_good<-all_pf[order(unlist(sapply(all_pf$var, function(x) which(env_vars == x)))),]
all_pf_good<-na.omit(all_pf_good)

ggsave("dotplot_H.pet.fallax_nnumeber_top_candidate_windows_all_vars_analyses.pdf", width = 18, height = 8)
p<-ggplot(data= all_pf_good, mapping = aes(x = reorder(var, top_candidate_count), top_candidate_count)) +
geom_point(aes(color = analysis), size = 2.7,alpha=0.4)+scale_x_discrete(guide = guide_axis(angle = 60))+theme(text = element_text(size=10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 scale_color_manual(values=c("#D95F02", "#1B9E77"))
p + aes(x = fct_inorder(var))
print(p)
dev.off()


ggsave("dotplot_H.pet.fallax_asending_Ordered_nnumeber_top_candidate_windows_all_vars_analyses.pdf", width = 18, height = 8)
p<-ggplot(data= all_pf_good, mapping = aes(x = reorder(var, top_candidate_count), top_candidate_count)) +
geom_point(aes(color = analysis), size = 2.7,alpha=0.4)+scale_x_discrete(guide = guide_axis(angle = 60))+theme(text = element_text(size=10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 scale_color_manual(values=c("#D95F02", "#1B9E77"))
#p + aes(x = fct_inorder(var))
print(p)
dev.off()

##################################                                                                                                                                 
###### H.petiolaris.petiolaris #######
##################################
rm(list = ls())

library(ggplot2)
library(forcats)
setwd("~/Documents/input_data/tests_analysis_for_paper/top_candidate_windows/top_cadidates_count_per_variable")

## oreder of vars
env_vars<-c("latitude","longitude","elevation","MAT","MWMT","MCMT","TD","MAP","MSP","AHM","SHM","DD_0","DD5","DD_18","DD18","NFFD","bFFP","eFFP","FFP","PAS","EMT","EXT","Eref","CMD","MAR","RH","BICARB","CA","CEC","K","MG","OM","P1","P2","PERCENT_CA","PERCENT_K","PERCENT_MG","PERCENT_NA","PH","Sodium","SOL_SALTS","Days_to_budding","Disk_diameter","Distance_of_first_branching_from_ground","DTF","Flower_FHDD_ratio","Guides_3_petals","Guides_individual_petal","Inflorescence_diameter","Internode_length","Leaf_area","Leaf_circular","Leaf_C_N_ratio","Leaf_curved_height","Leaf_curvedHeight_maxWidth","Leaf_curved_shape_index","Leaf_distal_eccentricity","Leaf_eccentricity_area_index","Leaf_ellipsoid","Leaf_height_mid-width","Leaf_maximum_height","Leaf_maximum_width","Leaf_obovoid","Leaf_perimeter","Leaf_proximal_eccentricity","Leaf_rectangular","Leaf_shape_index_external_I","Leaf_shape_index_external_II","Leaf_shape_index_internal","Leaf_total_C","Leaf_total_N","Leaf_width_mid-height","Leaf_width_widest_pos","Ligule_length","Ligule_LW_ratio","Ligules_number","Ligu le_width","LIR","Plant_height_at_flowering","Primary_branches","RGB_proportion_blue","RGB_proportion_green","RGB_proportion_red","Seed_area","Seed_circular","Seed_curved_height","Seed_curved_shape_index","Seed_distal_eccentricity","Seed_eccentricity","Seed_eccentricity_area_index","Seed_ellipsoid","Seed_height_mid_width","Seed_HW_ratio","Seed_maximum_height","Seed_maximum_width","Seed_ovoid","Seed_perimeter","Seed_proximal_eccentricity","Seed_rectangular","Seed_shape_index_external_I","Seed_shape_index_external_II","Seed_shape_index_internal","Seed_width_mid_height","Seed_width_widest_pos","SLA","Stem_diameter_at_flowering","Stem_diameter_final_after_5th_node","Stem_diameter_final_before_1st_node","TLN","Total_RGB")


## annus
pp_env_correc<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.petiolaris_BayPass_Corrected.table", header = TRUE)
pp_env_correc$analysis<-rep("corrected")
pp_env_correc$species<-rep("H.pet.petiolaris")
pp_env_correc$var<-gsub("llevation","elevation",pp_env_correc$var)


pp_env_uncorrec<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.petiolaris_Spearman.table", header = TRUE)
pp_env_uncorrec$analysis<-rep("uncorrected")
pp_env_uncorrec$species<-rep("H.pet.petiolaris")
pp_env_uncorrec$var<-gsub("llevation","elevation",pp_env_correc$var)

all_pp_env<-do.call("rbind", list(pp_env_correc,pp_env_uncorrec)) 
#all_annuus_env_ordered<-all_annuus_env[order(all_annuus_env$top_candidate_count),]
#row.names(all_annuus_env)<- all_annuus_env$var
#all_annuus_env$var<-NULL


pp_pheno_correc<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.petiolaris_Phenotype_Corrected.table", header = TRUE)
pp_pheno_correc$analysis<-rep("corrected")
pp_pheno_correc$species<-rep("H.pet.petiolaris")
#pp_pheno_correc<-pp_pheno_correc[!grepl("LIR", pp_pheno_correc$var),]


pp_pheno_uncorrec<-read.table(file ="total_top_candidate_winows_count_per_variable_Ordered_H.petiolaris.petiolaris_Phenotype_Uncorrected.table", header = TRUE)
pp_pheno_uncorrec$analysis<-rep("uncorrected")
pp_pheno_uncorrec$species<-rep("H.pet.petiolaris")
#pp_pheno_uncorrec<-pp_pheno_uncorrec[!grepl("LIR", pp_pheno_uncorrec$var),]


all_pp <- do.call("rbind", list(pp_env_correc, pp_env_uncorrec,pp_pheno_correc,pp_pheno_uncorrec)) 

all_pp_good<-all_pp[order(unlist(sapply(all_pp$var, function(x) which(env_vars == x)))),]
all_pp_good<-na.omit(all_pp_good)

ggsave("dotplot_H.pet.petiolaris_nnumeber_top_candidate_windows_all_vars_analyses.pdf", width = 18, height = 8)
p<-ggplot(data= all_pp_good, mapping = aes(x = reorder(var, top_candidate_count), top_candidate_count)) +
geom_point(aes(color = analysis), size = 2.7,alpha=0.4)+scale_x_discrete(guide = guide_axis(angle = 60))+theme(text = element_text(size=10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 scale_color_manual(values=c("#D95F02", "#1B9E77"))

p + aes(x = fct_inorder(var))
print(p)
dev.off()


ggsave("dotplot_H.pet.petiolaris_asending_Ordered_nnumeber_top_candidate_windows_all_vars_analyses.pdf", width = 18, height = 8)
p<-ggplot(data= all_pp_good, mapping = aes(x = reorder(var, top_candidate_count), top_candidate_count)) +
geom_point(aes(color = analysis), size = 2.7,alpha=0.4)+scale_x_discrete(guide = guide_axis(angle = 60))+theme(text = element_text(size=10))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 scale_color_manual(values=c("#D95F02", "#1B9E77"))

#p + aes(x = fct_inorder(var))
print(p)
dev.off()





