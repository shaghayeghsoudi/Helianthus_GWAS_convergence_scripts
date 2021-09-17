

#### analyzing convergent wondows ####
#rm(list = ls())
setwd("/data/users/soudi/paper_sunflowers/null_W")

### load climate and soil vars Corrected #####
## load corrected table
null_corr<-read.table(file= "convergent_2binning_soil_climate_BayPass.table", header = FALSE)
colnames(null_corr)<-c("data","type","analysis","Taxa1","Taxa2","comparison","variable","window_name","window5","Emp_p")

### break the main corercted file into each comparison type
null_corr_annarg<-null_corr[null_corr$comparison=="Annuus_Argophyllus",]
null_corr_annpf<-null_corr[null_corr$comparison=="Annuus_petfal",]
null_corr_annpp<-null_corr[null_corr$comparison=="Annuus_petpet",]
null_corr_argpf<-null_corr[null_corr$comparison=="Argophyllus_petfal",]
null_corr_argpp<-null_corr[null_corr$comparison=="Argophyllus_petpet",]
null_corr_pfpp<-null_corr[null_corr$comparison=="petfal_petpet",]

### find unique windows
null_corr_annarg_UNIQ<-null_corr_annarg[!duplicated(null_corr_annarg$window_name), ]
null_corr_annarg_UNIQ$status<-rep("corrected")

null_corr_annpf_UNIQ<-null_corr_annpf[!duplicated(null_corr_annpf$window_name), ]
null_corr_annpf_UNIQ$status<-rep("corrected")

null_corr_annpp_UNIQ<-null_corr_annpp[!duplicated(null_corr_annpp$window_name), ]
null_corr_annpp_UNIQ$status<-rep("corrected")

null_corr_argpf_UNIQ<-null_corr_argpf[!duplicated(null_corr_argpf$window_name), ]
null_corr_argpf_UNIQ$status<-rep("corrected")

null_corr_argpp_UNIQ<-null_corr_argpp[!duplicated(null_corr_argpp$window_name), ]
null_corr_argpp_UNIQ$status<-rep("corrected")

null_corr_pfpp_UNIQ<-null_corr_pfpp[!duplicated(null_corr_pfpp$window_name), ]
null_corr_pfpp_UNIQ$status<-rep("corrected")

all_climate_corrected<-do.call("rbind", list(null_corr_annarg_UNIQ,null_corr_annpf_UNIQ,null_corr_annpp_UNIQ,null_corr_argpf_UNIQ,null_corr_argpp_UNIQ,null_corr_pfpp_UNIQ))

all_climate_corrected<-all_climate_corrected[,c("window_name","comparison","status")]


aa<-merge(null_corr_annarg_UNIQ,null_corr_annpf_UNIQ, by.x = "window_name", by.y = "window_name")
bb<-merge(null_corr_annarg_UNIQ,null_corr_annpp_UNIQ, by.x = "window_name", by.y = "window_name")
cc<-merge(null_corr_annarg_UNIQ,null_corr_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
dd<-merge(null_corr_annarg_UNIQ,null_corr_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ee<-merge(null_corr_annarg_UNIQ,null_corr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


ff<-merge(null_corr_annpf_UNIQ,null_corr_annpp_UNIQ, by.x = "window_name", by.y = "window_name")
gg<-merge(null_corr_annpf_UNIQ,null_corr_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
hh<-merge(null_corr_annpf_UNIQ,null_corr_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ii<-merge(null_corr_annpf_UNIQ,null_corr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


jj<-merge(null_corr_annpp_UNIQ,null_corr_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
kk<-merge(null_corr_annpp_UNIQ,null_corr_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ll<-merge(null_corr_annpp_UNIQ,null_corr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")

mm<-merge(null_corr_argpf_UNIQ,null_corr_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
nn<-merge(null_corr_argpf_UNIQ,null_corr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")

oo<-merge(null_corr_argpp_UNIQ,null_corr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


##################################
##### load uncorected table ######
##################################
null_uncorr<-read.table(file = "convergent_2binning_soil_climate_Spearman.table", header = FALSE)
colnames(null_uncorr)<-c("data","type","analysis","Taxa1","Taxa2","comparison","variable","window_name","window5","Emp_p")


### break the main Uncorercted file into each comparison type
null_Uncorr_annarg<-null_uncorr[null_uncorr$comparison=="Annuus_Argophyllus",]
null_Uncorr_annpf<-null_uncorr[null_uncorr$comparison=="Annuus_petfal",]
null_Uncorr_annpp<-null_uncorr[null_uncorr$comparison=="Annuus_petpet",]
null_Uncorr_argpf<-null_uncorr[null_uncorr$comparison=="Argophyllus_petfal",]
null_Uncorr_argpp<-null_uncorr[null_uncorr$comparison=="Argophyllus_petpet",]
null_Uncorr_pfpp<-null_uncorr[null_uncorr$comparison=="petfal_petpet",]

### find unique windows
null_Uncorr_annarg_UNIQ<-null_Uncorr_annarg[!duplicated(null_Uncorr_annarg$window_name), ]
null_Uncorr_annarg_UNIQ$status<-rep("Uncorrected")

null_Uncorr_annpf_UNIQ<-null_Uncorr_annpf[!duplicated(null_Uncorr_annpf$window_name), ]
null_Uncorr_annpf_UNIQ$status<-rep("Uncorrected")

null_Uncorr_annpp_UNIQ<-null_Uncorr_annpp[!duplicated(null_Uncorr_annpp$window_name), ]
null_Uncorr_annpp_UNIQ$status<-rep("Uncorrected")

null_Uncorr_argpf_UNIQ<-null_Uncorr_argpf[!duplicated(null_Uncorr_argpf$window_name), ]
null_Uncorr_argpf_UNIQ$status<-rep("Uncorrected")

null_Uncorr_argpp_UNIQ<-null_Uncorr_argpp[!duplicated(null_Uncorr_argpp$window_name), ]
null_Uncorr_argpp_UNIQ$status<-rep("Uncorrected")

null_Uncorr_pfpp_UNIQ<-null_Uncorr_pfpp[!duplicated(null_Uncorr_pfpp$window_name), ]
null_Uncorr_pfpp_UNIQ$status<-rep("Uncorrected")

all_climate_Uncorrected<-do.call("rbind", list(null_Uncorr_annarg_UNIQ,null_Uncorr_annpf_UNIQ,null_Uncorr_annpp_UNIQ,null_Uncorr_argpf_UNIQ,null_Uncorr_argpp_UNIQ,null_Uncorr_pfpp_UNIQ))

all_climate_Uncorrected<-all_climate_Uncorrected[,c("window_name","comparison","status")]


aa<-merge(null_Uncorr_annarg_UNIQ,null_Uncorr_annpf_UNIQ, by.x = "window_name", by.y = "window_name")
bb<-merge(null_Uncorr_annarg_UNIQ,null_Uncorr_annpp_UNIQ, by.x = "window_name", by.y = "window_name")
cc<-merge(null_Uncorr_annarg_UNIQ,null_Uncorr_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
dd<-merge(null_Uncorr_annarg_UNIQ,null_Uncorr_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ee<-merge(null_Uncorr_annarg_UNIQ,null_Uncorr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


ff<-merge(null_Uncorr_annpf_UNIQ,null_Uncorr_annpp_UNIQ, by.x = "window_name", by.y = "window_name")
gg<-merge(null_Uncorr_annpf_UNIQ,null_Uncorr_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
hh<-merge(null_Uncorr_annpf_UNIQ,null_Uncorr_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ii<-merge(null_Uncorr_annpf_UNIQ,null_Uncorr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


jj<-merge(null_Uncorr_annpp_UNIQ,null_Uncorr_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
kk<-merge(null_Uncorr_annpp_UNIQ,null_Uncorr_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ll<-merge(null_Uncorr_annpp_UNIQ,null_Uncorr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")

mm<-merge(null_Uncorr_argpf_UNIQ,null_Uncorr_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
nn<-merge(null_Uncorr_argpf_UNIQ,null_Uncorr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")

oo<-merge(null_Uncorr_argpp_UNIQ,null_Uncorr_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")



########################
## find intersects for each species corrected vs. uncorrected (climate) for each comparison type ####
########################

ann_arg_intersect <- data.frame("windows"=(intersect(null_corr_annarg_UNIQ$window_name, null_Uncorr_annarg_UNIQ$window_name)))
ann_arg_intersect$comparison<-rep("Annuus_Argophyllus")
ann_arg_intersect$status<-rep("intersect")


ann_pf_intersect <- data.frame("windows"=(intersect(null_corr_annpf_UNIQ$window_name, null_Uncorr_annpf_UNIQ$window_name)))
ann_pf_intersect$comparison<-rep("Annuus_petfal")
ann_pf_intersect$status<-rep("intersect")


ann_pp_intersect <- data.frame("windows"=(intersect(null_corr_annpp_UNIQ$window_name, null_Uncorr_annpp_UNIQ$window_name)))
ann_pp_intersect$comparison<-rep("Annuus_petpet")
ann_pp_intersect$status<-rep("intersect")


arg_pf_intersect <- data.frame("windows"=(intersect(null_corr_argpf_UNIQ$window_name, null_Uncorr_argpf_UNIQ$window_name)))
arg_pf_intersect $comparison<-rep("Argophyllus_petfal")
arg_pf_intersect$status<-rep("intersect")


arg_pp_intersect <- data.frame("windows"=(intersect(null_corr_argpp_UNIQ$window_name, null_Uncorr_argpp_UNIQ$window_name)))
arg_pp_intersect $comparison<-rep("Argophyllus_petpet")
arg_pp_intersect$status<-rep("intersect")


pf_pp_intersect <- data.frame("windows"=(intersect(null_corr_pfpp_UNIQ$window_name, null_Uncorr_pfpp_UNIQ$window_name)))
pf_pp_intersect $comparison<-rep("petfal_petpet")
pf_pp_intersect$status<-rep("intersect")


all_intersect_input<-do.call("rbind", list(ann_arg_intersect,ann_pf_intersect,ann_pp_intersect,arg_pf_intersect,arg_pp_intersect,pf_pp_intersect))
colnames(all_intersect_input)<-c("window_name", "comparison","status")


### combine corrected, uncorrected and both 
all_climate_CorrUncorrCBoth<-do.call("rbind", list(all_climate_corrected,all_climate_Uncorrected,all_intersect_input))
write.table(all_climate_CorrUncorrCBoth, file = "convergent_sig_windows_climate_Corrected_Uncorercted_Intersected.table", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



#####################################
########### phenotype data ##########
#####################################
## load corrected table
null_corr_ph<-read.table(file= "convergent_2binning_phenotype_GWAS_corrected", header = FALSE)
colnames(null_corr_ph)<-c("data","type","analysis","Taxa1","Taxa2","comparison","variable","window_name","window5","Emp_p")

### break the main corercted file into each comparison type
null_corrP_annarg<-null_corr_ph[null_corr_ph$comparison=="Annuus_Argophyllus",]
null_corrP_annpf<-null_corr_ph[null_corr_ph$comparison=="Annuus_petfal",]
null_corrP_annpp<-null_corr_ph[null_corr_ph$comparison=="Annuus_petpet",]
null_corrP_argpf<-null_corr_ph[null_corr_ph$comparison=="Argophyllus_petfal",]
null_corrP_argpp<-null_corr_ph[null_corr_ph$comparison=="Argophyllus_petpet",]
null_corrP_pfpp<-null_corr_ph[null_corr_ph$comparison=="petfal_petpet",]

### find unique windows
null_corrP_annarg_UNIQ<-null_corrP_annarg[!duplicated(null_corrP_annarg$window_name), ]
null_corrP_annarg_UNIQ$status<-rep("corrected")

null_corrP_annpf_UNIQ<-null_corrP_annpf[!duplicated(null_corrP_annpf$window_name), ]
null_corrP_annpf_UNIQ$status<-rep("corrected")

null_corrP_annpp_UNIQ<-null_corrP_annpp[!duplicated(null_corrP_annpp$window_name), ]
null_corrP_annpp_UNIQ$status<-rep("corrected")

null_corrP_argpf_UNIQ<-null_corrP_argpf[!duplicated(null_corrP_argpf$window_name), ]
null_corrP_argpf_UNIQ$status<-rep("corrected")

null_corrP_argpp_UNIQ<-null_corrP_argpp[!duplicated(null_corrP_argpp$window_name), ]
null_corrP_argpp_UNIQ$status<-rep("corrected")

null_corrP_pfpp_UNIQ<-null_corrP_pfpp[!duplicated(null_corrP_pfpp$window_name), ]
null_corrP_pfpp_UNIQ$status<-rep("corrected")

all_phenotype_corrected<-do.call("rbind", list(null_corrP_annarg_UNIQ,null_corrP_annpf_UNIQ,null_corrP_annpp_UNIQ,null_corrP_argpf_UNIQ,null_corrP_argpp_UNIQ,null_corrP_pfpp_UNIQ))

all_phenotype_corrected<-all_phenotype_corrected[,c("window_name","comparison","status")]



aa<-merge(null_corrP_annarg_UNIQ,null_corrP_annpf_UNIQ, by.x = "window_name", by.y = "window_name")
bb<-merge(null_corrP_annarg_UNIQ,null_corrP_annpp_UNIQ, by.x = "window_name", by.y = "window_name")
cc<-merge(null_corrP_annarg_UNIQ,null_corrP_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
dd<-merge(null_corrP_annarg_UNIQ,null_corrP_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ee<-merge(null_corrP_annarg_UNIQ,null_corrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


ff<-merge(null_corrP_annpf_UNIQ,null_corrP_annpp_UNIQ, by.x = "window_name", by.y = "window_name")
gg<-merge(null_corrP_annpf_UNIQ,null_corrP_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
hh<-merge(null_corrP_annpf_UNIQ,null_corrP_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ii<-merge(null_corrP_annpf_UNIQ,null_corrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


jj<-merge(null_corrP_annpp_UNIQ,null_corrP_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
kk<-merge(null_corrP_annpp_UNIQ,null_corrP_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ll<-merge(null_corrP_annpp_UNIQ,null_corrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")

mm<-merge(null_corrP_argpf_UNIQ,null_corrP_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
nn<-merge(null_corrP_argpf_UNIQ,null_corrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")

oo<-merge(null_corrP_argpp_UNIQ,null_corrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")




##################################
##### load uncorected table ######
##################################
null_uncorrP<-read.table(file = "convergent_2binning_phenotype_GWAS_UNcorrected", header = FALSE)
colnames(null_uncorrP)<-c("data","type","analysis","Taxa1","Taxa2","comparison","variable","window_name","window5","Emp_p")


### break the main Uncorercted file into each comparison type
null_UncorrP_annarg<-null_uncorrP[null_uncorrP$comparison=="Annuus_Argophyllus",]
null_UncorrP_annpf<-null_uncorrP[null_uncorrP$comparison=="Annuus_petfal",]
null_UncorrP_annpp<-null_uncorrP[null_uncorrP$comparison=="Annuus_petpet",]
null_UncorrP_argpf<-null_uncorrP[null_uncorrP$comparison=="Argophyllus_petfal",]
null_UncorrP_argpp<-null_uncorrP[null_uncorrP$comparison=="Argophyllus_petpet",]
null_UncorrP_pfpp<-null_uncorrP[null_uncorrP$comparison=="petfal_petpet",]

### find unique windows
null_UncorrP_annarg_UNIQ<-null_UncorrP_annarg[!duplicated(null_UncorrP_annarg$window_name), ]
null_UncorrP_annarg_UNIQ$status<-rep("Uncorrected")

null_UncorrP_annpf_UNIQ<-null_UncorrP_annpf[!duplicated(null_UncorrP_annpf$window_name), ]
null_UncorrP_annpf_UNIQ$status<-rep("Uncorrected")

null_UncorrP_annpp_UNIQ<-null_UncorrP_annpp[!duplicated(null_UncorrP_annpp$window_name), ]
null_UncorrP_annpp_UNIQ$status<-rep("Uncorrected")

null_UncorrP_argpf_UNIQ<-null_UncorrP_argpf[!duplicated(null_UncorrP_argpf$window_name), ]
null_UncorrP_argpf_UNIQ$status<-rep("Uncorrected")

null_UncorrP_argpp_UNIQ<-null_UncorrP_argpp[!duplicated(null_UncorrP_argpp$window_name), ]
null_UncorrP_argpp_UNIQ$status<-rep("Uncorrected")

null_UncorrP_pfpp_UNIQ<-null_UncorrP_pfpp[!duplicated(null_UncorrP_pfpp$window_name), ]
null_UncorrP_pfpp_UNIQ$status<-rep("Uncorrected")

all_phenotype_Uncorrected<-do.call("rbind", list(null_UncorrP_annarg_UNIQ,null_UncorrP_annpf_UNIQ,null_UncorrP_annpp_UNIQ,null_UncorrP_argpf_UNIQ,null_UncorrP_argpp_UNIQ,null_UncorrP_pfpp_UNIQ))

all_phenotype_Uncorrected<-all_phenotype_Uncorrected[,c("window_name","comparison","status")]



aa<-merge(null_UncorrP_annarg_UNIQ,null_UncorrP_annpf_UNIQ, by.x = "window_name", by.y = "window_name")
bb<-merge(null_UncorrP_annarg_UNIQ,null_UncorrP_annpp_UNIQ, by.x = "window_name", by.y = "window_name")
cc<-merge(null_UncorrP_annarg_UNIQ,null_UncorrP_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
dd<-merge(null_UncorrP_annarg_UNIQ,null_UncorrP_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ee<-merge(null_UncorrP_annarg_UNIQ,null_UncorrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


ff<-merge(null_UncorrP_annpf_UNIQ,null_UncorrP_annpp_UNIQ, by.x = "window_name", by.y = "window_name")
gg<-merge(null_UncorrP_annpf_UNIQ,null_UncorrP_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
hh<-merge(null_UncorrP_annpf_UNIQ,null_UncorrP_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ii<-merge(null_UncorrP_annpf_UNIQ,null_UncorrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")


jj<-merge(null_UncorrP_annpp_UNIQ,null_UncorrP_argpf_UNIQ, by.x = "window_name", by.y = "window_name")
kk<-merge(null_UncorrP_annpp_UNIQ,null_UncorrP_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
ll<-merge(null_UncorrP_annpp_UNIQ,null_UncorrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")

mm<-merge(null_UncorrP_argpf_UNIQ,null_UncorrP_argpp_UNIQ, by.x = "window_name", by.y = "window_name")
nn<-merge(null_UncorrP_argpf_UNIQ,null_UncorrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")

oo<-merge(null_UncorrP_argpp_UNIQ,null_UncorrP_pfpp_UNIQ, by.x = "window_name", by.y = "window_name")



########################
## find intersects for each species corrected vs. uncorrected (climate) for each comparison type ####
########################

ann_arg_intersect <- data.frame("windows"=(intersect(null_corrP_annarg_UNIQ$window_name, null_UncorrP_annarg_UNIQ$window_name)))
ann_arg_intersect$comparison<-rep("Annuus_Argophyllus")
ann_arg_intersect$status<-rep("intersect")


ann_pf_intersect <- data.frame("windows"=(intersect(null_corrP_annpf_UNIQ$window_name, null_UncorrP_annpf_UNIQ$window_name)))
ann_pf_intersect$comparison<-rep("Annuus_petfal")
ann_pf_intersect$status<-rep("intersect")


ann_pp_intersect <- data.frame("windows"=(intersect(null_corrP_annpp_UNIQ$window_name, null_UncorrP_annpp_UNIQ$window_name)))
ann_pp_intersect$comparison<-rep("Annuus_petpet")
ann_pp_intersect$status<-rep("intersect")


arg_pf_intersect <- data.frame("windows"=(intersect(null_corrP_argpf_UNIQ$window_name, null_UncorrP_argpf_UNIQ$window_name)))
arg_pf_intersect $comparison<-rep("Argophyllus_petfal")
arg_pf_intersect$status<-rep("intersect")


arg_pp_intersect <- data.frame("windows"=(intersect(null_corrP_argpp_UNIQ$window_name, null_UncorrP_argpp_UNIQ$window_name)))
arg_pp_intersect $comparison<-rep("Argophyllus_petpet")
arg_pp_intersect$status<-rep("intersect")


pf_pp_intersect <- data.frame("windows"=(intersect(null_corrP_pfpp_UNIQ$window_name, null_UncorrP_pfpp_UNIQ$window_name)))
pf_pp_intersect $comparison<-rep("petfal_petpet")
pf_pp_intersect$status<-rep("intersect")


all_intersect_input<-do.call("rbind", list(ann_arg_intersect,ann_pf_intersect,ann_pp_intersect,arg_pf_intersect,arg_pp_intersect,pf_pp_intersect))
colnames(all_intersect_input)<-c("window_name", "comparison","status")


### combine corrected, uncorrected and both 
all_climate_CorrUncorrCBoth<-do.call("rbind", list(all_phenotype_corrected,all_phenotype_Uncorrected,all_intersect_input))
write.table(all_climate_CorrUncorrCBoth, file = "convergent_sig_windows_GWAS_Corrected_Uncorercted_Intersected.table", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


################################################################
#### assemble  both (climate and phenotype) plots together ####
# Grouped Bar Plot
setwd("~/Documents/input_data/tests_analysis_for_paper/null_w")
all_climate<-read.table(file = "convergent_sig_windows_climate_Corrected_Uncorercted_Intersected.table", header = TRUE)

counts_cli <- table(all_climate$status, all_climate$comparison)
counts_cli <- counts_cli[ c(3,1,2),]
colnames(counts_cli)<-c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")

###### Final plot from here (two plots assembled) #######


pdf (file = "barplots_convergent_windows_allchromstogether_recomadjusted_climate_phenotype_Intersection.pdf",width = 8, height= 8)
par(mfrow=c(1,2))
par(mar = c(12,4.8,5.1,4.1))
barplot(counts_cli, 
         col=c("#1B9E77","#D95F02","#7570B3"),
        beside=TRUE, ylim=c(0,2000),ylab = "number of convergent windows", cex.lab=1.5, cex.axis=1.5,width=c(0.03,0.03,0.03,0.03),cex = 0.9,las= 3)
#legend("topright", inset=c(-0.2,-0.01),
#       legend = c("both", "uncorrected","corrected"), cex = 0.75,y.intersp=1,
#       fill = c("#1B9E77","#D95F02","#7570B3"),bty="n")


all_pheno<-read.table(file = "convergent_sig_windows_GWAS_Corrected_Uncorercted_Intersected.table", header = TRUE)

counts_pheno <- table(all_pheno$status, all_pheno$comparison)
counts_pheno <- counts_pheno[c(3,1,2),]
colnames(counts_pheno)<-c("H.annuus-H.argophyllus","H.annuus_H.pet.falax","H.annuus_H.pet.petiolaris","H.argophyllus_H.pet.fallax","H.argophyllus_H.pet.petiolaris","H.pet.fallax-H.pet.petiolaris")

barplot(counts_pheno, 
         col=c("#1B9E77","#D95F02","#7570B3"),
        beside=TRUE, ylim=c(0,2000),ylab = "number of convergent windows", cex.lab=1.5, cex.axis=1.5,width=c(0.03,0.03,0.03,0.03),cex = 0.9,las=3)


dev.off()
#############################################





