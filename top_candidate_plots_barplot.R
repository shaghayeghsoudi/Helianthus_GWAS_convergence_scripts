
#####
## set working directory

setwd("/data/users/soudi/paper_sunflowers")

###################################################
## load corrected files (climate and soil)
###################################################
## annuus
ann_corr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_baypass_H.annuus", header = FALSE)
colnames(ann_corr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ann_corr<-ann_corr[!(ann_corr$outlier_count==0),]
ann_corr_top<-ann_corr[(ann_corr$outlier_count>ann_corr$p4),]
ann_corr_topUnique<-ann_corr_top[!duplicated(ann_corr_top$windows), ]
ann_corr_topUnique$status<-rep("corrected")

## argophyllus
arg_corr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_baypass_argophyllus", header = FALSE)
colnames(arg_corr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
arg_corr<-arg_corr[!(arg_corr$outlier_count==0),]
arg_corr_top<-arg_corr[(arg_corr$outlier_count>arg_corr$p4),]
arg_corr_topUnique<-arg_corr_top[!duplicated(arg_corr_top$windows), ]
arg_corr_topUnique$status<-rep("corrected")

## petfall
ppfal_corr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_baypass_H.petiolaris.fallax", header = FALSE)
colnames(ppfal_corr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppfal_corr<-ppfal_corr[!(ppfal_corr$outlier_count==0),]
ppfal_corr_top<-ppfal_corr[(ppfal_corr$outlier_count>ppfal_corr$p4),]
ppfal_corr_topUnique<-ppfal_corr_top[!duplicated(ppfal_corr_top$windows), ]
ppfal_corr_topUnique$status<-rep("corrected")

## pet.pet
ppet_corr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_baypass_H.petiolaris.Pet", header = FALSE)
colnames(ppet_corr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppet_corr<-ppet_corr[!(ppet_corr$outlier_count==0),]
ppet_corr_top<-ppet_corr[(ppet_corr$outlier_count>ppet_corr$p4),]
ppet_corr_topUnique<-ppet_corr_top[!duplicated(ppet_corr_top$windows), ]
ppet_corr_topUnique$status<-rep("corrected")


###################################################
### load uncorrected files (climate and soil) #####
###################################################
ann_Uncorr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_spearman_H.annuus", header = FALSE)
colnames(ann_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ann_Uncorr<-ann_Uncorr[!(ann_Uncorr$outlier_count==0),]
ann_Uncorr_top<-ann_Uncorr[(ann_Uncorr$outlier_count>ann_Uncorr$p4),]
ann_UncorrUnique<-ann_Uncorr_top[!duplicated(ann_Uncorr_top$windows), ]
ann_UncorrUnique$status<-rep("Uncorrected")


arg_Uncorr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_spearman_H.argophyllus", header = FALSE)
colnames(arg_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
arg_Uncorr<-arg_Uncorr[!(arg_Uncorr$outlier_count==0),]
arg_Uncorr_top<-arg_Uncorr[(arg_Uncorr$outlier_count>arg_Uncorr$p4),]
arg_Uncorr_topUnique<-arg_Uncorr_top[!duplicated(arg_Uncorr_top$windows), ]
arg_Uncorr_topUnique$status<-rep("Uncorrected")


ppfal_Uncorr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_spearman_H.petiolaris.fallax", header = FALSE)
colnames(ppfal_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppfal_Uncorr<-ppfal_Uncorr[!(ppfal_Uncorr$outlier_count==0),]
ppfal_Uncorr_top<-ppfal_Uncorr[(ppfal_Uncorr$outlier_count>ppfal_Uncorr$p4),]
ppfal_Uncorr_topUnique<-ppfal_Uncorr_top[!duplicated(ppfal_Uncorr_top$windows), ]
ppfal_Uncorr_topUnique$status<-rep("Uncorrected")



ppet_Uncorr<-read.table(file = "top_candidates_Mojtaba_climate/top_candidate_spearman_H.petiolaris.petiolaris", header = FALSE)
colnames(ppet_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppet_Uncorr<-ppet_Uncorr[!(ppet_Uncorr$outlier_count==0),]
ppet_Uncorr_top<-ppet_Uncorr[(ppet_Uncorr$outlier_count>ppet_Uncorr$p4),]
ppet_Uncorr_topUnique<-ppet_Uncorr_top[!duplicated(ppet_Uncorr_top$windows), ]
ppet_Uncorr_topUnique$status<-rep("Uncorrected")



## combine corrcetd and uncorected windows

ann_bothSpep<-do.call("rbind", list(ann_corr_topUnique, ann_UncorrUnique))
arg_bothSpep<- do.call("rbind", list(arg_corr_topUnique, arg_Uncorr_topUnique))
ppfal_bothSpep<-do.call("rbind", list(ppfal_corr_topUnique, ppfal_Uncorr_topUnique))
ppet_bothSpep<-do.call("rbind", list(ppet_corr_topUnique, ppet_Uncorr_topUnique))

all_sep<-do.call("rbind", list(ann_bothSpep,arg_bothSpep,ppfal_bothSpep,ppet_bothSpep))
all_sep_input<-all_sep[,c(2,8,9)]
colnames(all_sep_input)<-c("windows","species","status")

########
## find intersects for each species corrected vs. uncorrected
ann_intersect <- data.frame("windows"=(intersect(ann_corr_topUnique$windows, ann_UncorrUnique$windows)))
ann_intersect$species<-rep("H.annuus")
ann_intersect$status<-rep("intersect")


arg_intersect <- data.frame("windows"=(intersect(arg_corr_topUnique$windows, arg_Uncorr_topUnique$windows)))
arg_intersect$species<-rep("H.argophyllus")
arg_intersect$status<-rep("intersect")


pf_intersect <- data.frame("windows"=(intersect(ppfal_corr_topUnique$windows, ppfal_Uncorr_topUnique$windows)))
pf_intersect$species<-rep("H.petiolaris.fallax")
pf_intersect$status<-rep("intersect")


pp_intersect <- data.frame("windows"=(intersect(ppet_corr_topUnique$windows, ppet_Uncorr_topUnique$windows)))
pp_intersect $species<-rep("H.petiolaris.petiolaris")
pp_intersect$status<-rep("intersect")

all_intersect_input<-do.call("rbind", list(ann_intersect,arg_intersect,pf_intersect,pp_intersect))

#### merge all and intersect combinations
all<-do.call("rbind", list(all_sep_input,all_intersect_input))
write.table(all, file = "all_windows_corrceted_uncorrcted_GEA_with_intersection.table", col.names = TRUE, row.names= FALSE, quote= FALSE, sep = "\t")


#########################################################################################################
###### load phenotype data ######
#### load corrected files (phenotype) ######
############################################
#rm(list = ls())
setwd("/data/users/soudi/paper_sunflowers")
ann_corr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_corrected_H.annuus", header = FALSE)
colnames(ann_corr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ann_corr<-ann_corr[!(ann_corr$outlier_count==0),]
ann_corr_top<-ann_corr[(ann_corr$outlier_count>ann_corr$p4),]
ann_corr_topUnique<-ann_corr_top[!duplicated(ann_corr_top$windows), ]
ann_corr_topUnique$status<-rep("corrected")


arg_corr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_corrected_H.argophyllus", header = FALSE)
colnames(arg_corr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
arg_corr<-arg_corr[!(arg_corr$outlier_count==0),]
arg_corr_top<-arg_corr[(arg_corr$outlier_count>arg_corr$p4),]
arg_corr_topUnique<-arg_corr_top[!duplicated(arg_corr_top$windows), ]
arg_corr_topUnique$status<-rep("corrected")


ppfal_corr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_corrected_H.petiolaris.fallax", header = FALSE)
colnames(ppfal_corr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppfal_corr<-ppfal_corr[!(ppfal_corr$outlier_count==0),]
ppfal_corr_top<-ppfal_corr[(ppfal_corr$outlier_count>ppfal_corr$p4),]
ppfal_corr_topUnique<-ppfal_corr_top[!duplicated(ppfal_corr_top$windows), ]
ppfal_corr_topUnique$status<-rep("corrected")


ppet_corr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_corrected_H.petiolaris.petiolaris", header = FALSE)
colnames(ppet_corr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppet_corr<-ppet_corr[!(ppet_corr$outlier_count==0),]
ppet_corr_top<-ppet_corr[(ppet_corr$outlier_count>ppet_corr$p4),]
ppet_corr_topUnique<-ppet_corr_top[!duplicated(ppet_corr_top$windows), ]
ppet_corr_topUnique$status<-rep("corrected")



###################################################
### load uncorrected files (phenotype) #####
###################################################
ann_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_uncorrected_H.annus", header = FALSE)
colnames(ann_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ann_Uncorr<-ann_Uncorr[!(ann_Uncorr$outlier_count==0),]
ann_Uncorr_top<-ann_Uncorr[(ann_Uncorr$outlier_count>ann_Uncorr$p4),]
ann_UncorrUnique<-ann_Uncorr_top[!duplicated(ann_Uncorr_top$windows), ]
ann_UncorrUnique$status<-rep("Uncorrected")


arg_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_uncorrected_H.argophyllus", header = FALSE)
colnames(arg_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
arg_Uncorr<-arg_Uncorr[!(arg_Uncorr$outlier_count==0),]
arg_Uncorr_top<-arg_Uncorr[(arg_Uncorr$outlier_count>arg_Uncorr$p4),]
arg_Uncorr_topUnique<-arg_Uncorr_top[!duplicated(arg_Uncorr_top$windows), ]
arg_Uncorr_topUnique$status<-rep("Uncorrected")


ppfal_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_uncorrected_H.petiolaris.fallax", header = FALSE)
colnames(ppfal_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppfal_Uncorr<-ppfal_Uncorr[!(ppfal_Uncorr$outlier_count==0),]
ppfal_Uncorr_top<-ppfal_Uncorr[(ppfal_Uncorr$outlier_count>ppfal_Uncorr$p4),]
ppfal_Uncorr_topUnique<-ppfal_Uncorr_top[!duplicated(ppfal_Uncorr_top$windows), ]
ppfal_Uncorr_topUnique$status<-rep("Uncorrected")



ppet_Uncorr<-read.table(file = "top_candidates_Mojtaba_phenotype/top_candidate_GWAS_uncorrected_H.petiolaris.petiolaris", header = FALSE)
colnames(ppet_Uncorr)<-c("test_name","windows", "snp_count","outlier_count", "expected","p4","p8","pecies")
ppet_Uncorr<-ppet_Uncorr[!(ppet_Uncorr$outlier_count==0),]
ppet_Uncorr_top<-ppet_Uncorr[(ppet_Uncorr$outlier_count>ppet_Uncorr$p4),]
ppet_Uncorr_topUnique<-ppet_Uncorr_top[!duplicated(ppet_Uncorr_top$windows), ]
ppet_Uncorr_topUnique$status<-rep("Uncorrected")


## combine corrcetd and uncorected windows

ann_bothSpep<-do.call("rbind", list(ann_corr_topUnique, ann_UncorrUnique))
arg_bothSpep<- do.call("rbind", list(arg_corr_topUnique, arg_Uncorr_topUnique))
ppfal_bothSpep<-do.call("rbind", list(ppfal_corr_topUnique, ppfal_Uncorr_topUnique))
ppet_bothSpep<-do.call("rbind", list(ppet_corr_topUnique, ppet_Uncorr_topUnique))

all_sep<-do.call("rbind", list(ann_bothSpep,arg_bothSpep,ppfal_bothSpep,ppet_bothSpep))
all_sep_input<-all_sep[,c(2,8,9)]
colnames(all_sep_input)<-c("windows","species","status")


########
##### find intersects for each species corrected vs. uncorrected
ann_intersect <- data.frame("windows"=(intersect(ann_corr_topUnique$windows, ann_UncorrUnique$windows)))
ann_intersect$species<-rep("H.annuus")
ann_intersect$status<-rep("intersect")


arg_intersect <- data.frame("windows"=(intersect(arg_corr_topUnique$windows, arg_Uncorr_topUnique$windows)))
arg_intersect$species<-rep("H.argophyllus")
arg_intersect$status<-rep("intersect")


pf_intersect <- data.frame("windows"=(intersect(ppfal_corr_topUnique$windows, ppfal_Uncorr_topUnique$windows)))
pf_intersect$species<-rep("H.petiolaris.fallax")
pf_intersect$status<-rep("intersect")


pp_intersect <- data.frame("windows"=(intersect(ppet_corr_topUnique$windows, ppet_Uncorr_topUnique$windows)))
pp_intersect $species<-rep("H.petiolaris.petiolaris")
pp_intersect$status<-rep("intersect")

all_intersect_input<-do.call("rbind", list(ann_intersect,arg_intersect,pf_intersect,pp_intersect))

#### merge all and intersect combinations
all<-do.call("rbind", list(all_sep_input,all_intersect_input))
write.table(all, file = "all_windows_corrceted_uncorrcted_GWAS_with_intersection.table", col.names = TRUE, row.names= FALSE, quote= FALSE, sep = "\t")

###############################
###### plot the results #######

#### assemble  both (climate and phenotype) plots together ####
# Grouped Bar Plot
setwd("~/Documents/input_data/tests_analysis_for_paper/top_candidate_windows")
all_climate<-read.table(file = "all_windows_corrceted_uncorrcted_GEA_with_intersection.table", header = TRUE)

counts_cli <- table(all_climate$status, all_climate$species)
counts_cli <- counts_cli[ c(3,1,2),]
colnames(counts_cli)<-c("H.annuus","H.argophyllus","H.pet.fallax","H.pet.petiolaris")

###### Final plot from here (two plots assembled) #######
pdf (file = "barplots_top_candidate_windows_recomadjusted_climate_phenotype_with_Intersection.pdf",width = 10, height= 7)
par(mfrow=c(1,2))
par(mar = c(5.1,4.8,5.1,4.1))
barplot(counts_cli, 
        xlab="species", col=c("#1B9E77","#D95F02","#7570B3"),
        beside=TRUE, ylim=c(0,40000),ylab = "number of top candidate windows", cex.lab=1.5, cex.axis=1.5,width=c(0.03,0.03,0.03,0.03),cex = 0.65)
#legend("topright", inset=0.00009,
#       legend = c("both", "uncorrected","corrected"), 
#       fill = c("#1B9E77","#D95F02","#7570B3"),bty="n")

##
all_pheno<-read.table(file = "all_windows_corrceted_uncorrcted_GWAS_with_intersection.table", header = TRUE)

counts_pheno <- table(all_pheno$status, all_pheno$species)
counts_pheno <- counts_pheno[ c(3,1,2),]
colnames(counts_pheno)<-c("H.annuus","H.argophyllus","H.pet.fallax","H.pet.petiolaris")

barplot(counts_pheno, 
        xlab="species", col=c("#1B9E77","#D95F02","#7570B3"),
        beside=TRUE, ylim=c(0,40000),ylab = "number of top candidate windows", cex.lab=1.5, cex.axis=1.5,width=c(0.03,0.03,0.03,0.03),horiz = F,cex = 0.65)
#legend("topright", inset=c(-0.6,-0.01),
#       legend = c("both", "uncorrected","corrected"), cex = 0.75,
#       fill = c("#1B9E77","#D95F02","#7570B3"),bty="n")
dev.off()