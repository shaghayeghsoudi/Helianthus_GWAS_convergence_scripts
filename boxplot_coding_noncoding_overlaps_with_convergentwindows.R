

### this script plot (boxplot) 
setwd("~/Documents/input_data/tests_analysis_for_paper/coding_noncoding/out_put_results")

cds_reg<-read.table(file = "test_cds_reg_script.txt", header = TRUE)
cds_reg<-cds_reg[cds_reg$status=="convergent",]

cds_annarg<-cds_reg[cds_reg$type=="Annuus_Argophyllus" & cds_reg$genomic_pos=="CDS",]
reg_annarg<-cds_reg[cds_reg$type=="Annuus_Argophyllus" & cds_reg$genomic_pos=="regulatory",]


cds_annpf<-cds_reg[cds_reg$type=="Annuus_petfal"& cds_reg$genomic_pos=="CDS",]
reg_annpf<-cds_reg[cds_reg$type=="Annuus_petfal"& cds_reg$genomic_pos=="regulatory",]


cds_annpp<-cds_reg[cds_reg$type=="Annuus_petpet"& cds_reg$genomic_pos=="CDS",]
reg_annpp<-cds_reg[cds_reg$type=="Annuus_petpet"& cds_reg$genomic_pos=="regulatory",]


cds_argpf<-cds_reg[cds_reg$type=="Argophyllus_petfal"& cds_reg$genomic_pos=="CDS",]
reg_argpf<-cds_reg[cds_reg$type=="Argophyllus_petfal"& cds_reg$genomic_pos=="regulatory",]


cds_argpp<-cds_reg[cds_reg$type=="Argophyllus_petpet"& cds_reg$genomic_pos=="CDS",]
reg_argpp<-cds_reg[cds_reg$type=="Argophyllus_petpet"& cds_reg$genomic_pos=="regulatory",]


cds_pfpp<-cds_reg[cds_reg$type=="petfal_petpet"& cds_reg$genomic_pos=="CDS",]
reg_pfpp<-cds_reg[cds_reg$type=="petfal_petpet"& cds_reg$genomic_pos=="regulatory",]



cds_reg$type=factor(cds_reg$type,levels=c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","Argophyllus_petpet","petfal_petpet"))
cds_reg$genomic_pos=factor(cds_reg$genomic_pos,levels=c("CDS","regulatory"))
    
    
### Making the plot ####

lvl1=c("CDS","regulatory")
lvl2=c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","Argophyllus_petpet","petfal_petpet")
    factor1=as.factor(c(rep("CDS",nrow(cds_annarg)),rep("regulatory",nrow(reg_annarg)),rep("CDS",nrow(cds_annpf)),rep("regulatory",nrow(reg_annpf)),rep("CDS",nrow(cds_annpp)),rep("regulatory",nrow(reg_annpp)),rep("CDS",nrow(cds_argpf)),rep("regulatory",nrow(reg_argpf)), rep("CDS",nrow(cds_argpp)),rep("regulatory",nrow(reg_argpp)),rep("CDS",nrow(cds_pfpp)),rep("regulatory",nrow(reg_pfpp))))
    
    factor2=as.factor(c(rep("Annuus_Argophyllus",nrow(reg_annarg)),rep("Annuus_petfal",nrow(reg_annpf)),rep("Annuus_petpet",nrow(reg_annpp)),rep("Argophyllus_petfal",nrow(reg_argpf)),rep("Argophyllus_petpet",nrow(reg_argpp)),rep("petfal_petpet",nrow(reg_pfpp))))
plotgrp=factor(paste(factor2,factor1),levels=c(sapply(lvl2,paste,lvl1)))
    
pi=c(cds_annarg$proporlength,   reg_annarg$proporlength,   cds_annpf$proporlength,   reg_annpf$proporlength,   cds_annpp$proporlength,  reg_annpp$proporlength, 
cds_argpf$proporlength, reg_argpf$proporlength, cds_argpp$proporlength, reg_argpp$proporlength,cds_pfpp$proporlength,reg_pfpp$proporlength )
    
    
 par(las=1)
    par(mar=c(3,3,1,1))
    #cols=c("firebrick1","dodgerblue2","dodgerblue4","darkblue", "brown","brown4")

    at_1=c(1:2,4:5,7:8,10:11,13:14,16:17)
    at=c(1.5,4.5,7.5,10.5,13.5,16.5)
    boxplot(pi~plotgrp,notch=FALSE,at=at_1,xaxt="n",cex.axis=1.5,frame.plot =        FALSE,xlab="",ylab="",cex.lab=2.5,labs=2,col=c("firebrick","dodgerblue4","firebrick","dodgerblue4","firebrick","dodgerblue4","firebrick","dodgerblue4","firebrick","dodgerblue4","firebrick","dodgerblue4"),outline=FALSE,ylim=c(0,0.8))
    axis(1,at=at,labels=lvl2,tick=T,cex.axis=1)
    mtext(expression(Fst),las=1,side=2,line=6,cex=2)   
    
    
    
    
### test significance
#A.aegy1
median(cds_annarg$proporlength,na.rm=T)
mean(cds_annarg$proporlength,na.rm=T)
#quantile(nonregion_G1$pi,na.rm=T,c(0.025,0.975))
mean(reg_annarg$proporlength,na.rm=T)
#quantile(region_G1$pi,na.rm=T,c(0.025,0.975))    
mean(cds_annpf$proporlength,na.rm=T)
mean(reg_annpf$proporlength,na.rm=T)

mean(cds_annpp$proporlength,na.rm=T)
mean(reg_annpp$proporlength,na.rm=T)

mean(cds_argpf$proporlength,na.rm=T)
mean(reg_argpf$proporlength,na.rm=T)

mean(cds_argpp$proporlength,na.rm=T)
mean(reg_argpp$proporlength,na.rm=T)


mean(cds_pfpp$proporlength,na.rm=T)
mean(reg_pfpp$proporlength,na.rm=T)


wilcox.test(cds_annarg$proporlength,reg_annarg$proporlength)   
wilcox.test(cds_annpf$proporlength,reg_annpf$proporlength)    
wilcox.test(cds_annpp$proporlength,reg_annpp$proporlength)
wilcox.test(cds_argpf$proporlength,reg_argpf$proporlength)
wilcox.test(cds_argpp$proporlength,reg_argpp$proporlength)
wilcox.test(cds_pfpp$proporlength,reg_pfpp$proporlength)




    
        