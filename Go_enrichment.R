
library("biomaRt")
library("topGO")
library("enrichplot")
library("ggplot2")
library(ggpubr)

## http://yulab-smu.top/clusterProfiler-book/chapter12.html

setwd("~/Documents/input_data/tests_analysis_for_paper/go_enrichment")

windows<-read.table(file = "CONV_WIN_GENE_overlap_Ath_HanXRQr1-2_Ha412HOv2_NA", header = TRUE)
windows$arabidopsis<-sub("\\..*", "", windows$arabidopsis)


comparisons<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","petfal_petpet","Argophyllus_petpet")
analysis<-c("spearman","GWAS_corrected")

#windows_spearman<-windows[(windows$comparison == "Annuus_Argophyllus") & (windows$analysis=="spearman"), ]
#windows_spearman<-na.omit(windows_spearman)


#collect gene names from biomart
mart <- biomaRt::useMart(biomart = "plants_mart",
                        dataset = "athaliana_eg_gene",
                         host = 'plants.ensembl.org')
                         

#listAttributes(mart) ## ensembl_gene_id and go_id selected from getBM                        
GTOGO <- getBM(attributes = c( "ensembl_gene_id",
                                        "go_id"), mart = mart)
                                        
                                        
#examine result
head (GTOGO)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))
#examine result
head (geneID2GO)


#Make topGO data object
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
#int.genes <- sample(x = all.genes, size = 200) # some random genes 

#### int.genes are focal convergent genes ####

out_res<-NULL
for (i in 1:length(comparisons)){
	
	for (j in 1:length(analysis)){
		
		windows_focal<-windows[(windows$comparison == comparisons[i]) & (windows$analysis==analysis[j]), ]
        windows_focal<-windows_focal[!is.na(windows_focal$arabidopsis),]
        int.genes<-unique(windows_focal$arabidopsis)
        int.genes <- factor(as.integer(all.genes %in% int.genes))
        names(int.genes) = all.genes
        go.obj <- new("topGOdata", ontology='CC'      ## BP,MF, CC
              , allGenes = int.genes
              , annot = annFUN.gene2GO
              , gene2GO = geneID2GO              
              )

        results_elim <- runTest(go.obj, algorithm = "elim", statistic = "fisher")
        #results_classic <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
        
        results.tab <- GenTable(object = go.obj, elimFisher = results_elim)
        #results.tab<-GenTable(object = go.obj, classicFisher = results_classic)
        # Plot results using the GO hierarchical DAG:
        #showSigOfNodes(go.obj, score(results), firstSigNode=15, useInfo ='all')

        ### perform BH correction for P-values
        p.adj= round(p.adjust(results.tab$elimFisher, method= "BH"), digits = 4)
        results.tab$p.adj<-p.adj
        results.tab$comparison<-rep(comparisons[i])
        results.tab$analysis<-rep(analysis[j])
        out_res<-rbind (out_res,results.tab)


	}
}

write.csv(out_res, file = "out_res_GOEnrichment_TopGo_arabidopsis_homologs_convergent_windowos_ElimFisher_CC.csv",  quote = FALSE)

    
##### plot the results in a facet format

comparisons<-c("Annuus_Argophyllus","Annuus_petfal","Annuus_petpet","Argophyllus_petfal","petfal_petpet","Argophyllus_petpet")
analysis<-c("spearman","GWAS_corrected")


BP<-read.csv(file = "out_res_GOEnrichment_TopGo_arabidopsis_homologs_convergent_windowos_ElimFisher_BP.csv")
BP$gofunction<-rep("BP")



MF<-read.csv(file = "out_res_GOEnrichment_TopGo_arabidopsis_homologs_convergent_windowos_ElimFisher_MF.csv")
MF$gofunction<-rep("MF")


CC<-read.csv(file = "out_res_GOEnrichment_TopGo_arabidopsis_homologs_convergent_windowos_ElimFisher_CC.csv")
CC$gofunction<-rep("CC")

#for (j in 1:length(analysis)){
     
#     for (i in 1:length(comparisons)){
     	
    
	
BP<-BP[(BP$comparison==comparisons[i]) & (BP$analysis==analysis[j]), ]
BP<-BP[order(BP$p.adj),]
BP$Term <- factor(BP$Term,                                    # Factor 
                   levels = BP$Term[order(BP$Significant, decreasing = TRUE)])

 
MF<-MF[(MF$comparison==comparisons[i]) & (MF$analysis==analysis[j]), ]
MF<-MF[order(MF$p.adj),]
MF$Term <- factor(MF$Term,                                    # Factor 
                   levels = MF$Term[order(MF$Significant, decreasing = TRUE)])
 
 
 
CC<-CC[(CC$comparison==comparisons[i]) & (CC$analysis==analysis[j]), ]
CC<-CC[order(CC$p.adj),]
CC$Term <- factor(CC$Term,                                    # Factor 
                   levels = CC$Term[order(CC$Significant, decreasing = TRUE)])
 
 
all_functions<-do.call("rbind", list(BP,MF,CC))  
 
 
 gg<-ggplot(all_functions, aes (x = Term, y = Significant, fill = p.adj))+
     geom_bar(stat="identity")+
     facet_grid(gofunction ~ ., scales = "free_y",space = "free_y")  +
    coord_flip() + theme(panel.background = element_blank(),axis.text=element_text(size=8))+
     labs(y = "Go terms", x = "Number of genes")     
 
 
 pdf(file = paste("Go_enrichment_",comparisons[i],"_bars_",analysis[j],".pdf",sep = "" ),height= 5.7, width = 6)
 print(gg)
 dev.off()
 
#     }
#}
    
                 
