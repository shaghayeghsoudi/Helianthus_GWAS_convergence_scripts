rm(list = ls())
#library("GenomicRanges")
library("stringr")
library("biomaRt")
library("topGO")
library("enrichplot")
library("ggplot2")
#library("ggpubr")
#library("dplyr")

setwd("~/Dropbox/input_data/GO_annotaion_analysis/go_enrichment_picmin")

## load picmin significant windows
picmin_unclust<-read.delim(file = "picmin_fdr20_results_unclustered_4way_sh.txt", header = FALSE)
picmin_unclust$V9<-gsub("llevation" ,"Elevation",picmin_unclust$V9 )

## find list of significant variables in the picmin table
significat_vars<-unique(picmin_unclust$V9)
sig_soil<-c("K","PERCENT_K","MG","PH","PERCENT_NA","Sodium","PERCENT_MG","P1","CA","OM","BICARB","P2","SOL_SALTS","CEC")
sig_pheno<-c("RGB_proportion_blue","Leaf_curvedHeight_maxWidth","Leaf_maximum_height","Leaf_total_N","Leaf_proximal_eccentricity","DTF","Leaf_total_C","SLA")
sig_temp<-c("EXT","NFFD","FFP","eFFP","bFFP","EMT","MWMT","DD_18","latitude","DD5","MAT","longitude","DD_0","DD18")
sig_precip<-c("MSP","AHM","CMD","Eref","RH","MAP","SHM","Elevation","MAR")


## assign four variable type to the picmin table
picmin_unclust$var_type<-ifelse (picmin_unclust$V9%in%sig_soil, "soil",
                         ifelse (picmin_unclust$V9%in%sig_pheno, "phenotype",
                         ifelse (picmin_unclust$V9%in%sig_temp, "temprature",
                         ifelse (picmin_unclust$V9%in%sig_precip, "precipitation","-")
                         )))


colnames(picmin_unclust)<-c("w1","w2","w3","v4","chrom","start","end","window_id","variable","q","inversion_overlap","recomb_rate","recomb_quantile","over","unique_id","var_type")
picmin_unclust_good<-picmin_unclust[,c("chrom","start","end","window_id","var_type")]
write.table(picmin_unclust_good, file = "picmin_fdr20_results_unclustered_4way_with_variable_type_id.bed", col.names = FALSE, sep = "\t", quote= FALSE, row.names = FALSE)


###########################################
### prepare gff file ###
gff<-read.table(file = "HAN412_Eugene_curated_v1_1.just.GENES.gff3", header= FALSE)[,c(1:5,9)]
gff_good<-gff[!grepl("Ha412HOChr00",gff$V1),]
gff_good$cities <- gsub(";.*$", "", gff_good$V9)
gff_good$ensemble_gene_id <- str_split_fixed(gff_good$cities, "=", 2)[,2]

gff_good<-gff_good[,c("V1","V4","V5","ensemble_gene_id")]
colnames(gff_good)<-c("chrom","start","end","ensemble_gene_id")

gff_good$start<-(gff_good$start)-500
gff_good$start<-abs(gff_good$start)
gff_good$end<-(gff_good$end)+500

write.table(gff_good, file = "HAN412_Eugene_curated_v1_1_GENES.bed", col.names = FALSE, sep = "\t", quote= FALSE, row.names = FALSE)

###########################################
### perform GO-enrichment analysis ###
windows_picmin<-read.table(file = "picmin_fdr20_results_unclustered_4way_overlapped_with_HAN412_Eugene_with_variable_id.table", header = FALSE)

var_types<-unique(windows_picmin$V5)

out_res<-NULL
for (i in 1:length(var_types)){
    
    windows_picmin_focal<-windows_picmin[windows_picmin$V5%in%var_types[i],]
    windows_picmin_focal_unique<-windows_picmin_focal[!duplicated(windows_picmin_focal$V9),]   ## V9 is gene ID
    out_res<-rbind(windows_picmin_focal_unique,out_res)

}


joes<-read.delim(file = "FINAL_PANGENOME_df_V4.table", header = TRUE)

overlaps_with_xrq<-merge(out_res,joes, by.x = "V9",by.y ="HAN412")
overlaps_with_xrq$m1<-str_split_fixed(overlaps_with_xrq$XRQv2_match, ":", 2)[,1]
overlaps_with_xrq$m2<-sub(".*HanXRQ","",overlaps_with_xrq$m1)
overlaps_with_xrq$m3<-rep("HanXRQr2")

overlaps_with_xrq$ensembl_gene_id<-paste(overlaps_with_xrq$m3,overlaps_with_xrq$m2, sep = "_")

write.table(overlaps_with_xrq, file= "CONV_WIN_GENE_picmin_overlap_HanXRQv2_Ha412HOv2_with_variable_id.table", col.names= TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

##### run Go-enrichment analysis ##
## run Go-enrichment analysis
mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "hannuus_eg_gene",
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
overlaps_with_xrq<-read.delim(file = "CONV_WIN_GENE_picmin_overlap_HanXRQv2_Ha412HOv2_with_variable_id.table", header = TRUE)
var_types<-unique(overlaps_with_xrq$V5)


for (i in 1:length(var_types)){
overlaps_with_xrq_focal<-overlaps_with_xrq[overlaps_with_xrq$V5%in%var_types[i],]
int.genes<-unique(overlaps_with_xrq_focal$ensembl_gene_id)
int.genes <- factor(as.integer(all.genes %in% int.genes))

go_functions<-c("BP","MF","CC")
  for (j in 1:length(go_functions)){

   names(int.genes) = all.genes
   #go.obj <- new("topGOdata", ontology='MF'      ## BP,MF, CC
   go.obj <- new("topGOdata", ontology= go_functions[j]     ## BP,MF, CC
             
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
   p.adj= round(p.adjust(results.tab$elimFisher, method= "BH"), digits = 4)
   results.tab$p.adj<-p.adj
   results.tab$go_function<-go_functions[j]
   results.tab$var_type<-var_types[i]
   #-log(results.tab$p.adj)
   #results.tab$comparison<-rep(comparisons[i])

  #CC<-results.tab[order(results.tab$p.adj),]
  #results.tab$Term <- factor(CC$Term,                                    # Factor 
  #                        levels = CC$Term[order(CC$Significant, decreasing = TRUE)])
  write.table(results.tab, file =paste("out_res_GOEnrichment_TopGo_convergent_windowos_",var_types[i],"picmin_ElimFisher_",go_functions[j],".table", sep = ""),  quote = FALSE, sep = "\t", row.names = FALSE, col.names =TRUE )

   }
}


#### plot the results for each variable type (i.e., soil, temperature, precipitation, phenotype)
### load all tables 
all_tables<-list.files("out_res_GOEnrichment_TopGo_var_types", pattern = "out_res_GOEnrichment_TopGo_convergent_windowos_*",full.names = TRUE)

attackStats <- lapply(all_tables,function(x) {
    read.delim(x,  header=TRUE, sep = "\t")
    })

aa <- do.call("rbind", attackStats) 


var_types<-unique(aa$var_type)
go_functions<-unique(aa$go_function)


aa$Term<-gsub("histone H3-K4 demethylation, trimethyl-H...","histone H3-K4 demethylation, trimethyl-Histone",aa$Term)
aa$Term<-gsub("inorganic ion import across plasma membr...","inorganic ion import across plasma membrane",aa$Term)
aa$Term<-gsub("histone demethylase activity (H3-trimeth...","histone demethylase activity H3-trimethylase activity",aa$Term)
aa$Term<-gsub("ammonium transmembrane transporter activ...","ammonium transmembrane transporter activity",aa$Term)
aa$Term<-gsub("cyclic threonylcarbamoyladenosine biosyn...","cyclic threonylcarbamoyladenosine biosynthetic process",aa$Term)
aa$Term<-gsub("DNA synthesis involved in UV-damage exci...","DNA synthesis involved in UV-damage excision repair",aa$Term)
aa$Term<-gsub("nucleotide-excision repair, DNA gap fill...","nucleotide-excision repair, DNA gap filling",aa$Term)
aa$Term<-gsub("tetrahydrofolylpolyglutamate metabolic p...","tetrahydrofolylpolyglutamate metabolic process",aa$Term)
aa$Term<-gsub("tRNA threonylcarbamoyladenosine dehydrat...","tRNA threonylcarbamoyladenosine dehydratase",aa$Term)
aa$Term<-gsub("purine ribonucleoside triphosphate bindi...","purine ribonucleoside triphosphate binding",aa$Term)
aa$Term<-gsub("copper transmembrane transporter activit...","copper transmembrane transporter activity",aa$Term)
aa$Term<-gsub("deadenylation-independent decapping of n...","deadenylation-independent decapping of nuclear-transcribed mRNA",aa$Term)
aa$Term<-gsub("NAD+ nucleotidase, cyclic ADP-ribose gen...","NAD+ nucleotidase, cyclic ADP-ribose generating",aa$Term)
aa$Term<-gsub("DNA synthesis involved in UV-damage exci...","DNA synthesis involved in UV-damage excision repair",aa$Term)
aa$Term<-gsub("brassinosteroid mediated signaling pathw...","brassinosteroid mediated signaling pathway",aa$Term)
aa$Term<-gsub("putrescine biosynthetic process from arg...","putrescine biosynthetic process from arginine, using agmatinase",aa$Term)
aa$Term<-gsub("nucleotide-excision repair, DNA gap fill...","nucleotide-excision repair, DNA gap filling",aa$Term)
aa$Term<-gsub("tetrahydrofolylpolyglutamate metabolic p...","tetrahydrofolylpolyglutamate metabolic process",aa$Term)
aa$Term<-gsub("nuclear-transcribed mRNA catabolic proce...","nuclear-transcribed mRNA catabolic process, deadenylation-independent decay",aa$Term)
aa$Term<-gsub("deadenylation-independent decapping of n...","deadenylation-independent decapping of nuclear-transcribed mRNA",aa$Term)
aa$Term<-gsub("citrate transmembrane transporter activi...","citrate transmembrane transporter activity",aa$Term)
aa$Term<-gsub("DNA (cytosine-5-)-methyltransferase acti...","DNA (cytosine-5-)-methyltransferase activity",aa$Term)
aa$Term<-gsub("translation release factor activity, cod...","translation release factor activity, codon specific",aa$Term)



for ( in in 1:length(var_types)){

  aa_focal<-aa[aa$var_type==var_types[i],]
  
  
  our_res2<-NULL
  for (j in 1:length(go_functions)){

    aa_focal_term<-aa_focal[aa_focal$go_function==go_functions[j],]
    
    
    aa_focal_term<-aa_focal_term[order(aa_focal_term$p.adj),]
    aa_focal_term$Term <- factor(aa_focal_var$Term,                                    # Factor 
                  levels = aa_focal_var$Term[order(aa_focal_var$Significant, decreasing = TRUE)])
   
   
    our_res2<-rbind(aa_focal_term,our_res2)
   
   
   
   gg<-ggplot(our_res2, aes (x = Term, y = Significant, fill = p.adj))+
   geom_bar(stat="identity")+ scale_fill_gradient(low="blue", high="red") + 
   facet_grid(go_function ~ ., scales = "free_y",space = "free_y")  +
   coord_flip() + theme(panel.background = element_blank(),axis.text=element_text(size=12))+
   labs(y = "Number of significant genes", x = "GO terms") 

   pdf(file = paste("Go_enrichment_picmin_4way_unclustered_ElimFisher_per_variable_type_",var_types[i],".pdf",sep = "" ),height= 12, width = 12)
   print(gg)
   dev.off()

  }

}