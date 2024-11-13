########################################################
#   Go enrichment on groups of tissue-specific genes   #
########################################################

library(clusterProfiler)
library(dplyr)
library(readr)
library(tidyr)
library(reshape2) # to use decast

tau_scores <- read.table("tissue_specificity_median.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

Counts.tpm.flt <- read.table("counts_tpm_filtered.txt", sep="\t", stringsAsFactors=FALSE,header=TRUE)

go_anno <- read.table("GO_terms.csv", sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)
colnames(go_anno)[1] <- "GO_term"

gogene <- read.table('duo_GO_genes_LW_interpro.txt', sep="\t", stringsAsFactors=FALSE, header=TRUE)

go_molfunc <- go_anno[go_anno$Ontology=="molecular_function", ]
colnames(go_molfunc)[1] <- "GO_term"
go_biolproc <- go_anno[go_anno$Ontology=="biological_process", ]
colnames(go_biolproc)[1] <- "GO_term"
go_cellcomp <- go_anno[go_anno$Ontology=="cellular_component", ]
colnames(go_cellcomp)[1] <- "GO_term"

tissue_list <- unique(tau_scores$tissue)

#remove obsolete terms
go_anno <- go_anno[!grepl("OBSOLETE", go_anno$Definition),]
gogene <- semi_join(gogene, go_anno, by="GO_term") #from 54,501 associations to 53,976

# Creating vector on ontology terms and the tissues in which they appear
# (if an ontology term is absent in a tissue, the function will bug when this tissue is called)
gene_ontology <- left_join(gogene, go_anno, by="GO_term")
tissue_ont <- left_join(gene_ontology, tau_scores, by="Geneid")
tissue_ont <- tissue_ont[, -c(1, 2, 4:6, 8:9)]
tissue_ontology <- tissue_ont %>% drop_na(tissue) %>% unique() 
vec_molfunc <- tissue_ontology$tissue[tissue_ontology$Ontology=="molecular_function"]
vec_biolproc <- tissue_ontology$tissue[tissue_ontology$Ontology=="biological_process"]
vec_cellcomp <- tissue_ontology$tissue[tissue_ontology$Ontology=="cellular_component"]

#### For Molecular function ----

# extract association GO term-gene for GO terms corresponding to molecular functions
gogene_molfunc <- inner_join(go_molfunc, gogene, by="GO_term")
gogene_molfunc_red <- gogene_molfunc[, -c(2:5)]

go_enrichment <- function(tissue_name){
  tissue_vector <- tau_scores[tau_scores$tissue==tissue_name, ]
  
  if (tissue_name %in% vec_molfunc){
    result <- enricher(as.character(tissue_vector[,1]), 
                       universe=Counts.tpm.flt$geneid,
                       TERM2GENE=gogene_molfunc_red,
                       pvalueCutoff = 0.05)
    return(result)
  }
}

res_data <- lapply(tissue_list, go_enrichment)

res_table <- data.frame(tissue=character(), go_id=character(), gene_ratio=character(), bg_ratio=character(), p_adj=character(), 
                        stringsAsFactors=FALSE)

for (x in 1:length(tissue_list)) {
  tissue <- tissue_list[[x]]
  if (length(res_data[[x]]$ID) != 0 ){
    for (y in 1:length(res_data[[x]]$ID)){
      goid <- res_data[[x]]$ID[[y]]
      GeneRatio <- res_data[[x]]$GeneRatio[[y]]
      BgRatio <- res_data[[x]]$BgRatio[[y]]
      p.adj <- res_data[[x]]$p.adjust[[y]]
      res_table[nrow(res_table) + 1,] = c(tissue, goid, GeneRatio, BgRatio, p.adj)
    }  
  }
}

go_funct <- go_anno[, c(1, 3)]
colnames(go_funct)[1] <- "go_id"
final_table <- left_join(res_table, go_funct, by="go_id")
write_delim(final_table, file = "D:/GEA-data/Results_GO_enricher/GO_enrichment_tissue_median_molfunct.txt", delim = "\t")

#### For biological process ----

# extract association GO term-gene for GO terms corresponding to biological process
gogene_biolproc <- inner_join(go_biolproc, gogene, by="GO_term")
gogene_biolproc_red <- gogene_biolproc[, -c(2:5)]

go_enrichment_biolproc <- function(tissue_name){
  tissue_vector <- tau_scores[tau_scores$tissue==tissue_name, ]
  
  if (tissue_name %in% vec_biolproc){
    result <- enricher(as.character(tissue_vector[,1]), 
                       universe=Counts.tpm.flt$geneid,
                       TERM2GENE=gogene_biolproc_red,
                       pvalueCutoff = 0.05)
    return(result)
  }
}

res_data_biolproc <- lapply(tissue_list, go_enrichment_biolproc)

res_table_biolproc <- data.frame(tissue=character(), go_id=character(), gene_ratio=character(), bg_ratio=character(), 
                                 p_adj=character(), stringsAsFactors=FALSE)

for (x in 1:length(tissue_list)) {
  tissue <- tissue_list[[x]]
  if (length(res_data_biolproc[[x]]$ID) != 0 ){
    for (y in 1:length(res_data_biolproc[[x]]$ID)){
      goid <- res_data_biolproc[[x]]$ID[[y]]
      GeneRatio <- res_data_biolproc[[x]]$GeneRatio[[y]]
      BgRatio <- res_data_biolproc[[x]]$BgRatio[[y]]
      p.adj <- res_data_biolproc[[x]]$p.adjust[[y]]
      res_table_biolproc[nrow(res_table_biolproc) + 1,] = c(tissue, goid, GeneRatio, BgRatio, p.adj)
    }  
  }
}

go_funct <- go_anno[, c(1, 3)]
colnames(go_funct)[1] <- "go_id"
final_table_biolproc <- left_join(res_table_biolproc, go_funct, by="go_id")
write_delim(final_table_biolproc, file = "D:/GEA-data/Results_GO_enricher/GO_enrichment_tissue_median_biolproc.txt", delim = "\t")

#### For cellular component ----

# extract association GO term-gene for GO terms corresponding to biological process
gogene_cellcomp <- inner_join(go_cellcomp, gogene, by="GO_term")
gogene_cellcomp_red <- gogene_cellcomp[, -c(2:5)]

go_enrichment_cellcomp <- function(tissue_name){
  tissue_vector <- tau_scores[tau_scores$tissue==tissue_name, ]
  
  if (tissue_name %in% vec_cellcomp){
    result <- enricher(as.character(tissue_vector[,1]), 
                       universe=Counts.tpm.flt$geneid,
                       TERM2GENE=gogene_cellcomp_red,
                       pvalueCutoff = 0.05)
    return(result)
  }
}

res_data_cellcomp <- lapply(tissue_list, go_enrichment_cellcomp)

res_table_cellcomp <- data.frame(tissue=character(), go_id=character(), gene_ratio=character(), bg_ratio=character(), 
                                 p_adj=character(), stringsAsFactors=FALSE)

for (x in 1:length(tissue_list)) {
  tissue <- tissue_list[[x]]
  if (length(res_data_cellcomp[[x]]$ID) != 0 ){
    for (y in 1:length(res_data_cellcomp[[x]]$ID)){
      goid <- res_data_cellcomp[[x]]$ID[[y]]
      GeneRatio <- res_data_cellcompc[[x]]$GeneRatio[[y]]
      BgRatio <- res_data_cellcomp[[x]]$BgRatio[[y]]
      p.adj <- res_data_cellcomp[[x]]$p.adjust[[y]]
      res_table_cellcomp[nrow(res_table_cellcomp) + 1,] = c(tissue, goid, GeneRatio, BgRatio, p.adj)
    }  
  }
}

go_funct <- go_anno[, c(1, 3)]
colnames(go_funct)[1] <- "go_id"
final_table_cellcomp <- left_join(res_table_cellcomp, go_funct, by="go_id")
write_delim(final_table_cellcomp, file = "D:/GEA-data/Results_GO_enricher/GO_enrichment_tissue_median_cellcomp.txt", delim = "\t")
