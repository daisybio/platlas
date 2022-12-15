##### functional enrichment program
####panther
#@author: Leonora Raka
library(rbioapi)

test_dataset <-fread("www/all_diseased_TPM-0-1_disease-phenotype-all.CSV")

get_significant_genes <- function(dataset,pval_cutoff, fc_cutoff){ #momentan ensembl gene IDs
  #list of two vectors, upregulated and downregulated
  filtered_up <- dataset[dataset$log2FoldChange >= fc_cutoff & dataset$padj <= pval_cutoff,]
  filtered_down <- dataset[dataset$log2FoldChange <= -fc_cutoff & dataset$padj <= pval_cutoff,]
  filtered_up <- filtered_up[filtered_up$hgnc_symbol !="",]
  filtered_down <- filtered_down[filtered_down$hgnc_symbol !="",]
  up_genes <- filtered_up$hgnc_symbol
  down_genes <- filtered_down$hgnc_symbol
  resulting_list <- list("upregulated" = up_genes, "downregulated" = down_genes)
  return(resulting_list)
}


test_sig_list <- get_significant_genes(test_dataset,0.05, 1)

upregulated_genes <- test_sig_list[["upregulated"]]

annots <- rba_panther_info(what = "datasets")
organisms <- rba_panther_info(what = "organisms") #9606 = human

#annotation kann anders sein : molecular_function, biological_process, cellular_component, PANTHER GO Slim Molecular Function, PANTHER GO Slim Biological Process, PANTHER GO Slim Cellular Location
#genes <- c(p53,BRCA1,cdk2,Q99835,CDC42,CDK1,KIF23,PLK1,RAC2,RACGAP1,RHOA,RHOB,PHF14,RBM3,MSL1)

test_enriched_panther <- rba_panther_enrich(genes = upregulated_genes,
                               organism = 9606,
                               annot_dataset = "GO:0008150",
                               cutoff = 0.05, verbose = TRUE,test_type = "FISHER",
                               correction = "FDR")


test_mapping <- rba_panther_mapping(genes = c("Cd40", 7124, "ENSG00000203747", "P33681"),
                                    organism = 9606)




enrichr_libs <- rba_enrichr_libs()
enrichr_enrich <- rba_enrichr(gene_list = upregulated_genes,
                              gene_set_library = "GO_Molecular_Function_2015")



############intermission making a sample_annotation file 
library(stringr)
DBO <- "www/DBO-samples.tsv"
dbo <- read.table(DBO, header = FALSE, sep = "\t", quote = "")
colnames(dbo) <- c("sample","condition","platelet")

new_sample_annotation_fp <- "/Users/leonoraka/Desktop/Projekte/final_plateletAnalysis_CCSonly/data/sample_annotation/sample_names.tsv"
new_sample_annotation <- fread(new_sample_annotation_fp)

dbo$sample <- str_replace(dbo$sample,"[.]","-")

merged_annotation <- left_join(dbo,new_sample_annotation, by = "sample")
merged_annotation$group[merged_annotation$sample == "DBO-03"] <- "healthy_P1"
merged_annotation$group[merged_annotation$sample == "DBO-04"] <- "healthy_P2"
merged_annotation$group[merged_annotation$sample == "DBO-05"] <- "healthy_P3"
merged_annotation$group[merged_annotation$sample == "DBO-06"] <- "healthy_P4"
merged_annotation$group[merged_annotation$sample == "DBO-07"] <- "healthy_P5"
merged_annotation$group[merged_annotation$sample == "DBO-08"] <- "healthy_P6"
merged_annotation$group[merged_annotation$sample == "DBO-09"] <- "healthy_P7"
merged_annotation$group[merged_annotation$sample == "DBO-10"] <- "healthy_P8"


merged_annotation$original_sample_name <- merged_annotation$sample
merged_annotation$original_sample_name <- str_replace(merged_annotation$original_sample_name,"-",".")

#path 
other_new_annotation <- fread("/Users/leonoraka/Desktop/Projekte/final_plateletAnalysis_CCSonly/data/sample_annotation/processed_metadata_ccs_only_noSample6_filtered.csv")
other_new_annotation$sample_name_huge <- str_replace(other_new_annotation$sample_name_huge,"[_]","-")
colnames(other_new_annotation)[4] <- "sample"
new_merged_annot <- full_join(merged_annotation,other_new_annotation,by = "sample")

#fill NAs
new_merged_annot$RPs_MPs[is.na(new_merged_annot$RPs_MPs)] <- paste0(new_merged_annot$platelet,"s")
new_merged_annot$platelet[is.na(new_merged_annot$platelet)] <- str_replace(new_merged_annot$RPs_MPs,"s","")

new_merged_annot$condition[is.na(new_merged_annot$condition)|new_merged_annot$condition == "stable CAD"] <- new_merged_annot$disease
new_merged_annot$condition[new_merged_annot$sample %in% c("DBO-33","DBO-34","DBO-35", "DBO-36","DBO-43","DBO-44","DBO-47", "DBO-48","DBO-53", "DBO-54","DBO-63", "DBO-64", "DBO-111", "DBO-112","DBO-115", "DBO-116")] <- new_merged_annot$disease[new_merged_annot$sample %in% c("DBO-33","DBO-34","DBO-35", "DBO-36","DBO-43","DBO-44","DBO-47", "DBO-48","DBO-53", "DBO-54","DBO-63", "DBO-64", "DBO-111", "DBO-112","DBO-115", "DBO-116")]
new_merged_annot$disease[is.na(new_merged_annot$disease)& new_merged_annot$condition == "MI"] <- "ACS"
new_merged_annot$condition[new_merged_annot$condition == "MI"] <- "ACS"
new_merged_annot$disease[is.na(new_merged_annot$disease)] <- "healthy"
write.table(new_merged_annot,file = "www/incomplete_sample_names_all.tsv",sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


##simulate new sample names
simulated_annot <- copy(new_merged_annot)
simulated_annot$sample_name[is.na(simulated_annot$sample_name)] <- simulated_annot$RPs_MPs

simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-03","DBO-04")] <- paste0("21_",simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-03","DBO-04")])
simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-05","DBO-06")] <- paste0("22_",simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-05","DBO-06")])
simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-07","DBO-08")] <- paste0("23_",simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-07","DBO-08")])
simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-09","DBO-10")] <- paste0("24_",simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-09","DBO-10")])

simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-51","DBO-52")] <- paste0("25_",simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-51","DBO-52")])
simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-55","DBO-56")] <- paste0("26_",simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-55","DBO-56")])
simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-61","DBO-62")] <- paste0("27_",simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-61","DBO-62")])
simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-73","DBO-74")] <- paste0("28_",simulated_annot$sample_name[simulated_annot$sample %in% c("DBO-73","DBO-74")])


write.table(simulated_annot,file = "www/simulated_sample_names_all.tsv",sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#E-Mail an Kilian

write.table(merged_annotation,file = "www/complete_sample_names.tsv",sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)





test_merged_annotation_file <- fread("www/complete_sample_names.tsv")


fill_in_DBOs <- fread("/Users/leonoraka/Desktop/Projekte/final_plateletAnalysis_CCSonly/data/sample_annotation/sample_names.tsv")

list_of_DBOs <- simulated_annot$sample[simulated_annot$sample %in% fill_in_DBOs$sample]
fill_in_DBOs <- as.data.frame(fill_in_DBOs)
rownames(fill_in_DBOs) <- fill_in_DBOs$sample

simulated_annot$group[simulated_annot$sample %in% fill_in_DBOs$sample] <-  fill_in_DBOs[list_of_DBOs, 1]
simulated_annot$original_sample_name[is.na(simulated_annot$original_sample_name)] <- str_replace(simulated_annot$sample[is.na(simulated_annot$original_sample_name)],"-",".")
simulated_annot$condition <- simulated_annot$disease

simulated_annot[55,] <- c("DBO-65","CCS","MP","CCS_R28","DBO.65","RPs039","MPs", "CCS", NA, "m", "ASS, Ticagrelor","40","6_MPs")
simulated_annot[56,] <- c("DBO-66","CCS","RP","CCS_R29","DBO.66","RPs039","RPs", "CCS", NA, "m", "ASS, Ticagrelor","40","6_RPs")#"DBO-66"




#test the files for CCS
difex_file_path <- "/Users/leonoraka/Desktop/Projekte/final_plateletAnalysis_CCSonly/gene_analysis/results/ccs_rp_vs_mp_diff_gene_counts_noSample6.xlsx"
heat_CSS_path <- "/Users/leonoraka/Desktop/Projekte/final_plateletAnalysis_CCSonly/data/gene_counts/salmon_merged_gene_counts.tsv"
heat_file_CCS <- fread(heat_CSS_path)
library(openxlsx)
difex_file_CCS <- read.xlsx(difex_file_path) #eig muss read.csv2 werden -> rewrite 


#anpassung

difex_file_CCS$log10Pval <- NA
difex_file_CCS$effect_size <- NA
difex_file_CCS$significant <- NA

colnames(difex_file_CCS)[1] <- "geneID"

difex_correct <- difex_file_CCS[,c(1,3,4,5,6,7,8,10,11,12,2,9)]


write.csv2(difex_correct,file = "/Users/leonoraka/Desktop/Projekte/platlas/www/all_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV", quote = FALSE, row.names = FALSE, col.names = TRUE)
#write.csv2()

test_difex_file <- read.csv2("/Users/leonoraka/Desktop/Projekte/platlas/www/all_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV")


#fix colnames
heat_colnames <- colnames(heat_file_CCS)[2:ncol(heat_file_CCS)]
rownames(fill_in_DBOs) <-fill_in_DBOs$group 
translate <- fill_in_DBOs[heat_colnames,2]
colnames(heat_file_CCS)[2:ncol(heat_file_CCS)] <- translate

#
difex_to_map <- difex_correct[, c(1,11,12)]
rownames(difex_to_map) <- difex_to_map$geneID
heat_file_CCS$hgnc_symbol <- NA
heat_file_CCS$hgnc_symbol <-difex_to_map[heat_file_CCS$gene_id,2]
heat_file_CCS$gene_biotype <- difex_to_map[heat_file_CCS$gene_id,3]

heat_file_end <- heat_file_CCS[,c(1,50,51,2:49)]
colnames(heat_file_end)[1] <- "ensembl_gene_id"
write.csv2(heat_file_end,file = "/Users/leonoraka/Desktop/Projekte/platlas/www/heat_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV", quote = FALSE, row.names = FALSE, col.names = TRUE)


test_heat_file <- read.csv2("/Users/leonoraka/Desktop/Projekte/platlas/www/heat_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV")

















