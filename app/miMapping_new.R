#miMapping neu + weniger files
library(openxlsx)


tar_db <- read.xlsx("/Users/leonoraka/Downloads/hsa_MTI.xlsx")
rdb_db <- fread("/Users/leonoraka/Downloads/miRDB_v6.0_prediction_result.txt")
rdb_db_human <- rdb_db[grepl("hsa",rdb_db$V1),]

tar_db <- tar_db[,c(2,1,3,4,5,6,7,8,9)]
tar_db$exp_type <- NA
tar_db$exp_type[!grepl("(Weak)",tar_db$Support.Type)& !is.na(tar_db$Support.Type)] <- "experimentally_retrieved"
tar_db$exp_type[grepl("(Weak)",tar_db$Support.Type) | is.na(tar_db$Support.Type)] <- "not experimentally retrieved"

##miRNAs aus Datensätzen filtern:
import_css_miRNA_de <- read.xlsx("/Users/leonoraka/Desktop/Projekte/final_plateletAnalysis_CCSonly/miRNA_analysis/results/ccs_rp_vs_mp_diff_mir_counts_noSample6_filtered_pCutoff0.05.xlsx")
import_dis_miRNA_de <- read.csv2("/Users/leonoraka/Desktop/Projekte/platlas/www/disease_diff_mir_counts.CSV")
import_rp_mp_miRNA_de <- read.csv2("/Users/leonoraka/Desktop/Projekte/platlas/www/rp_vs_mp_diff_mir_counts.CSV")

css_tar_db <- tar_db[tar_db$miRNA %in% import_css_miRNA_de$X1,]
css_rdb_db <- rdb_db_human[rdb_db_human$V1 %in% import_css_miRNA_de$X1,]

dis_tar_db <- tar_db[tar_db$miRNA %in% import_dis_miRNA_de$X,]
dis_rdb_db <- rdb_db_human[rdb_db_human$V1 %in% import_dis_miRNA_de$X,]

rp_mp_tar_db <- tar_db[tar_db$miRNA %in% import_rp_mp_miRNA_de$X,]
rp_mp_rdb_db <- rdb_db_human[rdb_db_human$V1 %in% import_rp_mp_miRNA_de$X,]








#rdb sachen translaten
css_rdb_db <- rdb_translate(css_rdb_db)
css_tar_db <- tar_translate(css_tar_db)
fp <- "/Users/leonoraka/Desktop/Projekte/platlas/www/miRNA/new/"
write.table(css_rdb_db, paste0(fp,"ccs_rdb_db.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(css_tar_db, paste0(fp,"ccs_tar_db.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
dis_rdb_db <- rdb_translate(dis_rdb_db)
dis_tar_db <- tar_translate(dis_tar_db)
write.table(dis_rdb_db, paste0(fp,"dis_rdb_db.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(dis_tar_db, paste0(fp,"dis_tar_db.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
rp_mp_rdb_db <- rdb_translate(rp_mp_rdb_db)
rp_mp_tar_db <- tar_translate(rp_mp_tar_db)

write.table(rp_mp_rdb_db, paste0(fp,"rp_mp_rdb_db.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(rp_mp_tar_db, paste0(fp,"rp_mp_tar_db.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")


#ccs file in .CSV umwandeln

write.csv2(import_css_miRNA_de, "/Users/leonoraka/Desktop/Projekte/platlas/www/ccs_rp_vs_mp_diff_mir_counts_noSample6_filtered_pCutoff0.05.CSV", row.names = FALSE)
test_css <- read.csv2("/Users/leonoraka/Desktop/Projekte/platlas/www/ccs_rp_vs_mp_diff_mir_counts_noSample6_filtered_pCutoff0.05.CSV")




#functions:
rdb_translate <- function(filtered_tab){
  #translate rdb gene format
  attributes <- c("refseq_mrna", "ensembl_gene_id","hgnc_symbol")
  filters <- "refseq_mrna"
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  #attributes <- listAttributes(mart)
  #filters <- listFilters(mart)
  genes <- getBM(attributes = attributes,
                 filters = filters,
                 values = filtered_tab$V2,
                 mart = mart)
  
  
  
  colnames(filtered_tab)[2] = "refseq_mrna"
  complete_filtered_rdb <- left_join(filtered_tab,genes)
  complete_filtered_rdb <- unique(complete_filtered_rdb)
  complete_filtered_rdb <- complete_filtered_rdb[,c(2,5,1,3,4)]
  return(complete_filtered_rdb)
}
tar_translate <- function(filtered_tab){
  #translate rdb gene format
  attributes <- c("hgnc_symbol", "ensembl_gene_id")
  filters <- "hgnc_symbol"
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  #attributes <- listAttributes(mart)
  #filters <- listFilters(mart)
  genes <- getBM(attributes = attributes,
                 filters = filters,
                 values = filtered_tab$Target.Gene,
                 mart = mart)
  
  
  #View(genes)
  colnames(filtered_tab)[4] = "hgnc_symbol"
  #View(filtered_tab)
  complete_filtered_rdb <- left_join(filtered_tab,genes)
  complete_filtered_rdb <- unique(complete_filtered_rdb)
  #complete_filtered_rdb <- complete_filtered_rdb[,c(2,5,1,3,4)]
  return(complete_filtered_rdb)
}





















