######## get miRNA tar & rbs mappings 

library(openxlsx)

import_css_miRNA_de <- read.xlsx("/Users/leonoraka/Desktop/Projekte/final_plateletAnalysis_CCSonly/miRNA_analysis/results/ccs_rp_vs_mp_diff_mir_counts_noSample6_filtered_pCutoff0.05.xlsx")
#vergiss nich abzuspeichern mit .CSV & X statt X1 colname



old_css_miRNA_de <- read.csv2("/Users/leonoraka/Desktop/Projekte/platlas/www/cad_rp_vs_mp_diff_mir_counts.CSV")

new_miRNA <- import_css_miRNA_de$X1
old_miRNA <- old_css_miRNA_de$X

#new_miRNA[!(new_miRNA %in% old_miRNA)]

#old_miRNA[!(old_miRNA %in% new_miRNA)]

tar_db <- read.xlsx("/Users/leonoraka/Downloads/hsa_MTI.xlsx")
rdb_db <- fread("/Users/leonoraka/Downloads/miRDB_v6.0_prediction_result.txt")
rdb_db_human <- rdb_db[grepl("hsa",rdb_db$V1),]


tar_db <- tar_db[,c(2,1,3,4,5,6,7,8,9)]
filtered_tar_db <- tar_db[tar_db$miRNA %in% new_miRNA,]
experiment_tar_db <- filtered_tar_db[!grepl("(Weak)",filtered_tar_db$Support.Type)& !is.na(filtered_tar_db$Support.Type),]
all_tar_db_filtered <- filtered_tar_db
only_weak_tar <- filtered_tar_db[grepl("(Weak)",filtered_tar_db$Support.Type) | is.na(filtered_tar_db$Support.Type),]

filtered_rdb_db_human <- rdb_db_human[rdb_db_human$V1 %in% new_miRNA]


#translate rdb gene format
attributes <- c("refseq_mrna", "ensembl_gene_id")
filters <- "refseq_mrna"
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#attributes <- listAttributes(mart)
#filters <- listFilters(mart)
genes <- getBM(attributes = attributes,
               filters = filters,
               values = filtered_rdb_db_human$V2,
               mart = mart)



colnames(filtered_rdb_db_human)[2] = "refseq_mrna"
complete_filtered_rdb <- left_join(filtered_rdb_db_human,genes)
complete_filtered_rdb <- unique(complete_filtered_rdb)
complete_filtered_rdb <- complete_filtered_rdb[,c(2,4,1,3)]


#list up & downregulated & sig
upreg_miRNA <- import_css_miRNA_de$X1[import_css_miRNA_de$log2FoldChange > 0]
downreg_miRNA <- import_css_miRNA_de$X1[import_css_miRNA_de$log2FoldChange < 0]
sig_upreg <- import_css_miRNA_de$X1[import_css_miRNA_de$log2FoldChange > 0.5 & import_css_miRNA_de$padj <= 0.05 & !is.na(import_css_miRNA_de$padj)] #zur sicherheit weil eig nicht p<= 0.05 & log2FC 0.5
sig_downreg <- import_css_miRNA_de$X1[import_css_miRNA_de$log2FoldChange < -0.5 &import_css_miRNA_de$padj < 0.05& !is.na(import_css_miRNA_de$padj)]


up_tar_exp <- experiment_tar_db[experiment_tar_db$miRNA %in% upreg_miRNA, ]
up_tar_all <- all_tar_db_filtered[all_tar_db_filtered$miRNA %in% upreg_miRNA, ]
up_tar_weak <- only_weak_tar[only_weak_tar$miRNA %in% upreg_miRNA, ]
up_rdb <- complete_filtered_rdb[complete_filtered_rdb$V1 %in% upreg_miRNA, ]

down_tar_exp <- experiment_tar_db[experiment_tar_db$miRNA %in% downreg_miRNA, ]
down_tar_all <- all_tar_db_filtered[all_tar_db_filtered$miRNA %in% downreg_miRNA, ]
down_tar_weak <- only_weak_tar[only_weak_tar$miRNA %in% downreg_miRNA, ]
down_rdb <- complete_filtered_rdb[complete_filtered_rdb$V1 %in% downreg_miRNA, ]

sig_up_tar_exp <- experiment_tar_db[experiment_tar_db$miRNA %in% sig_upreg, ]
sig_up_tar_all <- all_tar_db_filtered[all_tar_db_filtered$miRNA %in% sig_upreg, ]
sig_up_tar_weak <- only_weak_tar[only_weak_tar$miRNA %in% sig_upreg, ]
sig_up_rdb <- complete_filtered_rdb[complete_filtered_rdb$V1 %in% sig_upreg, ]

sig_down_tar_exp <- experiment_tar_db[experiment_tar_db$miRNA %in% sig_downreg, ]
sig_down_tar_all <- all_tar_db_filtered[all_tar_db_filtered$miRNA %in% sig_downreg, ]
sig_down_tar_weak <- only_weak_tar[only_weak_tar$miRNA %in% sig_downreg, ]
sig_down_rdb <- complete_filtered_rdb[complete_filtered_rdb$V1 %in% sig_downreg, ]




fp <- "/Users/leonoraka/Desktop/Projekte/platlas/www/miRNA/new_results/"
write.table(up_tar_exp, paste0(fp,"E_up_tar_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(up_tar_weak, paste0(fp,"up_tar_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(up_rdb, paste0(fp,"up_rdb_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")

write.table(down_tar_exp, paste0(fp,"E_down_tar_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(down_tar_weak, paste0(fp,"down_tar_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(down_rdb, paste0(fp,"down_rdb_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")

write.table(sig_up_tar_exp, paste0(fp,"E_up_tar_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(sig_up_tar_weak, paste0(fp,"up_tar_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(sig_up_rdb, paste0(fp,"up_rdb_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")

write.table(sig_down_tar_exp, paste0(fp,"E_down_tar_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(sig_down_tar_weak, paste0(fp,"down_tar_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(sig_down_rdb, paste0(fp,"down_rdb_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")

###Nodes

N_rdb_tr <- translate_to_ensembl(complete_filtered_rdb$ensembl_gene_id)
N_tar_exp_tr <- translate_to_hgnc(experiment_tar_db$Target.Gene)
N_tar_weak_tr <- translate_to_hgnc(only_weak_tar$Target.Gene)

translate_to_hgnc <- function(to_translate){
  attributes <- c("hgnc_symbol", "ensembl_gene_id")
  filters <- "hgnc_symbol"
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  #attributes <- listAttributes(mart)
  #filters <- listFilters(mart)
  genes <- getBM(attributes = attributes,
               filters = filters,
               values = to_translate,
               mart = mart)
  return(genes)
}
translate_to_ensembl <- function(to_translate){
  attributes <- c("ensembl_gene_id", "hgnc_symbol")
  filters <- "ensembl_gene_id"
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  #attributes <- listAttributes(mart)
  #filters <- listFilters(mart)
  genes <- getBM(attributes = attributes,
                 filters = filters,
                 values = to_translate,
                 mart = mart)
  return(genes)
}

#complete_filtered_rdb
rdb_N <- left_join(complete_filtered_rdb,N_rdb_tr)
for_man_ex_tar <-experiment_tar_db
for_man_ex_tar <- for_man_ex_tar[,c(1,4)]
colnames(for_man_ex_tar)[2] <- "hgnc_symbol"
for_man_weak_tar <- only_weak_tar
for_man_weak_tar <- for_man_weak_tar[,c(1,4)]
colnames(for_man_weak_tar)[2] <- "hgnc_symbol"

tar_N_expr <- left_join(for_man_ex_tar, N_tar_exp_tr)
tar_N_weak <- left_join(for_man_weak_tar, N_tar_weak_tr)

rdb_N <- unique(rdb_N)
tar_N_expr <- unique(tar_N_expr)
tar_N_weak <- unique(tar_N_weak)

rdb_N <- rdb_N[,c(3,5,2)]

###Nodes 


N_up_tar_exp <- tar_N_expr[tar_N_expr$miRNA %in% upreg_miRNA, ]
#up_tar_all <- all_tar_db_filtered[all_tar_db_filtered$miRNA %in% upreg_miRNA, ]
N_up_tar_weak <- tar_N_weak[tar_N_weak$miRNA %in% upreg_miRNA, ]
N_up_rdb <- rdb_N[rdb_N$V1 %in% upreg_miRNA, ]

N_down_tar_exp <- tar_N_expr[tar_N_expr$miRNA %in% downreg_miRNA, ]
#down_tar_all <- all_tar_db_filtered[all_tar_db_filtered$miRNA %in% downreg_miRNA, ]
N_down_tar_weak <- tar_N_weak[tar_N_weak$miRNA %in% downreg_miRNA, ]
N_down_rdb <- rdb_N[rdb_N$V1 %in% downreg_miRNA, ]

N_sig_up_tar_exp <- tar_N_expr[tar_N_expr$miRNA %in% sig_upreg, ]
#sig_up_tar_all <- all_tar_db_filtered[all_tar_db_filtered$miRNA %in% sig_upreg, ]
N_sig_up_tar_weak <- tar_N_weak[tar_N_weak$miRNA %in% sig_upreg, ]
N_sig_up_rdb <- rdb_N[rdb_N$V1 %in% sig_upreg, ]

N_sig_down_tar_exp <- tar_N_expr[tar_N_expr$miRNA %in% sig_downreg, ]
#sig_down_tar_all <- all_tar_db_filtered[all_tar_db_filtered$miRNA %in% sig_downreg, ]
N_sig_down_tar_weak <- tar_N_weak[tar_N_weak$miRNA %in% sig_downreg, ]
N_sig_down_rdb <- rdb_N[rdb_N$V1 %in% sig_downreg, ]



write.table(N_up_tar_exp, paste0(fp,"E_N_up_tar_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(N_up_tar_weak, paste0(fp,"N_up_tar_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(N_up_rdb, paste0(fp,"N_up_rdb_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")

write.table(N_down_tar_exp, paste0(fp,"E_N_down_tar_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(N_down_tar_weak, paste0(fp,"N_down_tar_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(N_down_rdb, paste0(fp,"N_down_rdb_cad.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")

write.table(N_sig_up_tar_exp, paste0(fp,"E_N_up_tar_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(N_sig_up_tar_weak, paste0(fp,"N_up_tar_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(N_sig_up_rdb, paste0(fp,"N_up_rdb_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")

write.table(N_sig_down_tar_exp, paste0(fp,"E_N_down_tar_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(N_sig_down_tar_weak, paste0(fp,"N_down_tar_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(N_sig_down_rdb, paste0(fp,"N_down_rdb_cad_sig.tsv"), col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")


nodes_path <- "/Users/leonoraka/Desktop/Projekte/platlas/www/miRNA/new_results/Nodes/"
all_nodes <- list.files(nodes_path)
for (i in all_nodes) {
  read_path <- paste0(nodes_path,i)
  tab <- fread(read_path)
  write_path <- paste0(nodes_path,"Zgenes/",i)
  write.table(tab[,3], write_path , col.names = FALSE, row.names = FALSE,quote = FALSE, sep = "\t")
}


