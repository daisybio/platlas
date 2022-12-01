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


































