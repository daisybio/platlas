####panther script 2
library(topGO)
library(GO.db) #gibt es auch Panther?
library(biomaRt)
library(Rgraphviz)
library(PANTHER.db)

de_file <- read.csv2("/Users/leonoraka/Desktop/Projekte/platlas/www/all_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV")
genes_universe <- de_file$geneID
# testgenes <- c("ENSG00000143797","ENSG00000142611","ENSG00000142046","ENSG00000162814","ENSG00000185551","ENSG00000072364","ENSG00000166822","ENSG00000084453","ENSG00000182732",
# "ENSG00000049323","ENSG00000163374","ENSG00000141161","ENSG00000157741","ENSG00000148516","ENSG00000186115","ENSG00000106031","ENSG00000165730","ENSG00000112414",
# "ENSG00000070444","ENSG00000140992","ENSG00000274672","ENSG00000273896","ENSG00000274796","ENSG00000278712","ENSG00000274078","ENSG00000275165","ENSG00000278622","ENSG00000274129",
# "ENSG00000278605","ENSG00000170892","ENSG00000163520","ENSG00000112319","ENSG00000078043","ENSG00000109458","ENSG00000126070","ENSG00000137601","ENSG00000117650",
# "ENSG00000187741","ENSG00000159208","ENSG00000147421","ENSG00000187735","ENSG00000086848","ENSG00000163873","ENSG00000196632","ENSG00000115159","ENSG00000128872",
# "ENSG00000175497","ENSG00000146469","ENSG00000135074","ENSG00000007047","ENSG00000115827","ENSG00000175336","ENSG00000108443","ENSG00000102753","ENSG00000155545",
# "ENSG00000087088","ENSG00000078674","ENSG00000140090","ENSG00000285479","ENSG00000151067","ENSG00000177663","ENSG00000182158","ENSG00000174574","ENSG00000134874",
# "ENSG00000166913","ENSG00000173517","ENSG00000182247","ENSG00000100784","ENSG00000175348","ENSG00000099899","ENSG00000136158","ENSG00000081026","ENSG00000109184",
# "ENSG00000153531","ENSG00000167034","ENSG00000128510","ENSG00000128510","ENSG00000164258","ENSG00000090975","ENSG00000117834","ENSG00000115548","ENSG00000154380",
# "ENSG00000107796","ENSG00000285479","ENSG00000151067","ENSG00000155849","ENSG00000175470","ENSG00000106952","ENSG00000198944","ENSG00000139679","ENSG00000109670",
# "ENSG00000285479","ENSG00000151067","ENSG00000173218","ENSG00000213626","ENSG00000113368","ENSG00000064419","ENSG00000103150","ENSG00000163626","ENSG00000164949",
# "ENSG00000143494","ENSG00000143603","ENSG00000077942","ENSG00000168283","ENSG00000198908","ENSG00000168802","ENSG00000154760","ENSG00000276495","ENSG00000034053",
# "ENSG00000198040","ENSG00000155974","ENSG00000106799","ENSG00000101384","ENSG00000023171","ENSG00000129084"
# )

significant_df <- de_file[abs(de_file$log2FoldChange) > 0.5 & de_file$padj < 0.05,]
significant_df <- significant_df[!is.na(significant_df$geneID),]
testgenes <- significant_df$geneID

# topGO -------------------------------------------------------------------
# gene universe file (complete de_file)
# de_file <- read.csv2("/Users/leonoraka/Desktop/Projekte/platlas/www/all_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV")
# genes_universe <- de_file$geneID
testgenes <- testgenes[(testgenes %in% genes_universe)]

# create GO db for genes to be used using biomaRt - please note that this takes a while
db= useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="https://www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'ensembl_gene_id', 'namespace_1003'), filters='ensembl_gene_id', values=genes_universe, mart=db)
gene_2_GO=unstack(go_ids[,c(1,2)])

# remove any candidate genes without GO annotation
keep = testgenes %in% go_ids[,2]
keep =which(keep==TRUE)
testgenes=testgenes[keep]

# make named factor showing which genes are of interest
geneList=factor(as.integer(genes_universe %in% testgenes))
names(geneList)= genes_universe

#build object
GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)

# define test using the classic algorithm with fisher (refer to [1] if you want to understand how the different algorithms work)
classic_fisher_result = runTest(GOdata, algorithm='classic', statistic='fisher') #hierarchy not taken into consideration

# define test using the weight01 algorithm (default) with fisher
weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher')  # better

# generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
allGO=usedGO(GOdata)
all_res=GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))

#performing BH correction on our p values
p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_final=cbind(all_res,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]

#get list of significant GO before multiple testing correction
results.table.p= all_res_final[which(all_res_final$weightFisher<=0.001),]

#get list of significant GO after multiple testing correction
results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]

# PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
showSigOfNodes(GOdata, score(weight_fisher_result), useInfo = "none", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)


myterms =results.table.p$GO.ID # change it to results.table.bh$GO.ID if working with BH corrected values
mygenes = genesInTerm(GOdata, myterms)

var=c()
for (i in 1:length(myterms))
{
  myterm=myterms[i]
  mygenesforterm= mygenes[myterm][[1]]
  mygenesforterm=paste(mygenesforterm, collapse=',')
  var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}


# mit clusterProfiler -----------------------------------------------------

# install and load the clusterProfiler package
#install.packages("clusterProfiler")
library(clusterProfiler)

# define the list of genes
#genes <- c("TP53", "BRCA1", "BRCA2", "PTEN")

# run the FEA analysis using the GO database
res <- enrichGO(gene=testgenes,universe = genes_universe,OrgDb="org.Hs.eg.db", ont="BP",pAdjustMethod = "BH",keyType = 'ENSEMBL', readable =FALSE,pvalueCutoff=0.05)#,pvalueCutoff=0.05
#test2_BH <- enrichGO(gene=testgenes,universe = genes_universe,OrgDb="org.Hs.eg.db",pAdjustMethod = "BH",ont="BP",keyType = 'ENSEMBL',pvalueCutoff=0.05)
#test3_fdr <- enrichGO(gene=testgenes,universe = genes_universe,OrgDb="org.Hs.eg.db",pAdjustMethod = "fdr",ont="BP",keyType = 'ENSEMBL',pvalueCutoff=0.05)
#test4_bonferroni <- enrichGO(gene=testgenes,universe = genes_universe,OrgDb="org.Hs.eg.db",pAdjustMethod = "bonferroni",ont="BP",keyType = 'ENSEMBL',pvalueCutoff=0.05)
  #enrichGO(gene=genes, OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.05)
#GSEAres<- gseGO(gene = testgenes, OrgDb = "org.Hs.eg.db",ont = "BP",keyType = 'ENSEMBL',pvalueCutoff = 0.05,verbose = FALSE)
# create a table with the ID, Description, Count, GeneRatio, pvalue, GeneId, List, Pop.GeneRatio, Pop.GeneId, OddsRatio, FDR column
erg <- res@result#[, c("ID", "Description", "Count", "GeneRatio", "pvalue", "GeneId", "List", "Pop.GeneRatio", "Pop.GeneId", "OddsRatio", "FDR")]#"List", "Pop.GeneRatio", "Pop.GeneId", "OddsRatio", "FDR"
erg <- mutate(erg, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))
#library(dplyr)
erg <- tidyr::separate(erg,col=GeneRatio, into=c("Count2", "List Total"), sep="/")
erg <- tidyr::separate(erg,col=BgRatio, into=c("Pop Hits", "Pop Total"), sep="/")
erg$percent <- as.numeric(erg$Count2)/as.numeric(erg$`List Total`)
erg <- erg[,c("ID","Description","Count","percent","pvalue","geneID","List Total","Pop Hits","Pop Total","richFactor","qvalue")]

# rename the columns
names(erg) <- c("Category", "Term", "Count", "%", "PValue", "Genes", "List Total", "Pop Hits", "Pop Total", "Fold Enrichment", "q-value (BH)") #eig letztes sollte FDR sein

# view the table
print(table)



# Erstellung einer Funktion  ----------------------------------------------

functional_enrichment <- function(significant_genes, universal_genes, ontology, annotation_type){
  res <- enrichGO(gene=significant_genes,universe = universal_genes,OrgDb="org.Hs.eg.db", ont=ontology,pAdjustMethod = "BH",keyType = annotation_types)#'ENSEMBL'
  erg <- res@result#[, c("ID", "Description", "Count", "GeneRatio", "pvalue", "GeneId", "List", "Pop.GeneRatio", "Pop.GeneId", "OddsRatio", "FDR")]#"List", "Pop.GeneRatio", "Pop.GeneId", "OddsRatio", "FDR"
  erg <- mutate(erg, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))
  #library(dplyr)
  erg <- tidyr::separate(erg,col=GeneRatio, into=c("Count2", "List Total"), sep="/")
  erg <- tidyr::separate(erg,col=BgRatio, into=c("Pop Hits", "Pop Total"), sep="/")
  erg$percent <- as.numeric(erg$Count2)/as.numeric(erg$`List Total`)
  erg <- erg[,c("ID","Description","Count","percent","pvalue","geneID","List Total","Pop Hits","Pop Total","richFactor","qvalue")]
  
  # rename the columns
  names(erg) <- c("Category", "Term", "Count", "%", "PValue", "Genes", "List Total", "Pop Hits", "Pop Total", "Fold Enrichment", "q-value (BH)") #eig letztes sollte FDR sein
  return(erg)
}

functional_enrichment_input <- function(dataset, pval, foldchange, annot){
  upregulated <- dataset[dataset$log2FoldChange > foldchange & dataset$padj < pval,]
  downregulated <- dataset[dataset$log2FoldChange < -foldchange & dataset$padj < pval,]
  upregulated_genes <- upregulated[!is.na(upregulated),annot]
  downregulated_genes <- downregulated[!is.na(downregulated),annot]
  upregulated_genes <- as.vector(upregulated_genes)
  downregulated_genes <- as.vector(downregulated_genes)
  all_significant <- dataset[abs(dataset$log2FoldChange) > foldchange & dataset$padj < pval,]
  all_significant_genes <- all_significant[!is.na(all_significant),annot]
  all_significant_genes <- as.vector(all_significant_genes)
  universe <- as.vector(dataset[,annot])
  return(list("upregulated" = upregulated_genes, "downregulated" = downregulated_genes,"all" = all_significant_genes,"universe" = universe))
}



