# FER <- eventReactive(c(input$clickGO, input$clickGO2),ignoreInit = T,{
#   tabBox(
#     width = 12,
#     height = "1000px",
#     title = "Get the results",
library(enrichplot)
library(clusterProfiler)
library(DescTools)
#library(qvalue)
library(GOplot)
library(stats)
library(wesanderson)
library("RColorBrewer")
# GO ----------------------------------------------------------------------
#path3 <- "C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/GOE_"
#C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/
functional_enrichment <- function(significant_genes_up,significant_genes_down,universal_genes, ontology, annotation_type){
  res_up <- enrichGO(gene=significant_genes_up,universe = universal_genes,OrgDb="org.Hs.eg.db", ont=ontology,pAdjustMethod = "BH",keyType = annotation_type)#'ENSEMBL', pvalueCutoff = 0.05
  res_down <- enrichGO(gene=significant_genes_down,universe = universal_genes,OrgDb="org.Hs.eg.db", ont=ontology,pAdjustMethod = "BH",keyType = annotation_type)#'ENSEMBL', pvalueCutoff = 0.05
  
  erg_up <- res_up@result#[, c("ID", "Description", "Count", "GeneRatio", "pvalue", "GeneId", "List", "Pop.GeneRatio", "Pop.GeneId", "OddsRatio", "FDR")]#"List", "Pop.GeneRatio", "Pop.GeneId", "OddsRatio", "FDR"
  erg_down <- res_down@result
  erg_up$significantly_expressed <- "upregulated"
  erg_down$significantly_expressed <- "downregulated"
  erg <- rbind(erg_up,erg_down)
  
  erg <- mutate(erg, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))
  #library(dplyr)
  erg <- tidyr::separate(erg,col=GeneRatio, into=c("Count2", "List Total"), sep="/")
  erg <- tidyr::separate(erg,col=BgRatio, into=c("Pop Hits", "Pop Total"), sep="/")
  erg$percent <- as.numeric(erg$Count2)/as.numeric(erg$`List Total`)
  # print(erg$pvalue)
  # print("is infinite")
  # print(is.infinite(erg$pvalue))
  # print("is na")
  # print(is.na(erg$pvalue))
  # print("is null")
  # print(is.null(erg$pvalue))
  
  erg$fdr <-erg$p.adjust#p.adjust(erg$pvalue, method="BH")#[!is.na(erg$pvalue) & !is.infinite(erg$pvalue) & !is.null(erg$pvalue)]#
  erg$GO <- ontology
  erg <- erg[,c("GO","ID","Description","Count","percent","pvalue","geneID","p.adjust","List Total","Pop Hits","Pop Total","richFactor","qvalue","fdr", "significantly_expressed")]
  
  
  # rename the columns
  names(erg) <- c("GO","Category", "Term", "Count", "%", "PValue", "Genes", "p.adjust","List Total", "Pop Hits", "Pop Total", "Fold Enrichment", "q-value","FDR","significantly_expressed" ) #eig letztes sollte FDR sein
  return(list("tab" = erg, "obj_up" = res_up, "obj_down" = res_down))
}

gse_analysis <- function(genes, ontology, annotation_type){
  res <- gseGO(gene=genes,OrgDb="org.Hs.eg.db", ont=ontology,minGSSize = 1,maxGSSize = 500,pvalueCutoff =0.05, keyType = annotation_type)#'ENSEMBL', pvalueCutoff = 0.05
  #View(res)
  #print(res@result[["core_enrichment"]][1])
  #print(res@result[["ID"]][1])
  #print(length(unlist(res@result[["core_enrichment"]][1],"/")))
  #print(res@geneSets[["GO:0016236"]])
  #print(length(res@geneSets[["GO:0016236"]]))
  #res_down <- gseGO(gene=significant_genes_down_list,OrgDb="org.Hs.eg.db", ont=ontology,minGSSize = 1,maxGSSize = 500,pvalueCutoff = 0.05,keyType = annotation_type, by = "fgsea")
  return(list("obj" = res))# "obj_down" = res_down
}


functional_enrichment_input <- function(dataset, pval, foldchange, annot){
  upregulated <- dataset[dataset$log2FoldChange > foldchange & dataset$padj < pval,]
  downregulated <- dataset[dataset$log2FoldChange < -foldchange & dataset$padj < pval,]
  
  ##for gsea
 
  annot_name <- colnames(dataset)[annot]
 # upreg <- dataset[dataset$log2FoldChange > 0,]
#  downreg <- dataset[dataset$log2FoldChange < 0,]
 # View(upreg)
  gsea_list <- dataset[!is.na(dataset),c("log2FoldChange",annot_name)]
  colnames(gsea_list) = c("log2FoldChange","gene")
  gsea_list <- setNames(gsea_list$log2FoldChange, gsea_list$gene)
  gsea_list <- sort(gsea_list, decreasing = TRUE)
  #up_list <- upreg[!is.na(upreg),c("log2FoldChange",annot_name)]
  #colnames(up_list) = c("log2FoldChange","gene")
  #uplist <- setNames(up_list$log2FoldChange, up_list$gene)
  #uplist <- sort(uplist, decreasing = TRUE)
 # View(uplist)
  #down_list <- downreg[!is.na(downreg),c("log2FoldChange",annot_name)]
  #colnames(down_list) = c("log2FoldChange","gene")
  #downlist <- setNames(down_list$log2FoldChange, down_list$gene)
  #downlist <- sort(downlist, decreasing = TRUE)
  ###
  
  upregulated_genes <- upregulated[!is.na(upregulated),annot]
  downregulated_genes <- downregulated[!is.na(downregulated),annot]
  upregulated_genes <- as.vector(upregulated_genes)
  downregulated_genes <- as.vector(downregulated_genes)
  all_significant <- dataset[abs(dataset$log2FoldChange) > foldchange & dataset$padj < pval,]
  all_significant_genes <- all_significant[!is.na(all_significant),annot]
  all_significant_genes <- as.vector(all_significant_genes)
  universe <- as.vector(dataset[,annot])
  return(list("upregulated" = upregulated_genes, "downregulated" = downregulated_genes,"all" = all_significant_genes,"universe" = universe, "gsea_input" = gsea_list))#"up_list" = uplist, "down_list"= downlist
}

pathGO <- "www/"
#pGO2 <- "./www/"

clickedUDA1 = reactiveVal(isolate(input$clickGO))
clickedUDA2 = reactiveVal(isolate(input$clickGO2))


saved_GO = reactiveVal()
saved_circ_mat = reactiveVal()
saved_ontology = reactiveVal()
saved_GO_obj_up = reactiveVal()
saved_GO_obj_down = reactiveVal()
current_UDA = reactiveVal()
other_UDA = reactiveVal()

saved_GSEA_up = reactiveVal()
saved_GSEA_down = reactiveVal()
saved_GSEA = reactiveVal()
saved_annot4heatmap = reactiveVal()
#sav_cir =reactiveVal()
#sav_go_h_p = reactiveVal()
#sav_edges = reactiveVal()
GO <-  eventReactive(input$clickGO,{#c(input$clickGO,input$clickGO2)
  finalpath <- ""
  ontology <- ""
  UDA <- ""
  pval <- 0
  fc <- 0
  annot <- 0
  annotation_type <- ""
  HoDGO <- "disease - CAD"
  TPMGO <-  "TPM > 0.2"
  #adjustment_method <- ""
  #if(clickedUDA1() < input$clickGO){
    finalpath <-paths2frames("1",path,TPMGO,HoDGO,".CSV")
    ontology <- switch(input$GOSS,"biological process"="BP","cellular component"= "CC","molecular function" = "MF")
    #UDA <- switch(input$UDA, "upregulated" = "upregulated","downregulated" = "downregulated","all" = "all")
    annot <- switch(input$DVP1GO, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
    annot4heatmap <- switch(input$DVP1GO, "ENSEMBL ID"= 2, "HGNC symbol" = 1)
    saved_annot4heatmap(annot4heatmap)
    annotation_type <- switch(input$DVP1GO, "ENSEMBL ID"= "ENSEMBL", "HGNC symbol" = "SYMBOL")
    #adjustment_method <- switch(input$AdjMethod, "Benjamini-Hochberg"= "BH", "FDR" = "fdr")
    pval <- input$pCO1GO
    fc <- input$fcO1GO
    clickedUDA1(input$clickGO)
  # }else if(clickedUDA2() < input$clickGO2){
  #   finalpath <-paths2frames("1",path,input$TPM2GO,input$HoD2GO,".CSV")
  #   ontology <- switch(input$GOSS2,"biological process"="BP","cellular component"= "CC","molecular function" = "MF")
  #   #UDA <- switch(input$UDA2, "upregulated" = "upregulated","downregulated" = "downregulated","all" = "all")
  #   annot <- switch(input$DVP2GO, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
  #   annotation_type <- switch(input$DVP2GO, "ENSEMBL ID"= "ENSEMBL", "HGNC symbol" = "SYMBOL")
  #   #adjustment_method <- switch(input$AdjMethod2, "Benjamini-Hochberg"= "BH", "FDR" = "fdr")
  #   pval <- input$pCO2GO
  #   fc <- input$fcO2GO
  #   clickedUDA2(input$clickGO2)
  # }
  print(finalpath)
  if(finalpath != ""){
  matGO <- read.csv2(finalpath)
  annot_colname <- colnames(matGO)[annot]
  saved_circ_mat(matGO[,c(annot_colname,"log2FoldChange")])
  saved_ontology(ontology)
  input_list <- functional_enrichment_input(matGO, pval, fc, annot)
  
  sig_genes_up <- input_list[["upregulated"]]
  sig_genes_down <- input_list[["downregulated"]]
  
  universe_genes <- input_list[["universe"]]
  functional_enrichment_result <- functional_enrichment(sig_genes_up,sig_genes_down, universe_genes, ontology, annotation_type)
  #TODO
  gse_analysis_result <- gse_analysis(input_list[["gsea_input"]],ontology, annotation_type)
  
  
  #go_bar_result_upreg <- functional_enrichment(input_list[[]], universe_genes, ontology, annotation_type)
  significant_tab<- functional_enrichment_result$tab[functional_enrichment_result$tab$p.adjust <= 0.05,]
  saved_GO(significant_tab)
  saved_GO_obj_up(functional_enrichment_result$obj_up)
  saved_GO_obj_down(functional_enrichment_result$obj_down)
  saved_GSEA(gse_analysis_result$obj)
  saved_GSEA_result(gse_analysis_result$obj@result)
 # saved_GSEA_down(gse_analysis_result$obj_down)
  #print(saved_GO_obj_down())
  
  }
})

FEannot <- eventReactive(input$clickFeA,{
  res <- switch(input$FeA, "ENSEMBL ID" = 6, "ENTREZ ID" = 7, "HGNC symbol" = 8)
  res
})

output$select_GSEA <- renderUI({
  data <- NULL #"upregulated","downregulated"
  if(!is.null(saved_GSEA())){
    data <- saved_GSEA()
  }
  index_tab <- which(colnames(data@result) == "Description")
  colnames(data@result)[index_tab] = "Term"
  data@result <- filterGO(data@result,input$gos_of_interest)
  #View(data@result)
  selectInput("select_GSEA_term",
              label = "Select gene set of interest",
              choices = as.list(data@result$Term)
              )
})

output$select_GSEA2 <- renderUI({
  data <- NULL #"upregulated","downregulated"
  if(!is.null(saved_GSEA())){
    data <- saved_GSEA()
  }
  index_tab <- which(colnames(data@result) == "Description")
  colnames(data@result)[index_tab] = "Term"
  data@result <- filterGO(data@result,input$gos_of_interest)
  #View(data@result)
  selectInput("select_GSEA_term2",
              label = "Select gene set of interest",
              choices = as.list(data@result$Term)
  )
})

saved_GSEA_plot = reactiveVal()
saved_GSEA_result = reactiveVal()
saved_GSEA_plot_index = reactiveVal()
output$GO_GSEA <- renderPlot({
  data <- NULL #"upregulated","downregulated"
  if(!is.null(saved_GSEA())){
    data <- saved_GSEA()
    go_descr <- input$select_GSEA_term
  index_gsea <- which(data$Description == go_descr)
  saved_GSEA_plot_index(index_gsea)
  #View(data)
  
  saved_GSEA_plot(data)
  #temptab <- data@result
  #trace("gseaplot2", edit = TRUE)
  #View(data)
  
  try(my_gseaplot2(data, geneSetID = index_gsea , title = data$Description[index_gsea], pvalue_table = TRUE),silent = TRUE)
 
  }
  

})
saved_GSEA_heat_index = reactiveVal()
saved_GSEA_heat = reactiveVal()
saved_GSEA_heatmap = reactiveVal()
output$GSEA_heat <- renderPlot({
  data <- NULL #"upregulated","downregulated"
  if(!is.null(saved_GSEA())){
    data <- saved_GSEA()
    go_descr <- input$select_GSEA_term2
    index_gsea <- which(data$Description == go_descr)
    saved_GSEA_heat_index(index_gsea)
    #View(data)
    
    saved_GSEA_heat(data)
    #View(data)
    #print(unlist(str_split(data$core_enrichment[index_gsea],"/")))
    
    extract_genes <- as.vector(unlist(str_split(data$core_enrichment[index_gsea],"/")))
    #print(extract_genes)
    #fcs <- as.data.frame(data@geneList)
    #print(fcs)
    logFC <- data@geneList[extract_genes]
    lfc_df <- as.data.frame(logFC)
    #View(as.data.frame(list_with_fcs))
    #temptab <- data@result
    #trace("gseaplot2", edit = TRUE)
    #View(data)
    
    #try(my_gseaplot2(data, geneSetID = index_gsea , title = data$Description[index_gsea], pvalue_table = TRUE),silent = TRUE)
    #find countdf
    #map to sample names
    #extract genes
    ###hard coded bc no other data at the moment => diseased CCS
    merged_annotation <- data.table::fread("www/simulated_sample_names_all.tsv")
    countpath <- "www/heat_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV"
    merged_annotation <- as.data.frame(merged_annotation)
    rownames(merged_annotation) <- merged_annotation$original_sample_name
    merged_annotation$platelet_type <- gsub("s","",merged_annotation$RPs_MPs)
    #View(merged_annotation)
    ##heatmap table
    datei <- read.csv2(countpath)
    h <- data.frame(datei)
    h <- h[, -c(saved_annot4heatmap(),3)]
    names <-h[,1]
    h <- h[,-1]
    rownames(h) = make.names(names, unique=TRUE)
    
   # h <- h %>% select(sort(names(h),))
    neuer <- data.matrix(h)
    #zahlen manipulieren/normalisieren
    neuer<-log(neuer,2)
    neuer[neuer < 0] <- 0
    to_translate <- colnames(neuer)
    new_names <- merged_annotation[to_translate,]
    colnames(neuer) <- new_names$sample_name
    neuer <- neuer[ , !(colnames(neuer) %in% c("6_MPs","6_RPs"))]
    #View(neuer)
    
    
    ##filter leading edge
    filtered_heat <- neuer[rownames(lfc_df),]
    #order columns
    # Extract the suffixes
    suffixes <- sapply(strsplit(colnames(filtered_heat), "_"), function(x) tail(x, 1))
    prefixes <- sapply(strsplit(colnames(filtered_heat), "_"), function(x) head(x, 1))
    # Create a data frame of column names and their suffixes
    col_suffix_df <- data.frame(
      col_name = colnames(filtered_heat),
      suffix = suffixes,
      prefix = as.numeric(prefixes),
      stringsAsFactors = FALSE
    )
    # Order the columns by suffix and then by the numeric part of the column name
    col_suffix_df <- col_suffix_df[order(col_suffix_df$suffix, col_suffix_df$prefix),]
    
    # Get the new order of column names
    ordered_col_names <- col_suffix_df$col_name
    # Reorder the data frame columns
    filtered_heat<- filtered_heat[, ordered_col_names]
    #View(filtered_heat)
    #View(filtered_heat_new)

    acs <- c("25_MPs","25_RPs","26_MPs","26_RPs","27_MPs","27_RPs","28_MPs","28_RPs")
    filtered_heat <- filtered_heat[, !(colnames(filtered_heat) %in% acs)]
    
    #filtered_heat <- filtered_heat %>% select(sort(names(filtered_heat)))
    #threshold <- sort(rowVars(filtered_heat), decreasing = T)
    #filtered_heat = filtered_heat[which(rowVars(filtered_heat) >= threshold),]
    #counts <- counts[ , !(colnames(counts) %in% c("6_MPs","6_RPs"))]
    #View(filtered_heat)
    #print(type(filtered_heat))
    ###heatmap plot
    dbo <- data.frame("names" = merged_annotation$sample_name, "Condition" = merged_annotation$condition, "Platelet type" = merged_annotation$platelet_type)
    
    rownames(dbo) = dbo$names
    dbo <- dbo[,-1]
    
    #colnames(dbo) = c("Condition", "Platelet type")
    ##break
    #View(expression_table)
    ###
    try({
      pl_names = c("red","blue")#c(RPs = cbbPalette[1], MPs = cbbPalette[4])
    names(pl_names) = c("RP", "MP")
    d_names = c("lightblue","orange","darkgreen")
    names(d_names) = c("healthy","ACS", "CCS")
    ann_col = list(`Platelet type` = pl_names,Condition = d_names)
    dbo <- dbo[ colnames(filtered_heat), ]
    p<- pheatmap(filtered_heat, 
                 color = colorRampPalette(c("blue", "red"))(50),
                 border_color = NA,
                 annotation_colors = ann_col,
                 annotation_col = dbo,#dbo %>% select('Platelet type' = dbo$`Platelet type`, 'Disease' = dbo$Disease),
                 show_rownames = input$show_gene_names, show_colnames = TRUE, scale="row", cluster_cols = FALSE,
                 main = paste0("Expression of leading edge genes in: ", as.character(input$select_GSEA_term2))
    )
    saved_GSEA_heatmap(p)
    
    p
    })
    
    
    
  }
  
  
})
output$GSEATab <- renderDataTable({
  gsea <- saved_GSEA_result()
  index_tab <- which(colnames(gsea) == "Description")
  colnames(gsea)[index_tab] = "Term"
  gsea <- filterGO(gsea,input$gos_of_interest)
  gsea$core_enrichment <- str_replace_all(gsea$core_enrichment,"/", "; ")
  rownames(gsea) <- 1:nrow(gsea)
  datatable(
    cbind(' ' = '&oplus;', gsea), escape = -2,
    extensions = c("Buttons","Scroller"),
    options = list(
      columnDefs = list(
        list(visible = FALSE, targets = c(12)),
        list(orderable = FALSE, className = 'details-control', targets = 1)
      ),
      lengthMenu = c(5, 10, 15),
      width = 500,
      scrollX = TRUE,
      dom = 'Bfrtip',
      bPaginate = FALSE,
      buttons =list(list(
        extend = 'collection',
        buttons = c('csv', 'excel'),
        text = 'Download'
      )) 
    ),
    callback = JS("
  table.column(1).nodes().to$().css({cursor: 'pointer'});
  var format = function(d) {
    return '<div style=\"background-color:#eee; padding: .5em;\">  ' +
            'Genes: ' + d[12] + '</div>';
  };
  table.on('click', 'td.details-control', function() {
    var td = $(this), row = table.row(td.closest('tr'));
    if (row.child.isShown()) {
      row.child.hide();
      td.html('&oplus;');
    } else {
      row.child(format(row.data())).show();
      td.html('&CircleMinus;');
    }
  });"
    ))
                                         
})
output$GOTab <- renderDataTable({
  GO()
  gotab<- saved_GO()
  
  if(!is.null(gotab)){
  gotab<- as.data.table(gotab)
 
  gotab <- gotab[gotab$significantly_expressed == input$UDA, ]
  gotab <- filterGO(gotab,input$gos_of_interest)
  
  gotab <- gotab[, -"significantly_expressed"]

  
  gotab$Genes <- str_replace_all(gotab$Genes,"/", "; ")
  #annot <- FEannot()
  #gotab <- gotab[,c(1,2,3,4,5,annot,9)]
  datatable(
    cbind(' ' = '&oplus;', gotab), escape = -2,
    extensions = c("Buttons","Scroller"),
    options = list(
      columnDefs = list(
        list(visible = FALSE, targets = c(8)),
        list(orderable = FALSE, className = 'details-control', targets = 1)
      ),
      lengthMenu = c(5, 10, 15),
      width = 500,
      scrollX = TRUE,
      dom = 'Bfrtip',
      bPaginate = FALSE,
      buttons =list(list(
        extend = 'collection',
        buttons = c('csv', 'excel'),
        text = 'Download'
      )) 
    ),
    callback = JS("
  table.column(1).nodes().to$().css({cursor: 'pointer'});
  var format = function(d) {
    return '<div style=\"background-color:#eee; padding: .5em;\">  ' +
            'Genes: ' + d[8] + '</div>';
  };
  table.on('click', 'td.details-control', function() {
    var td = $(this), row = table.row(td.closest('tr'));
    if (row.child.isShown()) {
      row.child.hide();
      td.html('&oplus;');
    } else {
      row.child(format(row.data())).show();
      td.html('&CircleMinus;');
    }
  });"
    ))
  }
  }#escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 15, width = 700)
)
# output$GOplot <- renderPlot({
#   GOdf <- GO()
#   #b <- subset(GOdf,GOdf[,7] > 2)
#   b <- b[order(-b[,7]), ]
#   nbl <- unique(head(b[, c(6,7)],10),incomparables=FALSE, fromLast=FALSE, by = nbl[,1])
#   ggplot(nbl, aes(x = Functional_Category, y = Enrichment_FDR)) +        # Create barchart with ggplot2
#     geom_bar(stat = "identity",fill = "blue")+
#     coord_flip()
# 
# })
output$GOpG <- renderPlotly({
  tabi <- saved_GO()
  tabi <- tabi[tabi$significantly_expressed == input$UDA, ]
  tabi <- filterGO(tabi,input$gos_of_interest)
  tabi$percentage <- tabi$'%' *100
  #print(tabi$percentage)
  tabi <- as.data.table(tabi)
  tabi$percentage <- as.numeric(tabi$percentage)
  tabi <- tabi[order(tabi$percentage,decreasing = FALSE),]

  # form <- list(categoryorder = "array",
  #             categoryarray = tabi$id,
  #             title = "category")
  fig <- plot_ly(tabi, x = ~Category, y = ~percentage, type = 'bar',color = ~p.adjust,colors =c("#ba3131","#fff68f"),name = 'percentage of differentially expressed genes in category',text = ~Term)
  #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))
  fig <- fig %>% layout(title = "Percentage of differentially expressed genes in GO categories")
  fig
  #barmode = 'group'
 # )
  # fig
  
})

fthreshold <- eventReactive(input$clickGOFE,{
  num <- as.numeric(input$GOFE)
  end <- num
})
output$GOpFE <- renderPlotly({
  # n <- fthreshold()
  testtab <- saved_GO()#f_foldenrich(saved_GO(),n)
  testtab <- testtab[testtab$significantly_expressed == input$UDA, ]
  testtab <- filterGO(testtab,input$gos_of_interest)
  
  testtab <- testtab[order(testtab$FDR, decreasing = FALSE),]
  testtab <- testtab %>% mutate(FDR = -log10(FDR))
  form <- list(categoryorder = "array",
               categoryarray = testtab$Category,
               title = "category")
  print(colnames(testtab))
  fig2 <- plot_ly(testtab, x =~FDR, y = ~Category, type = 'bar', name = 'enrichment FDR in categories',color = ~p.adjust,colors = c("#ba3131","#fff68f"),text = ~Term, orientation = 'h' )
  #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))
  fig2 <- fig2 %>% layout(title = "Enrichment FDR in categories", xaxis = list(title = "enrichment FDR"),
                          yaxis = form,
                          margin = list(b = 100))
  fig2
  
})
output$GOp <- renderPlotly({
  gotab <- saved_GO()
  gotab <- gotab[gotab$significantly_expressed == input$UDA & gotab$p.adjust <= 0.05, ]
  gotab <- filterGO(gotab,input$gos_of_interest)
  #circ_data <- make_circ_data(gotab,saved_circ_mat(),saved_ontology())
  #circ_data_sorted <- circ_data[,-4] %>% arrange(adj_pval)

  gotab <- gotab%>% arrange(p.adjust)
  circ_end <- gotab[1:10,]
  form <- list(categoryorder = "array",
               categoryarray = circ_end$Category,
               title = "Category")
  fig2 <- plot_ly(circ_end, x =~Category, y = ~p.adjust, type = 'bar', name = 'Top 10 categories with lowest adjusted p-values',text = ~Term, marker = list(color = '#ba3131'))
  fig2 <- fig2 %>% layout(title = "Top 10 most significant GO categories", yaxis = list(title = "adjusted p-value"),
                          xaxis = form,
                          margin = list(b = 100))
  # fig2
  
})

output$GO_perc_term <- renderPlotly({
  gotab <- saved_GO()
  gotab <- gotab[gotab$significantly_expressed == input$UDA, ]
  gotab <- filterGO(gotab,input$gos_of_interest)
  gotab$percentage <- gotab$'%' *100
  #print(tabi$percentage)
  gotab <- as.data.table(gotab)
  gotab$percentage <- as.numeric(gotab$percentage)
  form <- list(categoryorder = "array",
               categoryarray = gotab$percentage,
               title = "Percentage")
  fig2 <- plot_ly(gotab, x =~percentage, y = ~p.adjust, color = ~p.adjust,colors = c("#ba3131","#fff68f"),name = "Scatterplot: percentage of genes in a GO Term and its adjusted p-value",mode = "markers", size = ~percentage, text = ~Term)
  fig2 <- fig2 %>% layout(title = "Scatterplot: percentage of genes in a GO Term and its adjusted p-value", yaxis = list(title = "adjusted p-value"),
                          xaxis = form,
                          margin = list(b = 100))
  # fig2
  
})

output$GO_padj_logFC <- renderPlotly({
  gotab <- saved_GO()
  gotab <- gotab[gotab$significantly_expressed == input$UDA, ]
  gotab <- filterGO(gotab,input$gos_of_interest)
  
  circ_data <- make_circ_data(gotab,saved_circ_mat(),saved_ontology())
  
  #gotab$percentage <- gotab$'%' *100
  #print(tabi$percentage)
  gotab <- as.data.table(circ_data)
  #gotab$percentage <- as.numeric(gotab$percentage)

  form <- list(categoryorder = "array",
               categoryarray = gotab$logFC,
               title = "logFC")
  fig2 <- plot_ly(gotab, x =~logFC, y = ~adj_pval, name = ~term,mode = "markers",color = ~term,customdata = ~term,size = ~logFC, text = ~genes, hovertemplate = paste('<br>gene: %{text}<br>', '<br>term: %{customdata}<br>'))
  fig2 <- fig2 %>% layout(title = "Scatterplot: logFC of genes in a GO Term and its adjusted p-value", yaxis = list(title = "adjusted p-value"),
                          xaxis = form,
                          margin = list(b = 100))
  # fig2
  
})

#observe({
#  if(!is.null(input$HoDGO)){
#    if(input$HoDGO == "disease - CAD"){
#      updateSelectInput(session, "TPMGO", choices = c("TPM > 0.2"))
#    }else{
#      updateSelectInput(session, "TPMGO", choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"))
#    }
#  }
#})

#TODO: noch einen Plot mit sortierung nach p.adjusted terms (dh die make_circ_data methode noch mal umschreiben)
make_circ_data <- function(erg_table, difex_tab,categ){

  erg_table$Genes <- str_replace_all(erg_table$Genes,"/", ",")
  erg_table <- separate_rows(erg_table,Genes, sep = ",", convert = FALSE)
  go_tab <- erg_table[,c("Category","Term","Genes","p.adjust")]

  colnames(go_tab) <- c("ID","term","genes","adj_pval")
  go_tab$category <- categ
 # go_tab <- go_tab[,c("category","ID","term","genes","adj_pval")]
  colnames(difex_tab) = c("ID","logFC")
  
  #difex_tab <- difex_tab[difex_tab$ID %in% go_tab$genes,]

  circ <- circle_dat(go_tab, difex_tab)
  circ$regulation <- ""
  circ$regulation[circ$logFC < 0] <- "downregulated"
  circ$regulation[circ$logFC > 0] <- "upregulated"
  circ <- circ %>% arrange(desc(abs(logFC)))

  return(circ)
}



output$GOC <- renderPlotly({
  gotab <- saved_GO()
  gotab <- gotab[gotab$significantly_expressed == input$UDA, ]
  gotab <- filterGO(gotab,input$gos_of_interest)
  circ_data <- make_circ_data(gotab,saved_circ_mat(),saved_ontology())
  chosen_terms <- circ_data[1:10,3] #top 10 highest logFC
  filtered_tab <- circ_data[circ_data$term %in% chosen_terms,]
  p <- plot_ly(filtered_tab,x = ~term, y = ~logFC, type = "violin",text = ~genes,color= ~regulation,hovertemplate = paste('%{x}', '<br>gene: %{text}<br>')) %>%
    add_markers(size = 5)
  p <- p %>% layout(title = "GO categories with the 10 highest DE genes")
  p
})

saved_GO_graph = reactiveVal()
output$GOZ <- renderPlot({
  obj <- NULL
  if(input$UDA == "upregulated"){
    obj <- saved_GO_obj_up()
  }else if(input$UDA == "downregulated"){
    obj <- saved_GO_obj_down()

  }
  
  if(!is.null(obj)){
    term_enriched <- any(obj@result[["p.adjust"]]<= 0.05)
    if(term_enriched){
      saved_GO_graph(obj)
      goplot(obj)
    }else{
      validate(
        need(term_enriched == TRUE, "0 enriched terms found (p-value < 0.05)")
      )
    }
    
  }
})

filterGO <- function(tabs, gos_of_interest){
  gos <- c("%platelet%", "%hemostasis%", "%coagulation%", "%blood%", "%thrombosis%", "%atherosclerosis%", "%coronary artery disease%", "%clotting%", "%fibrin clot formation%", "%prothrombin activation%", "%inflammatory response to injury%")
  res_tab <- NULL
  if(gos_of_interest){
    res_tab <- tabs[tabs$Term %like any% gos,]
  }else{
    res_tab <- tabs
  }
  return(res_tab)
}

saved_chord = reactiveVal()
output$circGOplot <- renderPlot({
  gotab <- saved_GO()
  gotab <- gotab[gotab$significantly_expressed == input$UDA, ]
  gotab <- filterGO(gotab,input$gos_of_interest)
  gotab <- gotab%>% arrange(p.adjust)
  
  top5 <- gotab[1:5,"Term"]
  
  circ_data <- make_circ_data(gotab,saved_circ_mat(),saved_ontology())

  fil_circ <- circ_data[circ_data$term %in% top5,]

  chord <- chord_dat(fil_circ, unique(fil_circ[, c("genes","logFC")]), unique(fil_circ$term))
  #View(chord)
  saved_chord(chord)

  GOChord(chord, space = 0.02, limit = c(3, 0), gene.order = 'logFC', gene.space = 0.25, gene.size = 3, ribbon.col = brewer.pal(n = 5, name = "Set2"), border.size = 0.01)
})
# output$GObar <- renderPlot({ #
#   tab <- saved_GO()
#   #todo anpasseds
#   panther_result_up_ref <- tab[tab$significantly_expressed == "upregulated",]
#   panther_result_down_ref <- tab[tab$significantly_expressed == "downregulated",]
#   
#   print("in output GO bar")
#   p <- myGObarplot_new(panther_result_up = panther_result_up_ref, 
#                     panther_result_down = panther_result_down_ref, 
#                     count_threshold = input$countGO, 
#                     fdr_threshold = input$fdrGO,
#                     padj_threshold = input$padjGO,
#                     gos_of_interest = c("%platelet%", "%hemostasis%", "%coagulation%"),
#                     background = "Reference genes")
#   p
# })


# Buttons -----------------------------------------------------------------





