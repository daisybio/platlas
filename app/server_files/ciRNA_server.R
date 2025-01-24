# CircRNA -----------------------------------------------------------------

complete_circRNA_df = reactiveVal()
upreg_circRNA = reactiveVal()
downreg_circRNA = reactiveVal()
ds_circRNA = reactiveVal()
cor_ciRNA_df = reactiveVal()
ciRNA_Network = reactiveVal()
sig_ciRNA = reactiveVal()
CRNA_dt <- eventReactive(input$clickCirc,{
  #print(input$ciDS)
  ciDS <- "MP vs. RP in CAD patients"
  enddf <- circRNA_path(ciDS)
  #View(enddf)
  enddf <- as.data.table(enddf)
  colnames(enddf)[1] = "X"
  ###View(enddf)
  ds_circRNA(ciDS)
  complete_circRNA_df(enddf)
  
  enddf <- enddf[(enddf$padj <= input$pCOcirc |(is.na(enddf$padj) & enddf$pvalue <= input$pCOcirc)),]
  up <- enddf[(enddf$log2FoldChange >= input$fCOcirc),]#& enddf$padj <= input$pCOcirc)
  down <- enddf[(enddf$log2FoldChange <= input$fCOcirc),]# & enddf$padj <= input$pCOcirc
  allci <- rbind(up,down)
  sig_ciRNA(allci)
  upreg_circRNA(up)
  downreg_circRNA(down)
  #circRNA_correlation <- read.csv2("www/circRNA_correlation.CSV")#circRNA_correlation.CSV
  #circRNA_genes <- read.table("www/circRNA_genes.tsv", sep = "\t", header = TRUE)#circRNA_genes.tsv
  #circ_RNA = merge(circRNA_correlation, circRNA_genes, by.y=c("ensembl_gene_id"), by.x=c("geneA"))
  #net <- readRDS("www/circRNA_network.rds")
  #ciRNA_Network(net)
  
  #cor_ciRNA_df(circ_RNA)
  ###View(cor_ciRNA_df())
})

ciRNA_nodes = reactiveVal()
ciRNA_up = reactiveVal()
ciRNA_down = reactiveVal()

cRNA_BP = reactiveVal()
cRNA_MF = reactiveVal()
cRNA_CC = reactiveVal()

CRNA_dt_2 <- eventReactive(input$clickCirc,{
  circRNA_correlation <- sig_ciRNA()#read.csv2("www/circRNA/circRNA_CCS.csv")#circRNA_correlation.CSV
  circRNA_correlation <- circRNA_correlation[,c("X","host_symbols")]

  circRNA_genes <- read.table("www/circRNA_genes_new.tsv", sep = "\t", header = TRUE)#circRNA_genes.tsv
  circ_RNA = merge(circRNA_correlation, circRNA_genes, by.y=c("hgnc_symbol"), by.x=c("host_symbols"))
  ##View(circ_RNA)
  #net <- readRDS("www/circRNA_network.rds")
  ciDS2 <- "MP vs. RP in CAD patients"
  expr_filepath <- switch(ciDS2,"MP vs. RP" = "www/all_healthy_TPM-0-1.CSV", "MP vs. RP in MI patients" = "www/all_diseased_TPM-0-1_disease-phenotype-all.CSV", "MP vs. RP in CAD patients" = "www/all_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV")
  matr <- read.csv2(expr_filepath)
  
  #full merge 
  circ_RNA_complete = merge(circ_RNA, matr, by.y=c("geneID"), by.x=c("ensembl_gene_id"))
  colnames(circ_RNA_complete)[3] = "circRNA"
  sig_circ_RNA_complete = circ_RNA_complete[((circ_RNA_complete$padj <= input$pCOcirc2 & abs(circ_RNA_complete$log2FoldChange) >= input$fCOcirc2) | (is.na(circ_RNA_complete$padj) & circ_RNA_complete$pvalue <= input$pCOcirc2 & abs(circ_RNA_complete$log2FoldChange) >= input$fCOcirc2)),]
  sig_circ_RNA_complete = unique(sig_circ_RNA_complete)
 
  
  #make nodes
  nodes = circ_RNA_complete[((circ_RNA_complete$padj <= input$pCOcirc2 & abs(circ_RNA_complete$log2FoldChange) >= input$fCOcirc2) | (is.na(circ_RNA_complete$padj) & circ_RNA_complete$pvalue <= input$pCOcirc2 & abs(circ_RNA_complete$log2FoldChange) >= input$fCOcirc2)),c("circRNA","ensembl_gene_id","host_symbols")]
  nodes = unique(nodes)
  ciRNA_nodes(nodes)
  ##View(nodes)
  ##View(sig_circ_RNA_complete)
  ##up und down
  up <- filter_genes_up_down(sig_circ_RNA_complete,input$pCOcirc2,input$fCOcirc2,input$ciAnnot,"upregulated")
  down <- filter_genes_up_down(sig_circ_RNA_complete,input$pCOcirc2,input$fCOcirc2,input$ciAnnot,"downregulated")
  ##View(down)
  all_cirna <- rbind(up,down)
  ##View(all_cirna)

  ciRNA_up(up)
  ciRNA_down(down)
  
  #functional enrichment
 # circ_RNA_genes <- unique(circ_RNA$geneA)
#  ##View(circ_RNA_genes)
#  filtered_matr <- matr[matr$geneID %in% circ_RNA_genes,]
#  ##View(filtered_matr)
  #rownames(matr)<- matr$geneID
  #filtered_m2 <- matr[circ_RNA_genes,]
  ###View(filtered_m2)
  #annotation <- switch(input$miAnnot, "ENSEMBL ID"= "geneID", "HGNC symbol" = "hgnc_symbol")
  #input_list <- functional_enrichment_input(filtered_matr,input$pCOcirc2,input$fCOcirc2,annotation)
  
#  sig_genes_up <- filtered_matr[filtered_matr$log2FoldChange > 0 & filtered_matr$padj < 0.05,annotation]#input_list[["upregulated"]]
#  sig_genes_down <- filtered_matr[filtered_matr$log2FoldChange < 0 & filtered_matr$padj < 0.05,annotation]
#  ##View(sig_genes_up)
#  ##View(sig_genes_down)
  
  
#  universe_genes <-filtered_matr[,annotation]
  #functional_enrichment(significant_genes_up,significant_genes_down,universal_genes, ontology, annotation_type)
#  annotation <- switch(input$miAnnot, "ENSEMBL ID"= "ENSEMBL", "HGNC symbol" = "SYMBOL")
  
#  circ_fe_BP <- functional_enrichment(sig_genes_up,sig_genes_down,universe_genes, "BP", annotation)#miRNA_path(input$miOnt, input$miEXP,"G",input$UDMI, input$miDB,input$miF2,input$miSign,"1")
#  circ_fe_MF <- functional_enrichment(sig_genes_up,sig_genes_down,universe_genes, "MF", annotation)
 # circ_fe_CC <- functional_enrichment(sig_genes_up,sig_genes_down,universe_genes, "CC", annotation)
  
#  ##View(circ_fe_BP$tab)
#  cRNA_BP(circ_fe_BP)
#  cRNA_MF(circ_fe_MF)
#  cRNA_CC(circ_fe_CC)
  
  
  #ciRNA_Network(net)
  cor_ciRNA_df(sig_circ_RNA_complete)
})


filter_genes_up_down <- function(dt,pa,fc,annotation,ud){
  annot <- switch(annotation,"ENSEMBL ID" = 1, "HGNC symbol" = 2)
  if(ud == "upregulated"){
    filtered <- dt[dt$padj <= pa & dt$log2FoldChange >= fc, ..annot]
  }else{
    filtered <- dt[dt$padj <= pa & dt$log2FoldChange <= -fc, ..annot]
  }
  return(filtered)
}


cRNA_ud_dt <- eventReactive(input$clickCUD,{
  use_dt <- NULL
  CRNA_dt()
  if(input$circUD == "upregulated"){
    use_dt = upreg_circRNA()
  }else if(input$circUD == "downregulated"){
    use_dt = downreg_circRNA()
  }
  use_dt<- use_dt %>% separate(colnames(use_dt)[1], c("location", "strand"), "_")%>% separate(location, c("chromosome", "coordinates"), ":")%>% separate(coordinates, c("start", "stop"), "-")
  use_dt <- use_dt[!(is.na(use_dt$strand)),]
  use_dt
  
})
saved_CIVP = reactiveVal()
output$CIVP <- renderPlot({
  CRNA_dt()
  ###View(upreg_circRNA())
  ###View(downreg_circRNA())
  upreg_circ <- upreg_circRNA()
  downreg_circ <- downreg_circRNA()
 # ##View(upreg_circ)
  top10_upreg_circ <-  upreg_circ[order(upreg_circ$padj, decreasing = FALSE),1] [1:10]
  ##View(upreg_circ)
  #print(top10_upreg_circ)#
  top10_downreg_circ <-  downreg_circ[order(downreg_circ$padj, decreasing = FALSE),1] [1:10]
  
  p <- myVolcanoPlot(complete_circRNA_df(),
                     data_type = "miRNA",
                     pCutoff = as.numeric(input$pCOcirc),
                     FCcutoff = as.numeric(input$fCOcirc),
                     xlim = c(-6, 6),
                     lab_selection = c(top10_upreg_circ,top10_downreg_circ),
                     symbol_type = "X", titel = paste0("Volcanoplot on the ", ds_circRNA(), " dataset"), subtitel = ""
  )
  # p <- EnhancedVolcano(complete_circRNA_df(),
  #                      lab = complete_circRNA_df()[,1],
  #                      x = 'log2FoldChange',
  #                      y = 'padj',
  #                      xlim = c(-5, 8),
  #                      #title = 'N061011 versus N61311',
  #                      pCutoff = as.numeric(input$pCOcirc),
  #                      FCcutoff = as.numeric(input$fCOcirc),
  #                      xlab = bquote(~Log[2]~ 'fold change'),
  #                      #title = ' ',
  #                      subtitle = ' ',##,
  #                      col=c('lightgrey', 'lightgrey', 'lightgrey', 'red'),
  #                      #shape = c(0, 0, 0, 0),
  #                      colAlpha = 4/5,
  #                      legendPosition = 'right',
  #                      legendLabSize = 14,
  #                      legendIconSize = 4.0)
  saved_CIVP(p)
  p
})


output$CIMA <- renderPlot({
  l <- CRNA_dt()
  if(!is.null(l)){
    d <- complete_circRNA_df()
   # ##View(d)
    d <- d[!(is.na(d$padj)),]
    d <- d %>% mutate(padj = str_replace(padj,"E", "e"))
    d <- d %>% mutate(padj = str_replace_all(padj,",", "."))
    d[,7] = as.numeric(d[,7])
    ##View(d)

   # annot <- "X"
    MA_title <- paste0("MA plot on the",ds_circRNA() , " dataset")
    p <- myMAplot_lim(d, fdr = as.numeric(input$pCOcirc), fc = as.numeric(input$fCOcirc), ylim = c(-4, 4), MA_title) #fdr ... + make ylim scalable
    p
  }
})



output$CTO <- renderDataTable(cRNA_ud_dt(),escape = FALSE, extensions = c("Buttons","Scroller"),options = list(lengthMenu = c(10, 15, 20), pageLength = 10, width = 700, scrollX = TRUE, dom = 'Bfrtip',
                                                                                                               buttons =list(list(
                                                                                                                 extend = 'collection',
                                                                                                                 buttons = c('csv', 'excel'),
                                                                                                                 text = 'Download'
                                                                                                               )))
)


output$ciNetwork <- renderVisNetwork({
  CRNA_dt_2()
  circ_RNA_obj <- ciRNA_Network()
  circ_RNA_obj %>% 
    visExport()
})


output$CircCor <- renderDataTable({
  CRNA_dt_2()
  cor_ciRNA_df()
},escape = FALSE,extensions = c("Buttons","Scroller"), options = list(lengthMenu = c(10, 15, 20), pageLength = 20, width = 700,scrollX = TRUE, dom = 'Bfrtip',bPaginate = FALSE,
                                 buttons =list(list(
                                   extend = 'collection',
                                   buttons = c('csv', 'excel'),
                                   text = 'Download'
                                 )))
)

output$CircBP <- renderPlotly({
  CRNA_dt_2()
  tab <-cor_ciRNA_df()
  ##View(tab)
  annot <- switch(input$ciAnnot,"ENSEMBL ID" = 1, "HGNC symbol" = 2)
  genes <- tab[ ,..annot]
  biotype <-tab[ ,"gene_biotype.x"]
  df <- data.frame(genes,biotype)
  colnames(df)[2]="biotype"
  ##View(df)
  #significant <- subset(df, pval< pcReactive() & abs(fc)> fcReactive())
  biotypes <- df %>% dplyr::count(biotype)
  p <- ggplot(biotypes, aes(x = n, y = biotype)) + geom_bar(stat = "identity", fill = "darkred") + xlab("Number of mappings based on correlation") + ylab("Biotypes") +ggtitle("Biotypes of genes, whose expression significantly correlates to circRNA expression")
  
  #saved_biotype_plot(p)
  p <- ggplotly(p)
  p
})


output$CInet <- renderForceNetwork({
  CRNA_dt_2()
  t <- ciRNA_nodes()
  
  if(!is.null(t)){
    u <- unlist(ciRNA_up())#miRNA_up()
    
    d <- unlist(ciRNA_down())#miRNA_down()
    ##View(d)
    ensemb <- data.frame(t[,1],t[,2], stringsAsFactors = FALSE)
    hgnc <-data.frame(t[,1],t[,3], stringsAsFactors = FALSE)
    colnames(ensemb) <- c("from","to")
    colnames(hgnc) <- c("from", "to")
    if(input$ciAnnot=="ENSEMBL ID"){
      edges <- ensemb
    }else{
      edges <- hgnc
    }
    nodes <- data.frame(name = unique(c(edges$from, edges$to)), stringsAsFactors = FALSE)
    ##View(t)
    #nur u & d in nodes
    nodes$id <- 0:(nrow(nodes) - 1)
    nodes$type <- "Gene"
    nodes$type[startsWith(nodes$name,"hsa")|startsWith(nodes$name,"chr")] <- "circ"
    nodes$type[nodes$name %in% u] <- "u"
    print(u)
    print(nodes$name %in% u)
    print(nodes$name)
    nodes$type[nodes$name %in% d] <- "d"
    ##View(nodes)
    miedges <- edges %>%
      left_join(nodes, by = c("from" = "name")) %>%
      select(-from) %>%
      rename(from = id) %>%
      left_join(nodes, by = c("to" = "name")) %>%
      select(-to) %>%
      rename(to = id)
    forceNetwork(Links = miedges, Nodes = nodes,
                 Source = "from", Target = "to",
                 #Value = "", 
                 NodeID = "name", Group = "type",opacity = 0.9,fontSize = 11, zoom = TRUE,colourScale = JS(ColourScale))}
})

output$Circ_mscore <- renderPlotly({
  CRNA_dt_2()
  corr_df <- cor_ciRNA_df()
  corr_df_filtered <- switch (input$circ_cor,
    "mscore" = corr_df[,c("geneA", "hgnc_symbol","geneB","mscor")],
    "partial correlation" = corr_df[,c("geneA", "hgnc_symbol","geneB","pcor")], 
    "correlation" = corr_df[,c("geneA", "hgnc_symbol","geneB","cor")]
  )
  corr_df_filtered$statistic <- switch (input$circ_cor,
                              "mscore" = corr_df$mscor,
                              "partial correlation" = corr_df$pcor, 
                              "correlation" = corr_df$cor
  )
  ###View(corr_df_filtered)
  #"ENSEMBL ID", "HGNC symbol"
  if(input$ciAnnot == "ENSEMBL ID"){
    corr_df_filtered$pairing <- paste0(corr_df_filtered$geneA," & ",corr_df_filtered$geneB)
    corr_df_filtered$geneID <- corr_df_filtered$geneA
    #corr_df_filtered <- corr_df_filtered[,-"geneA"]
  }else{
    corr_df_filtered$pairing <- paste0(corr_df_filtered$hgnc_symbol," & ",corr_df_filtered$geneB)
    corr_df_filtered$geneID <- corr_df_filtered$hgnc_symbol
    #corr_df_filtered <- corr_df_filtered[,-"hgnc_symbol"]
  }
  ###View(corr_df_filtered)
  corr_df_filtered <- corr_df_filtered[abs(corr_df_filtered$statistic) >= abs(input$circ_cutoff), ]
   form <- list(categoryorder = "array",
                categoryarray = corr_df_filtered$pairing,
                title = "Gene-ciRNA pairing")
   fig2 <- plot_ly(corr_df_filtered, x =~statistic, y = ~pairing, type = 'bar',  orientation = 'h' ,marker = list(color = '#ba3131'),customdata = ~geneB,text = ~geneID, hovertemplate = paste('<br>gene: %{text}<br>', '<br>ciRNA: %{customdata}<br>'))
   #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))name = paste0(input$circ_cor, " per pairing"),
   fig2 <- fig2 %>% layout(title = paste0(input$circ_cor, " per pairing"), xaxis = list(title = input$circ_cor),
                           yaxis = form,
                           margin = list(b = 100))
   fig2
  # 
})


output$Circ_gene <- renderPlotly({
  CRNA_dt_2()
  corr_df <- cor_ciRNA_df()
  if(input$ciAnnot == "ENSEMBL ID"){
    corr_df$geneID <- corr_df$ensembl_gene_id
  }else{
    corr_df$geneID <- corr_df$host_symbols
  }
  count_corr <- corr_df %>% dplyr::count(geneID)
  ###View(count_corr)
  count_corr <- count_corr[count_corr$n >= input$circ_cutoff_gene,]
  form <- list(categoryorder = "array",
               categoryarray = count_corr$geneID,
               title = "genes")
  fig2 <- plot_ly(count_corr, x =~geneID, y = ~n, type = 'bar',  marker = list(color = '#ba3131'),customdata = ~geneID,text = ~n, hovertemplate = paste('<br>number of circRNAs: %{text}<br>', '<br>gene ID: %{customdata}<br>'))
  #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))
  fig2 <- fig2 %>% layout(title = paste0("number of circRNAs per genes"), xaxis = form,
                          yaxis = list(title = "number of circRNAs"),
                          margin = list(b = 100))
  fig2
})


output$CircGO <- renderDataTable({
  CRNA_dt_2()
  #cor_ciRNA_df()c("Biological Process", "Cellular Component", "Molecular Function")
  GO_obj <- switch(input$circGO_ont, "Biological Process" = cRNA_BP(), "Cellular Component" = cRNA_CC(), "Molecular Function" = cRNA_MF())
  gotab <- GO_obj$tab[GO_obj$tab$p.adjust <= 0.05,]
  gotab <- as.data.table(gotab)
  
  gotab <- gotab[gotab$significantly_expressed == input$circGO_UDA, ]
  #gotab <- filterGO(gotab,input$gos_of_interest)
  
  gotab <- gotab[, -"significantly_expressed"]
  
  
  gotab$Genes <- str_replace_all(gotab$Genes,"/", "; ")
  #annot <- FEannot()
  #gotab <- gotab[,c(1,2,3,4,5,annot,9)]
  datatable(
    cbind(' ' = '&oplus;', gotab), escape = -2,
    options = list(
      columnDefs = list(
        list(visible = FALSE, targets = c(8)),
        list(orderable = FALSE, className = 'details-control', targets = 1)
      ),
      lengthMenu = c(5, 10, 15),
      width = 500,
      scrollX = TRUE
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
},
)







