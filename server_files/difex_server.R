library(plotly)
library(DT)
library(matrixStats)

#___________matrix erstellen___________
clicked1 = reactiveVal(isolate(input$click1))
clicked2 = reactiveVal(isolate(input$click2))
Titl <- reactiveVal("")
nameHod <- reactiveVal("")
nameTPM1 <- reactiveVal("")
pcReactive <- reactiveVal("")
fcReactive <- reactiveVal("")
gene_name_column <- reactiveVal()
difex_title <- reactiveVal()

Mat <- eventReactive(c(input$click1,input$click2),{#ignoreInit = T,
  ma <- data.frame()
  if(clicked1() < input$click1){#####komisch
    print("in mat1")
    finalpath <-paths2frames("1",path,input$TPM1,input$HoD,".CSV")
    print(finalpath)
    matriks2 <- read.csv2(finalpath)
    ma <- matriks2
    Titl('Volcano Plot: MP vs. RP')
    difex_title("MP vs. RP")
    nameHod(as.character(input$HoD))
    nameTPM1(as.character(input$TPM1))
    pcReactive(as.numeric(input$pCO1))
    fcReactive(as.numeric(input$fcO1))
    gncol <- switch (input$DVP1, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
    gene_name_column(gncol)
    clicked1(input$click1)
  }else if(clicked2() < input$click2){
    finalpath <-paths2frames("2",path,input$TPM2,input$HoD02,".CSV")
    matriks3 <- read.csv2(finalpath)
    ma <- matriks3
    Titl('Volcano Plot: CCS vs. MI')
    difex_title("CCS vs. MI")
    nameHod(as.character(input$HoD02))
    nameTPM1(as.character(input$TPM2))
    pcReactive(as.numeric(input$pCO2))
    fcReactive(as.numeric(input$fcO2))
    gncol <- switch (input$DVP2, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
    gene_name_column(gncol)
    clicked2(input$click2)
  }
  if(nrow(ma)>0){
     neu <- ma[!(ma[,7]=="#NV"),]
     neu
  
  }
 
})

observe({
  if(!is.null(input$HoD)){
    if(input$HoD == "disease - CCS"){
      updateSelectInput(session, "TPM1", choices = c("TPM > 0.2"))
    }else{
      updateSelectInput(session, "TPM1", choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"))
    }
  }
})



#___________heatmap___________ # ich glaub code-wiederholung
clicked1_H = reactiveVal(isolate(input$click1))
clicked2_H = reactiveVal(isolate(input$click2))
hm = reactiveVal()
sav_path = reactiveVal()
gene_name_column_heat = reactiveVal()
heat <-  eventReactive(c(input$click1,input$click2),{ 
  finalpath <- ""
  if(clicked1_H() < input$click1){
    finalpath <-paths2frames("1",path1,input$TPM1,input$HoD,".CSV")
    clicked1_H(input$click1)
    nameHod(as.character(input$HoD))
    nameTPM1(as.character(input$TPM1))
    pcReactive(as.numeric(input$pCO1))
    fcReactive(as.numeric(input$fcO1))
    gncol <- switch (input$DVP1, "ENSEMBL ID"= 2, "HGNC symbol" = 1) #which to delete
    gene_name_column_heat(gncol)
    sav_path(finalpath) #warum beim zweiten nicht
    
    
  }else if(clicked2_H() < input$click2){
    finalpath <-paths2frames("2",path1,input$TPM2,input$HoD02,".CSV")
    clicked2_H(input$click2)
    sav_path(finalpath)
    nameHod(as.character(input$HoD02))
    nameTPM1(as.character(input$TPM2))
    pcReactive(as.numeric(input$pCO2))
    fcReactive(as.numeric(input$fcO2))
    gncol <- switch (input$DVP1, "ENSEMBL ID"= 2, "HGNC symbol" = 1) #which to delete
    gene_name_column_heat(gncol)
    
  }
  if(finalpath != ""){
    datei <- read.csv2(finalpath)
  h <- data.frame(datei)
  #View(h)
  h <- h[, -c(gene_name_column_heat(),3)]
  names <-h[,1]
  h <- h[,-1]
  rownames(h) = make.names(names, unique=TRUE)
  neuer <- data.matrix(h)
  neuer<-log(neuer,2)
  neuer[neuer < 0] <- 0
  merged_annotation <- fread("www/simulated_sample_names_all.tsv")
  merged_annotation <- as.data.frame(merged_annotation)
  rownames(merged_annotation) <- merged_annotation$original_sample_name
  to_translate <- colnames(neuer)
  new_names <- merged_annotation[to_translate,]
  colnames(neuer) <- new_names$sample_name
  if(nameHod()=="disease - CCS"){
    #threshold <- sort(rowVars(counts), decreasing = T)[150]
    #counts = counts[which(rowVars(counts) >= threshold),]
    neuer <- neuer[ , !(colnames(neuer) %in% c("6_MPs","6_RPs"))]
  }
  hm(neuer)
  }
  
})


# output$heatmap <- renderPlot({
# bg <- colorRampPalette(c('blue','green','yellow'))
# heatmap.2(heat(), main = paste("Heatmap (", as.character(nameTPM1()) ,as.character(nameHod()), ")"), trace = "none", margins = c(10,12),col = bg )
saved_heatmap_plot = reactiveVal()
output$heatmap <- renderPlot({
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  heat()
  if(!is.null(hm())){
  #DBO <- "www/complete_sample_names.tsv"
  counts <- hm()
  #View(counts)
  if(nameHod()=="disease - CCS"){
    threshold <- sort(rowVars(counts), decreasing = T)[150]
    counts = counts[which(rowVars(counts) >= threshold),]
    #counts <- counts[ , !(colnames(counts) %in% c("6_MPs","6_RPs"))]
    acs <- c("25_MPs","25_RPs","26_MPs","26_RPs","27_MPs","27_RPs","28_MPs","28_RPs")
    counts <- counts[ , !(colnames(counts) %in% acs)]
    # counts <- counts[, c("8_MPs","8_RPs","11_RPs","12_MPs","13_MPs","4_MPs","5_MPs","1_MPs","9_MPs","2_MPs","7_MPs","10_MPs","11_MPs","14_MPs","10_RPs",
    #                      "19_RPs","9_RPs","1_RPs","5_RPs","15_RPs","14_RPs","12_RPs","13_RPs","7_RPs","16_MPs","16_RPs","20_MPs","20_RPs","3_MPs",
    #                      "3_RPs","19_MPs","18_MPs","18_RPs","17_MPs","17_RPs","15_MPs","2_RPs","4_RPs")]
    #View(counts)
    
  }
  heatmap_title <- paste0("Heatmap: ", difex_title(), " in ", nameHod(), " dataset")
  dbo_file <- fread("www/simulated_sample_names_all.tsv")
  dbo <- data.frame("names" = dbo_file$sample_name, "Condition" = dbo_file$condition, "Platelet type" = dbo_file$platelet )
  rownames(dbo) = dbo$names
  dbo <- dbo[,-1]
  #colnames(dbo) = c("Condition", "Platelet type")
  pl_names = c(cbbPalette[1],cbbPalette[4])#c(RPs = cbbPalette[1], MPs = cbbPalette[4])
  names(pl_names) = c("RP", "MP")
  d_names = c("lightblue","orange","darkgreen")
  names(d_names) = c("healthy","ACS", "CCS")
  ann_col = list(`Platelet type` = pl_names,Condition = d_names)
  p<- pheatmap(counts, 
               color = colorRampPalette(c("blue", "red"))(50),
               border_color = NA,
               annotation_colors = ann_col,
               annotation_col = dbo,#dbo %>% select('Platelet type' = dbo$`Platelet type`, 'Disease' = dbo$Disease),
               show_rownames = FALSE, show_colnames = TRUE, scale="row",
               main = heatmap_title
  )
  saved_heatmap_plot(p)
  p
  }
})
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
saved_pcabi = reactiveVal()
output$pcabi<- renderPlot({
  heat()
  if(!is.null(sav_path())){
  cmeta <- hm()
  #View(cmeta)
  dis_metatab <- metatab
  dis_metatab <-dis_metatab[row.names(dis_metatab) %in% colnames(cmeta),]
  #View(dis_metatab)
  p <- pca(cmeta, metadata = dis_metatab, removeVar = 0.1)
  b<- biplot(p,colby = 'platelet', colkey = c('MP' = 'blue', 'RP' = 'red', "0" = 'black'),
             shape = 'condition', shapekey = c('CCS'=15, 'ACS'=17, 'healthy' = 15 ),#'Grade 3'=8
             hline = 0, vline = c(-25, 0, 25),
             legendPosition = 'right', legendLabSize = 16, legendIconSize = 8.0,
             title = "Principal component analysis bi-plot",
             
  )
  #View(cmeta)
  #View(b)
  saved_pcabi(b)
  
  b
  }
})
saved_pcah = reactiveVal()
output$pcah <- renderPlot({ #pse veq diseased
  heat()
  if(nameTPM1()!=""){
    tp <- nameTPM1()
 # nt <- switch(tp,"TPM > 0.1"= "TPM-0-1","TPM > 0.3"="TPM-0-3","TPM > 1"="TPM-1","TPM > 2"="TPM-2")
#  all_dis <- c("www/heat_diseased_",nt,"_phenotype-disease-all.CSV")
  #pall_dis <- paste(all_dis, collapse = "")
  #CAD_MI <- read.csv2(pall_dis)
    
  #dis_metatab <- metatab[metatab[,1]!= "healthy", ]
  #cmeta <- pickGeneAnnot("e",CAD_MI)
  cmeta <- hm()
  dis_metatab <- metatab
  dis_metatab <-dis_metatab[row.names(dis_metatab) %in% colnames(cmeta),]
  # dis_metatab <- dis_metatab[]
  #View(cmeta)
  #dis_metatab <- dis_metatab[,-7]
  #View(dis_metatab)
  dis_metatab <- dis_metatab[,-c(1:12)]
  p <- pca(cmeta, metadata = dis_metatab, removeVar = 0.1)
  #View(p)
  b <-eigencorplot(p,metavars = colnames(dis_metatab),posColKey = 'top',main = 'PC1-27 correlations')
  #b <- ggplotly(b)
  saved_pcah(b)
  #View(b)
  b
  }
})

# }
# )
# VPannot <- eventReactive(c(input$clickDVP1,input$clickDVP2),{
#   res <- switch (input$DVP, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
#   res
# })

# ___________Volcano Plot___________
savedVP_plot = reactiveVal()
output$VP <- renderPlot({
  d <- Mat()
  if(!is.null(d)){
  d <- d[!(is.na(d$padj)),]
  d <- d %>% mutate(padj = str_replace(padj,"E", "e"))
  d <- d %>% mutate(padj = str_replace_all(padj,",", "."))
  #View(d)
  
  d[,7] = as.numeric(d[,7])
  annot <- gene_name_column()
  #order 
  upreg_genes <- d[d$log2FoldChange > fcReactive() & d$padj < pcReactive(), ]
  upreg_genes <- upreg_genes[upreg_genes[,annot]!= "",]
  top10_upreg_genes <- upreg_genes[order(upreg_genes$padj, decreasing = FALSE),annot] [1:10]
  
  downreg_genes <- d[d$log2FoldChange < -fcReactive() & d$padj < pcReactive(), ]
  downreg_genes <- downreg_genes[downreg_genes[,annot]!= "",]
  top10_downreg_genes <- downreg_genes[order(downreg_genes$padj, decreasing = FALSE), annot][1:10]
  
  #print(top10_upreg_genes)
  #print(top10_downreg_genes)
  #upreg gene labels
  #downreg gene labels
  #print(paste0("lab_variable: ", annot))
  p <-myVolcanoPlot(d,
                    data_type = "gene",
                    pCutoff = pcReactive(),
                    FCcutoff = fcReactive(),
                    xlim = c(-6, 6),
                    lab_selection = c(top10_upreg_genes,top10_downreg_genes),
                    symbol_type = colnames(d)[annot] , titel = Titl(), subtitel = paste(nameHod(),"patients dataset considering", nameTPM1())
                    )
  # p <- EnhancedVolcano(d,
  #                      lab = d[,annot],
  #                      x = 'log2FoldChange',
  #                      y = 'padj',
  #                      xlim = c(-5, 8),
  #                      title = Titl(),
  #                      subtitle = paste(nameHod(),"patients dataset considering", nameTPM1()),
  #                      # tryCatch(
  #                      pCutoff = pcReactive(),
  #                      #  error = function(e) {e$message <- paste("log10(p-value) too small, choose a smaller p-value than: ",pcReactive())}),
  #                      FCcutoff = fcReactive(),
  #                      # pointSize = 3.0,
  #                      # labSize = 3.0,
  #                      xlab = bquote(~Log[2]~ 'fold change'),
  #                      # title = ' ',
  #                      #subtitle = ' ',##,
  #                      col=c('lightgrey', 'lightgrey', 'lightgrey', 'red'),
  #                      #shape = c(0, 0, 0, 0),
  #                      colAlpha = 4/5,
  #                      legendPosition = 'right',
  #                      legendLabSize = 14,
  #                      legendIconSize = 4.0,
  # )#),
  savedVP_plot(p)
  # p<- p +
  #   ggplot2::coord_cartesian(xlim=c(-6, 9)) +
  #   ggplot2::scale_x_continuous(
  #     breaks=seq(-6,9, 1))
  # p<- ggplotly(p +aes(x= log2FoldChange, y= -log10(padj)) )#+ aes(x= log2FoldChange, y= -log10(pvalue)
   p
  }
})


#___________biotypesPlot___________

matt<-eventReactive(c(input$click1,input$click2),{
  print("in matt ")
  if(!is.null(Mat())){
  #View(Mat())
  tab <- Mat()
  if(gene_name_column() == 11){
    #print("is hgnc")
    tab <- tab[!is.na(tab$hgnc_symbol),]
  }
  genes <- tab[ ,gene_name_column()]
  fc <- tab[ ,3]
  pval<- tab[ ,7]
  biotype <-tab[ ,12]
  df <- data.frame(genes,fc,pval,biotype)
  significant <- subset(df, pval< pcReactive() & abs(fc)> fcReactive())
  #View(significant)
  biotypes <- significant %>% dplyr::count(biotype)
  biotypes
  }
})
saved_biotype_plot = reactiveVal()
output$biotype1<- renderPlotly({
  d<- matt()
  if(!is.null(d)){
  #   View(d)
  #   p <- barplot(d,
  #              main = paste("Biotypes of the significant Genes",nameTPM1()),
  #              xlab = "Number of Genes",
  #              ylab = "Biotypes",
  #              names.arg = c("protein coding", "misc_RNA","Mt_rRNA", "proc. pseudogene", "Mt_tRNA", "transcr_unproc_pseudog"),
  #              col = "darkred",
  #              horiz = FALSE,
  #              cex.names=0.8,
  #              las=1
  # )
    p <- ggplot(d, aes(x = n, y = biotype)) + geom_bar(stat = "identity", fill = "darkred") + xlab("Number of Genes") + ylab("Biotypes") +ggtitle("Biotypes of significant Genes")
  
  saved_biotype_plot(p)
  p <- ggplotly(p)
  p
  }
})



#___________tabelle___________
saved_ud_tab =reactiveVal()
UD<-eventReactive(c(input$click1,input$click2,input$click4),ignoreInit = T,{
  # genes <- Mat()[ ,1]
  # fc <- Mat()[ ,3]
  # pval<- Mat()[ ,7]
  # bt <- Mat()[ ,12]
  # df <- data.frame(genes,fc,pval,bt)
  tab1 <- Mat()[,-10]
  
  ## tab1 %>% drop_na(padj)
  # tab1 <- Mat
  # )
  tab1 <- tab1[!(is.na(tab1$padj)),]
  tab1 <- tab1 %>% mutate(padj = str_replace(padj,"E", "e"))
  tab1 <- tab1 %>% mutate(padj = str_replace_all(padj,",", "."))
  tab1[,7] = as.numeric(tab1[,7])
  #View(tab1)
  
  significant <- subset(tab1, tab1[ ,7]<=pcReactive() & abs(tab1[ ,3])>=fcReactive())
  #significant %>% drop_na()
  upreg <- subset(significant,significant[ ,3] > 0)
  #View(upreg)
  #upreg %>% drop_na()
  ordUp <- upreg[order(+upreg[,3]), ]
  
  #ordUp %>% drop_na()
  endtab <- data.frame()
  downreg <- subset(significant,significant[ ,3] < 0)
  #downreg %>% drop_na()
  ordDown<- downreg[order(+downreg[,3]), ]
  if(input$UpDown == "upregulated"){
    endtab <- ordUp
    print(NROW(endtab))
  }else if(input$UpDown == "downregulated"){
    endtab <- ordDown
    print(NROW(endtab))
  }
  #  gen <- endtab
  saved_ud_tab(endtab)
  neu <- names2links(endtab)
  
})

output$TopGenesTab <- DT::renderDataTable(UD(),
                                      escape = FALSE, options = list(lengthMenu = c(5, 10, 15), pageLength = 15, width = 500)
)


##MA plot
output$DMA_plot <- renderPlot({
  d <- Mat()
  if(!is.null(d)){
    d <- d[!(is.na(d$padj)),]
    d <- d %>% mutate(padj = str_replace(padj,"E", "e"))
    d <- d %>% mutate(padj = str_replace_all(padj,",", "."))
    #View(d)
    
    d[,7] = as.numeric(d[,7])
    annot <- gene_name_column()
    MA_title <- paste0("MA plot: ", difex_title(), " in ", nameHod(), " dataset")
    p <- myMAplot_lim(d, fdr = pcReactive(), fc = fcReactive(), ylim = c(-4, 4), MA_title) #fdr ... + make ylim scalable

    p
  }
})

output$Dbiotype_plot <- renderPlot({ ##passt nicht ganz mit barplot
  heat()
  d <- hm()
  if(!is.null(d)){
    #View(d)
    if(gene_name_column_heat() == 2){
      id_type = "ensembl_gene_id"
    }else{
      id_type = "hgnc_symbol"
    }
    Bt_title <- paste0("Biotype plot: ", difex_title(), " in ", nameHod(), " dataset")
    p<- myBiotypeBarplot(as.data.frame(d), Bt_title, id_type = id_type)#hier noch namenparameter
    
    #p$biotype_plot
    p$biotype_plot
  }
})

output$Dexpression_plot <- renderPlot({
  heat()
  d <- hm()
  if(!is.null(d)){
    if(nameHod()=="disease - CCS"){
      EA_title <-  paste0("Expression analysis plot: ", difex_title(), " in ", nameHod(), " dataset before filtering")
    }
    EA_title <- paste0("Expression analysis plot: ", difex_title(), " in ", nameHod(), " dataset after filtering")
    #View(d)
    p <- myCountAnalysisPlot(as.data.frame(d), "gene", EA_title)#fdr ... + make ylim scalable
    
    p
  }
})

