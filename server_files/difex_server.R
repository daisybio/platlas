library(plotly)
library(crosstalk)
# DIFEX <- eventReactive(c(input$click1, input$click2),ignoreInit = T,{ 
#  
# })
# output$DIFEXtab <- renderUI({
#   DIFEX()
# })


#___________matrix erstellen___________
clicked1 = reactiveVal(isolate(input$click1))
clicked2 = reactiveVal(isolate(input$click2))
Titl <- reactiveVal("")
nameHod <- reactiveVal("")
nameTPM1 <- reactiveVal("")
pcReactive <- reactiveVal("")
fcReactive <- reactiveVal("")
gene_name_column <- reactiveVal()

Mat <- eventReactive(c(input$click1,input$click2),{#ignoreInit = T,
  ma <- data.frame()
  if(clicked1() < input$click1){#####komisch
    print("in mat1")
    #Click1<- TRUE
    finalpath <-paths2frames("1",path,input$TPM1,input$HoD,".CSV")
    print(finalpath)
    matriks2 <- read.csv2(finalpath)
    ma <- matriks2
    Titl('Volcano Plot: MP vs. RP')
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


#___________heatmap___________ # ich glaub code-wiederholung
clicked1_H = reactiveVal(isolate(input$click1))
clicked2_H = reactiveVal(isolate(input$click2))
hm = reactiveVal()
sav_path = reactiveVal()
heat <-  eventReactive(c(input$click1,input$click2),{ 
  finalpath <- ""
  if(clicked1_H() < input$click1){
    finalpath <-paths2frames("1",path1,input$TPM1,input$HoD,".CSV")
    clicked1_H(input$click1)
    nameHod(as.character(input$HoD))
    nameTPM1(as.character(input$TPM1))
    pcReactive(as.numeric(input$pCO1))
    fcReactive(as.numeric(input$fcO1))
    sav_path(finalpath)
    
    
  }else if(clicked2_H() < input$click2){
    finalpath <-paths2frames("2",path1,input$TPM2,input$HoD02,".CSV")
    clicked2_H(input$click2)
    sav_path(finalpath)
    nameHod(as.character(input$HoD02))
    nameTPM1(as.character(input$TPM2))
    pcReactive(as.numeric(input$pCO2))
    fcReactive(as.numeric(input$fcO2))
    
    
  }
  if(finalpath != ""){
    datei <- read.csv2(finalpath)
  h <- data.frame(datei)
  h <- h[, -c(1,3)]
  names <-h[,1]
  h <- h[,-1]
  rownames(h) = make.names(names, unique=TRUE)
  neuer <- data.matrix(h)
  neuer<-log(neuer,2)
  neuer[neuer < 0] <- 0
  endtab<- neuer
  #print(neuer)
  #View(neuer)
  hm(endtab)
  }
  
})


# output$heatmap <- renderPlot({
# bg <- colorRampPalette(c('blue','green','yellow'))
# heatmap.2(heat(), main = paste("Heatmap (", as.character(nameTPM1()) ,as.character(nameHod()), ")"), trace = "none", margins = c(10,12),col = bg )
saved_heatmap_plot = reactiveVal()
output$heatmap <- renderPlot({
  heat()
  if(!is.null(hm())){
  DBO <- "www/DBO-samples.tsv"
  dbo <- read.table(DBO, header = FALSE, sep = "\t", quote = "")
  rownames(dbo) = dbo[,1]
  dbo <- dbo[,-1]
  colnames(dbo) = c("Condition", "Platelet type")
  pl_names = c("red","blue")
  names(pl_names) = c("RP", "MP")
  d_names = c("lightblue","orange","darkgreen")
  names(d_names) = c("healthy","MI", "stable CAD")
  ann_col = list(`Platelet type` = pl_names,Condition = d_names)
  p<- pheatmap(hm(), 
               color = inferno(10),
               annotation_colors = ann_col,
               annotation_col = dbo,#dbo %>% select('Platelet type' = dbo$`Platelet type`, 'Disease' = dbo$Disease),
               show_rownames = FALSE, show_colnames = TRUE, scale="row"
  )
  saved_heatmap_plot(p)
  p
  }
  # plot_ly(x= colnames(hm()),y = rownames(hm()),z = hm(), colors = colorRamp(c("blue","green" ,"yellow")), type = "heatmap")
  
})
saved_pcabi = reactiveVal()
output$pcabi<- renderPlot({
  heat()
  if(!is.null(sav_path())){
    CAD_MI <- read.csv2(sav_path())
  cmeta <- pickGeneAnnot("e",CAD_MI)
  dis_metatab <- metatab
  dis_metatab <-dis_metatab[row.names(dis_metatab) %in% colnames(cmeta),]
  p <- pca(cmeta, metadata = dis_metatab, removeVar = 0.1)
  b<- biplot(p,colby = 'platelet_condition', colkey = c('MP' = 'blue', 'RP' = 'red', "0" = 'black'),
             # colLegendTitle = 'types of platelets',
             # encircle config
             #encircle = TRUE,
             #encircleFill = TRUE,
             # ellipse = TRUE,
             # ellipseConf = 0.95,
             # ellipseFill = TRUE,
             # ellipseAlpha = 1/4,
             # ellipseLineSize = 1.0,
             # xlim = c(-125,125), ylim = c(-50, 80),
             shape = 'medical_condition', shapekey = c('stable CAD'=15, 'MI'=17, 'healthy' = 15 ),#'Grade 3'=8
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
output$pcah <- renderPlot({
  heat()
  if(nameTPM1()!=""){
    tp <- nameTPM1()
  nt <- switch(tp,"TPM > 0.1"= "TPM-0-1","TPM > 0.3"="TPM-0-3","TPM > 1"="TPM-1","TPM > 2"="TPM-2")
  all_dis <- c("www/heat_diseased_",nt,"_phenotype-disease-all.CSV")
  pall_dis <- paste(all_dis, collapse = "")
  CAD_MI <- read.csv2(pall_dis)
  dis_metatab <- metatab[metatab[,1]!= "healthy", ]
  cmeta <- pickGeneAnnot("e",CAD_MI)
  # dis_metatab <- metatab
  dis_metatab <-dis_metatab[row.names(dis_metatab) %in% colnames(cmeta),]
  # dis_metatab <- dis_metatab[]
  dis_metatab <- dis_metatab[,-7]
  p <- pca(cmeta, metadata = dis_metatab, removeVar = 0.1)
  b<-eigencorplot(p,metavars = c("stable_CAD","MI","MP","RP"),posColKey = 'top',main = 'PC1-27 correlations')
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
  d[,7] = as.numeric(d[,7])
  annot <- gene_name_column()
  p <- EnhancedVolcano(d,
                       lab = d[,annot],
                       x = 'log2FoldChange',
                       y = 'padj',
                       xlim = c(-5, 8),
                       title = Titl(),
                       subtitle = paste(nameHod(),"patients dataset considering", nameTPM1()),
                       # tryCatch(
                       pCutoff = pcReactive(),
                       #  error = function(e) {e$message <- paste("log10(p-value) too small, choose a smaller p-value than: ",pcReactive())}),
                       FCcutoff = fcReactive(),
                       # pointSize = 3.0,
                       # labSize = 3.0,
                       xlab = bquote(~Log[2]~ 'fold change'),
                       # title = ' ', 
                       #subtitle = ' ',##,
                       col=c('lightgrey', 'lightgrey', 'lightgrey', 'red'),
                       #shape = c(0, 0, 0, 0),
                       colAlpha = 4/5,
                       legendPosition = 'right',
                       legendLabSize = 14,
                       legendIconSize = 4.0,
  )#),
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
  genes <- Mat()[ ,1]
  fc <- Mat()[ ,3]
  pval<- Mat()[ ,7]
  biotype <-Mat()[ ,12]
  df <- data.frame(genes,fc,pval,biotype)
  significant <- subset(df, pval< pcReactive() & abs(fc)> fcReactive())
  #View(significant)
  biotypes <- significant %>% count(biotype)
  biotypes
  #View(biotypes)
  # protein_coding<- as.numeric(lengths(subset(significant,biotype == "protein_coding"))[1])
  # misc_RNA <- as.numeric(lengths(subset(significant,biotype == "misc_RNA"))[1])
  # Mt_rRNA<- as.numeric(lengths(subset(significant,biotype == "Mt_rRNA"))[1])
  # processed_pseudogene <- as.numeric(lengths(subset(significant,biotype == "processed_pseudogene"))[1])
  # Mt_tRNA<- as.numeric(lengths(subset(significant,biotype == "Mt_tRNA"))[1])
  # transcribed_unprocessed_pseudogene <- as.numeric(lengths(subset(significant,biotype == "transcribed_unprocessed_pseudogene"))[1])
  # biot2 <- c(protein_coding,misc_RNA,Mt_rRNA,processed_pseudogene,Mt_tRNA,transcribed_unprocessed_pseudogene)
  # 
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
  # tab1 %>% drop_na(padj)
  # tab1 <- Mat
  # )
  significant <- subset(tab1, tab1[ ,7]<=pcReactive() & abs(tab1[ ,3])>=fcReactive())
  #significant %>% drop_na()
  upreg <- subset(significant,significant[ ,3] > 0)
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

output$TopGenesTab <- renderDataTable(UD(),
                                      escape = FALSE, options = list(lengthMenu = c(5, 10, 15), pageLength = 15, width = 500)
)




