# UIMI <- eventReactive(input$clickMIF,{
#   tagList( 
#     tabBox(
#       tabPanel(
#         id="t5","DE miRNA", 
#         span("Significant differentially expressed miRNA",style = "color: black; font-size: 16px" ),
#         br(),
#         downloadButton("dlDTMI", "Download table"),
#         br(),
#         dataTableOutput("DTMI"),
#         tags$head(tags$style("#DTMI table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
#         
#         tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
#                                border-top: 1px solid #000000;}",media="screen", type="text/css"))
#         
#       ),
#       tabPanel(
#         id="t6","Volcano Plot", 
#         
#         plotOutput("VPMI"),
#         downloadButton("dlVPMI", "Download plot in PNG format")
#         
#       ),
#       tabPanel(
#         id="t7","HeatMap", 
#         
#         plotOutput("HMMI"),
#         downloadButton("dlHMMI", "Download plot in PNG format")
#       )
#     )
#   )
# }
# )
# MITarg <- eventReactive(input$clickMIF2,{
#   tagList(
#     tabBox(
#       width = 12,
#       height = "1000px",
#       tabPanel(
#         id= "tE","miRNA Targets",
#         downloadButton("dlMIDT", "Download table"),
#         br(),br(),
#         dataTableOutput("MIDT"),
#         tags$head(tags$style("#DTMIT table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
# 
#         tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
#                                border-top: 1px solid #000000;}",media="screen", type="text/css"))
#       ),
#       tabPanel(
#         id="t10","miRNA Target Enrichment",
#         #selectInput(inputId = "MSIG" , label = "Choose significant or all mRNA" , choices = c("significant","all"), multiple = FALSE),
#         #actionButton(inputId = "clickSIG",label = "Get mRNA ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
#         span("The 30 Most significant GO terms" ,style = "color: black; font-size: 18px"),br(),
#         downloadButton("dlDTMIT", "Download table"),
#         br(),br(),
#         dataTableOutput("DTMIT"),
#         tags$head(tags$style("#DTMIT table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
# 
#         tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
#                                border-top: 1px solid #000000;}",media="screen", type="text/css"))
# 
#       ),
#       tabPanel(
#         id="t11","miRNA Target Network Graph",
#         span("Search if miRNA Targets are differentially expressed, by first defining differentially expressed genes with the filters below. The results show miRNAs (darkblue nodes) and their mapping genes (light blue nodes), while green nodes are upregulated genes and red nodes are downregulated genes.",
#              style = "color: black; font-size: 18px"),
#         br(),
#         selectInput(inputId = "mdTPM" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
#         numericInput(inputId = "mdpc",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
#         numericInput(inputId = "mdfc",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
#         actionButton(inputId = "md_click",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000"),
#         forceNetworkOutput("MInet", height = "850px")
# 
#       ),
#       tabPanel(
#         id="t12","miRNA Target Enrichment Graph",
#         span("miRNA Target Enrichment Graph",style = "color: black; font-size: 24px"),br(),br(),
#         forceNetworkOutput("MIGOnet", height = "850px")
# 
#       ),
#       tabPanel(
#         title = "percentage of target genes in GO categories", solidHeader = TRUE,# width = 12,
#         plotlyOutput("MIpG",height = "800px"),
#         #downloadButton("dlPNGGO", "Download plot in PNG format")
#       ),
#       tabPanel(
#         title = "enrichment FDR per GO categories", solidHeader = TRUE,# width = 12,
#         #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
#         #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
#         #br(),
#         plotlyOutput("MIpFE",height = "800px"),
#         #downloadButton("dlPNGGO", "Download plot in PNG format")
#       )
#     )
#   )
# })


# output$UIMI <-renderUI({
#   UIMI()
#   
# })   
# output$MITarg <-renderUI({
#   MITarg()
#   
# })   


# MiRNA Section -----------------------------------------------------------

getDTVP <- eventReactive(input$clickMIF,{
  df <- data.frame()
  mipa <- switch(input$miF,"Mature Platelets vs. Reticulated - all data" = "www/rp_vs_mp_diff_mir_counts.CSV", "in MI vs. CAD patients" = "www/disease_diff_mir_counts.CSV" , "only CAD patients"= "www/cad_rp_vs_mp_diff_mir_counts.CSV")
  df <- read.csv2(mipa,header = TRUE, na.strings = "_")
  pval <- input$mpCO
  fc <- input$mfCO
  df <- df[abs(df$log2FoldChange)>= fc & df$padj < pval, ]
  df <- df[!(is.na(df$baseMean)),]
  names2milinks(df)
})


getVP <- eventReactive(input$clickMIF,{
  df <- data.frame()
  if(input$miF == "Mature Platelets vs. Reticulated - all data"){
    #
    rp_vs_mp <-read.csv2("www/rp_vs_mp_diff_mir_counts.CSV",header = TRUE, na.strings = "_")
    df <- rp_vs_mp
  }else if(input$miF == "in MI vs. CAD patients"){
    dis_diff <-read.csv2("www/disease_diff_mir_counts.CSV",header = TRUE, na.strings = "_")
    df<- dis_diff
  }else if(input$miF == "only CAD patients"){
    CAD_rp_vs_mp <- read.csv2("www/cad_rp_vs_mp_diff_mir_counts.CSV",header = TRUE, na.strings = "_")
    df <- CAD_rp_vs_mp
  }
}
)
getHM <- eventReactive(input$clickMIF,{
  end <- matrix()
  df <- data.frame()
  cdf <- data.frame()
  if(input$miF == "Mature Platelets vs. Reticulated - all data"){
    sig_rp_vs_mp <-read.csv2("www/rp_vs_mp_sign_diff_mir_counts.CSV",header = TRUE, na.strings = "_")
    df <- sig_rp_vs_mp
    counts <- read.csv2("www/counts_normalized.CSV",header = TRUE, na.strings = "_")
    cdf <- counts
  }else if(input$miF == "in MI vs. CAD patients"){
    sig_dis_diff<-read.csv2("www/disease_sign_diff_mir_counts.CSV",header = TRUE, na.strings = "_")
    df<- sig_dis_diff
    counts_disease <- read.csv2("www/counts_normalized_disease.CSV",header = TRUE, na.strings = "_")
    cdf<- counts_disease
  }else if(input$miF == "only CAD patients"){
    sig_CAD_rp_vs_mp <- read.csv2("www/cad_rp_vs_mp_sign_diff_mir_counts.CSV",header = TRUE, na.strings = "_") # ccs_rp_vs_mp_diff_mir_counts_noSample6_filtered_pCutoff0.05.CSV bzw das signifikante
    df <- sig_CAD_rp_vs_mp
    counts_disease <- read.csv2("www/counts_normalized_disease.CSV",header = TRUE, na.strings = "_")
    cdf <- counts_disease
  }
  neu<- cdf[grep(paste(df[,1],collapse="|"),cdf[,1]),]
  names <-neu[,1]
  neu <-neu[,-1]
  rownames(neu) = make.names(names, unique=TRUE)
  #end<- neu
  merged_annotation <- fread("www/simulated_sample_names_all.tsv")
  #View(merged_annotation)
  merged_annotation <- as.data.frame(merged_annotation)
  rownames(merged_annotation) <- merged_annotation$original_sample_name
  to_translate <- colnames(neu)
  
  print(to_translate)
  new_names <- merged_annotation[to_translate,]
  #View(new_names)
  #print(to_translate)
  #print(new_names)
  #View(new_names)
  colnames(neu) <- new_names$sample_name
  #View(neuer)
  if(input$miF == "only CAD patients"){
    #threshold <- sort(rowVars(counts), decreasing = T)[150]
    #counts = counts[which(rowVars(counts) >= threshold),]
    neu <- neu[ , !(colnames(neu) %in% c("6_MPs","6_RPs"))]
  }
  neu
  # counts_disease <- data.matrix(counts_disease)
  #sample_annotation <- colnames(neu)
})

 m2g = reactiveVal("")
 g2g = reactiveVal("")
# nod = reactiveVal("")
# edg = reactiveVal("")
 minod = reactiveVal("")
 mi_name = reactiveVal("")
 miRNA_up = reactiveVal()
 miRNA_down = reactiveVal()
 miRNA_GO_obj = reactiveVal()
 miRNA_saved_circ = reactiveVal()
 miRNA_saved_ontology = reactiveVal()
 
getMIPaths <- eventReactive(input$clickMIF2,{
  m2Genes <- get_miRNA_datasets(input$miF2, input$miDB, input$mTpCO, input$mTfCO, input$miEXP, input$miSign)
  m2Genes_tab <- m2Genes[[input$UDMI]]
  mi_name(input$miF2)
  #View(m2Genes_tab)
  sav_circ <- getGeneMat(input$miF2,input$miAnnot)
  # if(input$UDMI == "upregulated"){
  #   sav_circ <- m2Genes$up_expr
  # }else if(input$UDMI == "downregulated"){
  #   sav_circ <- m2Genes$down_expr
  # }
  # annot <- switch(input$miAnnot, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
  # annot_colname <- colnames(sav_circ)[annot]
  if(input$UDMI == "upregulated"){
    sav_circ <- sav_circ[sav_circ$log2FoldChange > 0,]
  }else if(input$UDMI == "downregulated"){
    sav_circ <- sav_circ[sav_circ$log2FoldChange < 0,]
  }
  miRNA_saved_circ(sav_circ)
  #View(m2Genes)
  ontology <- switch(input$miOnt,"Biological Process"="BP","Cellular Component"= "CC","Molecular Function" = "MF")
  miRNA_saved_ontology(ontology)
  #View(m2Genes[["upregulated"]])
  #View(m2Genes[["downregulated"]])
  up <- extract_genes_from_m2g(m2Genes[["upregulated"]],input$miDB,input$miAnnot)
  down <- extract_genes_from_m2g(m2Genes[["downregulated"]],input$miDB,input$miAnnot)
  universe <- getUniverse(input$miF2,input$miAnnot)
  
  annotation <- switch(input$miAnnot, "ENSEMBL ID"= "ENSEMBL", "HGNC symbol" = "SYMBOL")
  ge2GO <- functional_enrichment(up,down,universe, ontology, annotation)#miRNA_path(input$miOnt, input$miEXP,"G",input$UDMI, input$miDB,input$miF2,input$miSign,"1")
  mir_Nodes <- makeNodes(m2Genes_tab,input$miDB)
  miRNA_up(up)
  miRNA_down(down)
  #print(mir_Nodes)
  #mG <- read.table(m2Genes,header = FALSE)
  m2g(m2Genes_tab)
 # gG <- read.table(ge2GO,quote = "",sep = ",", header = TRUE)
  #mi_tab <- read.table(mir_Nodes, sep = "\t", header = FALSE)
  #View(ge2GO$tab)
  g2g(ge2GO$tab)
  go_obj <- NULL
  if(input$UDMI== "upregulated"){
    go_obj <- ge2GO$obj_up
  }else if(input$UDMI== "downregulated"){
    go_obj <- ge2GO$obj_down
  }
  miRNA_GO_obj(go_obj)
  #nod(GO_Nodes)
  #edg(GO_Edges)
  minod(mir_Nodes)
  #finalpath <- m2Genes
  #mG <- read.delim(m2Genes, header=FALSE, na.strings = NA,comment.char = "",quote = "\"", sep = "\t",dec = ".")
  #if(input$miDB == "miRDB"){ 
  #  colnames(mG)=c("Target Gene (Ref-Seq ID)", "Target Gene (Ensembl ID)", "miRNA", "prediction score")
  #}else if(input$miDB == "miRTarBase"){
  #  colnames(mG)=c("miRNA","miRTarBase ID", "Species (mRNA)", "Target Gene", "Target Gene (Entrez ID)", "Species (Target Gene)","Experiments","Support Type","References (PMID)")
  #}
  #end <- 
  m2Genes_tab
})
output$MIDT <- renderDataTable(getMIPaths(),escape = FALSE, options = list(lengthMenu = c(10,15,20), pageLength = 10, scrollX = TRUE,width = 500))
output$DTMI <- renderDataTable(getDTVP(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 5,scrollX = TRUE,width = 200)
)
saved_VPMI = reactiveVal()
output$VPMI <- renderPlot({
  tab <- getVP()
 if(!is.null(tab)){ 
   tab <- tab[!(is.na(tab$padj)),]
   tab <- tab %>% mutate(padj = str_replace(padj,"E", "e"))
   tab <- tab %>% mutate(padj = str_replace_all(padj,",", "."))
   #View(d)
   
   tab[,7] = as.numeric(tab[,7])
   #annot <- gene_name_column()
   #order 
   upreg_genes <- tab[tab$log2FoldChange > as.numeric(input$mfCO) & tab$padj < as.numeric(input$mpCO), ]
   upreg_genes <- upreg_genes[upreg_genes[,1]!= "",]
   top10_upreg_genes <- upreg_genes[order(upreg_genes$padj, decreasing = FALSE),1] [1:10]
   
   downreg_genes <- tab[tab$log2FoldChange < -as.numeric(input$mfCO) & tab$padj < as.numeric(input$mpCO), ]
   downreg_genes <- downreg_genes[downreg_genes[,1]!= "",]
   top10_downreg_genes <- downreg_genes[order(downreg_genes$padj, decreasing = FALSE), 1][1:10]
   
  #View(getVP())
  p <- myVolcanoPlot(tab,
                data_type = "miRNA",
                pCutoff = as.numeric(input$mpCO),
                FCcutoff = as.numeric(input$mfCO),
                xlim = c(-6, 6),
                lab_selection = c(top10_upreg_genes,top10_downreg_genes),
                symbol_type = colnames(getVP())[1] , titel = "", subtitel = "")

  saved_VPMI(p)
  p}
})
saved_HMMI = reactiveVal()
output$HMMI <- renderPlot({
  counts <- getHM()
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  #View(counts)
  if(!is.null(counts)){
    if(input$miF == "only CAD patients"){
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
  dbo_file <- fread("www/simulated_sample_names_all.tsv")
  dbo <- data.frame("names" = dbo_file$sample_name, "Condition" = dbo_file$condition, "Platelet type" = dbo_file$platelet )
  rownames(dbo) = dbo$names
  dbo <- dbo[,-1]
  pl_names = c(cbbPalette[1],cbbPalette[4])#c(RPs = cbbPalette[1], MPs = cbbPalette[4])
  names(pl_names) = c("RP", "MP")
  d_names = c("lightblue","orange","darkgreen")
  names(d_names) = c("healthy","ACS", "CCS")
  #ann_col = list(`Platelet type` = pl_names,Condition = d_names)
  # colnames(dbo) = c("Disease", "Platelet type")
  # pl_names = c("red","blue")
  # names(pl_names) = c("RP", "MP")
  # d_names = c("orange","darkgreen")
  # names(d_names) = c("MI", "stable CAD")
  # 
  #View(counts)
  ann_col = list(`Platelet type` = pl_names,Disease = d_names)
  p <- pheatmap(counts,
                color = colorRampPalette(c("blue", "red"))(50),
                annotation_colors = ann_col,
                annotation_col = dbo,#dbo %>% select('Platelet type' = dbo$`Platelet type`, 'Disease' = dbo$Disease),
                show_rownames = FALSE, show_colnames = TRUE, scale="row")
  saved_HMMI(p)
  p
  }
  #saved_HMMI()
  
}
)

output$MAMI <- renderPlot({
  tab <- getVP()
  if(!is.null(tab)){
    tab <- tab[!(is.na(tab$padj)),]
    tab <- tab %>% mutate(padj = str_replace(padj,"E", "e"))
    tab <- tab %>% mutate(padj = str_replace_all(padj,",", "."))
    #View(d)
    
    tab[,7] = as.numeric(tab[,7])
    annot <- gene_name_column()
    p <- myMAplot_lim(tab, fdr = as.numeric(input$mpCO), fc = as.numeric(input$mfCO), ylim = c(-4, 4)) #fdr ... + make ylim scalable
    
    p
  }
  
})

#myPCAplot(vst(dds_ccs), intgroup = "RPs_MPs") + 
#geom_label_repel(aes(label=colnames(dds_ccs)), size = 2, show.legend = FALSE)
#

output$PCAMI <- renderPlot({
  tab <- getHM()
  if(!is.null(tab)){
    dis_metatab <- metatab
    dis_metatab <-dis_metatab[row.names(dis_metatab) %in% colnames(tab),]
    p <- pca(tab, metadata = dis_metatab, removeVar = 0.1)
    b<- biplot(p,colby = 'platelet', colkey = c('MP' = 'blue', 'RP' = 'red', "0" = 'black'),
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
               shape = 'condition', shapekey = c('CCS'=15, 'ACS'=17, 'healthy' = 15 ),#'Grade 3'=8
               hline = 0, vline = c(-25, 0, 25),
               legendPosition = 'right', legendLabSize = 16, legendIconSize = 8.0,
               title = "Principal component analysis bi-plot",
               
    )
    b
  }
  
})


output$ExpMI <- renderPlot({
  d <- getHM()
  if(!is.null(d)){
    #View(d)
    p <- myCountAnalysisPlot(as.data.frame(d), "miRNA", "test title")#fdr ... + make ylim scalable
    
    p
  }
})

# miRNA Targets & GO ------------------------------------------------------

saved_DTMIT = reactiveVal()
output$DTMIT <- renderDataTable({
  g2g <- as.data.table(g2g())#[,c(1,2,3,4,5,6)]
  g2g <- g2g[g2g$significantly_expressed == input$UDMI, ]
 # g2g <- filterGO(gotab,input$gos_of_interest)
  
  g2g <- g2g[, -"significantly_expressed"]
  #View(gotab)
  
  g2g$Genes <- str_replace_all(g2g$Genes,"/", "; ")
  datatable(
    cbind(' ' = '&oplus;', g2g), escape = -2,
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
  #  ,escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 15)
}


)

output$MIGOnet <- renderPlot({
  obj <- miRNA_GO_obj()
  if(!is.null(obj)){
    term_enriched <- any(obj@result[["p.adjust"]]<= 0.05)
    if(term_enriched){
      goplot(obj)
    }else{
      validate(
        need(term_enriched == TRUE, "0 enriched terms found (p-value < 0.05)")
      )
    }
  
  }
})


output$MIpG<- renderPlotly({
  tabi <- g2g()
  tabi <- tabi[tabi$significantly_expressed == input$UDMI, ]
 #tabi <- filterGO(tabi,input$gos_of_interest)
  tabi$percentage <- tabi$'%' *100
  #print(tabi$percentage)
  tabi <- as.data.table(tabi)
  tabi$percentage <- as.numeric(tabi$percentage)
  tabi <- tabi[order(tabi$percentage,decreasing = FALSE),]
  #View(tabi)
  # form <- list(categoryorder = "array",
  #             categoryarray = tabi$id,
  #             title = "category")
  fig <- plot_ly(tabi, x = ~Category, y = ~percentage, type = 'bar', name = 'percentage of differentially expressed genes in category',text = ~Term,marker = list(color = '#ba3131'))
  # barmode = 'group'
  #)
  # fig
  
})

output$MIpFE <-renderPlotly({
  # n <- fthreshold()
  testtab <- g2g()#f_foldenrich(saved_GO(),n)
  testtab <- testtab[testtab$significantly_expressed == input$UDMI, ]
  #testtab <- filterGO(testtab,input$gos_of_interest)
  testtab <- testtab[order(testtab$FDR, decreasing = FALSE),]
  testtab <- testtab %>% mutate(FDR = -log10(FDR))
  form <- list(categoryorder = "array",
               categoryarray = testtab$Category,
               title = "category")
  fig2 <- plot_ly(testtab, x =~FDR, y = ~Category, type = 'bar', name = 'enrichment FDR in categories',text = ~Term, orientation = 'h' ,marker = list(color = '#ba3131'))
  #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))
  fig2 <- fig2 %>% layout(title = "enrichment FDR in categories", xaxis = list(title = "enrichment FDR"),
                          yaxis = form,
                          margin = list(b = 100))
  # fig2
  
})

output$MIGO_perc_term <- renderPlotly({
  gotab <- g2g()
  gotab <- gotab[gotab$significantly_expressed == input$UDMI, ]
  #gotab <- filterGO(gotab,input$gos_of_interest)
  gotab$percentage <- gotab$'%' *100
  #print(tabi$percentage)
  gotab <- as.data.table(gotab)
  gotab$percentage <- as.numeric(gotab$percentage)
  form <- list(categoryorder = "array",
               categoryarray = gotab$percentage,
               title = "Percentage")
  fig2 <- plot_ly(gotab, x =~percentage, y = ~p.adjust, name = 'Top 10 categories with lowest adjusted p-values',mode = "markers", size = ~percentage, text = ~Term)
  fig2 <- fig2 %>% layout(title = "adjusted p-values in categories", yaxis = list(title = "adjusted p-value"),
                          xaxis = form,
                          margin = list(b = 100))
  # fig2
  
})
output$MIGO_padj_logFC <- renderPlotly({
  gotab <- g2g()
  gotab <- gotab[gotab$significantly_expressed == input$UDMI, ]
  #gotab <- filterGO(gotab,input$gos_of_interest)
  
  circ_data <- make_circ_data(gotab,miRNA_saved_circ(),miRNA_saved_ontology())
  
  #gotab$percentage <- gotab$'%' *100
  #print(tabi$percentage)
  gotab <- as.data.table(circ_data)
  #gotab$percentage <- as.numeric(gotab$percentage)
  View(gotab)
  form <- list(categoryorder = "array",
               categoryarray = gotab$logFC,
               title = "logFC")
  fig2 <- plot_ly(gotab, x =~logFC, y = ~adj_pval, name = ~term,mode = "markers",color = ~term,customdata = ~term,size = ~logFC, text = ~genes, hovertemplate = paste('<br>gene: %{text}<br>', '<br>term: %{customdata}<br>'))
  fig2 <- fig2 %>% layout(title = "adjusted p-values in categories", yaxis = list(title = "adjusted p-value"),
                          xaxis = form,
                          margin = list(b = 100))
  # fig2
  
})



output$MIGOC <- renderPlotly({
  gotab <- g2g()
  gotab <- gotab[gotab$significantly_expressed == input$UDMI, ]
  #gotab <- filterGO(gotab,input$gos_of_interest)
  circ_data <- make_circ_data(gotab,miRNA_saved_circ(),miRNA_saved_ontology())
  chosen_terms <- circ_data[1:10,3] #top 10 highest logFC
  filtered_tab <- circ_data[circ_data$term %in% chosen_terms,]
  p<- plot_ly(filtered_tab,x = ~term, y = ~logFC, type = "violin",text = ~genes,color= ~regulation,hovertemplate = paste('%{x}', '<br>gene: %{text}<br>')) %>%
    add_markers(size = 5)
  p
})



ColourScale <- 'd3.scaleOrdinal()
            .domain(["mir", "Gene", "u", "d"])
           .range(["#0000ff","#a0b6fe","#008000","#d01212"]);'

output$MInet <- renderForceNetwork({
  #minet()
  t <- minod()
  if(!is.null(t)){
  #datasetname <- mi_name()
  #miud_liste<-getDIFmi(datasetname,input$mdTPM,input$mdfc,input$mdpc)
  u <- miRNA_up()
  d <- miRNA_down()
  ensemb <- data.frame(t[,1],t[,3])
  hgnc <-data.frame(t[,1],t[,2])
  colnames(ensemb) <- c("from","to")
  colnames(hgnc) <- c("from", "to")
  if(input$miAnnot=="ENSEMBL ID"){
    edges <- ensemb
  }else{
    edges <- hgnc
  }
  nodes <- data.frame(name = unique(c(edges$from, edges$to)), stringsAsFactors = FALSE)
  nodes$id <- 0:(nrow(nodes) - 1)
  nodes$type <- "Gene"
  nodes$type[startsWith(nodes$name,"hsa-")] <- "mir"
  nodes$type[nodes[,1] %in% u] <- "u"
  nodes$type[nodes[,1] %in% d] <- "d"
  # View(nodes)
  miedges <- edges %>%
    left_join(nodes, by = c("from" = "name")) %>%
    select(-from) %>%
    rename(from = id) %>%
    left_join(nodes, by = c("to" = "name")) %>%
    select(-to) %>%
    rename(to = id)
  # View(miedges)
  forceNetwork(Links = miedges, Nodes = nodes,
               Source = "from", Target = "to",
               #Value = "", 
               NodeID = "name", Group = "type",opacity = 0.9,fontSize = 11, zoom = TRUE,colourScale = JS(ColourScale))}
})

# output$MIGOC <- renderPlotly({
#   gotab <- saved_GO()
#   gotab <- gotab[gotab$significantly_expressed == input$UDA, ]
#   gotab <- filterGO(gotab,input$gos_of_interest)
#   circ_data <- make_circ_data(gotab,saved_circ_mat(),saved_ontology())
#   chosen_terms <- circ_data[1:10,3] #top 10 highest logFC
#   filtered_tab <- circ_data[circ_data$term %in% chosen_terms,]
#   p<- plot_ly(filtered_tab,x = ~term, y = ~logFC, type = "violin",text = ~genes,color= ~regulation,hovertemplate = paste('%{x}', '<br>gene: %{text}<br>')) %>%
#     add_markers(size = 5)
#   p
# })

output$MIGOp <- renderPlotly({
  gotab <- g2g()
  gotab <- gotab[gotab$significantly_expressed == input$UDMI & gotab$p.adjust <= 0.05, ]
  #gotab <- filterGO(gotab,input$gos_of_interest)
  #circ_data <- make_circ_data(gotab,saved_circ_mat(),saved_ontology())
  #circ_data_sorted <- circ_data[,-4] %>% arrange(adj_pval)
  #View(gotab)
  gotab <- gotab%>% arrange(p.adjust)
  circ_end <- gotab[1:10,]
  form <- list(categoryorder = "array",
               categoryarray = circ_end$Category,
               title = "Category")
  fig2 <- plot_ly(circ_end, x =~Category, y = ~p.adjust, type = 'bar', name = 'Top 10 categories with lowest adjusted p-values',text = ~Term, marker = list(color = '#ba3131'))
  fig2 <- fig2 %>% layout(title = "adjusted p-values in categories", yaxis = list(title = "adjusted p-value"),
                          xaxis = form,
                          margin = list(b = 100))
  # fig2
  
})

miRNA_fp <- "/Users/leonoraka/Desktop/Projekte/platlas/www/miRNA/new/"
get_miRNA_datasets <- function(dataset, db, pC, fC, experimentally, significant){
  data_name <- switch(dataset, "Mature Platelets vs. Reticulated - all data" = "rp_mp","in diseased patients" = "dis","only CAD patients" = "ccs"  )
  database_name <- switch(db,"miRTarBase"= "tar", "miRDB" = "rdb")
  miRNA_targets_filepath <- paste0(miRNA_fp,data_name,"_",database_name,"_db.tsv")
  miRNA_expr_filepath <- switch(dataset,"Mature Platelets vs. Reticulated - all data" = "www/rp_vs_mp_diff_mir_counts.CSV", "in MI vs. CAD patients" = "www/disease_diff_mir_counts.CSV", "only CAD patients" = "www/ccs_rp_vs_mp_diff_mir_counts_noSample6_filtered_pCutoff0.05.CSV")
  miRNA_targets <- fread(miRNA_targets_filepath)
  miRNA_expr <- read.csv2(miRNA_expr_filepath,header = TRUE, na.strings = "_")
  target_colnames <- switch(db, "miRDB" = c("Target Gene (Ref-Seq ID)", "Target Gene (HGNC Symbol)", "miRNA", "prediction score","Target Gene (ENSEMBL ID)"),"miRTarBase" = c("miRNA","miRTarBase ID", "Species (mRNA)", "Target Gene", "Target Gene (Entrez ID)", "Species (Target Gene)","Experiments","Support Type","References (PMID)","exp_type", "Target Gene (Ensembl ID)"))
  colnames(miRNA_targets)= target_colnames
  if(experimentally && db == "miRTarBase"){
    miRNA_targets <- miRNA_targets[miRNA_targets$exp_type == "experimentally_retrieved",]
  }
  if(significant){
    miRNA_expr_filtered <- miRNA_expr[!is.na(miRNA_expr$padj) & miRNA_expr$padj <= pC & abs(miRNA_expr$log2FoldChange) >= fC,1]
    miRNA_targets <- miRNA_targets[miRNA_targets$miRNA %in% miRNA_expr_filtered,]
  }
  #up_miRNA <- miRNA_expr[miRNA_expr$log2FoldChange > 0 ,]
  #down_miRNA <- miRNA_expr[miRNA_expr$log2FoldChange < 0 ,]
  up_miRNA_expr <- miRNA_expr[miRNA_expr$log2FoldChange > 0 ,1]
  down_miRNA_expr <- miRNA_expr[miRNA_expr$log2FoldChange < 0 ,1]
  #universe <- getUniverse()
  #View(universe)
  up_miRNA_targets <- miRNA_targets[miRNA_targets$miRNA %in% up_miRNA_expr,]
  down_miRNA_targets <- miRNA_targets[miRNA_targets$miRNA %in% down_miRNA_expr,]
  
  
  return(list("upregulated" = up_miRNA_targets,"downregulated" = down_miRNA_targets))
}


makeNodes <- function(miRNATargets, db){
  #View(miRNATargets)
  nodes <- NULL
  if(db == "miRTarBase"){
    nodes <- miRNATargets[,c(1,4,11)]
  }else if(db == "miRDB"){
    nodes <- miRNATargets[,c(3,2,5)]
  }
  colnames(nodes) = c("miRNA", "hgnc_symbol", "ensembl_id")
  # nodes$regulation <- "no regulation"
  # if(annot == "ENSEMBL ID"){
  #     nodes$regulation[nodes$ensembl_id %in% up] <- "upregulated"
  #     nodes$regulation[nodes$ensembl_id %in% down] <- "downregulated"
  # }else{
  #     nodes$regulation[nodes$hgnc_symbol %in% up] <- "upregulated"
  #     nodes$regulation[nodes$hgnc_symbol %in% down] <- "downregulated"
  # }
  return(nodes)
}

library(purrr)
extract_genes_from_m2g <- function(miRNA_Targets, db, annot_type){
  extracted <- NULL
  if(db == "miRTarBase"){
    annotate <- switch(annot_type, "ENSEMBL ID" = 11,"HGNC symbol" = 4)
    extracted <- miRNA_Targets[,..annotate]
  }else if(db == "miRDB"){
    annotate <- switch(annot_type, "ENSEMBL ID" = 5,"HGNC symbol" = 2)
    extracted <- miRNA_Targets[,..annotate]
  }
  extracted <- unlist(extracted)
  extracted_new <- extracted[!is.na(extracted)] 
  return(extracted_new)
}

getUniverse <- function(input_datapath,annot_type){
  miRNA_expr_filepath <- switch(input_datapath,"Mature Platelets vs. Reticulated - all data" = "www/all_healthy_TPM-0-1.CSV", "in MI vs. CAD patients" = "www/all_diseased_TPM-0-1_disease-phenotype-all.CSV", "only CAD patients" = "www/all_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV")
  miRNA2 <-  switch(input_datapath,"Mature Platelets vs. Reticulated - all data" = "www/all_diseased_TPM-0-1_phenotype-disease-all.CSV", "in MI vs. CAD patients" = "", "only CAD patients" = "")
  extra_list <- c()
  gncol <- switch(annot_type, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
  if(miRNA2 != ""){
    matr2<- read.csv2(miRNA2)
    extra_list <- matr2[,gncol]
  }
  matr <- read.csv2(miRNA_expr_filepath)
  genes <- matr[,gncol]
  all_genes <- c(genes,extra_list)
  return(all_genes)
}

getGeneMat <- function(input_datapath,annot_type){
  miRNA_expr_filepath <- switch(input_datapath,"Mature Platelets vs. Reticulated - all data" = "www/all_healthy_TPM-0-1.CSV", "in MI vs. CAD patients" = "www/all_diseased_TPM-0-1_disease-phenotype-all.CSV", "only CAD patients" = "www/all_diseased_TPM-0-2_phenotype-disease-stable-CAD.CSV")
  miRNA2 <-  switch(input_datapath,"Mature Platelets vs. Reticulated - all data" = "www/all_diseased_TPM-0-1_phenotype-disease-all.CSV", "in MI vs. CAD patients" = "", "only CAD patients" = "")
  extra_df <- NULL
  gncol <- switch(annot_type, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
  if(miRNA2 != ""){
    matr2<- read.csv2(miRNA2)
    annot_colname <- colnames(matr2)[gncol]
    extra_df <- matr2[,c(annot_colname,"log2FoldChange")]
  }
  matr <- read.csv2(miRNA_expr_filepath)
  annot_colname <- colnames(matr)[gncol]
  genes <- matr[,c(annot_colname,"log2FoldChange")]
  all_genes <- rbind(genes,extra_df)
  return(all_genes)
}




