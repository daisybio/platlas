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
    sig_CAD_rp_vs_mp <- read.csv2("www/cad_rp_vs_mp_sign_diff_mir_counts.CSV",header = TRUE, na.strings = "_")
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
nod = reactiveVal("")
edg = reactiveVal("")
minod = reactiveVal("")
mi_name = reactiveVal("")

getMIPaths <- eventReactive(input$clickMIF2,{
  m2Genes <- miRNA_path("", input$miEXP,"",input$UDMI, input$miDB,input$miF2,input$miSign,"0")
  #print(m2Genes)
  mi_name(input$miF2)
  ge2GO <- miRNA_path(input$miOnt, input$miEXP,"G",input$UDMI, input$miDB,input$miF2,input$miSign,"1")
  #print(ge2GO)
  GO_Nodes <- miRNA_path(input$miOnt, input$miEXP,"G",input$UDMI, input$miDB,input$miF2,input$miSign,"2")
  #print(GO_Nodes)
  GO_Edges <- miRNA_path(input$miOnt, input$miEXP,"G",input$UDMI, input$miDB,input$miF2,input$miSign,"3")
  #print(GO_Edges)
  mir_Nodes <- miRNA_path("", input$miEXP,"N",input$UDMI, input$miDB,input$miF2,input$miSign,"0")
  #print(mir_Nodes)
  #mG <- read.table(m2Genes,header = FALSE)
  m2g(m2Genes)
  gG <- read.table(ge2GO,quote = "",sep = ",", header = TRUE)
  mi_tab <- read.table(mir_Nodes, sep = "\t", header = FALSE)
  g2g(gG)
  nod(GO_Nodes)
  edg(GO_Edges)
  minod(mi_tab)
  finalpath <- m2Genes
  mG <- read.delim(m2Genes, header=FALSE, na.strings = NA,comment.char = "",quote = "\"", sep = "\t",dec = ".")
  if(input$miDB == "miRDB"){ 
    colnames(mG)=c("Target Gene (Ref-Seq ID)", "Target Gene (Ensembl ID)", "miRNA", "prediction score")
  }else if(input$miDB == "miRTarBase"){
    colnames(mG)=c("miRNA","miRTarBase ID", "Species (mRNA)", "Target Gene", "Target Gene (Entrez ID)", "Species (Target Gene)","Experiments","Support Type","References (PMID)")
  }
  end <- mG
})
output$MIDT <- renderDataTable(getMIPaths(),escape = FALSE, options = list(lengthMenu = c(10,15,20), pageLength = 15, width = 700))
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

  # p <- EnhancedVolcano(getVP(),
  #                      lab = getVP()[,1],
  #                      x = 'log2FoldChange',
  #                      y = 'padj',####gucken
  #                      xlim = c(-5, 8),
  #                      #title = 'N061011 versus N61311',
  #                      pCutoff = as.numeric(input$mpCO),
  #                      FCcutoff = as.numeric(input$mfCO),
  #                      xlab = bquote(~Log[2]~ 'fold change'),
  #                      #title = ' ',
  #                      subtitle = ' ',##,
  #                      col=c('lightgrey', 'lightgrey', 'lightgrey', 'red'),
  #                      #shape = c(0, 0, 0, 0),
  #                      colAlpha = 4/5,
  #                      legendPosition = 'right',
  #                      legendLabSize = 14,
  #                      legendIconSize = 4.0
  # )
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


saved_DTMIT = reactiveVal()
output$DTMIT <- renderDataTable({
  g2g <- g2g()[,c(1,2,3,4,5,6)]
  saved_DTMIT(g2g)
  g2g
  datatable(
    cbind(' ' = '&oplus;', g2g), escape = -2,
    options = list(
      columnDefs = list(
        list(visible = FALSE, targets = c(7)),
        list(orderable = FALSE, className = 'details-control', targets = 1)
      )
    ),
    callback = JS("
  table.column(1).nodes().to$().css({cursor: 'pointer'});
  var format = function(d) {
    return '<div style=\"background-color:#eee; padding: .5em;\">  ' +
            'Genes: ' + d[7] + '</div>';
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

output$MIpG<- renderPlotly({
  tabi <- g2g()
  tabi$percentage <- (tabi$Genes_in_list / tabi$Total_genes)*100
  #print(tabi$percentage)
  tabi <- tabi[order(tabi$percentage),]
  form <- list(categoryorder = "array",
               categoryarray = tabi$id,
               title = "category")
  fig <- plot_ly(tabi, x = ~id, y = ~percentage, type = 'bar', name = 'percentage of target genes in category',text = ~Functional_Category ,marker = list(color = '#ba3131'))
  #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))
  fig <- fig %>% layout(title = "percentage of target genes in GO categories", xaxis = form,
                        yaxis = list(title = "percentage of target genes in category"),
                        margin = list(b = 100)
                        # barmode = 'group'
  )
  # fig
  
})

output$MIpFE <-renderPlotly({
  # n <- fthreshold()
  testtab <- g2g()#f_foldenrich(saved_GO(),n)
  testtab <- testtab[order(testtab$Enrichment_FDR, decreasing = FALSE),]
  testtab <- testtab %>% mutate(Enrichment_FDR = -log10(Enrichment_FDR))
  form <- list(categoryorder = "array",
               categoryarray = testtab$id,
               title = "category")
  fig2 <- plot_ly(testtab, x =~Enrichment_FDR, y = ~id, type = 'bar', name = 'enrichment FDR in categories',text = ~Functional_Category, orientation = 'h' ,marker = list(color = '#ba3131'))
  #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))
  fig2 <- fig2 %>% layout(title = "enrichment FDR in categories", xaxis = list(title = "enrichment FDR"),
                          yaxis = form,
                          margin = list(b = 100))
  # fig2
  
})



minet <- eventReactive(input$md_click,{
  
})
ColourScale <- 'd3.scaleOrdinal()
            .domain(["mir", "Gene", "u", "d"])
           .range(["#0000ff","#a0b6fe","#008000","#d01212"]);'

output$MInet <- renderForceNetwork({
  minet()
  t <- minod()
  datasetname <- mi_name()
  miud_liste<-getDIFmi(datasetname,input$mdTPM,input$mdfc,input$mdpc)
  u <- miud_liste[[1]]
  d <- miud_liste[[2]]
  #View(u)
  #View(d)
  ensemb <- data.frame(t[,1],t[,3])
  hgnc <-data.frame(t[,1],t[,2])
  colnames(ensemb) <- c("from","to")
  colnames(hgnc) <- c("from", "to")
  edges <- ensemb
  nodes <- data.frame(name = unique(c(edges$from, edges$to)), stringsAsFactors = FALSE)
  nodes$id <- 0:(nrow(nodes) - 1)
  nodes$type <- "Gene"
  nodes$type[startsWith(nodes$name,"hsa-")] <- "mir"
  nodes$type[nodes[,1] %in% u[,1]] <- "u"
  nodes$type[nodes[,1] %in% d[,1]] <- "d"
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
               NodeID = "name", Group = "type",opacity = 0.9,fontSize = 11, zoom = TRUE,colourScale = JS(ColourScale))
})

output$MIGOnet <- renderForceNetwork({
  edgestab<-read.table(edg(),quote = "",sep = ",", header = TRUE)
  nodestab<-read.table(nod(), quote = "",comment.char = "", dec = ".",sep = ",", header = TRUE ,na.strings = "NA", encoding = "automatic")
  nodestab$tabid <- 0:(NROW(nodestab)-1)
  edgestab$from <- trimws(edgestab$from, which = c("both"))
  edgestab$to <- trimws(edgestab$to, which = c("both"))
  nodestab$id<-trimws(nodestab$id, which = c("both"))
  nodestab$type <- "term"
  #View(nodestab)
  edges <- edgestab %>%
    left_join(nodestab, by = c("from" = "id")) %>%
    select(-from) %>%
    rename(from = tabid) %>%
    left_join(nodestab, by = c("to" = "id")) %>%
    select(-to) %>%
    rename(to = tabid)
  forceNetwork(Links = edges, Nodes = nodestab[, c(1,2,5,7)],
               Source = "from", Target = "to",
               Value = "width", NodeID = "id", Group = "type",opacity = 0.9,fontSize = 15 ,zoom = TRUE)
  
})
