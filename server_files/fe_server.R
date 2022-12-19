# FER <- eventReactive(c(input$clickGO, input$clickGO2),ignoreInit = T,{
#   tabBox(
#     width = 12,
#     height = "1000px",
#     title = "Get the results",
#     tabPanel(
#       id="GOtab",title = "Top 500 GO Terms", solidHeader = TRUE,
#       "",
#       selectInput(inputId = "FeA" , label = "Choose gene annotation in GO categories" , choices = c("ENSEMBL ID","ENTREZ ID","HGNC symbol"), multiple = FALSE),
#       actionButton(inputId = "clickFeA",label = "Get GOs",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
#       br(), 
#       br(),
#       downloadButton("downloadGO", "Download table"),
#       br(),
#       br(),
#       dataTableOutput(outputId = "GOTab"),#width = 12,
#       
#       # .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th
#       tags$head(tags$style("#GOTab table {background-color: #DCDCDC; color : #000000}", media="screen", type="text/css")),
#       tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th
#                              {border-top: 1px solid #000000}"))
#     ),
#     tabPanel(
#       title = "percentage of DE genes in GO categories", solidHeader = TRUE,# width = 12,
#       br(),
#       plotlyOutput("GOpG",height = "800px"),
#       #downloadButton("dlPNGGO", "Download plot in PNG format")
#     ),
#     tabPanel(
#       title = "enrichment FDR per GO categories", solidHeader = TRUE,# width = 12,
#       #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
#       #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
#       #br(),
#       br(),
#       plotlyOutput("GOpFE",height = "800px"),
#       #downloadButton("dlPNGGO", "Download plot in PNG format")
#     ),
#     tabPanel(
#       title = "increased/decreased categories and their terms & their genes", solidHeader = TRUE,# width = 12,
#       span("Increased/decreased categories and their terms in relation to upregulated and downregulated genes in the corresponding category",style = "color: black; font-size: 19px"),
#       br(),
#       downloadButton("dlGOC", "Download plot in PNG format"),
#       br(),
#       plotOutput("GOC",height = "800px"),
#       #downloadButton("dlPNGGO", "Download plot in PNG format")
#     ),
#     tabPanel(
#       title = "GO categories and their z-score", solidHeader = TRUE,# width = 12,
#       span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
#       br(),
#       withMathJax(),
#       helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
#       br(),
#       br(),
#       br(),
#       plotlyOutput("GOZ",height = "800px"),
#       #downloadButton("dlPNGGO", "Download plot in PNG format")
#     ),
#     
#     tabPanel(
#       id= "tN", 
#       #uiOutput("HMO")
#       
#       title = "GO Network", background = "red" , solidHeader = TRUE,
#       #plotOutput(""),
#       forceNetworkOutput("GOnet", height = "850px")
#       #  br(),
#       # downloadButton("dlKPMN", "Download plot in PNG format")
#       
#     ),
#     tabPanel(
#       id= "tH", 
#       #uiOutput("HMO")
#       
#       title = "GO Hierarchy Tree", background = "red" , solidHeader = TRUE,
#       #plotOutput(""),
#       imageOutput("GOH", height = "850px"),
#       br(),
#       downloadButton("dlGOH", "Download plot in PNG format")
#       
#     )
#   )
# })
# 
# output$DTGO <- renderUI({
#   FER()
# })
library(clusterProfiler)

# GO ----------------------------------------------------------------------
#path3 <- "C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/GOE_"
#C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/
functional_enrichment <- function(significant_genes, universal_genes, ontology, annotation_type){
  res <- enrichGO(gene=significant_genes,universe = universal_genes,OrgDb="org.Hs.eg.db", ont=ontology,pAdjustMethod = "BH",keyType = annotation_type)#'ENSEMBL'
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

pathGO <- "www/"
#pGO2 <- "./www/"

clickedUDA1 = reactiveVal(isolate(input$clickGO))
clickedUDA2 = reactiveVal(isolate(input$clickGO2))


saved_GO = reactiveVal()
sav_cir =reactiveVal()
sav_go_h_p = reactiveVal()
sav_edges = reactiveVal()
GO <-  eventReactive(c(input$clickGO,input$clickGO2),{
  finalpath <- ""
  ontology <- ""
  UDA <- ""
  pval <- 0
  fc <- 0
  annot <- 0
  annotation_type <- ""
  if(clickedUDA1() < input$clickGO){
    finalpath <-paths2frames("1",path,input$TPMGO,input$HoDGO,".CSV")
    ontology <- switch(input$GOSS,"biological process"="BP","cellular component"= "CC","molecular function" = "MF")
    UDA <- switch(input$UDA, "upregulated" = "upregulated","downregulated" = "downregulated","all" = "all")
    annot <- switch(input$DVP1GO, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
    annotation_type <- switch(input$DVP1GO, "ENSEMBL ID"= "ENSEMBL", "HGNC symbol" = "SYMBOL")
    pval <- input$pCO1GO
    fc <- input$fcO1GO
    clickedUDA1(input$clickGO)
  }else if(clickedUDA2() < input$clickGO2){
    finalpath <-paths2frames("1",path,input$TPM2GO,input$HoD2GO,".CSV")
    ontology <- switch(input$GOSS2,"biological process"="BP","cellular component"= "CC","molecular function" = "MF")
    UDA <- switch(input$UDA2, "upregulated" = "upregulated","downregulated" = "downregulated","all" = "all")
    annot <- switch(input$DVP2GO, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
    annotation_type <- switch(input$DVP2GO, "ENSEMBL ID"= "ENSEMBL", "HGNC symbol" = "SYMBOL")
    pval <- input$pCO2GO
    fc <- input$fcO2GO
    clickedUDA2(input$clickGO2)
  }
  print(finalpath)
  if(finalpath != ""){
  matGO <- read.csv2(finalpath)
  input_list <- functional_enrichment_input(matGO, pval, fc, annot)
  sig_genes <- input_list[[UDA]]
  universe_genes <- input_list[["universe"]]
  functional_enrichment_result <- functional_enrichment(sig_genes, universe_genes, ontology, annotation_type)
  saved_GO(functional_enrichment_result)
  }
})

FEannot <- eventReactive(input$clickFeA,{
  res <- switch(input$FeA, "ENSEMBL ID" = 6, "ENTREZ ID" = 7, "HGNC symbol" = 8)
  res
})

output$GOTab <- renderDataTable({
  GO()
  gotab<- saved_GO()
  
  if(!is.null(gotab)){
  gotab<- as.data.table(gotab)
  View(gotab)
  gotab$Genes <- str_replace_all(gotab$Genes,"/", "; ")
  #annot <- FEannot()
  #gotab <- gotab[,c(1,2,3,4,5,annot,9)]
  datatable(
    cbind(' ' = '&oplus;', gotab), escape = -2,
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
  tabi$percentage <- tabi$'%' *100
  #print(tabi$percentage)
  tabi <- tabi[order(tabi$percentage),]
  form <- list(categoryorder = "array",
               categoryarray = tabi$id,
               title = "category")
  fig <- plot_ly(tabi, x = ~Category, y = ~percentage, type = 'bar', name = 'percentage of differentially expressed genes in category',text = ~Term ,marker = list(color = '#ba3131'))
  #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))
  fig <- fig %>% layout(title = "percentage of differentially expressed genes in GO categories", xaxis = form,
                        yaxis = list(title = "percentage of DE genes in category"),
                        margin = list(b = 100)
                        # barmode = 'group'
  )
  # fig
  
})

fthreshold <- eventReactive(input$clickGOFE,{
  num <- as.numeric(input$GOFE)
  end <- num
})
output$GOpFE <- renderPlotly({
  # n <- fthreshold()
  testtab <- saved_GO()#f_foldenrich(saved_GO(),n)
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

output$GOC <- renderPlot({
  GOCircle(sav_cir())
})


output$GOZ <- renderPlotly({
  zdat <- sav_cir()[,c(2,3,8)]
  zdat <- zdat[!(duplicated(zdat)),]
  zdat <- zdat[order(abs(zdat$zscore),decreasing = TRUE),]
  form <- list(categoryorder = "array",
               categoryarray = zdat$ID,
               title = "category")
  fig2 <- plot_ly(zdat, x =~zscore, y = ~ID, type = 'bar', name = 'z-score in categories',text = ~term, orientation = 'h' ,color = ~zscore,
                  colors = c("#0066b2","#ba3131"))#,marker = list(color = '#ba3131')
  #fig <- fig %>% add_trace(y = ~numInCat, name = 'number of genes in category', marker = list(color = '#ffe5ea'))
  fig2 <- fig2 %>% layout(title = "z-scores in categories", xaxis = list(title = "z-score"),
                          yaxis = form,
                          margin = list(b = 100))
})
output$GOnet <- renderForceNetwork({
  GOEtab<-sav_edges()#read.table(sav_edges(),quote = "",sep = " ", header = FALSE)
  GOE <- GOEtab[,c(1,2)]
  colnames(GOE)<- c("from","to")
  nodes <- data.frame(name = unique(c(GOE$from, GOE$to)), stringsAsFactors = FALSE)
  nodes$id <- 0:(nrow(nodes) - 1)
  nodes$type <- "type"
  GOE_edges <- GOE %>%
    left_join(nodes, by = c("from" = "name")) %>%
    select(-from) %>%
    rename(from = id) %>%
    left_join(nodes, by = c("to" = "name")) %>%
    select(-to) %>%
    rename(to = id)
  forceNetwork(Links = GOE_edges, Nodes = nodes,
               Source = "from", Target = "to",
               #Value = "", 
               NodeID = "name", Group = "type",opacity = 0.9,fontSize = 10 ,zoom = TRUE)
  
})

output$GOH <- renderImage({
  GO()
  
  #neu <- str_split(sav_go_h_p(), "www/")[[1]]
  
  filename <- normalizePath(sav_go_h_p())
  #print(filename)
  #tags$img(src = filename)
  #  print(filename)
  list(src = filename,
       width = "900px",
       height = "600px"
  )
}, deleteFile = FALSE)
# })

