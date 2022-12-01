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

# GO ----------------------------------------------------------------------
#path3 <- "C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/GOE_"
#C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/
pathGO <- "www/"
#pGO2 <- "./www/"

clickedUDA1 = reactiveVal(isolate(input$clickGO))
clickedUDA2 = reactiveVal(isolate(input$clickGO2))


saved_GO = reactiveVal()
sav_cir =reactiveVal()
sav_go_h_p = reactiveVal()
sav_edges = reactiveVal()
GO <-  eventReactive(c(input$clickGO,input$clickGO2),{
  # midtab <- data.frame()
  #udi <- UDA()
  #p <- udi[[1]]
  #c <- udi[[2]]
  if(clickedUDA1() < input$clickGO){
    fe_result<- fe_shiny_paths("M",input$TPMGO,input$HoDGO,input$GOSS,input$UDA)
    entab_file <- paste(pathGO,fe_result[[1]], sep= "")
    circ_file <- paste(pathGO,fe_result[[2]], sep= "")
    edg_file <- paste(pathGO,fe_result[[3]], sep= "")
    im_file <- paste(pathGO,fe_result[[4]], sep= "")
    entab <- read.table(entab_file, header = TRUE, sep = ",", quote = "")
    cir <- read.table(circ_file, sep = "\t", header= TRUE, quote = "")
    edge <- read.table(edg_file, header = TRUE, sep = ",", quote = "")
    sav_go_h_p(im_file)
    sav_edges(edge)
    on_fi <- switch(input$GOSS,"biological process"="BP","cellular component"= "CC","molecular function" = "MF")
    fi_on_entab <- entab[grep(on_fi,entab[,9]),c(1,2,3,4,5,6,7,8,9)] ##filter für GenNamen
    saved_GO(fi_on_entab)
    fi_cir <- cir[cir[,1]==on_fi, ]
    sav_cir(fi_cir)
    clickedUDA1(input$clickGO)
  }else if(clickedUDA2() < input$clickGO2){
    fe_result<- fe_shiny_paths("C",input$TPM2GO,input$HoD2GO,input$GOSS2,input$UDA2)
    entab_file <- paste(pathGO,fe_result[[1]], sep= "")
    circ_file <- paste(pathGO,fe_result[[2]], sep= "")
    edg_file <- paste(pathGO,fe_result[[3]], sep= "")
    im_file <- paste(pathGO,fe_result[[4]], sep= "")
    entab <- read.table(entab_file, header = TRUE, sep = ",", quote = "")
    cir <- read.table(circ_file, sep = "\t", header= TRUE, quote = "")
    edge <- read.table(edg_file, header = TRUE, sep = ",", quote = "")
    sav_go_h_p(im_file)
    sav_edges(edge)
    on_fi <- switch(input$GOSS2,"biological process"="BP","cellular component"= "CC","molecular function" = "MF")
    #
    fi_on_entab <- entab[grep(on_fi,entab[,9]),c(1,2,3,4,5,6,7,8,9)] ##filter fuer GenNamen
    saved_GO(fi_on_entab)
    fi_cir <- cir[cir[,1]==on_fi, ]
    sav_cir(fi_cir)
    clickedUDA2(input$clickGO2)
  }
  #saved_GO(midtab)
  #endtab<-names2GOlinks(midtab)
  saved_GO()
})

FEannot <- eventReactive(input$clickFeA,{
  res <- switch(input$FeA, "ENSEMBL ID" = 6, "ENTREZ ID" = 7, "HGNC symbol" = 8)
  res
})

output$GOTab <- renderDataTable({
  gotab<- GO()
  annot <- FEannot()
  gotab <- gotab[,c(1,2,3,4,5,annot,9)]
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
    ))}#escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 15, width = 700)
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
  tabi$percentage <- (tabi$Genes_in_list / tabi$Total_genes)*100
  #print(tabi$percentage)
  tabi <- tabi[order(tabi$percentage),]
  form <- list(categoryorder = "array",
               categoryarray = tabi$id,
               title = "category")
  fig <- plot_ly(tabi, x = ~id, y = ~percentage, type = 'bar', name = 'percentage of differentially expressed genes in category',text = ~Functional_Category ,marker = list(color = '#ba3131'))
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

