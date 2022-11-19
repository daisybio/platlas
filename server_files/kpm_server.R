KPM<- eventReactive(c(input$clickK1, input$clickK2),ignoreInit = T,{ 
  tabBox(
    width =12,
    height = "1000px",
    title = "Get the results",
    tabPanel(
      id= "tT", "Gene Ontologies (Top 30 most significant)",
      # uiOutput("DTO")
      #selectInput(inputId = "UpDown" , label = "Choose up or downregulated genes" , choices = c("upregulated","downregulated"), multiple = FALSE),
      #actionButton(inputId = "click4",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
      
      downloadButton("downloadKPMD", "Download table"),
      br(),
      br(),
      dataTableOutput(outputId = "KPMD"),
      tags$head(tags$style("#KPMD table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
      
      tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                             border-top: 1px solid #000000;}",media="screen", type="text/css"))
      
    ),
    tabPanel(
      id= "tH", 
      #uiOutput("HMO")
      
      title = "KPM Network", background = "red" , solidHeader = TRUE,
      #plotOutput(""),
      forceNetworkOutput("KPMnet", height = "850px")
      
      
      
    ),
    tabPanel(
      id= "tH", 
      #uiOutput("HMO")
      
      title = "GO Hierarchy Tree", background = "red" , solidHeader = TRUE,
      #plotOutput(""),
      span("Top 30 most significant GO-Terms",style = "color: black; font-size: 18px"),
      imageOutput("KPMH", height = "500px"),
      br(),
      downloadButton("dlKPMH", "Download plot in PNG format")
      
    ),
    tabPanel(
      id= "tN", 
      #uiOutput("HMO")
      
      title = "GO Network", background = "red" , solidHeader = TRUE,
      #plotOutput(""),
      forceNetworkOutput("KPMGOnet", height = "850px")
      #  br(),
      # downloadButton("dlKPMN", "Download plot in PNG format")
      
    ),
    tabPanel(
      title = "percentage of DE genes in GO categories", solidHeader = TRUE,# width = 12,
      span("Top 30 most significant GO-Terms",style = "color: black; font-size: 18px"),
      br(),
      plotlyOutput("KPMpG",height = "800px"),
      downloadButton("dlPNGGO", "Download plot in PNG format")
    ),
    tabPanel(
      title = "enrichment FDR per GO categories", solidHeader = TRUE,# width = 12,
      #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
      #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
      #br(),
      span("Top 30 most significant GO-Terms",style = "color: black; font-size: 18px"),
      br(), 
      plotlyOutput("KPMFE",height = "800px"),
      #downloadButton("dlPNGGO", "Download plot in PNG format")
    )
    
  )
})

output$KPMO <- renderUI({
  KPM()
})

# KPM ---------------------------------------------------------------------


clickedKPM1 = reactiveVal(isolate(input$clickK1))
clickedKPM2 = reactiveVal(isolate(input$clickK2))
h_p = reactiveVal("")
n1_p = reactiveVal("")
n2_p = reactiveVal("")
kp_p = reactiveVal("")
Kpm <- eventReactive(c(input$clickK1,input$clickK2),{
  ma <- data.frame()
  if(clickedKPM1() < input$clickK1){
    finalpath <-KPMpaths2frames("1",input$KO1,input$KTPM1,input$KD1,input$KUD1,".txt","1","")
    h_path <- KPMpaths2frames("1",input$KO1,input$KTPM1,input$KD1,input$KUD1,".png","2","")
    n1_path <- KPMpaths2frames("1",input$KO1,input$KTPM1,input$KD1,input$KUD1,".txt","4","")
    n2_path <- KPMpaths2frames("1",input$KO1,input$KTPM1,input$KD1,input$KUD1,".txt","5","")
    kp_path <- KPMpaths2frames("1","no",input$KTPM1,input$KD1,input$KUD1,".sif","6",input$KES1)
    print(finalpath)
    h_p(h_path)
    n1_p(n1_path)
    n2_p(n2_path)
    kp_p(kp_path)
    matriks2 <- read.table(finalpath,quote = "",sep = ",", header = TRUE)
    #View(matriks2)
    #matriks2 <- getGOs2tab(matriks2)
    
    if(input$KES1 == "ENSEMBL"){#"ENSEMBL", "gene symbol"
      #gen_names <- matriks2[,5]
      matriks2 <- getEorS("e",matriks2)
    }else if(input$KES1 == "gene symbol"){
      matriks2 <- getEorS("s",matriks2)
    }
    ma <- matriks2
    clickedKPM1(input$clickK1)
  }else if(clickedKPM2() < input$clickK2){
    finalpath <-KPMpaths2frames("2",input$KO2,input$KTPM2,input$KD2,input$KUD2,".txt","1","")
    h_path <- KPMpaths2frames("2",input$KO2,input$KTPM2,input$KD2,input$KUD2,".png","2","")
    n1_path <- KPMpaths2frames("2",input$KO2,input$KTPM2,input$KD2,input$KUD2,".txt","4","")
    n2_path <- KPMpaths2frames("2",input$KO2,input$KTPM2,input$KD2,input$KUD2,".txt","5","")
    kp_path <- KPMpaths2frames("2","no",input$KTPM2,input$KD2,input$KUD2,".sif","6",input$KES2)
    # print(finalpath)
    h_p(h_path)
    n1_p(n1_path)
    n2_p(n2_path)
    kp_p(kp_path)
    matriks3 <- read.table(finalpath,quote = "",sep = ",", header = TRUE)
    #matriks3 <- getGOs2tab(matriks3)
    #View(matriks3)
    if(input$KES2 == "ENSEMBL"){
      matriks3 <- getEorS("e",matriks3)
    }else if(input$KES2 == "gene symbol"){
      matriks3 <- getEorS("s",matriks3)
    }
    ma <- matriks3
    clickedKPM2(input$clickK2)
  }
  #View(ma)
  neu <- ma
})


output$KPMFE <- renderPlotly({
  # n <- fthreshold()
  testtab <- Kpm()#f_foldenrich(saved_GO(),n)
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

output$KPMpG <- renderPlotly({
  tabi <- Kpm()
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


output$KPMD <-renderDataTable({
  kp<-Kpm()
  datatable(
    cbind(' ' = '&oplus;', kp), escape = -2,
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
output$KPMnet <- renderForceNetwork({
  KPMtab<-read.table(kp_p(),quote = "",sep = " ", header = FALSE)
  neuKPM <- KPMtab[,c(1,3)]
  colnames(neuKPM)<- c("from","to")
  nodes <- data.frame(name = unique(c(neuKPM$from, neuKPM$to)), stringsAsFactors = FALSE)
  nodes$id <- 0:(nrow(nodes) - 1)
  nodes$type <- "type"
  KPMedges <- neuKPM %>%
    left_join(nodes, by = c("from" = "name")) %>%
    select(-from) %>%
    rename(from = id) %>%
    left_join(nodes, by = c("to" = "name")) %>%
    select(-to) %>%
    rename(to = id)
  forceNetwork(Links = KPMedges, Nodes = nodes,
               Source = "from", Target = "to",
               #Value = "", 
               NodeID = "name", Group = "type",opacity = 0.9,fontSize = 10 ,zoom = TRUE)
  
})

output$KPMH <- renderImage({
  filename <- normalizePath(h_p())
  print(filename)
  list(src = filename,
       width = "550px",
       height = "500px"
  )
}, deleteFile = FALSE)

output$KPMGOnet <- renderForceNetwork({
  edgestab<-read.table(n2_p(),quote = "",sep = ",", header = TRUE)
  #nodestab<-read.table(n1_p(), quote = "",comment.char = "", dec = ".",sep = ",", header = TRUE ,na.strings = "NA", encoding = "automatic")
  #nodestab$tabid <- 0:(NROW(nodestab)-1)
  #edgestab$from <- trimws(edgestab$from, which = c("both"))
  #edgestab$to <- trimws(edgestab$to, which = c("both"))
  #nodestab$id<-trimws(nodestab$id, which = c("both"))
  #nodestab$type <- "type"
  
  nodes <- data.frame(name = unique(c(as.character(edgestab$from), as.character(edgestab$to))), stringsAsFactors = FALSE)
  #View(nodes)
  nodes$id <- 0:(nrow(nodes) - 1)
  nodes$type <- "type"
  #View(nodes)
  nodes$name <- as.character(nodes$name)
  edges <- edgestab %>%
    left_join(nodes, by = c("from" = "name")) %>%
    select(-from) %>%
    rename(from = id) %>%
    left_join(nodes, by = c("to" = "name")) %>%
    select(-to) %>%
    rename(to = id)
  #View(edges)
  forceNetwork(Links = edges, Nodes = nodes,
               Source = "from", Target = "to",
               Value = "width", NodeID = "name", Group = "type",opacity = 0.9,fontSize = 10 ,zoom = TRUE)
  
})


