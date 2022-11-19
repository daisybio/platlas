

IN1 <- eventReactive(input$Iclick1,{
  tabBox(
    width = 12,
    height = "600px",
    title = "Get the results",
    tabPanel(
      
      id="","List of Genes", 
      #"Volcano plot",
      #uiOutput("VPO")
      #actionButton(inputId = "Iclick1",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000"),
      
      dataTableOutput(outputId = "IListOfGenes"),
      tags$head(tags$style("#TopGenesTab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
      
      tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
      
    ),
    tabPanel(
      id="tV","Pie Plot", 
      #"Volcano plot",
      #uiOutput("VPO")
      plotlyOutput("IGLplot")
      
    )
  )
})
output$Integ1 <- renderUI({
  IN1()
})

IN2 <- eventReactive(input$Iclick2,{
  tabBox(
    width = 12,
    height = "800px",
    title = "Get the results",
    tabPanel(
      
      id="tV","List of gene - miRNA Target gene mappings", 
      #"Volcano plot",
      #uiOutput("VPO")
      #actionButton(inputId = "Iclick2",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000"),
      
      dataTableOutput(outputId = "IMIResult"),
      tags$head(tags$style("#TopGenesTab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
      
      tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
      
    ),
    tabPanel(
      id="tV","Network Plot", 
      #"Volcano plot",
      #uiOutput("VPO")
      forceNetworkOutput("IGLnet", height = "780px")
      
    )
  )
})
output$Integ2 <- renderUI({
  IN2()
})

IN3 <- eventReactive(input$Iclick3,{
  tabBox(
    width = 12,
    height = "1000px",
    title = "Get the results",
    tabPanel(
      
      id="tV","List of LSVs on given gene list / DEG", 
      #"Volcano plot",
      #uiOutput("VPO")
      #actionButton(inputId = "Iclick3",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000"),
      
      dataTableOutput(outputId = "IDASResult"),
      tags$head(tags$style("#TopGenesTab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
      
      tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
      
    ),
    tabPanel(
      id="tV","Plots", 
      br(),
      column(6,plotlyOutput("IDASpie",height = "400px")),
      column(6,plotlyOutput("IDAStoe",height = "400px")),
      br(),br(),br(),br(),br(),br(),br(),br(),br(), br(), br(), br(), br(),br(), br(), br(), br(), br(), br(), br(),br(),br(),
      
      #fluidRow(
      plotlyOutput("IDASALL",height = "400px"),
      #)
      
      
    )
  )
})

output$Integ3 <- renderUI({
  IN3()
})

# Integrate ---------------------------------------------------------------



getDASdifINT <- function(sigi,left,right,genesintb ,manuallististrue){
  if(manuallististrue == FALSE){
    kleiner0fc <- switch(left, "Reticulated Platelets in healthy patients" = "differentially expressed (MP)","Mature Platelets in patients with stable CCS" =  "differentially expressed (MP)", "healthy patients (Mature Platelets)" = "differentially expressed (healthy)","healthy patients (Reticulated Platelets)"= "differentially expressed (healthy)" )
    groesser0fc <- switch(left, "Reticulated Platelets in healthy patients" = "differentially expressed (RP)","Mature Platelets in patients with stable CCS" =  "differentially expressed (RP)", "healthy patients (Mature Platelets)" = "differentially expressed (stable CSS)","healthy patients (Reticulated Platelets)"= "differentially expressed (stable CSS)")
    genesintb$category[genesintb$log2FoldChange < 0] <- kleiner0fc
    genesintb$category[genesintb$log2FoldChange > 0] <- groesser0fc
    genesMP <- genesintb[genesintb$category == kleiner0fc,1]
    genesRP <- genesintb[genesintb$category == groesser0fc,1]
    sigi$category <- "not differentially expressed"
    sigi$category[sigi[,1] %in% genesMP] <- kleiner0fc
    sigi$category[sigi[,1] %in% genesRP] <- groesser0fc
    commontab<-inner_join(genesintb, sigi, by = c("geneID" = "V1"))
    commontab <- commontab[,-c(10,11,12,13,14)]
    colnames(commontab)[c(6,7,8,9)] = c("DE in","E(dPSI) per LSV junction","P(|dPSI|>=0.20) per LSV junction","Type of splicing event")
    fil <- sigi[,c(1,9)]
    fil<- fil[!duplicated(fil),]
    x<-fil%>%
      group_by(category)%>%
      count()
    events <- commontab%>%
      group_by(commontab[,9])%>%
      count()
  }else{
    sigi$category <- "not in list"
    sigi$category[sigi[,1] %in% genesintb] <- "in list"
    genesintb <- as.data.frame(genesintb)
    colnames(genesintb) = "geneID"
    commontab<-inner_join(genesintb, sigi, by = c("geneID" = "V1"))
    commontab <- commontab[,-c(6,7,8,9,10)]
    colnames(commontab)[c(2,3,4,5)] = c("DE in","E(dPSI) per LSV junction","P(|dPSI|>=0.20) per LSV junction","Type of splicing event")
    fil <- sigi[,c(1,5)]
    fil<- fil[!duplicated(fil),]
    x<-fil%>%
      group_by(category)%>%
      count()
    
    events <- commontab%>%
      group_by(commontab[,5])%>%
      count()
  }
  allevents <- sigi%>%
    group_by(V4)%>%
    count()
  
  colnames(events) = c("tose","n")
  
  return(list(commontab,x,events,allevents))
  
}




sav_intersectGeneList <- reactiveVal()
sav_IallGenes <- reactiveVal()
InterGene <- eventReactive(input$Iclick1,{
  List_of_Genes1 <- NULL#c()
  List_of_Genes2 <- NULL#c()
  print(input$INT)
  if(input$INT == "Insert Genelist 1"){
    inGenes <- getGeneText(input$IGL)
    List_of_Genes1 <- inGenes
  }
  if(input$INT == "Insert Genelist 2"){
    inGenes <- getGeneText(input$IGL2)
    List_of_Genes2 <- inGenes
  }
  if(input$INT == "Mature vs. Reticulated 1"){
    kppath <-paths2frames("1",path,input$ITPM1, input$IHoD1,".CSV") # paths2frames
    
    genesintab <- read.csv2(kppath)
    #
    genesintab <- getIDif(genesintab, input$IpCO1, input$IfCO1,input$IUDMI1,3)
    View(genesintab)
    List_of_Genes1 <- genesintab[,1]
  }else if(input$INT == "CCS vs. ACS 1"){
    kppath <-paths2frames("2",path,input$ITPM2, input$IHoD02,".CSV")
    genesintab <- read.csv2(kppath)
    genesintab <- getIDif(genesintab, input$IpCO2, input$IfCO2,input$IUDMI2,3)
    List_of_Genes1 <- genesintab[,1]
  }
  if(input$INTZ == "Mature vs. Reticulated 2" ){
    kppath <-paths2frames("1",path,input$ITPM3, input$IHoD3, ".CSV")
    genesintab <- read.csv2(kppath)
    genesintab <- getIDif(genesintab, input$IpCO3, input$IfCO3,input$IUDMI3,3)
    List_of_Genes2 <- genesintab[,1]
  }else if(input$INTZ =="CCS vs. ACS 2"){
    kppath <-paths2frames("2",path,input$ITPM4, input$IHoD04, ".CSV")
    genesintab <- read.csv2(kppath)
    genesintab <- getIDif(genesintab, input$IpCO4, input$IfCO4,input$IUDMI4,3)
    List_of_Genes2 <- genesintab[,1]
  }
  # print(List_of_Genes1)
  #  print(List_of_Genes2)
  endlist <- intsect(List_of_Genes1,List_of_Genes2)
  sav_intersectGeneList(endlist)
  eins = as.data.frame( List_of_Genes1)
  zw = as.data.frame( List_of_Genes2)
  colnames(eins) = "genes"
  colnames(zw) = "genes"
  combined <- rbind(eins,zw)
  #View(combined)
  sav_IallGenes(combined)
  endlist <- as.data.frame(endlist)
  colnames(endlist) = "genes"
  endlist
})

####für U & D Genes sorgen


#sav_integMIRNA = reactiveVal()
sav_NetMiRNA = reactiveVal()
IntegMiRNA <- eventReactive(input$Iclick2,{
  List_of_Genes1 <- c()
  if(input$INTMI == "Insert Genelist"){
    inGenes <- getGeneText(input$MGIGL)
    List_of_Genes1 <- inGenes
  }
  if(input$INTMI == "Mature vs. Reticulated" ){
    kppath <-paths2frames("1",path,input$MGITPM1, input$MGIHoD1,".CSV")
    genesintab <- read.csv2(kppath)
    genesintab <- getIDif(genesintab, input$MGIpCO1, input$MGIfCO1,input$MGIUDMI1,3)
    List_of_Genes1 <- genesintab[,1]
  }else if(input$INTMI == "CCS vs. ACS 2"){
    kppath <-paths2frames("2",path,input$MGITPM2, input$MGIHoD02,".CSV")
    genesintab <- read.csv2(kppath)
    genesintab <- getIDif(genesintab, input$MGIpCO2, input$MGIfCO2,input$MGIUDMI2,3)
    List_of_Genes1 <- genesintab[,1]
  }
  mir_Nodesp <- miRNA_path("", input$ImiEXP,"N",input$IUDMI, input$ImiDB,input$ImiF2,FALSE,"0")
  mir_Nod <- read.table(mir_Nodesp, sep = "\t", header = FALSE)
  # mirnap <- ImiRNA_path(input$ImiF2)
  #mirna_I_tab <- read.csv2(mirnap)
  MTAB <- getMTAB(input$ImiF2)
  
  sig_mirna_I <- getIDif(MTAB,input$MIpCO2,input$MIfCO2, input$IUDMI,2)
  
  endNodes <- getImirnagene(List_of_Genes1,sig_mirna_I[,1], mir_Nod)
  #sav_integMIRNA(endList)
  sav_NetMiRNA(endNodes)
  endNodes
})
output$IListOfGenes <-renderDataTable(InterGene(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 10))
output$IMIResult <- renderDataTable(IntegMiRNA(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 10))


output$IGLplot <- renderPlotly({
  some <- sav_intersectGeneList()
  all1  <- sav_IallGenes()
  df <- as.data.frame(all1)
  df$category <- "all"
  df$category[df[,1] %in% some] <- "intersect"
  
  x<-df%>%
    group_by(category)%>%
    count()
  fig <- plot_ly(x, labels = ~category, values = ~n, type = 'pie',textposition = 'inside',
                 textinfo = 'percent',
                 insidetextfont = list(color = '#FFFFFF'),
                 hoverinfo = 'percent + text',
                 #text = ~paste(n, 'genes'),
                 marker = list(colors = DAS_colors,
                               line = list(color = '#FFFFFF', width = 1)),
                 #The 'pull' attribute can also be used to create space between the sectors
                 showlegend = TRUE)
  fig <- fig %>% layout(title = 'Percentages of all genes and the intersection of input sets of genes',
                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  
  fig
})


output$IGLnet <- renderForceNetwork({
  t <- sav_NetMiRNA()
  #View(t)
  ensemb <- data.frame(t[,1],t[,3])
  hgnc <-data.frame(t[,1],t[,2])
  colnames(ensemb) <- c("from","to")
  colnames(hgnc) <- c("from", "to")
  edges <- ensemb
  nodes <- data.frame(name = unique(c(edges$from, edges$to)), stringsAsFactors = FALSE)
  nodes$id <- 0:(nrow(nodes) - 1)
  nodes$type <- "Gene"
  nodes$type[startsWith(nodes$name,"hsa-")] <- "mir"
  #nodes$type[nodes[,1] %in% u[,1]] <- "u"
  #nodes$type[nodes[,1] %in% d[,1]] <- "d"
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

ID_tab = reactiveVal()
ID_pie = reactiveVal()
ID_toe = reactiveVal()
ID_alltoe = reactiveVal()
IntegDIFDAS <- eventReactive(input$Iclick3,ignoreInit = T,{
  List_of_Genes1 <- c()
  genesintab <- ""
  manuallististrue <- FALSE
  if(input$INTD == "Insert Genelist"){
    inGenes <- getGeneText(input$MGIGL)
    List_of_Genes1 <- inGenes
    genesintab <- List_of_Genes1
    manuallististrue <- TRUE
  }
  if(input$INTD == "Mature vs. Reticulated" ){
    kppath <-paths2frames("1",path,input$DSITPM1, input$DSIHoD1,".CSV")
    genesintab <- read.csv2(kppath)
    #getDiffTab4DAS(tab,lfc,pval)
    genesintab <- getDiffTab4DAS(genesintab, input$DSIfCO1, input$DSIpCO1)
    List_of_Genes1 <- genesintab[,1]
  }else if(input$INTD == "CCS vs. ACS 2"){
    kppath <-paths2frames("2",path,input$DSITPM2, input$DSIHoD02,".CSV")
    genesintab <- read.csv2(kppath)
    genesintab <- getDiffTab4DAS(genesintab, input$DSIfCO2, input$DSIpCO2)
    List_of_Genes1 <- genesintab[,1]
  }
  
  DASlist <- getDASdata(as.character(input$IDAS1))
  left <- DASlist[[2]]
  right <-DASlist[[3]]
  tb <- DASlist[[1]]
  dassig <- getValidEasy(tb,input$IDASe,input$IDASp)
  IDD <-getDASdifINT(dassig,left,right,genesintab ,manuallististrue)
  ID_tab(IDD[[1]])
  ID_pie(IDD[[2]])
  ID_toe(IDD[[3]])
  ID_alltoe(IDD[[4]])
  ID_tab()
})


output$IDASResult <- renderDataTable(IntegDIFDAS(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 8))
output$IDASALL <- renderPlotly({
  #IntegDIFDAS()
  allevents <- ID_alltoe()
  fig <- plot_ly(allevents, labels = ~V4, values = ~n, type = 'pie',textposition = 'inside',
                 textinfo = 'percent',
                 insidetextfont = list(color = '#FFFFFF'),
                 hoverinfo = 'percent + text',
                 #text = ~paste(n, 'genes'),
                 # marker = list(colors = colors2,
                 #               line = list(color = '#FFFFFF', width = 1)),
                 # #The 'pull' attribute can also be used to create space between the sectors
                 showlegend = TRUE)
  fig <- fig %>% layout(title = 'Percentages of types of significant local splicing events in significant local splicing variations',
                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  
  fig
})
output$IDASpie <- renderPlotly({
  #DIFDAS()
  x<- ID_pie()
  fig <- plot_ly(x, labels = ~category, values = ~n, type = 'pie',textposition = 'inside',
                 textinfo = 'percent',
                 insidetextfont = list(color = '#FFFFFF'),
                 hoverinfo = 'percent + text',
                 #text = ~paste(n, 'genes'),
                 marker = list(colors = DAS_colors,
                               line = list(color = '#FFFFFF', width = 1)),
                 #The 'pull' attribute can also be used to create space between the sectors
                 showlegend = TRUE)
  fig <- fig %>% layout(title = 'Percentages of (differentially) expressed genes with significant local splicing events ',
                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  
  fig
  
})

output$IDAStoe <- renderPlotly({
  #DIFDAS()
  colors2 <- c('rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)')
  events <- ID_toe()
  fig <- plot_ly(events, labels = ~tose, values = ~n, type = 'pie',textposition = 'inside',
                 textinfo = 'percent',
                 insidetextfont = list(color = '#FFFFFF'),
                 hoverinfo = 'percent + text',
                 #text = ~paste(n, 'genes'),
                 marker = list(colors = colors2,
                               line = list(color = '#FFFFFF', width = 1)),
                 #The 'pull' attribute can also be used to create space between the sectors
                 showlegend = TRUE)
  fig <- fig %>% layout(title = 'Percentages of types of significant local splicing events in (differentially) expressed genes',
                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  
  fig
  
})


