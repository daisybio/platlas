# Das <- eventReactive(input$clickDAS1,{
#   tabBox(
#     width =12,
#     height = "1500px",
#     title = "Get the results",
#     tabPanel(
#       id= "tT", "E(PSI) per junction",
#       # uiOutput("DTO")
#       #selectInput(inputId = "UpDown" , label = "Choose up or downregulated genes" , choices = c("upregulated","downregulated"), multiple = FALSE),
#       #actionButton(inputId = "click4",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
#       br(),
#       br(),
#       plotlyOutput("DASBX", height = "800px"),
#       #downloadButton("downloadDASBX", "Download table"),
#     ),
#     tabPanel(
#       id= "tH", 
#       #uiOutput("HMO")
#       
#       title = "Differential Alternative Splicing Activity", background = "red" , solidHeader = TRUE,
#       br(),
#       downloadButton("dlDASVP", "Download plot in PNG format"),
#       br(),
#       plotOutput("DASVP", height = "800px"),
#       #forceNetworkOutput("KPMnet", height = "850px")
#       
#       
#       
#     ),
#     tabPanel(
#       id= "tH", 
#       #uiOutput("HMO")
#       
#       title = "Significantly alternatively spliced Genes", background = "red" , solidHeader = TRUE,
#       #textOutput("DASueber"),
#       selectInput(inputId = "DASUD", label = "Choose significantly alternatively spliced LSV list ", choices = c("Condition 1", "Condition 2"), multiple = FALSE),
#       actionButton(inputId = "clickDASDT",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
#       br(),
#       br(),
#       downloadButton("dlDASDT", "Download Table"),
#       br(),
#       dataTableOutput(outputId = "DASDT"),
#       tags$head(tags$style("#DASDT table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
#       
#       tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
#                              border-top: 1px solid #000000;}",media="screen", type="text/css"))
#       
#     ),
#     tabPanel(
#       id= "tN", 
#       #uiOutput("HMO")
#       
#       title = "Occurences of significant local splicing variations per Gene", background = "red" , solidHeader = TRUE,
#       plotlyOutput("DASBP",height = "800px"),
#       #forceNetworkOutput("KPMGOnet", height = "850px")
#       #  br(),
#       # downloadButton("dlKPMN", "Download plot in PNG format")
#       
#     ),
#     tabPanel(
#       id= "tN",
#       #uiOutput("HMO")
#       
#       title = "LSV in differentially expressed genes", background = "red" , solidHeader = TRUE,
#       span("Choose the filters, which define differentially expressed genes and get a list of the LSVs in those Genes, that also shows if they are upregulated/downregulated",style = "color: black; font-size: 18px"),
#       selectInput(inputId = "DASTPM" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
#       numericInput(inputId = "DASpc",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
#       numericInput(inputId = "DASfc",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
#       actionButton(inputId = "DD_click",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000"),
#       
#       br(),
#       br(),
#       # fluidRow({width =12
#       column(6,plotlyOutput("DASDIF",height = "400px")),
#       column(6,plotlyOutput("DASDIFTOE",height = "400px")),
#       #}),
#       
#       
#       br(),
#       br(),
#       downloadButton("dlDASDIFtab", "Download Table"),
#       br(),
#       dataTableOutput(outputId = "DASDIFtab"),
#       tags$head(tags$style("#DASDIFtab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
#       
#       tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
#                              border-top: 1px solid #000000;}",media="screen", type="text/css"))
#       
#       #forceNetworkOutput("KPMGOnet", height = "850px")
#       #  br(),
#       # downloadButton("dlKPMN", "Download plot in PNG format")
#       
#     ),
#     tabPanel(
#       id= "tN", 
#       #uiOutput("HMO")
#       
#       title = "Types of alternative splicing events in the filtered data", background = "red" , solidHeader = TRUE,
#       plotlyOutput("DASPIE",height = "800px"),
#       #forceNetworkOutput("KPMGOnet", height = "850px")
#       #  br(),
#       # downloadButton("dlKPMN", "Download plot in PNG format")
#       
#     )
#     
#   )
#})

# output$DASO <- renderUI({
#   Das()
# })



# DAS ---------------------------------------------------------------------


DAS_left = reactiveVal("")
DAS_right = reactiveVal("")
DAS_sigtab = reactiveVal()
DAS_tab = reactiveVal()
DAS_um_tab = reactiveVal()
DAS_left_sig_tab =  reactiveVal()
DAS_right_sig_tab =  reactiveVal()

getDiffTab4DAS <- function(tab,lfc,pval){
  # if(mprp ==TRUE){
  #   endt <- tab[tab[,3]<(-lfc) & tab[,7] < pval,]
  # }else{
  endt <- tab[abs(tab[,3])>lfc & tab[,7] < pval,]
  # }
  
  endt <- endt[, c(1,3,7,11,12)]
  return(endt)
}

getDiffTab4DAS4HS <- function(tab, lfc, pval){
  endt <- tab[abs(as.numeric(tab[,3]))>lfc & as.numeric(tab[,4]) < pval,]
  endt <- endt[, c(1,3,4)]
  colnames(endt) = c("geneID","log2FoldChange", "padj")
  return(endt)
}

getDASdif <- function(sigi,left,right,tpm,fc,pval,name_of_column){
  nt <- switch(tpm, "TPM > 0.1"="TPM-0-1","TPM > 0.3"="TPM-0-3","TPM > 1"= "TPM-1","TPM > 2" = "TPM-2","TPM > 0.2"="TPM-0-2" )
  pa <- ""
  kleiner0fc <- switch(left, "Reticulated Platelets in healthy patients" = "differentially expressed (MP)","Mature Platelets in patients with stable CCS" =  "differentially expressed (MP)", "healthy patients (Mature Platelets)" = "differentially expressed (healthy)","healthy patients (Reticulated Platelets)"= "differentially expressed (healthy)" )
  groesser0fc <- switch(left, "Reticulated Platelets in healthy patients" = "differentially expressed (RP)","Mature Platelets in patients with stable CCS" =  "differentially expressed (RP)", "healthy patients (Mature Platelets)" = "differentially expressed (stable CSS)","healthy patients (Reticulated Platelets)"= "differentially expressed (stable CSS)")
  pa <- switch(left, "Reticulated Platelets in healthy patients" = paste0("www/all_healthy_",nt,".CSV"),"Mature Platelets in patients with stable CCS" =  paste0("www/all_diseased_",nt,"_phenotype-disease-stable-CAD.CSV"), "healthy patients (Mature Platelets)" = paste0("www/heat_healthyCAD_",nt,"_MP.CSV"),"healthy patients (Reticulated Platelets)"= paste0("www/heat_healthyCAD_",nt,"_RP.CSV"))
  tab <- read.csv2(pa)
  trueifMPRP <- switch(left, "Reticulated Platelets in healthy patients" = TRUE,"Mature Platelets in patients with stable CCS" = TRUE, "healthy patients (Mature Platelets)" = FALSE,"healthy patients (Reticulated Platelets)"= FALSE)
  endtab <- ""
  c1 <- c()
  c2 <- c()
  colnew <- 0
  if(trueifMPRP == TRUE){
    endtab <- getDiffTab4DAS(tab,fc,pval)
    c1 <- c(10,11,12,13,14)
    c2 <- c(6,7,8,9)
    colnew <- 9
  }else{
    endtab <- getDiffTab4DAS4HS(tab,fc,pval)
    c1 <- c(8,9,10,11,12)
    c2 <- c(4,5,6,7)
    colnew <- 7
  }
  endtab$category[endtab$log2FoldChange < 0] <- kleiner0fc
  endtab$category[endtab$log2FoldChange > 0] <- groesser0fc
  genesMP <- endtab[endtab$category == kleiner0fc,1]
  genesRP <- endtab[endtab$category == groesser0fc,1]
  sigi$category <- "not differentially expressed"
  sigi$category[sigi[,1] %in% genesMP] <- kleiner0fc
  sigi$category[sigi[,1] %in% genesRP] <- groesser0fc
  commontab<-inner_join(endtab, sigi, by = c("geneID" = "V1"))
  commontab <- commontab[,-c1]
  
  colnames(commontab)[c2] = c("DE in","E(dPSI) per LSV junction","P(|dPSI|>=0.20) per LSV junction","Type of splicing event")
  
  fil <- sigi[,c(1,9)]
  fil<- fil[!duplicated(fil),]
  x<-fil%>%
    dplyr::group_by(category)%>%
    dplyr::count()
  
  events <- commontab%>%
    dplyr::group_by(commontab[,colnew])%>%
    dplyr::count()
  
  allevents <- sigi%>%
    dplyr::group_by(V4)%>%
    dplyr::count()
  
  colnames(events) = c("tose","n")
  
  
  return(list(commontab,x,events,allevents))
  
}


#  kp_p = reactiveVal("")
DD_tab = reactiveVal()
DD_pie = reactiveVal()
DD_toe = reactiveVal()
D_alltoe = reactiveVal()  
DAStable <- eventReactive(input$clickDAS1,{
  das1 <- "Mature Platelets vs Reticulated Platelets (stable CCS)"
  DASlist <- getDASdata(das1)
  #DAS_left(DASlist[[2]])
  #DAS_right(DASlist[[3]])
  das_tab <- DASlist[[1]]
  DAS_tab(das_tab)
  
  input_tab_das <- switch(input$DASga,"ENSEMBL ID"= das_tab[,c(1,3,4,5,6,7)],"HGNC symbol"=das_tab[,c(2,3,4,5,6,7)])
  #View(input_tab_das)
  dassig <- getValidEasy(input_tab_das,input$DASe,input$DASp)
  DAS_left("Mature Platelets in patients with stable CCS")
  print(paste0("left: ",DAS_left()))
  print(paste0("right: ",DAS_right()))
  #updateSelectInput(session,"DASUD",choices = c(DAS_left(),DAS_right()))
  DAS_sigtab(dassig)
  #View(DAS_sigtab())
  #DASDataTable()
  DAS_um_tab(input_tab_das)#umformen(input_tab_das)
})

#observe({
#  if(!is.null(input$DAS1)){
#    if(input$DAS1 == "Mature Platelets vs Reticulated Platelets (stable CCS)"){
#      updateSelectInput(session, "DASTPM", choices = c("TPM > 0.2"))
#    }else{
#      updateSelectInput(session, "DASTPM", choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"))
#    }
#  }
#})


DIFDAS <- eventReactive(input$DD_click,{
  das_tpm <- "TPM > 0.2"
  colname <- switch(input$DASga,"ENSEMBL ID"= das_tab[,c(1,3,4,5,6,7)],"HGNC symbol"=das_tab[,c(2,3,4,5,6,7)])
  Difdas <- getDASdif(DAS_sigtab(),DAS_left(),DAS_right(),das_tpm,input$DASfc, input$DASpc,colname)
  DD_tab(Difdas[[1]])
  DD_pie(Difdas[[2]])
  DD_toe(Difdas[[3]])
  D_alltoe(Difdas[[4]])
})


DAS_end_sig_tab = reactiveVal()#vllt reactiveValues()

DASDataTable <- eventReactive(input$clickDAS1,{
  # r <- NULL
  # l <- NULL
  # e <- NULL
  # if(as.character(input$DASUD) == DAS_left()){
  #   l <-  DAS_sigtab()[DAS_sigtab()$V2 < 0, ]
  #   l <-  l[order(l$V2, decreasing = FALSE),] 
  #   DAS_end_sig_tab(l)
  #   e <- l
  # }else if(as.character(input$DASUD) == DAS_right()){
  #   r <- DAS_sigtab()[DAS_sigtab()$V2 > 0, ]
  #   r <-  r[order(r$V2, decreasing = TRUE),]
  #   DAS_end_sig_tab(r)
  #   e <- r
  # }
  # e <- e[,c(1,2,3,4)]
  # e <- e[!duplicated(e),]
  # colnames(e) <- c('Gene ID','E(dPSI) per LSV junction','P(|dPSI|>=0.20) per LSV junction','Type of splicing event')
  # e
  DAS_sigtab()
})
saved_DASVP = reactiveVal()
output$DASVP <- renderPlot({
  ##View(DAS_um_tab())
  
  
  
  
  # p <- myVolcanoPlot(DAS_um_tab(),
  #                     data_type = "gene",
  #                     pCutoff = pcReactive(),
  #                     FCcutoff = fcReactive(),
  #                     xlim = c(-6, 6),
  #                     lab_selection = c(top10_upreg_genes,top10_downreg_genes),
  #                     symbol_type = colnames(d)[annot] , titel = Titl(), subtitel = paste(nameHod(),"patients dataset considering", nameTPM1())
  # )
  p <- EnhancedVolcano(DAS_um_tab(),
                       lab = DAS_um_tab()[,1],
                       x = 'V2',
                       y = 'V3',
                       title = 'Differential Alternative Splicing',
                       subtitle = paste0(DAS_left()," vs ",DAS_right()),
                       pCutoff = 0.05,
                       FCcutoff = 2,
                       pointSize = 3.0,
                       col=c('lightgrey', 'lightgrey', 'lightgrey', 'red'),
                       #shape = c(0, 0, 0, 0),
                       colAlpha = 4/5,
                       legendLabSize = 14,
                       
                       labSize = 3.0,
                       legendLabels=c('Not sig.','Log(base 2) E(dPSI)x10 per LSV junction,','P(|dPSI|>=0.20)-1 per LSV junction ',
                                      'Log(base 2) E(dPSI)x10 per LSV junction and P(|dPSI|>=0.20)-1 per LSV junction'),
                       legendPosition = 'right',
                       
                       legendIconSize = 4.0
  )
  saved_DASVP(p)
  p
})

output$DASBX <- renderPlotly({
  DAStable()
  if(!is.null(DAS_sigtab())){
    das_df<- DAS_sigtab()
    melt_data <- NULL
    if(input$DASga == "ENSEMBL ID"){
      melt_data <- reshape2::melt(das_df, id = c("Ensembl ID","Type","Confidence","DELTAPSI"), value.name = "PSI value") 
    }else{
      melt_data <- reshape2::melt(das_df, id = c("Gene Name","Type","Confidence","DELTAPSI"), value.name = "PSI value")
    }
    #View(melt_data)
    #con1 <- as.data.frame(as.numeric(DAS_sigtab()[,8]))
  #print(con1)
    #con2 <- as.data.frame(as.numeric(DAS_sigtab()[,7]))
    #colnames(con1) = c("E(PSI)")
    #colnames(con2) = c("E(PSI)")
    #con1[,2] <- DAS_left()
    #con2[,2] <- DAS_right()
    #combined <- rbind(con1,con2)
    fig <- plot_ly(data = melt_data,y = ~`PSI value`, type = "box", quartilemethod="linear", color = ~variable, colors = c("#c16868","#990012"))
    fig <- fig %>% layout(title = paste0("PSI values ", "MPs", " vs ", "RPs"), yaxis = list(title = "PSI values",titlefont = f))
}
  })

output$DASDT <- renderDataTable(DASDataTable(),escape = FALSE,extensions = c("Buttons","Scroller"),options = list(lengthMenu = c(10, 15, 20), pageLength = 15,scrollX=TRUE, dom = 'Bfrtip',bPaginate = FALSE,
                                                                                                                  buttons =list(list(
                                                                                                                    extend = 'collection',
                                                                                                                    buttons = c('csv', 'excel'),
                                                                                                                    text = 'Download'
                                                                                                                  ))))

output$DASSC <- renderPlotly({
  if(!is.null(DAS_left())){
  df = data.frame(as.numeric(DAS_tab()[,5]),as.numeric(DAS_tab()[,6]))
  sigdif = data.frame(as.numeric(DAS_sigtab()[,5]),as.numeric(DAS_sigtab()[,6]))
  scattab <- getScatter(df,sigdif)
  colnames(scattab) = c("start", "end" , "junct")
  fig <- plot_ly(scattab, x = ~start, y =  ~end, type = 'scatter', mode = 'markers', symbol = ~junct, symbols = c("circle","o"),marker = list(size = 14), color = ~junct, colors = c("#e6b6b0","#990012"))
  #fig <- fig %>% add_trace(x = start1,y = ~end1, name = "significantly alternatively spliced junctions", symbol = ~end1,symbols = "o",type = "scatter", mode = "markers")
}})

output$DASBP <- renderPlotly({
  if(!is.null(DAS_sigtab())){
  oc <- dplyr::count(DAS_sigtab(), var = DAS_sigtab()[,1])
  colnames(oc) = c('Gene', "occurrence")
  oc <- oc[order(oc$occurrence,decreasing = TRUE),]
  form <- list(categoryorder = "array",
               categoryarray = oc$Gene,
               title = "Gene")
  fig <- plot_ly(oc,x = ~Gene,y = ~occurrence,
                 name = "Occurrences of significant local splicing variations in Genes",
                 type = "bar",
                 color = ~occurrence,
                 colors = c("#fff68f","#ba3131")
  )
  fig <- fig %>% layout(title = "Occurrences of significant local splicing variations in Genes", xaxis = form)
}})

DAS_colors <- c('rgb(211,94,96)','rgb(114,147,203)','rgb(128,133,133)' )

output$DASDIF <- renderPlotly({
  if(!is.null(DAS_left())){  
  DIFDAS()
  x<- DD_pie()
  font_1 <- list(
    size = 12)
  fig <- plot_ly(x, labels = ~category, values = ~n, type = 'pie',textposition = 'inside',
                 textinfo = 'percent',
                 insidetextfont = list(color = '#FFFFFF'),
                 hoverinfo = 'percent + text',
                 #text = ~paste(n, 'genes'),
                 marker = list(colors = DAS_colors,
                               line = list(color = '#FFFFFF', width = 1)),
                 #The 'pull' attribute can also be used to create space between the sectors
                 showlegend = TRUE)
  fig <- fig %>% layout(title = list(text = 'Percentages of (differentially) expressed genes with significant local splicing events', font = font_1),
                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  
  fig
  
}})

output$DASDIFTOE <- renderPlotly({
  if(!is.null(DAS_left())){
  DIFDAS()
  colors2 <- c('rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)')
  font_1 <- list(
    size = 12)
  events <- DD_toe()
  fig <- plot_ly(events, labels = ~tose, values = ~n, type = 'pie',textposition = 'inside',
                 textinfo = 'percent',
                 insidetextfont = list(color = '#FFFFFF'),
                 hoverinfo = 'percent + text',
                 #text = ~paste(n, 'genes'),
                 marker = list(colors = colors2,
                               line = list(color = '#FFFFFF', width = 1)),
                 #The 'pull' attribute can also be used to create space between the sectors
                 showlegend = TRUE)
  fig <- fig %>% layout(title = list(text = 'Percentages of types of significant local splicing events in (differentially) expressed genes', font = font_1),
                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  
  fig
  
}})

output$DASDIFtab <-renderDataTable(DD_tab(),escape = FALSE,extensions = c("Buttons","Scroller"), options = list(lengthMenu = c(10, 15, 20), pageLength = 8, dom = 'Bfrtip',bPaginate = FALSE,
                                                                                                                buttons =list(list(
                                                                                                                  extend = 'collection',
                                                                                                                  buttons = c('csv', 'excel'),
                                                                                                                  text = 'Download'
                                                                                                                )) ))

output$DASPIE <- renderPlotly({
  if(!is.null(DAS_sigtab())){ 
  #DIFDAS()
  das_df<- DAS_sigtab()
  das_df <- das_df %>% 
    mutate(Type = strsplit(as.character(Type), ",")) %>% 
    unnest(Type)
  allevents <- das_df[,c(1,2)]%>%
    dplyr::group_by(Type)%>%
    dplyr::count()
  allevents <- allevents %>%
    mutate(Types = case_when(
      Type == "ES" ~ "ES: Exon Skipping",
      Type == "A5" ~ "A5: Alternative 5′ splice site",
      Type == "A3" ~ "A3: Alternative 3′ splice site",
      Type == "IR" ~ "IR: Intron retention"
    ))
  
  
  fig <- plot_ly(allevents, labels = ~Types, values = ~n, type = 'pie',textposition = 'inside',
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
  
  
}})

