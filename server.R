source("www/functions.R")

server <- function(input,output,session){

  UIMI <- eventReactive(input$clickMIF,{
    tagList( 
      tabBox(
        tabPanel(
          id="t5","DE miRNA", 
          span("Significant differentially expressed miRNA",style = "color: black; font-size: 16px" ),
          br(),
          downloadButton("dlDTMI", "Download table"),
          br(),
          dataTableOutput("DTMI"),
          tags$head(tags$style("#DTMI table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
          
          tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
          
        ),
        tabPanel(
          id="t6","Volcano Plot", 
          
          plotOutput("VPMI"),
          downloadButton("dlVPMI", "Download plot in PNG format")
          
        ),
        tabPanel(
          id="t7","HeatMap", 
  
          plotOutput("HMMI"),
          downloadButton("dlHMMI", "Download plot in PNG format")
        )
      )
    )
  }
  )
  MITarg <- eventReactive(input$clickMIF2,{
    tagList( 
      tabBox(
        width = 12, 
        height = "1000px",
        tabPanel(
          id= "tE","miRNA Targets",
          downloadButton("dlMIDT", "Download table"),
          br(),br(),
          dataTableOutput("MIDT"),
          tags$head(tags$style("#DTMIT table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
          
          tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
        ),
        tabPanel(
          id="t10","miRNA Target Enrichment", 
          #selectInput(inputId = "MSIG" , label = "Choose significant or all mRNA" , choices = c("significant","all"), multiple = FALSE),
          #actionButton(inputId = "clickSIG",label = "Get mRNA ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
          span("The 30 Most significant GO terms" ,style = "color: black; font-size: 18px"),br(),
          downloadButton("dlDTMIT", "Download table"),
          br(),br(),
          dataTableOutput("DTMIT"),
          tags$head(tags$style("#DTMIT table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
          
          tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
          
          ),
        tabPanel(
          id="t11","miRNA Target Network Graph", 
          span("Search if miRNA Targets are differentially expressed, by first defining differentially expressed genes with the filters below. The results show miRNAs (darkblue nodes) and their mapping genes (light blue nodes), while green nodes are upregulated genes and red nodes are downregulated genes.",
          style = "color: black; font-size: 18px"),
          br(),
          selectInput(inputId = "mdTPM" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
          numericInput(inputId = "mdpc",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
          numericInput(inputId = "mdfc",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
          actionButton(inputId = "md_click",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000"),
          forceNetworkOutput("MInet", height = "850px")
          
        ),
        tabPanel(
          id="t12","miRNA Target Enrichment Graph", 
          span("miRNA Target Enrichment Graph",style = "color: black; font-size: 24px"),br(),br(),
          forceNetworkOutput("MIGOnet", height = "850px")
          
        ),
        tabPanel(
          title = "percentage of target genes in GO categories", solidHeader = TRUE,# width = 12,
          plotlyOutput("MIpG",height = "800px"),
          #downloadButton("dlPNGGO", "Download plot in PNG format")
        ),
        tabPanel(
          title = "enrichment FDR per GO categories", solidHeader = TRUE,# width = 12,
          #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
          #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
          #br(),
          plotlyOutput("MIpFE",height = "800px"),
          #downloadButton("dlPNGGO", "Download plot in PNG format")
        )
      )
    )
  })

  CRUI <- eventReactive(input$clickCirc,{
    tagList(
    tabBox(
      width = 12, 
      height = "1000px",
      tabPanel(id = "CIVP1", "Differentially Expressed circRNA", "",
               br(),
               downloadButton("dlCIVP", "Download plot in PNG format"),
               br(),
               plotOutput("CIVP", height = "800px")
      ),
      tabPanel(id = "CR1", "List of differentially expressed circRNA", "",
               selectInput(inputId = "circUD" , label = "Choose up or downregulated ciRNAs" , choices = c("upregulated","downregulated"), multiple = FALSE),
               actionButton(inputId = "clickCUD",label = "Get ciRNAs ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
               br(),
               downloadButton("dlCTO", "Download Table"),
               br(),
               dataTableOutput("CTO"),
               tags$head(tags$style("#CTO table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
               
               tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                                    border-top: 1px solid #000000;}",media="screen", type="text/css"))
               
               ),
      # tabPanel(id = "CR2", "Heatmap", "",
      #          plotOutput("CHM")
      # ),
      tabPanel(id = "CR3", "Mapped out circRNA", "",
               plotOutput("CGviz")
      )
      
    )
    )
  })
  
  output$UIMI <-renderUI({
    UIMI()
    
  })   
  output$MITarg <-renderUI({
    MITarg()
    
  })   
  
  output$Crna <- renderUI({
    CRUI()
  })
  
  output$DASO <- renderUI({
    Das()
  })
   
  DIFEX <- eventReactive(c(input$click1, input$click2),ignoreInit = T,{ 
      tabBox(
        width = 12,
        height = "1000px",
        title = "Get the results",
        tabPanel(
          id="tV","Volcano plot", 
          #"Volcano plot",
          br(), 
          selectInput(inputId = "DVP" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
          actionButton(inputId = "clickDVP",label = "Get plot ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
          br(),
          br(),
          #uiOutput("VPO")
          plotOutput("VP",height = "750px"),
          br(),
          
          downloadButton("dlPNG", "Download plot in PNG format")
        ),
        tabPanel(
          id= "tT", "Differentially Expressed Genes",
          # uiOutput("DTO")
          selectInput(inputId = "UpDown" , label = "Choose up or downregulated genes" , choices = c("upregulated","downregulated"), multiple = FALSE),
          actionButton(inputId = "click4",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
          downloadButton("downloadAll", "Download table"),
          br(),
          br(),
          dataTableOutput(outputId = "TopGenesTab"),
          tags$head(tags$style("#TopGenesTab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
          
          tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
          
          ),
        tabPanel(
          id= "tB", 
          # uiOutput("BPO")
          title = "Biotypes of the Detected Genes", solidHeader = TRUE, # background = "white" ,
          plotOutput("biotype1",height = "800px"),
          br(),
          downloadButton("dlPNG2", "Download plot in PNG format")
        ),
        tabPanel(
          id= "tH", 
          #uiOutput("HMO")
          title = "Counts per samples per genes", background = "red" , solidHeader = TRUE,
          br(),
          downloadButton("dlPNG3", "Download plot in PNG format"),
          br(),
          plotOutput("heatmap",height = "800px"),
          br(),
         # 
          
        ),
        tabPanel(
          id= "tP", 
          #uiOutput("HMO")
          title = "PCA Analysis Biplot", background = "red" , solidHeader = TRUE,
          br(),
          downloadButton("dlpcabi", "Download plot in PNG format"),
          br(),
          plotOutput("pcabi",height = "800px"),
          br(),
          # downloadButton("dlPNG3", "Download plot in PNG format")
          
        ),
        tabPanel(
          id= "tP", 
          #uiOutput("HMO")
          title = "PCA Analysis: CCS vs. MI and MP vs. RP", background = "red" , solidHeader = TRUE,
          br(),
          downloadButton("dlpcah", "Download plot in PNG format"),
          br(),
          plotOutput("pcah",height = "800px"),
          br(),
          # downloadButton("dlPNG3", "Download plot in PNG format")
          
        )
      )
    })
  output$DIFEXtab <- renderUI({
   DIFEX()
  })
  
  FER <- eventReactive(c(input$clickGO, input$clickGO2),ignoreInit = T,{
    tabBox(
      width = 12,
      height = "1000px",
      title = "Get the results",
      tabPanel(
        id="GOtab",title = "Top 500 GO Terms", solidHeader = TRUE,
        "",
        selectInput(inputId = "FeA" , label = "Choose gene annotation in GO categories" , choices = c("ENSEMBL ID","ENTREZ ID","HGNC symbol"), multiple = FALSE),
        actionButton(inputId = "clickFeA",label = "Get GOs",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
       br(), 
       br(),
       downloadButton("downloadGO", "Download table"),
        br(),
        br(),
        dataTableOutput(outputId = "GOTab"),#width = 12,
        
        # .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th
        tags$head(tags$style("#GOTab table {background-color: #DCDCDC; color : #000000}", media="screen", type="text/css")),
        tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th
                             {border-top: 1px solid #000000}"))
        ),
      tabPanel(
        title = "percentage of DE genes in GO categories", solidHeader = TRUE,# width = 12,
        br(),
        plotlyOutput("GOpG",height = "800px"),
        #downloadButton("dlPNGGO", "Download plot in PNG format")
      ),
      tabPanel(
        title = "enrichment FDR per GO categories", solidHeader = TRUE,# width = 12,
        #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
        #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
        #br(),
        br(),
        plotlyOutput("GOpFE",height = "800px"),
        #downloadButton("dlPNGGO", "Download plot in PNG format")
      ),
      tabPanel(
        title = "increased/decreased categories and their terms & their genes", solidHeader = TRUE,# width = 12,
        span("Increased/decreased categories and their terms in relation to upregulated and downregulated genes in the corresponding category",style = "color: black; font-size: 19px"),
        br(),
        downloadButton("dlGOC", "Download plot in PNG format"),
        br(),
        plotOutput("GOC",height = "800px"),
        #downloadButton("dlPNGGO", "Download plot in PNG format")
      ),
      tabPanel(
        title = "GO categories and their z-score", solidHeader = TRUE,# width = 12,
        span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
        br(),
        withMathJax(),
        helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
        br(),
        br(),
        br(),
        plotlyOutput("GOZ",height = "800px"),
        #downloadButton("dlPNGGO", "Download plot in PNG format")
      ),
      
      tabPanel(
        id= "tN", 
        #uiOutput("HMO")
        
        title = "GO Network", background = "red" , solidHeader = TRUE,
        #plotOutput(""),
        forceNetworkOutput("GOnet", height = "850px")
        #  br(),
        # downloadButton("dlKPMN", "Download plot in PNG format")
        
      ),
      tabPanel(
        id= "tH", 
        #uiOutput("HMO")
        
        title = "GO Hierarchy Tree", background = "red" , solidHeader = TRUE,
        #plotOutput(""),
        imageOutput("GOH", height = "850px"),
        br(),
        downloadButton("dlGOH", "Download plot in PNG format")
        
      )
    )
  })
  
  output$DTGO <- renderUI({
    FER()
  })
  
  Das <- eventReactive(input$clickDAS1,{
    tabBox(
      width =12,
      height = "1500px",
      title = "Get the results",
      tabPanel(
        id= "tT", "E(PSI) per junction",
        # uiOutput("DTO")
        #selectInput(inputId = "UpDown" , label = "Choose up or downregulated genes" , choices = c("upregulated","downregulated"), multiple = FALSE),
        #actionButton(inputId = "click4",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
        br(),
        br(),
        plotlyOutput("DASBX", height = "800px"),
        #downloadButton("downloadDASBX", "Download table"),
      ),
      tabPanel(
        id= "tH", 
        #uiOutput("HMO")
        
        title = "Differential Alternative Splicing Activity", background = "red" , solidHeader = TRUE,
        br(),
        downloadButton("dlDASVP", "Download plot in PNG format"),
        br(),
        plotOutput("DASVP", height = "800px"),
        #forceNetworkOutput("KPMnet", height = "850px")
        
        
        
      ),
      tabPanel(
        id= "tH", 
        #uiOutput("HMO")
        
        title = "Significantly alternatively spliced Genes", background = "red" , solidHeader = TRUE,
        #textOutput("DASueber"),
        selectInput(inputId = "DASUD", label = "Choose significantly alternatively spliced LSV list ", choices = c("Condition 1", "Condition 2"), multiple = FALSE),
        actionButton(inputId = "clickDASDT",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
        br(),
        br(),
        downloadButton("dlDASDT", "Download Table"),
        br(),
        dataTableOutput(outputId = "DASDT"),
        tags$head(tags$style("#DASDT table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
        
        tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                             border-top: 1px solid #000000;}",media="screen", type="text/css"))
        
      ),
      tabPanel(
        id= "tN", 
        #uiOutput("HMO")
        
        title = "Occurences of significant local splicing variations per Gene", background = "red" , solidHeader = TRUE,
        plotlyOutput("DASBP",height = "800px"),
        #forceNetworkOutput("KPMGOnet", height = "850px")
        #  br(),
        # downloadButton("dlKPMN", "Download plot in PNG format")
        
      ),
      tabPanel(
        id= "tN",
        #uiOutput("HMO")

        title = "LSV in differentially expressed genes", background = "red" , solidHeader = TRUE,
        span("Choose the filters, which define differentially expressed genes and get a list of the LSVs in those Genes, that also shows if they are upregulated/downregulated",style = "color: black; font-size: 18px"),
        selectInput(inputId = "DASTPM" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
        numericInput(inputId = "DASpc",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
        numericInput(inputId = "DASfc",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
        actionButton(inputId = "DD_click",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000"),
        
        br(),
        br(),
       # fluidRow({width =12
         column(6,plotlyOutput("DASDIF",height = "400px")),
         column(6,plotlyOutput("DASDIFTOE",height = "400px")),
        #}),
        
        
        br(),
       br(),
       downloadButton("dlDASDIFtab", "Download Table"),
       br(),
        dataTableOutput(outputId = "DASDIFtab"),
        tags$head(tags$style("#DASDIFtab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
        
        tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                             border-top: 1px solid #000000;}",media="screen", type="text/css"))
        
        #forceNetworkOutput("KPMGOnet", height = "850px")
        #  br(),
        # downloadButton("dlKPMN", "Download plot in PNG format")

      ),
      tabPanel(
        id= "tN", 
        #uiOutput("HMO")
        
        title = "Types of alternative splicing events in the filtered data", background = "red" , solidHeader = TRUE,
        plotlyOutput("DASPIE",height = "800px"),
        #forceNetworkOutput("KPMGOnet", height = "850px")
        #  br(),
        # downloadButton("dlKPMN", "Download plot in PNG format")
        
      )
      
    )
  })
 
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
  

#___________matrix erstellen___________
  clicked1 = reactiveVal(isolate(input$click1))
  clicked2 = reactiveVal(isolate(input$click2))
  Titl <- reactiveVal("")
  nameHod <- reactiveVal("")
  nameTPM1 <- reactiveVal("")
  pcReactive <- reactiveVal("")
  fcReactive <- reactiveVal("")
  
  Mat <- eventReactive(c(input$click1,input$click2),{#ignoreInit = T,
    ma <- data.frame()
    if(clicked1() < input$click1){
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
      clicked2(input$click2)
    }
  
    neu <- ma[!(ma[,7]=="#NV"),]
   # neu %>% drop_na(padj)
    
    #neu <- ma
    #View(neu)
    neu
  
  })
  

#___________heatmap___________
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
  
})
  
  
  # output$heatmap <- renderPlot({
  # bg <- colorRampPalette(c('blue','green','yellow'))
  # heatmap.2(heat(), main = paste("Heatmap (", as.character(nameTPM1()) ,as.character(nameHod()), ")"), trace = "none", margins = c(10,12),col = bg )
  saved_heatmap_plot = reactiveVal()
  output$heatmap <- renderPlot({
    heat()
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
   # plot_ly(x= colnames(hm()),y = rownames(hm()),z = hm(), colors = colorRamp(c("blue","green" ,"yellow")), type = "heatmap")
    
  })
  saved_pcabi = reactiveVal()
  output$pcabi<- renderPlot({
    heat()
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
  })
  saved_pcah = reactiveVal()
  output$pcah <- renderPlot({
    heat()
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
  })
  
# }
# )
  VPannot <- eventReactive(input$clickDVP,{
    res <- switch (input$DVP, "ENSEMBL ID"= 1, "HGNC symbol" = 11)
    res
  })

# ___________Volcano Plot___________
  savedVP_plot = reactiveVal()
  output$VP <- renderPlot({
    d <- Mat()
    d <- d[!(is.na(d$padj)),]
    d <- d %>% mutate(padj = str_replace(padj,"E", "e"))
    d <- d %>% mutate(padj = str_replace_all(padj,",", "."))
    d[,7] = as.numeric(d[,7])
    annot <- VPannot()
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
     p
    
      })
  

#___________biotypesPlot___________
    
     matt<-eventReactive(c(input$click1,input$click2),{
    print("in matt ")
    genes <- Mat()[ ,1]
    fc <- Mat()[ ,3]
    pval<- Mat()[ ,7]
    biotype <-Mat()[ ,12]
    df <- data.frame(genes,fc,pval,biotype)
    significant <- subset(df, pval< pcReactive() & abs(fc)> fcReactive())
    protein_coding<- as.numeric(lengths(subset(significant,biotype == "protein_coding"))[1])
    misc_RNA <- as.numeric(lengths(subset(significant,biotype == "misc_RNA"))[1])
    Mt_rRNA<- as.numeric(lengths(subset(significant,biotype == "Mt_rRNA"))[1])
    processed_pseudogene <- as.numeric(lengths(subset(significant,biotype == "processed_pseudogene"))[1])
    Mt_tRNA<- as.numeric(lengths(subset(significant,biotype == "Mt_tRNA"))[1])
    transcribed_unprocessed_pseudogene <- as.numeric(lengths(subset(significant,biotype == "transcribed_unprocessed_pseudogene"))[1])
    biot2 <- c(protein_coding,misc_RNA,Mt_rRNA,processed_pseudogene,Mt_tRNA,transcribed_unprocessed_pseudogene)
   
    })
     saved_biotype_plot = reactiveVal()
     output$biotype1<- renderPlot({
      p <- barplot(matt(),
              main = paste("Biotypes of the significant Genes",nameTPM1()),
              xlab = "Number of Genes",
              ylab = "Biotypes",
              names.arg = c("protein coding", "misc_RNA","Mt_rRNA", "proc. pseudogene", "Mt_tRNA", "transcr_unproc_pseudog"),
              col = "darkred",
              horiz = FALSE,
              cex.names=0.8,
              las=1
              )
      saved_biotype_plot(p)
      p
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

    
#___________download_BUTTONS___________
     
     output$downloadAll <- downloadHandler(
       filename = function() {
         paste("Differentially_Expressed_Genes", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(saved_ud_tab(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     output$downloadGO <- downloadHandler(
       filename = function() {
         paste("Functional_Enrichment_Results", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(saved_GO(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     output$dlGOC <- downloadHandler(
       filename = function() {
         paste("circular_functional_enrichment_plot", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         sav_cir()
         dev.off()
       },
       contentType = 'image/png'
     )
     
     output$dlPNG <- downloadHandler(
       filename = function() {
         paste("VolcanoPlot", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         savedVP_plot()
         dev.off()
       },
       contentType = 'image/png'
     )
     output$dlPNG2 <- downloadHandler(
       
       filename = function() {
         paste("Biotypes", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         saved_biotype_plot()
         dev.off()
       },
       contentType = 'image/png'
     )
     output$dlPNG3 <- downloadHandler(
       filename = function() {
         paste("Heatmap", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         saved_heatmap_plot()
         dev.off()
         },
       contentType = 'image/png'
     )
     output$dlpcabi <- downloadHandler(
       filename = function() {
         paste("Heatmap", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         saved_pcabi()
         dev.off()
       },
       contentType = 'image/png'
     )
     output$dlpcah <- downloadHandler(
       filename = function() {
         paste("Heatmap", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         saved_pcah()
         dev.off()
       },
       contentType = 'image/png'
     )
     output$dlPNGGO <- downloadHandler(
       filename = function() {
         paste("GO_Barplot", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         GOdf <- GO()
         b <- subset(GOdf,GOdf[,7] > 2)
         b <- b[order(-b[,7]), ]
         nbl <- unique(head(b[, c(6,7)],10),incomparables=FALSE, fromLast=FALSE, by = nbl[,1])
         ggplot(nbl, aes(x = categorization, y = fold_enrichment)) +        # Create barchart with ggplot2
           geom_bar(stat = "identity",fill = "blue")+
           coord_flip()
         
          },
       contentType = 'image/png'
     )
     output$dlDTMI <- downloadHandler(
       filename = function() {
         paste("miRNA_Diff_Exp_Res", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(getDTVP(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     output$dlVPMI <- downloadHandler(
       filename = function() {
         paste("miRNA_Volcanoplot", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         saved_VPMI()
         dev.off()
       },
       contentType = 'image/png'
     )
     output$dlHMMI <- downloadHandler(
       filename = function() {
         paste("miRNA_Heatmap", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         saved_HMMI()
         dev.off()
       },
       contentType = 'image/png'
     )
     output$dlMIDT <- downloadHandler(
       filename = function() {
         paste("miRNA_Targets", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(getMIPaths(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     output$dlDTMIT <- downloadHandler(
       filename = function() {
         paste("miRNA_Target_Enrichment", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(saved_DTMIT(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     output$dlCIVP <- downloadHandler(
       filename = function() {
         paste("VolcanoPlot_circRNA", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         saved_CIVP()
         dev.off()
       },
       contentType = 'image/png'
     )
     output$dlCTO <- downloadHandler(
       filename = function() {
         paste("circRNAs", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(cRNA_ud_dt(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     output$downloadKPMD <- downloadHandler(
       filename = function() {
         paste("Functional_Enrichment_KPM", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(Kpm(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     output$dlDASVP <- downloadHandler(
       filename = function() {
         paste("VolcanoPlot_DAS", ".png", sep = "")
       },
       content = function(file) {
         png(file)
         saved_DASVP()
         dev.off()
       },
       contentType = 'image/png'
     )
     output$dlDASDT <- downloadHandler(
       filename = function() {
         paste("DiffAlternative_Splicing", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(DASDataTable(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     output$dlDASDIFtab <- downloadHandler(
       filename = function() {
         paste("DAS_and_DE", ".csv", sep = "")
       },
       content = function(file) {
         write.csv2(DD_tab(), file, row.names = FALSE)
       },
       contentType = "text/csv"
     )
     
#------------------------------------------GO
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
   # output$GOTab2 <- renderDataTable(GO2(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 5, width = 700)
   # )
#------------------------------miRNA
#sig_data = reactiveVal(0)
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
  # 
  
  
  rownames(neu) = make.names(names, unique=TRUE)
  end<- neu
  
  
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
output$DTMI <- renderDataTable(getDTVP(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 5, width = 700)
)
saved_VPMI = reactiveVal()
output$VPMI <- renderPlot({
 p <- EnhancedVolcano(getVP(),
                  lab = getVP()[,1],
                  x = 'log2FoldChange',
                  y = 'padj',####gucken
                  xlim = c(-5, 8),
                  #title = 'N061011 versus N61311',
                  pCutoff = as.numeric(input$mpCO),
                  FCcutoff = as.numeric(input$mfCO),
                  xlab = bquote(~Log[2]~ 'fold change'),
                  #title = ' ',
                  subtitle = ' ',##,
                  col=c('lightgrey', 'lightgrey', 'lightgrey', 'red'),
                  #shape = c(0, 0, 0, 0),
                  colAlpha = 4/5,
                  legendPosition = 'right',
                  legendLabSize = 14,
                  legendIconSize = 4.0
  )
 saved_VPMI(p)
 p
})
saved_HMMI = reactiveVal()
 output$HMMI <- renderPlot({
   DBO <- "www/DBO-samples.tsv"
   dbo <- read.table(DBO, header = FALSE, sep = "\t", quote = "")
   rownames(dbo) = dbo[,1]
   dbo <- dbo[,-1]
   colnames(dbo) = c("Disease", "Platelet type")
   pl_names = c("red","blue")
   names(pl_names) = c("RP", "MP")
   d_names = c("orange","darkgreen")
   names(d_names) = c("MI", "stable CAD")
   
   ann_col = list(`Platelet type` = pl_names,Disease = d_names)
   p <- pheatmap(getHM(),
            color = inferno(10),
            annotation_colors = ann_col,
            annotation_col = dbo,#dbo %>% select('Platelet type' = dbo$`Platelet type`, 'Disease' = dbo$Disease),
            show_rownames = TRUE, show_colnames = TRUE, scale="row")
   saved_HMMI(p)
   p
   
 }
 )
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

 #-----------------------------------------------CircRNA
  complete_circRNA_df = reactiveVal()
  upreg_circRNA = reactiveVal()
  downreg_circRNA = reactiveVal()
  CRNA_dt <- eventReactive(input$clickCirc,{
    #print(input$ciDS)
    enddf <- circRNA_path(input$ciDS)
    #View(enddf)
    complete_circRNA_df(enddf)
    up <- enddf[enddf$log2FoldChange >= input$fCOcirc & enddf$padj <= input$pCOcirc,]
    down <- enddf[enddf$log2FoldChange <= input$fCOcirc & enddf$padj <= input$pCOcirc,]
    upreg_circRNA(up)
    downreg_circRNA(down)
    
  })
  
  cRNA_ud_dt <- eventReactive(input$clickCUD,{
    use_dt <- NULL
    if(input$circUD == "upregulated"){
      use_dt = upreg_circRNA()
    }else if(input$circUD == "downregulated"){
      use_dt = downreg_circRNA()
    }
   # View(use_dt)
    use_dt<- use_dt %>% separate(colnames(use_dt)[1], c("location", "strand"), "_")%>% separate(location, c("chromosome", "coordinates"), ":")%>% separate(coordinates, c("start", "stop"), "-")
    use_dt <- use_dt[!(is.na(use_dt$strand)),]
    use_dt
    
  })
  saved_CIVP = reactiveVal()
  output$CIVP <- renderPlot({
    CRNA_dt()
    p <- EnhancedVolcano(complete_circRNA_df(),
                                 lab = complete_circRNA_df()[,1],
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 xlim = c(-5, 8),
                                 #title = 'N061011 versus N61311',
                                 pCutoff = as.numeric(input$pCOcirc),
                                 FCcutoff = as.numeric(input$fCOcirc),
                                 xlab = bquote(~Log[2]~ 'fold change'),
                                 #title = ' ',
                                 subtitle = ' ',##,
                                 col=c('lightgrey', 'lightgrey', 'lightgrey', 'red'),
                                 #shape = c(0, 0, 0, 0),
                                 colAlpha = 4/5,
                                 legendPosition = 'right',
                                 legendLabSize = 14,
                                 legendIconSize = 4.0)
    saved_CIVP(p)
    p
  })
  
  output$CTO <- renderDataTable(cRNA_ud_dt(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 10, width = 700)
  )
  
  
                
 
 
 #----------------------------------------------KPM
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
  #--------------------------------------------------------DAS
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
    #View(tab)
    #View(lfc)
    #View(pval)
    endt <- tab[abs(as.numeric(tab[,3]))>lfc & as.numeric(tab[,4]) < pval,]
    endt <- endt[, c(1,3,4)]
    colnames(endt) = c("geneID","log2FoldChange", "padj")
    return(endt)
  }
  
  getDASdif <- function(sigi,left,right,tpm,fc,pval){
    nt <- switch(tpm, "TPM > 0.1"="TPM-0-1","TPM > 0.3"="TPM-0-3","TPM > 1"= "TPM-1","TPM > 2" = "TPM-2" )
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
    #View(endtab)
    sigi$category <- "not differentially expressed"
    sigi$category[sigi[,1] %in% genesMP] <- kleiner0fc
    sigi$category[sigi[,1] %in% genesRP] <- groesser0fc
    commontab<-inner_join(endtab, sigi, by = c("geneID" = "V1"))
    commontab <- commontab[,-c1]
    colnames(commontab)[c2] = c("DE in","E(dPSI) per LSV junction","P(|dPSI|>=0.20) per LSV junction","Type of splicing event")
    
    fil <- sigi[,c(1,9)]
    fil<- fil[!duplicated(fil),]
    x<-fil%>%
      group_by(category)%>%
      count()
    
    events <- commontab%>%
      group_by(commontab[,colnew])%>%
      count()
    
    allevents <- sigi%>%
      group_by(V4)%>%
      count()
    
    colnames(events) = c("tose","n")
    
    
    return(list(commontab,x,events,allevents))
    
  }
  
  
#  kp_p = reactiveVal("")
  DD_tab = reactiveVal()
  DD_pie = reactiveVal()
  DD_toe = reactiveVal()
  D_alltoe = reactiveVal()  
  DAStable <- eventReactive(input$clickDAS1,{
    DASlist <- getDASdata(as.character(input$DAS1))
    DAS_left(DASlist[[2]])
    DAS_right(DASlist[[3]])
    DAS_tab(DASlist[[1]])
    dassig <- getValidEasy(DASlist[[1]],input$DASe,input$DASp)
    print(paste0("left: ",DAS_left()))
    print(paste0("right: ",DAS_right()))
    updateSelectInput(session,"DASUD",choices = c(DAS_left(),DAS_right()))
    DAS_sigtab(dassig)
    #print(DAS_sigtab())
    #DASDataTable()
    DAS_um_tab(umformen(DASlist[[1]]))
  })
  
  
  DIFDAS <- eventReactive(input$DD_click,{
    Difdas <- getDASdif(DAS_sigtab(),DAS_left(),DAS_right(),input$DASTPM,input$DASfc, input$DASpc)
    DD_tab(Difdas[[1]])
    DD_pie(Difdas[[2]])
    DD_toe(Difdas[[3]])
    D_alltoe(Difdas[[4]])
  })
  
  
  DAS_end_sig_tab = reactiveVal()#vllt reactiveValues()
  
  DASDataTable <- eventReactive(input$clickDASDT,{
    r <- NULL
    l <- NULL
    e <- NULL
    if(as.character(input$DASUD) == DAS_left()){
      l <-  DAS_sigtab()[DAS_sigtab()$V2 < 0, ]
      l <-  l[order(l$V2, decreasing = FALSE),] 
      DAS_end_sig_tab(l)
      e <- l
    }else if(as.character(input$DASUD) == DAS_right()){
      r <- DAS_sigtab()[DAS_sigtab()$V2 > 0, ]
      r <-  r[order(r$V2, decreasing = TRUE),]
      DAS_end_sig_tab(r)
      e <- r
    }
    e <- e[,c(1,2,3,4)]
    e <- e[!duplicated(e),]
    colnames(e) <- c('Gene ID','E(dPSI) per LSV junction','P(|dPSI|>=0.20) per LSV junction','Type of splicing event')
    e
    
  })
  saved_DASVP = reactiveVal()
  output$DASVP <- renderPlot({
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
    con1 <- as.data.frame(as.numeric(DAS_sigtab()[,8]))
    #print(con1)
    con2 <- as.data.frame(as.numeric(DAS_sigtab()[,7]))
    colnames(con1) = c("E(PSI)")
    colnames(con2) = c("E(PSI)")
    con1[,2] <- DAS_left()
    con2[,2] <- DAS_right()
    combined <- rbind(con1,con2)
    fig <- plot_ly(data = combined,y = ~combined[,1], type = "box", quartilemethod="linear", color = ~V2, colors = c("#c16868","#990012"))
    fig <- fig %>% layout(title = paste0("E(PSI) ", DAS_left(), " vs ", DAS_right()), yaxis = list(title = "E(PSI)",titlefont = f))
  })
  
  output$DASDT <- renderDataTable(DASDataTable(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 15))
  
  output$DASSC <- renderPlotly({
    df = data.frame(as.numeric(DAS_tab()[,5]),as.numeric(DAS_tab()[,6]))
    sigdif = data.frame(as.numeric(DAS_sigtab()[,5]),as.numeric(DAS_sigtab()[,6]))
    scattab <- getScatter(df,sigdif)
    colnames(scattab) = c("start", "end" , "junct")
    fig <- plot_ly(scattab, x = ~start, y =  ~end, type = 'scatter', mode = 'markers', symbol = ~junct, symbols = c("circle","o"),marker = list(size = 14), color = ~junct, colors = c("#e6b6b0","#990012"))
    #fig <- fig %>% add_trace(x = start1,y = ~end1, name = "significantly alternatively spliced junctions", symbol = ~end1,symbols = "o",type = "scatter", mode = "markers")
  })
  
  output$DASBP <- renderPlotly({
    oc <- count(DAS_sigtab(), var = DAS_sigtab()[,1])
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
  })
  
   DAS_colors <- c('rgb(211,94,96)','rgb(114,147,203)','rgb(128,133,133)' )
  
  output$DASDIF <- renderPlotly({
    DIFDAS()
    x<- DD_pie()
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
  
  output$DASDIFTOE <- renderPlotly({
    DIFDAS()
    colors2 <- c('rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)')
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
    fig <- fig %>% layout(title = 'Percentages of types of significant local splicing events in (differentially) expressed genes',
                          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    
    fig
    
  })
  
  output$DASDIFtab <-renderDataTable(DD_tab(),escape = FALSE, options = list(lengthMenu = c(10, 15, 20), pageLength = 8))
  
  output$DASPIE <- renderPlotly({
    DIFDAS()
    allevents <- D_alltoe()
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
  
  
  
 # --------------------homepage
  
  output$slickr <- renderSlickR({
    imgs <- list.files("www/home/", pattern=".png", full.names = TRUE)
    slickR(imgs)
  })
  output$start <- renderImage({
    list(src = "www/platlas-black.jpeg",
         contentType = 'image/jpeg',
         width = 650,
         height = 400,
         alt = "This is alternate text")
      }, deleteFile = FALSE)
    
  
}

#shinyApp(ui = ui, server = server)