library(shiny)
library(shinydashboard)
library(EnhancedVolcano)
library(gplots)
library(dplyr)    
library(viridis)
library(pheatmap)
library(biomaRt)
library(org.Hs.eg.db)
library("ggplot2")
library("networkD3")
library("plotly")
library("ggsci")
library("GOplot")
library(PCAtools)
library("DT")
library(janitor)
library(shinyWidgets)
library(slickR)
library(stringr)
library(tidyr)  ##neu installiert
#library(Gviz)   ##neu installiert
library(shinycssloaders)

ui<- dashboardPage(
  skin = "red",
  title = "Platlas",
  
  dashboardHeader(
    #title = "Platlas" ,
    tags$li(class = "dropdown",
            tags$style(".main-header {max-height: 100px}"),
            tags$style(".main-header .logo {height: 100px}"),
            tags$style(".main-header .sidebar-toggle {padding: 40px 40px}")
            
            
    ),
    title = tags$img( height = 100,width = 200,src = "NeuLogo14.PNG")#,#,#), #, "Platlas"
    #titleWidth = 300,
    
    
    
  ),
  
  dashboardSidebar(
    tags$style(".left-side, .main-sidebar {padding-top: 100px;}"),
    #tags$style(".skin-red .left-side, .skin-purple .main-sidebar, .skin-red .wrapper {
    #             background-color:#222D32CC; }"),#font-color: color:;#778899#222D32
    #tags$style(".skin-red .left-side, .skin-red .main-sidebar, .sidebar shiny-bound-input {color:#222D32;}"),#fff
    # width = 300,
    sidebarMenu(
      #tags$style(".skin-red .left-side, .skin-red .main-sidebar, .sidebar shiny-bound-input {color:#222D32;}"),
      br(), 
      menuItem(
        "Home",
        tabName = "Home",
        icon = icon("home")
        #font-color="#222D32"
        
      ),
      br(),
      menuItem(
        "Differential Expression Results",
        tabName = "DifExRes"
        #  # icon = icon("spinner")
        #   #font-color="#222D32"
      ),
      br(),
      menuItem(
        "Functional Enrichment Results",
        tabName = "GO"
        #icon = icon("address-card")#,
        
      ),
      
      br(),
      menuItem(
        "miRNA Analysis",
        tabName = "miA"
        #icon = icon("user-md")#,
        
      ),
      br(),
      menuItem(
        "circRNA Analysis",
        tabName = "cA"
        # icon = icon("ambulance")#,
        
      ),
      br(),
      menuItem(
        "De novo Network Enrichment",
        tabName = "KPM"
        #icon = icon("flask")#,
        
      ),
      br(),
      menuItem(
        "Differential Alternative Splicing",
        tabName = "DAS"
        #icon = icon("flask")#,
        
      ),
      br(),
      menuItem(
        "Integrative Analysis",
        tabName = "CA"
        #icon = icon("calendar")#,
        
      )
      
    )
  ),
  #)#)
  dashboardBody(
    #setBackgroundColor(color ="#222D32CC", shinydashboard = TRUE),
    setBackgroundColor(
      color = "ghostwhite",
      #gradient = c("linear", "radial"),
      #direction = c("bottom", "top", "right", "left"),
      shinydashboard = TRUE
    ),
    tabItems(
      tabItem(
        tabName = "Home",
        # "Home Page under construction",
        fluidRow(
          # column(4, offset = 4,imageOutput("start"))
          tabBox(
            tabPanel(
              title = span("Welcome to Platlas!", 
                           style = "color: black; font-size: 24px; font-weight: bold"),
              
              span("This platelet atlas presents the platelet transcriptome data in the subpages shown on the left.",style = "color: black; font-size: 18px"), br(), br(), br(),
              span("It especially focusses on visualizing and presenting platelet data which differentiate between reticulated platelets (RPs) and mature platelets (MPs).",style = "color: black; font-size: 18px"),br(), br(), br(),
              span("It also considers the transcriptomic differences between patients with coronary artery disease (CAD) and patients with myocardial infarction",style = "color: black; font-size: 18px"),br(), br(), br(),
            ),
            tabPanel(id = "id2", "image",
                     
                     imageOutput("patienten")
            )
            #slickROutput("slickr", width = "100%", height = "100%")
          ),
          #tabBox(
          box(
            background = "black",
            imageOutput("start")
          )
          #)
        ),
        br(),
        fluidRow(
          box(
            title = span("Platlas computes a variety of visualizations",style = "color: black; font-size: 18px;font-weight: bold"),
            span("By interactively browsing through our results, you can access different results. These results are presented by plots and data tables." ,style = "color: black; font-size: 18px"),
            
            slickROutput("slickr", width = "100%", height = "100%")
            
          )
          
          
        )
        
        
      ),

# Differential Expression Analysis ----------------------------------------

      tabItem( 
        tabName = "DifExRes",
        #tags$style('.nav-tabs-custom .nav-tabs li.active {
        #border-top-color: #d73925;
        #           }"'),
        fluidRow(
          width = 12,
          tabBox(
            title = "Filter the results",
            # The id lets us use input$tabset1 on the server to find the current tab
            
            id = "VPfilter", #height = "200px", 
            # tabPanel(id= "t1","Guide", "The Differential Expression results, browsable on this page by using different filters, are obtained from the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(), br() ,"By clicking on the first filter (Mature vs Reticulated) and filtering between the datasets of the patients, the different thresholds of read count normalizations (with TPM), p-value cutoffs and foldchange cutoffs, the differentially expressed genes of mature platelets versus reticulated platelets, which were computed by using DESeq2, visualized by downloadable plots and represented by downloadable data tables will be shown on the right side of the page.",br(),br() ,"Analogously by clicking on the second filter (CCS vs ACS), refinement of the search is possible by the same parameters as mentioned before and the differentially expressed genes of patients with Chronic coronary syndrome and the patients with Acute coronary syndrome, computed by using DESeq2, will be shown on the right side of the page, represented by downloadable data tables and plots as well."),
            tabPanel(id="t2","Mature vs. Reticulated", 
                     span("Mature Platelets vs. Reticulated Platelets",style = "color: black; font-size: 14px; font-weight: bold"),
                     selectInput(inputId = "HoD" , label = "Choose from patient datasets" , choices = c("healthy","diseased", "disease - CCS" , "disease - ACS"), multiple = FALSE),
                     selectInput(inputId = "TPM1" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c(), multiple = FALSE), #"TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2" 
                     numericInput(inputId = "pCO1",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
                     numericInput(inputId = "fcO1",value = 1.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "DVP1" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
                     #actionButton(inputId = "clickDVP1",label = "Get plot ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
                     
                     actionButton(inputId = "click1",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="t3","CCS vs. ACS", 
                     span("Diseased: Chronic coronary syndrome vs. Acute coronary syndrome",style = "color: black; font-size: 14px; font-weight: bold"),
                      selectInput(inputId = "HoD02" , label = "Choose from patient datasets" , choices = c("diseased","disease - mature platelets", "disease - reticulated platelets"), multiple = FALSE),
                     selectInput(inputId = "TPM2" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     numericInput(inputId = "pCO2",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "fcO2",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "DVP2" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
                     #actionButton(inputId = "clickDVP2",label = "Get plot ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
                     
                     actionButton(inputId = "click2",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000")
                     
                     
                     
            )
          ),
          box(
            title = span("Differential Expression Results - Guide", 
                         style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
            "The Differential Expression results, browsable on this page by using different filters, are obtained from the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(), br() ,"By clicking on the first filter (Mature vs Reticulated) and filtering between the datasets of the patients, the different thresholds of read count normalizations (with TPM), p-value cutoffs and foldchange cutoffs, the differentially expressed genes of mature platelets versus reticulated platelets, which were computed by using DESeq2, visualized by downloadable plots and represented by downloadable data tables will be shown on the bottom of the page.",br(),br() ,"Analogously by clicking on the second filter (CCS vs ACS), refinement of the search is possible by the same parameters as mentioned before and the differentially expressed genes of patients with Chronic coronary syndrome and the patients with Acute coronary syndrome, computed by using DESeq2, will be shown on the bottom of the page, represented by downloadable data tables and plots as well." )
          
          
        ),
        br(),
        fluidRow(
         # uiOutput("DIFEXtab")
          tabBox(
            width = 12,
            height = "1000px",
            title = "Get the results",
            tabPanel(
              id="tV","Volcano plot", 
              #"Volcano plot",
              br(), 
              #uiOutput("VPO")
              plotOutput("VP",height = "750px") %>% withSpinner(color= "#000000") ,
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
              dataTableOutput(outputId = "TopGenesTab") %>% withSpinner(color= "#000000"),
              tags$head(tags$style("#TopGenesTab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
              
              tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
              
            ),
            tabPanel(
              id= "tB",
              # uiOutput("BPO")
              title = "Biotypes of the Detected Genes", solidHeader = TRUE, # background = "white" ,
              plotlyOutput("biotype1",height = "800px")%>% withSpinner(color= "#000000"),
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
              plotOutput("heatmap",height = "800px") %>% withSpinner(color= "#000000"),
              br(),
              # 
              
            ),
            tabPanel(
              id= "tMA", 
              #uiOutput("HMO")
              title = "MA plot", background = "red" , solidHeader = TRUE,
              br(),
              #downloadButton("dlPNG3", "Download plot in PNG format"),
              br(),
              plotOutput("DMA_plot",height = "800px") %>% withSpinner(color= "#000000"),
              br(),
              # 
              
            ),
            tabPanel(
              id= "tBP", 
              #uiOutput("HMO")
              title = "Biotype plot", background = "red" , solidHeader = TRUE,
              br(),
              #downloadButton("dlPNG3", "Download plot in PNG format"),
              br(),
              plotOutput("Dbiotype_plot",height = "800px") %>% withSpinner(color= "#000000"),
              br(),
              # 
              
            ),
            tabPanel(
              id= "tEA", 
              #uiOutput("HMO")
              title = "Expression analysis plot", background = "red" , solidHeader = TRUE,
              br(),
              #downloadButton("dlPNG3", "Download plot in PNG format"),
              br(),
              plotOutput("Dexpression_plot",height = "800px") %>% withSpinner(color= "#000000"),
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
              plotOutput("pcabi",height = "800px") %>% withSpinner(color= "#000000"),
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
              plotOutput("pcah",height = "800px") %>% withSpinner(color= "#000000"),
              br(),
              # downloadButton("dlPNG3", "Download plot in PNG format")
              
            )
          )
        )
   
      ),
# Functional Enrichment ---------------------------------------------------
      tabItem(
        tabName = "GO",
        # "Functional Enrichment Results Page under construction",
        # tags$style('.nav-tabs-custom .nav-tabs li.active {
        # border-top-color: #d73925;
        #            }"'),
        fluidRow(
          width = 12,
          # height = "00px",
          tabBox(
            title = "Filter the results",
            # The id lets us use input$tabset1 on the server to find the current tab
            #height = "1500px",
            id = "VPfilter", #height = "200px", 
            tabPanel(id="t4","Mature vs. Reticulated", 
                     span("Mature Platelets vs. Reticulated Platelets",style = "color: black; font-size: 14px; font-weight: bold"),
                     selectInput(inputId = "HoDGO" , label = "Choose from patient datasets" , choices = c("healthy","diseased", "disease - CCS" , "disease - ACS"), multiple = FALSE),
                     selectInput(inputId = "TPMGO" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     numericInput(inputId = "pCO1GO",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "fcO1GO",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "DVP1GO" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
                     
                     selectInput(inputId = "GOSS" , label = "Choose ontology" , choices = c("biological process","cellular component", "molecular function"), multiple = FALSE),
                     #selectInput(inputId = "UDA" , label = "Choose Up/Down/All" , choices = c("upregulated","downregulated","all"), multiple = FALSE),
                     #selectInput(inputId = "AdjMethod" , label = "Choose the p-value adjustment method/correction method" , choices = c("Benjamini-Hochberg","FDR"), multiple = FALSE),
                     
                     
                     actionButton(inputId = "clickGO",label = "get Table ",icon = icon('arrow'), style = "color: #FFF; background-color:#000000; border-color: #000000")
                     
            ),
            tabPanel(id="t5","CCS vs. MI", 
                     span("Diseased: Chronic coronary syndrome vs. Myocardial infarction",style = "color: black; font-size: 14px; font-weight: bold"),
                     selectInput(inputId = "HoD2GO" , label = "Choose from patient datasets" , choices = c("diseased","disease - mature platelets", "disease - reticulated platelets"), multiple = FALSE),
                     selectInput(inputId = "TPM2GO" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     numericInput(inputId = "pCO2GO",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "fcO2GO",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "DVP2GO" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
                     
                     selectInput(inputId = "GOSS2" , label = "Choose ontology" ,choices = c("biological process","cellular component", "molecular_function"), multiple = FALSE),
                     #selectInput(inputId = "UDA2" , label = "Choose Up/Down/All" , choices = c("upregulated","downregulated","all"), multiple = FALSE),
                     #selectInput(inputId = "AdjMethod2" , label = "Choose the p-value adjustment method/correction method" , choices = c("Benjamini-Hochberg","FDR"), multiple = FALSE),
                     
                     
                     actionButton(inputId = "clickGO2",label = "get Table", style = "color: #FFF; background-color: #000000; border-color: #000000")
            )
          ),
          box(
            title = span("Functional Enrichment Results - Guide", 
                         style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
            "This page shows the Functional Enrichment Results, browsable by using two general filters that represent different approaches to the genetic analysis of platelets: The first being the differentiation of mature platelets versus reticulated platelets and the second being the comparative analysis between two diseased groups of patients - patients with Chronic coronary syndrome and patients with Acute coronary syndrome. The results are based on the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(),br(),
            "The shown results of the functional enrichment analysis were retrieved from the use of ShinyGO and are browsable by the TPM (transcript per million) read count normalization threshold of the original data, the patient datasets, the desired ontology (Biological Process, Cellular Component and Molecular function) and the GO-Entries discovered on upregulated genes or downregulated genes in both filters. We used a log2(Foldchange) cutoff of 0.8 and p-value cutoff of 0.05 to compute the differentially expressed genes",br(),br(),
            "The representation of these results was done with various plots and a data table of the discovered GO- Entries."
          )
          
        ),
        br(),
        fluidRow(
          #uiOutput("DTGO")
          tabBox(
            width = 12,
            height = "1000px",
            title = "Get the results",
            selectInput(inputId = "UDA" , label = "Choose Up/Downregulated" , choices = c("upregulated","downregulated"), multiple = FALSE, selected = "upregulated"),
            tabPanel(
              id="GOtab",title = "Top 500 GO Terms", solidHeader = TRUE,
              "",
              #selectInput(inputId = "FeA" , label = "Choose gene annotation in GO categories" , choices = c("ENSEMBL ID","ENTREZ ID","HGNC symbol"), multiple = FALSE),
              #actionButton(inputId = "clickFeA",label = "Get GOs",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
              #br(), 
              #br(),
              #selectInput(inputId = "UDA" , label = "Choose Up/Down/All" , choices = c("upregulated","downregulated"), multiple = FALSE),
              downloadButton("downloadGO", "Download table"),
              br(),
              br(),
              dataTableOutput(outputId = "GOTab")%>% withSpinner(color= "#000000"),#width = 12,
              
              # .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th
              tags$head(tags$style("#GOTab table {background-color: #DCDCDC; color : #000000}", media="screen", type="text/css")),
              tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th
                             {border-top: 1px solid #000000}"))
            ),
            tabPanel(
              title = "percentage of DE genes in GO categories", solidHeader = TRUE,# width = 12,
              br(),
              plotlyOutput("GOpG",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),
            tabPanel(
              title = "enrichment FDR per GO categories", solidHeader = TRUE,# width = 12,
              #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
              #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
              #br(),
              br(),
              plotlyOutput("GOpFE",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),
            tabPanel(
              title = "increased/decreased categories and their terms & their genes", solidHeader = TRUE,# width = 12,
              span("Increased/decreased categories and their terms in relation to upregulated and downregulated genes in the corresponding category",style = "color: black; font-size: 19px"),
              br(),
              downloadButton("dlGOC", "Download plot in PNG format"),
              br(),
              plotlyOutput("GOC",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),
            tabPanel(
              title = "GO categories and their z-score", solidHeader = TRUE,# width = 12,
              #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
              #br(),
              #withMathJax(),
              #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
              plotOutput("GOZ",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),#GOp
            tabPanel(
              title = "Top 10 most significant GO categories", solidHeader = TRUE,# width = 12,
              #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
              #br(),
              #withMathJax(),
              #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
              plotlyOutput("GOp",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),#GO_perc_term
            tabPanel(
              title = "Scatterplot- p.adjusted vs. identity percentage", solidHeader = TRUE,# width = 12,
              #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
              #br(),
              #withMathJax(),
              #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
              plotlyOutput("GO_perc_term",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),#GO_padj_logFC
            tabPanel(
              title = "Scatterplot- p.adjusted vs. gene logFC", solidHeader = TRUE,# width = 12,
              #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
              #br(),
              #withMathJax(),
              #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
              plotlyOutput("GO_padj_logFC",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            )
            # tabPanel(
            #   id= "tN", 
            #   #uiOutput("HMO")
            #   
            #   title = "GO Network", background = "red" , solidHeader = TRUE,
            #   #plotOutput(""),
            #   forceNetworkOutput("GOnet", height = "850px")%>% withSpinner(color= "#000000")
            #   #  br(),
            #   # downloadButton("dlKPMN", "Download plot in PNG format")
            #   
            # ),
            # tabPanel(
            #   id= "tH", 
            #   #uiOutput("HMO")
            #   
            #   title = "GO Hierarchy Tree", background = "red" , solidHeader = TRUE,
            #   #plotOutput(""),
            #   imageOutput("GOH", height = "850px")%>% withSpinner(color= "#000000"),
            #   br(),
            #   downloadButton("dlGOH", "Download plot in PNG format")
            #   
            # )
          )
        )
        
        
        
        
        
        
      ),
      
      tabItem(
        tabName = "miA",
        #"miRNA Analysis Page under construction"
        fluidRow(
          width = 12,
          tabBox(
            tabPanel(title ="miRNA analysis - Guide", 
                     span("miRNA analysis - Guide", style = "color: black; font-size: 18px; font-weight: bold"),br(),br(),
                     'The miRNA Analysis of platelets consists of two analytic parts, where the first one is the differential expression analysis of miRNA and the second part is the miRNA target enrichment, both only regarding the differences between mature platelets and reticulated platelets.' ,br(),br(),
                     'The miRNA differential expression analysis results were retrieved by using DESeq2 and are browsable by just choosing the data set on which the analysis should be done on, while the adjusted p-value threshold is  < 0.01 and the Log2(foldchange) threshold is set >1. This analysis is represented by a volcanoplot and a heatmap, and the significantly differentially expressed miRNAs are listed in a data table.', br(),br(),
                     'The miRNA target enrichment analysis consists of (significantly) differentially expressed miRNAs that are first mapped to their target genes using two different databases, chooseable in this browsing session, miRTarBase and miRDBase, while only genes that were found in the transcriptome analysis of platelets are chosen as target genes. Those target genes are then functionally enriched by using ShinyGO. This analysis is represented by different data tables showing the mappings of miRNAs to target genes and also target genes to functional category. The relationships of those mappings are then represented by network visualizations. The filtering to get the desired results consists of choosing the dataset in which you want to start the analysis, the gene ontology to which the desired functional categories correspond, the upregulated or downregulated miRNA of which the analysis is based, which corresponds to the miRNAs differentially expressed in mature versus reticulated platelets, and the database of which to choose the target genes, while when chosen miRTarBase one also has the option of choosing only the miRNA targets retrieved experimentally. Also the search can be filtered by choosing if the analysis should be done only on significantly differentially expressed miRNAs or not.'
                     
            ),
            tabPanel(
              id="tx",span("Differential Expression Results",style = "color: black; font-size: 14px; font-weight: bold"),
              title = "miRNA differential Expression in Mature Platelets vs. Reticulated Platelets" ,background = "black" ,solidHeader = TRUE,
              "",
              selectInput(inputId = "miF" , label = "Choose the dataset in which you want to start the analysis" , choices = c("Mature Platelets vs. Reticulated - all data", "in MI vs. CAD patients", "only CAD patients"), multiple = FALSE),
              numericInput(inputId = "mpCO",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
              numericInput(inputId = "mfCO",value = 1.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
              
              br(),
              br(),
              actionButton(inputId = "clickMIF",label = "Start Analysis", style = "color: #FFF; background-color: #000000; border-color: #FFF")
            ),
            
            tabPanel(
              id="ty",span("Target Enrichment",style = "color: black; font-size: 14px; font-weight: bold"),
              title = "miRNA target enrichment in Mature Platelets vs. Reticulated Platelets" ,background = "black" ,solidHeader = TRUE,
              "",
              selectInput(inputId = "miF2" , label = "Choose the dataset in which you want to start the analysis" , choices = c("Mature Platelets vs. Reticulated - all data", "in diseased patients", "only CAD patients"), multiple = FALSE),
              selectInput(inputId = "miOnt" , label = "Choose ontology" , choices = c("Biological Process", "Cellular Component", "Molecular Function"), multiple = FALSE),
              selectInput(inputId = "UDMI" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE),
              selectInput(inputId = "miDB" , label = "Choose database to map" , choices = c("miRTarBase", "miRDB"), multiple = FALSE),
              
              checkboxInput(inputId = "miEXP", label = "Choose miRNA Targets retrieved experimentally (miRTarBase)",value = FALSE),
              checkboxInput(inputId = "miSign", label = "Choose only significant miRNAs to map",value = FALSE),
              
              
              br(),
              br(),
              actionButton(inputId = "clickMIF2",label = "Start Analysis", style = "color: #FFF; background-color: #000000; border-color: #FFF")
              
            )
          ),
          #uiOutput("UIMI")
          tagList( 
            tabBox(
              tabPanel(
                id="t5","DE miRNA", 
                span("Significant differentially expressed miRNA",style = "color: black; font-size: 16px" ),
                br(),
                downloadButton("dlDTMI", "Download table"),
                br(),
                dataTableOutput("DTMI")%>% withSpinner(color= "#000000"),
                tags$head(tags$style("#DTMI table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
                
                tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
                
              ),
              tabPanel(
                id="t6","Volcano Plot", 
                
                plotOutput("VPMI")%>% withSpinner(color= "#000000"),
                downloadButton("dlVPMI", "Download plot in PNG format")
                
              ),
              tabPanel(
                id="t7","HeatMap", 
                
                plotOutput("HMMI")%>% withSpinner(color= "#000000"),
                downloadButton("dlHMMI", "Download plot in PNG format")
              ),
              tabPanel(
                id="t8","MA plot", 
                
                plotOutput("MAMI")%>% withSpinner(color= "#000000"),
                #downloadButton("dlHMMI", "Download plot in PNG format")
              ),
              tabPanel(
                id="t9","PCA plot", 
                
                plotOutput("PCAMI")%>% withSpinner(color= "#000000"),
                #downloadButton("dlHMMI", "Download plot in PNG format")
              ),
              tabPanel(
                id="t10","Expression plot", 
                
                plotOutput("ExpMI")%>% withSpinner(color= "#000000"),
                #downloadButton("dlHMMI", "Download plot in PNG format")
              )
            )
            
          )
          
          
        ),
        fluidRow(
          #uiOutput("MITarg")
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
        )
        
        
      ),
      tabItem(
        tabName = "cA",
        fluidRow(
          tabBox(
            title = "Filter the results",
            # The id lets us use input$tabset1 on the server to find the current tab
            
            id = "VPfilter", #height = "200px", 
            
            tabPanel(id="t2","Get differentially expressed circRNA", 
                     selectInput(inputId = "ciDS" , label = "Choose from patient datasets" , choices = c("MP vs. RP","MP vs. RP in CAD patients", "MP vs. RP in MI patients" , "CAD vs MI"), multiple = FALSE),
                     numericInput(inputId = "pCOcirc",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "fCOcirc",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
                     
                     actionButton(inputId = "clickCirc",label = "get Table ",icon = icon('arrow'), style = "color: #FFF; background-color:#000000; border-color: #000000")#,
                     
            )
          ),
          #fluidRow(
          box(
            title = span("circRNA Analysis - Guide", 
                         style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
            "This page shows the circular RNA differential expression results. The results are based on the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(),br(),
            "They are browsable by the patient datasets, while also being able to define the significance by using the desired log2(Foldchange) and p-value filters .",br(),br(),
            "The representation of these results was done with a downloadable data table and plots."
          )
        ),
        fluidRow(
          uiOutput("Crna")
        )
        #))
      ),
      tabItem(
        tabName = "KPM",
        #width =12,
        fluidRow(
          tabBox(
            # width =12, 
            tabPanel(id="t2","Mature vs. Reticulated", 
                     selectInput(inputId = "KTPM1" , label = "Choose TPM cutoff" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "KD1" , label = "Choose on which on dataset of patients the KPM-Algorithm will be performed on their significantly differentially expressed datasets" , choices = c("healthy","diseased", "diseased CCS", "diseased MI"), multiple = FALSE),
                     selectInput(inputId = "KUD1" , label = "Choose the upregulated/downregulated GOs" , choices = c("upregulated", "downregulated"), multiple = FALSE),
                     selectInput(inputId = "KO1" , label = "Choose Ontology" , choices = c("biological process", "cellular component", "molecular function"), multiple = FALSE),
                     selectInput(inputId = "KES1" , label = "Choose KPM network gene annotation" , choices = c("ENSEMBL", "gene symbol"), multiple = FALSE),
                     
                     actionButton(inputId = "clickK1",label = "get Results ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")#,
            ),     
            tabPanel(id="t3","CCS vs ACS", 
                     selectInput(inputId = "KTPM2" , label = "Choose TPM cutoff" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "KD2" , label = "Choose on which on dataset of patients the KPM-Algorithm will be performed on their significantly differentially expressed datasets" , choices = c("diseased", "diseased - mature platelets", "diseased - reticulated platelets"), multiple = FALSE),
                     selectInput(inputId = "KUD2" , label = "Choose the upregulated/downregulated GOs" , choices = c("upregulated", "downregulated"), multiple = FALSE),
                     selectInput(inputId = "KO2" , label = "Choose Ontology" , choices = c("biological process", "cellular component", "molecular function"), multiple = FALSE),
                     selectInput(inputId = "KES2" , label = "Choose KPM network gene annotation" , choices = c("ENSEMBL", "gene symbol"), multiple = FALSE),
                     
                     actionButton(inputId = "clickK2",label = "get Results ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")#,
            )
          ),
          box(
            title = span("De novo functional enrichment analysis using KeyPathwayMiner - Guide",style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
            "The de novo functional enrichment analysis is performed here by processing the dataset of each patient group's differentially expressed genes with the KeyPathwayMiner algorithm in order to identify new networks between genes that are then functionally categorized in Gene Ontologies by ShinyGO.",br(),br(),
            'The results are first categorized in two filters: The first being the differentiation of mature platelets versus reticulated platelets and the second being the comparative analysis between two diseased groups of patients - patients with Chronic coronary syndrome and patients with Acute coronary syndrome. The results are based on the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).',br(),br(),
            'The browsability of these results is ensured by being able to filter by the TPM (transcript per million) read count normalization threshold of the differentially expressed genes, the desired patient dataset and the up- or downregulated genes of this dataset, that correspond to the different sets in the before mentioned categories, the desired Ontology and the desired gene ID format.', br(), br(),
            'The representation of these results is then shown in the results by data tables that show the mapping of the genes to the functional categories, different network plots that also visualize those mappings and the relationships between the functional categories, and a Gene Ontology hierarchy plot to visualize the most important functional categories and the relationship between each other.' 
          )
        ),
        fluidRow(
          uiOutput("KPMO")
        )
      ),
      
      
      tabItem(
        tabName = "DAS",
        #width =12,
        fluidRow(
          #tabBox(
          # width =12, 
          box(id="t2",span("Configure your search",style = "color: black; font-size: 16px; font-weight: bold"), 
              selectInput(inputId = "DAS1" , label = "Choose differential Alternative Splicing Dataset" , choices = c("Reticulated Platelets vs Mature  Platelets (healthy)","Mature Platelets vs Reticulated Platelets (stable CCS)", "healthy vs stable CCS (Mature Platelets)", "healthy vs stable CCS (Reticulated Platelets)"), multiple = FALSE),
              numericInput(inputId = "DASe",value = 0.2 ,label = paste0('Choose E(','\U0394','PSI) per LSV junction cutoff'), min = 0 , max = 1, step = 0.05),
              numericInput(inputId = "DASp",value = 0.95 ,label = paste0('Choose P(|','\U0394','PSI|>=0.20) cutoff'), min = 0 , max = 1 , step = 0.05),
              
              actionButton(inputId = "clickDAS1",label = "get Results ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")#,
          ),     
          #),
          box(
            title = span("Differential Alternative Splicing - Guide",style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
            "The differential alternative splicing analysis results, that are shown on this page, were retrieved from the use of MAJIQ and Voila on the different datasets, of which the most important are browsable here. The dataset used throughout Platlas is based on the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS)." ,br(),br(),
            "This analysis shows what local splicing variations per gene are differentially spliced between mature platelets versus reticulated platelets or between healthy and diseased patients. To filter the degree of differentiality the delta E(PSI) is used, while 0,2 is the default value. This means that only local splicing variations per gene are chosen that are 20% more (less) spliced in the corresponding datasets. The P(|deltaPSI|>= 0.20) filter in this case functions as a p-value threshold, where the default value is 0.95, but the higher the number is the lower the actual p-value." ,br(),br(),
            "The results here are represented by various plots that show the differences in E(PSI) throughout the subsets of the chosen datasets, like a volcano plot adjusted to show the genes with the highest differences in percent spliced in values between mature platelets versus reticulated platelets or healthy patients versus diseased patients, a data table that lists those genes and other different general plots."
          )
        ),
        fluidRow(
          uiOutput("DASO")
        )
      ),
      
      
      
      tabItem(
        tabName = "CA",
        #"Integration Analysis Page under construction",
        fluidRow(
          width = 12,
          height = "50px",
          box(
            width = 12,
            title = span("Integrative Analysis - Guide",style = "color: white; font-size: 18px; font-weight: bold"),background = "black",solidHeader = TRUE,
            span("1.",style = "color: white; font-size: 20px; font-weight: bold"), "Get the intersection of gene sets, where you can choose between the previously shown datasets and your own list of genes (seperated by tab, a simple space or newline character)", br(),
            span("2.",style = "color: white; font-size: 20px; font-weight: bold"), "Get the intersection of a gene set, where you can choose between the previously shown datasets and your own list of genes (seperated by tab, a simple space or newline character) and miRNA Target Genes",br(),
            span("3.",style = "color: white; font-size: 20px; font-weight: bold"), "Get the intersection of a gene set and genes with significant local splicing variations",
            
          )
        ),
        # br(),
        fluidRow(
          width = 12,
          height = "20px",
          box(
            width = 3,
            column(1,span("1",style = "color: white; font-size: 24px; font-weight: bold")),
            column(2,actionButton(inputId = "Iclick1",label = "Click here to start analysis", style = "color: #fff; background-color: #000000; border-color: #fff")),
            background = "black",solidHeader = TRUE,
          )
        ),
        br(),
        fluidRow(
          width = 12,
          tabBox(
            title = "Intersection of gene sets",
            # The id lets us use input$tabset1 on the server to find the current tab
            id = "INT", #height = "200px", 
            # tabPanel(id= "t1","Guide", "The Differential Expression results, browsable on this page by using different filters, are obtained from the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(), br() ,"By clicking on the first filter (Mature vs Reticulated) and filtering between the datasets of the patients, the different thresholds of read count normalizations (with TPM), p-value cutoffs and foldchange cutoffs, the differentially expressed genes of mature platelets versus reticulated platelets, which were computed by using DESeq2, visualized by downloadable plots and represented by downloadable data tables will be shown on the right side of the page.",br(),br() ,"Analogously by clicking on the second filter (CCS vs ACS), refinement of the search is possible by the same parameters as mentioned before and the differentially expressed genes of patients with Chronic coronary syndrome and the patients with Acute coronary syndrome, computed by using DESeq2, will be shown on the right side of the page, represented by downloadable data tables and plots as well."),
            tabPanel(id="INT1","Insert Genelist 1", 
                     "Insert Genelist (EnsemblIDs separated either by space, new line or tab)",
                     textAreaInput(inputId = "IGL", label = "Your Genelist", value = "", width = "400px",height = "100px", placeholder = NULL) 
                     #actionButton(inputId = "Iclick0",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="INT2","Mature vs. Reticulated 1", 
                     "Mature Platelets vs. Reticulated Platelets",
                     selectInput(inputId = "ITPM1" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "IHoD1" , label = "Choose from patient datasets" , choices = c("healthy","diseased", "disease - stable CAD" , "disease - MI"), multiple = FALSE),
                     numericInput(inputId = "IpCO1",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
                     numericInput(inputId = "IfCO1",value = 1.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "IUDMI1" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE)#,
                     
                     #actionButton(inputId = "Iclick1",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="INT3","CCS vs. ACS 1", 
                     "Diseased: Chronic coronary syndrome vs. Myocardial infarction",
                     selectInput(inputId = "ITPM2" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "IHoD02" , label = "Choose from patient datasets" , choices = c("diseased","disease - mature platelets", "disease - reticulated platelets"), multiple = FALSE),
                     numericInput(inputId = "IpCO2",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "IfCO2",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "IUDMI2" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE)#,
                     
                     #actionButton(inputId = "Iclick0",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            )
          ),
          tabBox(
            title = "",
            # The id lets us use input$tabset1 on the server to find the current tab
            id = "INTZ", #height = "200px", 
            # tabPanel(id= "t1","Guide", "The Differential Expression results, browsable on this page by using different filters, are obtained from the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(), br() ,"By clicking on the first filter (Mature vs Reticulated) and filtering between the datasets of the patients, the different thresholds of read count normalizations (with TPM), p-value cutoffs and foldchange cutoffs, the differentially expressed genes of mature platelets versus reticulated platelets, which were computed by using DESeq2, visualized by downloadable plots and represented by downloadable data tables will be shown on the right side of the page.",br(),br() ,"Analogously by clicking on the second filter (CCS vs ACS), refinement of the search is possible by the same parameters as mentioned before and the differentially expressed genes of patients with Chronic coronary syndrome and the patients with Acute coronary syndrome, computed by using DESeq2, will be shown on the right side of the page, represented by downloadable data tables and plots as well."),
            tabPanel(id="INT4","Insert Genelist 2", 
                     "Insert Genelist (EnsemblIDs separated either by space, new line or tab)",
                     textAreaInput(inputId = "IGL2", label = "Your Genelist", value = "", width = "400px",height = "100px", placeholder = NULL) 
                     #actionButton(inputId = "Iclick0",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="INT5","Mature vs. Reticulated 2", 
                     "Mature Platelets vs. Reticulated Platelets",
                     selectInput(inputId = "ITPM3" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "IHoD3" , label = "Choose from patient datasets" , choices = c("healthy","diseased", "disease - stable CAD" , "disease - MI"), multiple = FALSE),
                     numericInput(inputId = "IpCO3",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
                     numericInput(inputId = "IfCO3",value = 1.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "IUDMI3" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE)#,
                     
                     #actionButton(inputId = "Iclick1",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="INT6","CCS vs. ACS 2", 
                     "Diseased: Chronic coronary syndrome vs. Myocardial infarction",
                     selectInput(inputId = "ITPM4" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "IHoD04" , label = "Choose from patient datasets" , choices = c("diseased","disease - mature platelets", "disease - reticulated platelets"), multiple = FALSE),
                     numericInput(inputId = "IpCO4",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "IfCO4",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "IUDMI4" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE)#,
                     
                     #actionButton(inputId = "Iclick1",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            )
          )
        ),
        br(),
        fluidRow(
          uiOutput("Integ1")
          
        ),
        br(),
        fluidRow(
          width = 12,
          height = "20px",
          box(
            width = 3,
            column(1,span("2",style = "color: white; font-size: 24px; font-weight: bold")),
            column(2,actionButton(inputId = "Iclick2",label = "Click here to start analysis", style = "color: #fff; background-color: #000000; border-color: #fff")),
            background = "black",solidHeader = TRUE,
          )
        ),
        br(),
        fluidRow(
          tabBox(
            title = "Filter miRNA gene targets",
            # The id lets us use input$tabset1 on the server to find the current tab
            id = "INTMI", #height = "200px", 
            # tabPanel(id= "t1","Guide", "The Differential Expression results, browsable on this page by using different filters, are obtained from the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(), br() ,"By clicking on the first filter (Mature vs Reticulated) and filtering between the datasets of the patients, the different thresholds of read count normalizations (with TPM), p-value cutoffs and foldchange cutoffs, the differentially expressed genes of mature platelets versus reticulated platelets, which were computed by using DESeq2, visualized by downloadable plots and represented by downloadable data tables will be shown on the right side of the page.",br(),br() ,"Analogously by clicking on the second filter (CCS vs ACS), refinement of the search is possible by the same parameters as mentioned before and the differentially expressed genes of patients with Chronic coronary syndrome and the patients with Acute coronary syndrome, computed by using DESeq2, will be shown on the right side of the page, represented by downloadable data tables and plots as well."),
            tabPanel(id="INTM1","Insert Genelist", 
                     "Insert Genelist (EnsemblIDs separated either by space, new line or tab)",
                     textAreaInput(inputId = "MGIGL", label = "Your Genelist", value = "", width = "400px",height = "100px", placeholder = NULL) 
                     #actionButton(inputId = "Iclick0",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="INTM2","Mature vs. Reticulated", 
                     "Mature Platelets vs. Reticulated Platelets",
                     selectInput(inputId = "MGITPM1" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "MGIHoD1" , label = "Choose from patient datasets" , choices = c("healthy","diseased", "disease - stable CAD" , "disease - MI"), multiple = FALSE),
                     numericInput(inputId = "MGIpCO1",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
                     numericInput(inputId = "MGIfCO1",value = 1.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "MGIUDMI1" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE)#,
                     
                     #actionButton(inputId = "Iclick1",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="INTM3","CCS vs. ACS", 
                     "Diseased: Chronic coronary syndrome vs. Myocardial infarction",
                     selectInput(inputId = "MGITPM2" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "MGIHoD02" , label = "Choose from patient datasets" , choices = c("diseased","disease - mature platelets", "disease - reticulated platelets"), multiple = FALSE),
                     numericInput(inputId = "MGIpCO2",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "MGIfCO2",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "MGIUDMI2" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE)#,
                     
                     #actionButton(inputId = "Iclick0",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            )
          ),
          tabBox(
            #title = "Filter miRNAs",
            tabPanel(
              id="ty",
              title = "miRNAs" ,background = "black" ,solidHeader = TRUE,
              "",
              selectInput(inputId = "ImiF2" , label = "Choose the dataset in which you want to start the analysis" , choices = c("Mature Platelets vs. Reticulated - all data", "in diseased patients", "only CAD patients"), multiple = FALSE),
              selectInput(inputId = "ImiOnt" , label = "Choose ontology" , choices = c("Biological Process", "Cellular Component", "Molecular Function"), multiple = FALSE),
              selectInput(inputId = "IUDMI" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE),
              selectInput(inputId = "ImiDB" , label = "Choose database to map" , choices = c("miRTarBase", "miRDB"), multiple = FALSE),
              
              numericInput(inputId = "MIpCO2",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
              numericInput(inputId = "MIfCO2",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
              
              checkboxInput(inputId = "ImiEXP", label = "Choose miRNA Targets retrieved experimentally (miRTarBase)",value = FALSE)#,
              # checkboxInput(inputId = "miSign", label = "Choose only significant miRNAs to map",value = FALSE),
              
              # br(),
              #  br(),
              #actionButton(inputId = "clickIMIF2",label = "Start Analysis", style = "color: #FFF; background-color: #000000; border-color: #FFF")
              
            )
          )
        ),
        br(),
        fluidRow(
          uiOutput("Integ2")
        ),
        br(),
        fluidRow(
          width = 12,
          height = "20px",
          box(
            width = 3,
            column(1,span("3",style = "color: white; font-size: 24px; font-weight: bold")),
            column(2, actionButton(inputId = "Iclick3",label = "Click here to start analysis", style = "color: #fff; background-color: #000000; border-color: #fff")),
            background = "black",solidHeader = TRUE
            
          )
        ),
        br(),
        fluidRow(
          tabBox(
            title = "LSVs in gene sets",
            # The id lets us use input$tabset1 on the server to find the current tab
            id = "INTD", #height = "200px", 
            # tabPanel(id= "t1","Guide", "The Differential Expression results, browsable on this page by using different filters, are obtained from the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(), br() ,"By clicking on the first filter (Mature vs Reticulated) and filtering between the datasets of the patients, the different thresholds of read count normalizations (with TPM), p-value cutoffs and foldchange cutoffs, the differentially expressed genes of mature platelets versus reticulated platelets, which were computed by using DESeq2, visualized by downloadable plots and represented by downloadable data tables will be shown on the right side of the page.",br(),br() ,"Analogously by clicking on the second filter (CCS vs ACS), refinement of the search is possible by the same parameters as mentioned before and the differentially expressed genes of patients with Chronic coronary syndrome and the patients with Acute coronary syndrome, computed by using DESeq2, will be shown on the right side of the page, represented by downloadable data tables and plots as well."),
            tabPanel(id="INTD1","Insert Genelist", 
                     "Insert Genelist (EnsemblIDs seperated either by space, new line or tab)",
                     textAreaInput(inputId = "DSIGL", label = "Your Genelist", value = "", width = "400px",height = "100px", placeholder = NULL) 
                     #actionButton(inputId = "Iclick0",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="INTD2","Mature vs. Reticulated", 
                     "Mature Platelets vs. Reticulated Platelets",
                     selectInput(inputId = "DSITPM1" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "DSIHoD1" , label = "Choose from patient datasets" , choices = c("healthy","diseased", "disease - stable CAD" , "disease - MI"), multiple = FALSE),
                     numericInput(inputId = "DSIpCO1",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
                     numericInput(inputId = "DSIfCO1",value = 1.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
                     #selectInput(inputId = "DSIUDMI1" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE)#,
                     
                     #actionButton(inputId = "Iclick1",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            ),
            tabPanel(id="INTD3","CCS vs. ACS", 
                     "Diseased: Chronic coronary syndrome vs. Myocardial infarction",
                     selectInput(inputId = "DSITPM2" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     selectInput(inputId = "DSIHoD02" , label = "Choose from patient datasets" , choices = c("diseased","disease - mature platelets", "disease - reticulated platelets"), multiple = FALSE),
                     numericInput(inputId = "DSIpCO2",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "DSIfCO2",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
                     #selectInput(inputId = "DSIUDMI2" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE)#,
                     
                     #actionButton(inputId = "Iclick0",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000")
                     
            )
          ),
          tabBox(
            #title = "LSVs in gene sets",
            tabPanel(
              id="t2",#span("Configure your search",style = "color: black; font-size: 16px; font-weight: bold"), 
              #"Filter LSVs",
              title = "Filter LSVs",
              selectInput(inputId = "IDAS1" , label = "Choose differential Alternative Splicing Dataset" , choices = c("Reticulated Platelets vs Mature  Platelets (healthy)","Mature Platelets vs Reticulated Platelets (stable CCS)", "healthy vs stable CCS (Mature Platelets)", "healthy vs stable CCS (Reticulated Platelets)"), multiple = FALSE),
              numericInput(inputId = "IDASe",value = 0.2 ,label = paste0('Choose E(','\U0394','PSI) per LSV junction cutoff'), min = 0 , max = 1, step = 0.05),
              numericInput(inputId = "IDASp",value = 0.95 ,label = paste0('Choose P(|','\U0394','PSI|>=0.20) cutoff'), min = 0 , max = 1 , step = 0.05),
              
            )
            
          )),
        br(),
        fluidRow(
          uiOutput("Integ3")
          
        ),
      )
    )
  )
)
