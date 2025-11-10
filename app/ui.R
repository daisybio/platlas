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
library(stringr)
library(tidyr)  ##neu installiert
library(shinyjs)
library(tidyverse)
#library(Gviz)   ##neu installiert
library(shinycssloaders)
library(visNetwork)
library(BiocManager)
library(BiocVersion)
#library(shinydlplot)
library(periscope)
library(orca)
library(readxl)
library(shinyjs)
library(bslib)
#library(httr)
#library(jsonlite)
#library(htmltools)

#all input data
#Difex -> HoD
#GO -> HoDGO
#miA -> miF, miF2
#cA -> ciDS, ciDS2
#DAS -> DAS1

css <- "
.nav li a.disabled {
background-color: #aaa !important;
color: #333 !important;
cursor: not-allowed !important;
border-color: #aaa !important;
}"
# ---- study link + local PDF (served from www/) ----
study_doi <- "10.1093/eurheartj/ehaf694"                 # <-- set your DOI
study_doi_url <- paste0("https://doi.org/", study_doi)
#study_pdf_rel <- "EHJ_paper.pdf"          # file at: www/papers/platlas_study.pdf

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
    shinyjs::useShinyjs(),
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
      # br(),
      # menuItem(
      #   "De novo Network Enrichment",
      #   tabName = "KPM"
      #   #icon = icon("flask")#,
      # 
      # ),
      br(),
      menuItem(
        "Differential Alternative Splicing",
        tabName = "DAS"
        #icon = icon("flask")#,
        
      )#,
      # br(),
      # menuItem(
      #   "Integrative Analysis",
      #   tabName = "CA"
      #   #icon = icon("calendar")#,
      #   
      # )
      
    )
  ),
  #)#)
  dashboardBody(
    #setBackgroundColor(color ="#222D32CC", shinydashboard = TRUE),
    shinyjs::useShinyjs(),
    shinyjs::inlineCSS(css),
    setBackgroundColor(
      color = "ghostwhite",
      #gradient = c("linear", "radial"),
      #direction = c("bottom", "top", "right", "left"),
      shinydashboard = TRUE
    ),
    tags$head(
      tags$style(
        HTML(
          ".modebar-btn[data-title='Download plot as a pdf'] {display: block !important;}"
        )
      )
    ),
    tabItems(
      tabItem(
        tabName = "Home",
        # "Home Page under construction",
        fluidRow(
          # column(4, offset = 4,imageOutput("start"))
          tabBox(
            width = 12,
            tabPanel(
              title = span("Welcome to Platlas!", 
                           style = "color: black; font-size: 24px; font-weight: bold"),
              
              span("Platlas is a comprehensive platelet atlas designed to provide access to platelet transcriptome data presented in our study. This platform focuses on visualizing and presenting data that differentiates between reticulated platelets (RPs) and mature platelets (MPs). It offers valuable insights into the transcriptomic differences observed in patients with coronary artery disease (CAD).",style = "color: black; font-size: 18px"),
              tags$a(href = study_doi_url, target = "_blank", rel = "noopener",
                     style = "display:inline-block;padding:.3rem .6rem;border:1px solid #ddd;border-radius:999px;color:#333;text-decoration:none;",
                     tagList(icon("link"), span(" DOI "), tags$code(study_doi))),
              #span("This platelet atlas presents the platelet transcriptome data in the subpages shown on the left.",style = "color: black; font-size: 18px"), br(), br(), br(),
              #span("It especially focusses on visualizing and presenting platelet data which differentiate between reticulated platelets (RPs) and mature platelets (MPs).",style = "color: black; font-size: 18px"),br(), br(), br(),
              #span("It also considers the transcriptomic differences between patients with coronary artery disease (CAD) and patients with myocardial infarction",style = "color: black; font-size: 18px"),br(), br(), br(),
            )#,
#            tabPanel(id = "id2", "image",
#                     imageOutput("patienten")
#            )
            #slickROutput("slickr", width = "100%", height = "100%")
          )#,
          #tabBox(

          #)
        ),
        br(),
        fluidRow(
          box( 
            width = 12,
            #background = "black",
           # column(width = 7, imageOutput("start")),
          #  column(width = 5,  span("Platlas computes a variety of visualizations",style = "color: white; font-size: 18px;font-weight: bold"),br(), br(), br(),
           #        span("By interactively browsing through our results, you can access different results. These results are presented by plots and data tables." ,style = "color: white; font-size: 18px"))
            #imageOutput("start")
          column(width = 12, imageOutput("summary")),#height = 100,width = 200,
          span("Created with " ,style = "color: black; font-size: 14px"),
          a("bioRender", href="https://www.biorender.com")
          
        )
        )
        ,
        br(),
         fluidRow(
           box(
             width = 12,
             title = span("Cite Our Study:", style = "color: black; font-size: 18px; font-weight: bold"),
             
             # Title
             div(style = "font-size:18px; font-weight:600; color:black;",
                 "Reticulated platelets in coronary artery disease: a multidimensional approach unveils prothrombotic signalling and novel therapeutic targets"
             ),
             
             # Journal meta + DOI + Published date (NO 'Free' badge)
             div(style = "margin-top:8px; color:#444;",
                 span("European Heart Journal, ehaf694 · "),
                 tags$a(href = study_doi_url, target = "_blank", rel = "noopener", study_doi_url),
                 span(" · Published: 31 August 2025")
             ),
             
             # Authors
             div(style = "margin-top:8px; color:#444;",
                 "Kilian Kirmes, Jiaying Han, Melissa Klug, Conor J Bloxham, Olena Babyak, Judith Bernett, ",
                 "Lis Arend, Quirin Manz, Leonora Raka, Leon Schwartz, Markus Hoffmann, Marc Rosenbaum, ",
                 "Jürgen Ruland, Octavia-Andreea Ciora, Zakaria Louadi, Olga Tsoy, Khalique Newaz, ",
                 "Jessica Modica, Carola Conca Dioguardi, Clelia Peano, Michaela Müller, Donato Santovito, ",
                 "Giacomo Viggiani, Stephanie Kühne, Moritz von Scheidt, Leo Nicolai, Tianjiao Wu, ",
                 "Jan Baumbach, Mauro Chiarito, Karl-Ludwig Laugwitz, Gianluigi Condorelli, ",
                 "Philip W J Raake, Markus List, Isabell Bernlochner, Dario Bongiovanni"
             ),
             
             tags$hr(),
             
             # How to cite (exact text you provided)
             
             # Actions: Quick preview (local PDF), Download PDF, Open DOI
             div(class = "btn-group",
                 #actionButton("home_preview_pdf",
                #              label = tagList(icon("file-pdf-o"), " Quick preview"),
                 #             class = "btn btn-danger"),
                 #tags$a(href = study_pdf_rel,
                #        class = "btn btn-default",
                #        download = "Platlas_study.pdf",
                #        "Download PDF"),
                 tags$a(href = study_doi_url,
                        target = "_blank", rel = "noopener",
                        class = "btn btn-primary",
                        tagList(icon("external-link"), " Open DOI"))
             ),
             tags$br(),
           )
           
           
         ),
  br(),
        
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
            
            id = "VPfilter1", #height = "200px", 
            # tabPanel(id= "t1","Guide", "The Differential Expression results, browsable on this page by using different filters, are obtained from the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).",br(), br() ,"By clicking on the first filter (Mature vs Reticulated) and filtering between the datasets of the patients, the different thresholds of read count normalizations (with TPM), p-value cutoffs and foldchange cutoffs, the differentially expressed genes of mature platelets versus reticulated platelets, which were computed by using DESeq2, visualized by downloadable plots and represented by downloadable data tables will be shown on the right side of the page.",br(),br() ,"Analogously by clicking on the second filter (CCS vs ACS), refinement of the search is possible by the same parameters as mentioned before and the differentially expressed genes of patients with Chronic coronary syndrome and the patients with Acute coronary syndrome, computed by using DESeq2, will be shown on the right side of the page, represented by downloadable data tables and plots as well."),
            tabPanel(id="DIFEX_t1","Reticulated Platelets vs. Mature Platelets", 
                     span("Reticulated Platelets vs. Mature Platelets",style = "color: black; font-size: 14px; font-weight: bold"),br(),
                     span("Read count normalization: TPM (Transcripts per million) threshhold > 0.2",style = "color: black; font-size: 14px"),
                    #splitLayout(
                     # pickerInput(inputId = "HoD" , label = "Choose from patient datasets" , choices = c("disease - CCS","healthy","diseased" , "disease - ACS"), multiple = FALSE,choicesOpt = list(
                      #  disabled = c("disease - CCS","healthy","diseased" , "disease - ACS") %in% c("healthy","diseased" , "disease - ACS")
                      #)),
                      #selectInput(inputId = "TPM1" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c(), multiple = FALSE), #"TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2" 
                     #),
                     #splitLayout(
                      numericInput(inputId = "pCO1",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
                      numericInput(inputId = "fcO1",value = 1.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
                    # ),
                     selectInput(inputId = "DVP1" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
                     #actionButton(inputId = "clickDVP1",label = "Get plot ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
                     
                     actionButton(inputId = "click1",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
                    # actionButton(inputId = "InfoDIFEX",label = "Start Analysis",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")
            )#,
            # tabPanel(id="DIFEX_t2","CCS vs. MI", disabled = TRUE,
            #          span("Diseased: Chronic coronary syndrome vs. Acute coronary syndrome",style = "color: black; font-size: 14px; font-weight: bold"),
            #           selectInput(inputId = "HoD02" , label = "Choose from patient datasets" , choices = c("diseased","disease - mature platelets", "disease - reticulated platelets"), multiple = FALSE),
            #          selectInput(inputId = "TPM2" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
            #          numericInput(inputId = "pCO2",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
            #          numericInput(inputId = "fcO2",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
            #          selectInput(inputId = "DVP2" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
            #          #actionButton(inputId = "clickDVP2",label = "Get plot ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
            #          
            #          actionButton(inputId = "click2",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000")
            #          
            #          
            #          
            # )
          ),
          box(
            title = span("Differential Expression Results - Guide", 
                         style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
            "The Differential Expression results, accessible on this page using various filters, are derived from transcriptome analyses of reticulated and mature platelets",br(), br() ,"By applying p-value and fold-change cutoffs, you can explore the differentially expressed genes between mature and reticulated platelets. These results, computed using DESeq2, are presented through a range of downloadable visualizations and data tables.",br(),br() ,"Additionally, the platform provides the ability to visualize volcano plots, a comprehensive list of differentially expressed genes, biotypes of detected genes, MA plots, and other detailed analyses to deepen your understanding of platelet transcriptomic differences." )
          
          
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
              #downloadablePlotUI("VP_obj",downloadtypes = c("png", "csv"),download_hovertext = "Download the plot and data here!",height = "500px",btn_halign = "left")
              downloadButton("dl_VP_DIFEX1", "Download plot in PNG format"),#dlPNG
              downloadButton("dl_VP_DIFEX2", "Download plot in PDF format")#dlPNG
              
            ),
            tabPanel(
              id= "tT", "Differentially Expressed Genes",
              # uiOutput("DTO")
              selectInput(inputId = "UpDown" , label = "Choose up or downregulated genes" , choices = c("upregulated","downregulated"), multiple = FALSE),
              actionButton(inputId = "click4",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
              #downloadButton("downloadAll", "Download table"),
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
              tags$head(
                tags$style(
                  HTML(
                    ".modebar-btn[data-title='Download plot as a pdf'] {display: block !important;}"
                  )
                )
              ),
              title = "Biotypes of the Detected Genes", solidHeader = TRUE, # background = "white" ,
              plotlyOutput("biotype1",height = "800px")%>% withSpinner(color= "#000000"),
              br(),
              #downloadButton("", "Download plot in PNG format")dlPNG2
            ),
            tabPanel(
              id= "tH", 
              #uiOutput("HMO")
              title = "Counts per samples per genes", background = "red" , solidHeader = TRUE,
              br(),
              br(),
              plotOutput("heatmap",height = "800px") %>% withSpinner(color= "#000000"),
              br(),
              # 
              downloadButton("dl_heatmap_DIFEX1", "Download plot in PNG format"),#dlPNG3
              downloadButton("dl_heatmap_DIFEX2", "Download plot in PDF format"),#dlPNG3
              
            ),
            tabPanel(
              id= "tMA", 
              #uiOutput("HMO")
              title = "MA plot", background = "red" , solidHeader = TRUE,
              br(),
              br(),
              plotOutput("DMA_plot",height = "800px") %>% withSpinner(color= "#000000"),
              br(),
              # 
              downloadButton("dl_MA_DIFEX1", "Download plot in PNG format"),
              downloadButton("dl_MA_DIFEX2", "Download plot in PDF format"),
            ),
            # tabPanel(
            #   id= "tBP", 
            #   #uiOutput("HMO")
            #   title = "Biotype plot", background = "red" , solidHeader = TRUE,
            #   br(),
            #   br(),
            #   plotOutput("Dbiotype_plot",height = "800px") %>% withSpinner(color= "#000000"),
            #   br(),
            #   # 
            #   downloadButton("dl_BT_DIFEX1", "Download plot in PNG format"),
            #   downloadButton("dl_BT_DIFEX2", "Download plot in PDF format"),
            #   
            # ),
            tabPanel(
              id= "tEA", 
              #uiOutput("HMO")
              title = "Expression analysis plot", background = "red" , solidHeader = TRUE,
              br(),
              br(),
              plotOutput("Dexpression_plot",height = "800px") %>% withSpinner(color= "#000000"),
              br(),
              # 
              downloadButton("dl_EXP_DIFEX1", "Download plot in PNG format"),
              downloadButton("dl_EXP_DIFEX2", "Download plot in PDF format"),
              
            ),
            tabPanel(
              id= "tP", 
              #uiOutput("HMO")
              title = "PCA Analysis Biplot", background = "red" , solidHeader = TRUE,
              br(),
              br(),
              plotOutput("pcabi",height = "800px") %>% withSpinner(color= "#000000"),
              br(),
              # downloadButton("dlPNG3", "Download plot in PNG format")
              downloadButton("dl_PCA_DIFEX1", "Download plot in PNG format"),
              downloadButton("dl_PCA_DIFEX2", "Download plot in PDF format"),
              
            ),
            tabPanel(
              # uiOutput("DTO")
              title = "Exclusive Genes",
              selectInput(inputId = "MPXRP" , label = "Choose genes only expressed in MPs or RPs (unfiltered)" , choices = c("expressed only in MPs","expressed only in RPs"), multiple = FALSE),
              br(),
              br(),
              dataTableOutput(outputId = "DIFXlist") %>% withSpinner(color= "#000000"),
              tags$head(tags$style("#DIFXlist table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
              
              tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
              # downloadButton("dlPNG3", "Download plot in PNG format")

            ),
            tabPanel(
              title = "Read counts",
              selectInput(inputId = "UpDown2" , label = "Choose up or downregulated genes" , choices = c("upregulated","downregulated"), multiple = FALSE),
              uiOutput("select_Genes")%>% withSpinner(color = "#222D32"),
              
              #tags$head(tags$style("#DIFVIOtab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
              
              #tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
              #                 border-top: 1px solid #000000;}",media="screen", type="text/css")),
              plotlyOutput("DIFVIO",height = "600px")%>% withSpinner(color= "#000000"),
              dataTableOutput(outputId = "DIFVIOtab")%>% withSpinner(color= "#000000")
              
              
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
            id = "VPfilter2", #height = "200px", 
            tabPanel(id="FEX_t1","Reticulated Platelets vs. Mature Platelets", 
                     span("Reticulated Platelets vs. Mature Platelets",style = "color: black; font-size: 14px; font-weight: bold"),br(),
                     span("Read count normalization: TPM (Transcripts per million) threshhold > 0.2",style = "color: black; font-size: 14px"),
            
                     #pickerInput(inputId = "HoDGO" , label = "Choose from patient datasets" , choices = c("disease - CCS","healthy","diseased",  "disease - ACS"), multiple = FALSE,choicesOpt = list(
                     #  disabled = c("disease - CCS","healthy","diseased" , "disease - ACS") %in% c("healthy","diseased" ,"disease - ACS")
                     #)),
                    # selectInput(inputId = "TPMGO" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
                     #selectInput(inputId = "TPMGO" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c(), multiple = FALSE), #"TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2" 
                     
                     selectInput(inputId = "DVP1GO" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
                     numericInput(inputId = "pCO1GO",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "fcO1GO",value =  0.5 , label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
                     
                     selectInput(inputId = "GOSS" , label = "Choose ontology" , choices = c("biological process","cellular component", "molecular function"), multiple = FALSE),
                     #selectInput(inputId = "UDA" , label = "Choose Up/Down/All" , choices = c("upregulated","downregulated","all"), multiple = FALSE),
                     #selectInput(inputId = "AdjMethod" , label = "Choose the p-value adjustment method/correction method" , choices = c("Benjamini-Hochberg","FDR"), multiple = FALSE),
                     
                     
                     actionButton(inputId = "clickGO",label = "get Table ",icon = icon('arrow'), style = "color: #FFF; background-color:#000000; border-color: #000000")
                     
            )#,
            # tabPanel(id="FEX_t2","CCS vs. MI", disabled = TRUE,
            #          span("Diseased: Chronic coronary syndrome vs. Myocardial infarction",style = "color: black; font-size: 14px; font-weight: bold"),
            #          selectInput(inputId = "HoD2GO" , label = "Choose from patient datasets" , choices = c("diseased","disease - mature platelets", "disease - reticulated platelets"), multiple = FALSE),
            #          selectInput(inputId = "TPM2GO" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
            #          selectInput(inputId = "DVP2GO" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
            # 
            #          numericInput(inputId = "pCO2GO",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
            #          numericInput(inputId = "fcO2GO",value =  1.5 , label = "Foldchange cutoff", min = 0 , max = 4.0, step = 0.5),
            #          
            #          selectInput(inputId = "GOSS2" , label = "Choose ontology" ,choices = c("cellular component", "biological process","molecular_function"), multiple = FALSE),
            #          #selectInput(inputId = "UDA2" , label = "Choose Up/Down/All" , choices = c("upregulated","downregulated","all"), multiple = FALSE),
            #          #selectInput(inputId = "AdjMethod2" , label = "Choose the p-value adjustment method/correction method" , choices = c("Benjamini-Hochberg","FDR"), multiple = FALSE),
            #          
            #          
            #          actionButton(inputId = "clickGO2",label = "get Table", style = "color: #FFF; background-color: #000000; border-color: #000000")
            # )
          ),
          box(
            title = span("Functional Enrichment Results - Guide", 
                         style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
            "This page presents the Functional Enrichment Results, which can be used to explore the functional enrichment differences between mature platelets and reticulated platelets.",br(),br(),
            "The functional enrichment analysis was performed using the enrichGO() function from the clusterProfiler package. The results are fully customizable, allowing users to browse patient datasets, adjust criteria for differential gene expression (p-value and fold-change cutoffs), and select the desired ontology (Biological Process, Cellular Component, or Molecular Function). Gene annotation outputs are also configurable (ENSEMBL ID or HGNC).",br(),br(),
            "The results are represented through a variety of visualizations, including detailed plots and a data table of the identified Gene Ontology (GO) entries, providing a comprehensive view of platelet-related functional enrichments as well as Gene Set Enrichment Analysis (GSEA)."
          )
          
        ),
        br(),
        fluidRow(
          #uiOutput("DTGO")
          tabBox(
            width = 12,
            height = "1000px",
            title = "Get the results",
            #maybe regulate this
            selectInput(inputId = "UDA" , label = "Choose Up/Downregulated" , choices = c("upregulated","downregulated"), multiple = FALSE, selected = "upregulated"),
            #slider only gos of interest
            conditionalPanel("input.GOSS == 'biological process' || input.GOSS2 == 'biological process'",
              prettySwitch(
              "gos_of_interest",
              label = "Plot only the GO Terms of interest (containing the words platelet, hemostasis and coagulation):",
              value = FALSE,
              status = "default",
              slim = FALSE,
              fill = TRUE,
              bigger = FALSE,
              inline = FALSE,
              width = NULL
            )            
            ),
           
            tabPanel(
              id="GOtab",title = "GO Terms", solidHeader = TRUE,
              "",
              #selectInput(inputId = "FeA" , label = "Choose gene annotation in GO categories" , choices = c("ENSEMBL ID","ENTREZ ID","HGNC symbol"), multiple = FALSE),
              #actionButton(inputId = "clickFeA",label = "Get GOs",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
              #br(), 
              #br(),
              #selectInput(inputId = "UDA" , label = "Choose Up/Down/All" , choices = c("upregulated","downregulated"), multiple = FALSE),
              #downloadButton("downloadGO", "Download table"),
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
              title = "GO categories with the 10 highest DE genes", solidHeader = TRUE,# width = 12,
              #span("Increased/decreased categories and their terms in relation to upregulated and downregulated genes in the corresponding category",style = "color: black; font-size: 19px"),
              br(),
              downloadButton("dlGOC", "Download plot in PNG format"),
              br(),
              plotlyOutput("GOC",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),
            tabPanel(
              title = "GO category graph", solidHeader = TRUE,# width = 12,
              #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
              #br(),
              #withMathJax(),
              #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
              plotOutput("GOZ",height = "800px")%>% withSpinner(color= "#000000"),
              downloadButton("dl_graph_FE1", "Download plot in PNG format"),
              downloadButton("dl_graph_FE2", "Download plot in PDF format"),
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
            ),
            tabPanel(
              title = "Gene Set Enrichment Analysis Plot", solidHeader = TRUE,# width = 12,
              #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
              #br(),
              #withMathJax(),
              #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
              uiOutput("select_GSEA") %>% withSpinner(color = "#222D32"),
              
              plotOutput("GO_GSEA",height = "600px")%>% withSpinner(color= "#000000"),
              downloadButton("dl_graph_GSEA1", "Download plot in PNG format"),
              downloadButton("dl_graph_GSEA2", "Download plot in PDF format"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),
            tabPanel(
              title = "Leading Edge Genes", solidHeader = TRUE,# width = 12,
              #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
              #br(),
              #withMathJax(),
              #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
              uiOutput("select_GSEA2") %>% withSpinner(color = "#222D32"),
              prettySwitch(
                "show_gene_names",
                label = "Show gene names on the y-axis",
                value = FALSE,
                status = "default",
                slim = FALSE,
                fill = TRUE,
                bigger = FALSE,
                inline = FALSE,
                width = NULL
              ),
              
              plotOutput("GSEA_heat",height = "600px")%>% withSpinner(color= "#000000"),
              downloadButton("dl_graph_GSEA_hm_1", "Download plot in PNG format"),
              downloadButton("dl_graph_GSEA_hm_2", "Download plot in PDF format"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            ),
            tabPanel(
              id="GSEAtab",title = "GSEA results", solidHeader = TRUE,
              "",
              #selectInput(inputId = "FeA" , label = "Choose gene annotation in GO categories" , choices = c("ENSEMBL ID","ENTREZ ID","HGNC symbol"), multiple = FALSE),
              #actionButton(inputId = "clickFeA",label = "Get GOs",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
              #br(), 
              #br(),
              #selectInput(inputId = "UDA" , label = "Choose Up/Down/All" , choices = c("upregulated","downregulated"), multiple = FALSE),
              #downloadButton("downloadGO", "Download table"),
              br(),
              br(),
              dataTableOutput(outputId = "GSEATab")%>% withSpinner(color= "#000000"),#width = 12,
              
              # .table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th
              tags$head(tags$style("#GSEATab table {background-color: #DCDCDC; color : #000000}", media="screen", type="text/css")),
              tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th
                             {border-top: 1px solid #000000}"))
            ),
            tabPanel(
              title = "circular GO plot", solidHeader = TRUE,# width = 12,
              #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
              #br(),
              #withMathJax(),
              #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
              #uiOutput("select_GSEA") %>% withSpinner(color = "#222D32"),
              
              plotOutput("circGOplot",height = "800px", width = "800px")%>% withSpinner(color= "#000000"),
              downloadButton("dl_graph_circGOplot1", "Download plot in PNG format"),
              downloadButton("dl_graph_circGOplot2", "Download plot in PDF format"),
              #downloadButton("dlPNGGO", "Download plot in PNG format")
            )
            # tabPanel(
            #   title = "GO barplot", solidHeader = TRUE,# width = 12,
            #   numericInput(inputId = "padjGO",value = 0.05 , label = "adjusted p-value cutoff", min = 1e-15 , max = 0.05, step = 0.01),
            #   numericInput(inputId = "fdrGO",value =  0.01 , label = "FDR cutoff", min = 0.00 , max = 1, step = 0.05),
            #   numericInput(inputId = "countGO",value =  2 , label = "gene count cutoff", min = 0 , max = 20.0, step = 5),
            #   
            #   #plotOutput("GObar",height = "800px")%>% withSpinner(color= "#000000"),
            #   #downloadButton("dlPNGGO", "Download plot in PNG format")
            # )

          )
        )
        
        
        
        
        
        
      ),
      

# miRNA Analysis ----------------------------------------------------------

      tabItem(
        tabName = "miA",
        #"miRNA Analysis Page under construction"
        fluidRow(
          width = 12,
          tabBox(
            tabPanel(title ="miRNA analysis - Guide", 
                     span("miRNA analysis - Guide", style = "color: black; font-size: 18px; font-weight: bold"),br(),br(),
                     'The platelet miRNA analysis consists of two analytical parts, the first is the miRNA differential expression analysis and the second is the miRNA target enrichment, both only regarding the differences between mature and reticulated platelets.' ,br(),br(),
                     'The results of the miRNA differential expression analysis were obtained using DESeq2 and can be browsed by simply selecting the adjusted p-value threshold and the Log2 fold-change threshold. This analysis is represented by a volcano plot and a heatmap, and the significantly differentially expressed miRNAs are listed in a data table.', br(),br(),
                     'The miRNA target enrichment analysis consists of (significantly) differentially expressed miRNAs that are first mapped to their target genes using two different databases, selectable in this browsing session, miRTarBase and miRDBase, where only genes found in the platelet transcriptome analysis are selected as target genes. These target genes are then functionally enriched using ShinyGO. This analysis is represented by different data tables showing the mappings of miRNAs to target genes and also target genes to functional category. The relationships of these mappings are then represented by network visualisations. Filtering to obtain the desired results consists of the gene ontology to which the desired functional categories correspond, the up- or downregulated miRNA on which the analysis is based, which corresponds to the miRNAs differentially expressed in mature versus reticulated platelets, and the database from which to select the target genes, although when miRTarBase is selected it is also possible to select only the experimentally retrieved miRNA targets. The search can also be filtered by selecting whether or not the analysis should only be performed on significantly differentially expressed miRNAs.'
                     
            ),
            tabPanel(
              id="tx",span("Differential Expression Results",style = "color: black; font-size: 14px; font-weight: bold"),
              title = "miRNA differential Expression in Reticulated Platelets vs. Mature Platelets" ,background = "black" ,solidHeader = TRUE,
              "",
             # selectInput(inputId = "miF" , label = "Choose the dataset in which you want to start the analysis" , choices = c("Mature Platelets vs. Reticulated - all data", "in MI vs. CAD patients", "only CAD patients"), multiple = FALSE),
              #pickerInput(inputId = "miF" , label = "Choose the dataset in which you want to start the analysis" , choices = c("only CAD patients","Mature Platelets vs. Reticulated - all data","in MI vs. CAD patients"), multiple = FALSE,choicesOpt = list(
              #  disabled = c("only CAD patients","Mature Platelets vs. Reticulated - all data","in MI vs. CAD patients") %in% c("Mature Platelets vs. Reticulated - all data","in MI vs. CAD patients")
              #)),
              numericInput(inputId = "mpCO",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
              numericInput(inputId = "mfCO",value = 1 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
              
              br(),
              br(),
              actionButton(inputId = "clickMIF",label = "Start Analysis", style = "color: #FFF; background-color: #000000; border-color: #FFF")
            ),
            
            tabPanel(
              id="ty",span("Target Enrichment",style = "color: black; font-size: 14px; font-weight: bold"),
              title = "miRNA target enrichment in Reticulated Platelets vs. Mature Platelets" ,background = "black" ,solidHeader = TRUE,
              "",
              #pickerInput(inputId = "miF2" , label = "Choose the dataset in which you want to start the analysis" , choices = c("only CAD patients","Mature Platelets vs. Reticulated - all data","in MI vs. CAD patients"), multiple = FALSE,choicesOpt = list(
              #  disabled = c("only CAD patients","Mature Platelets vs. Reticulated - all data","in MI vs. CAD patients") %in% c("Mature Platelets vs. Reticulated - all data","in MI vs. CAD patients")
              #)),
             # selectInput(inputId = "miF2" , label = "Choose the dataset in which you want to start the analysis" , choices = c("Mature Platelets vs. Reticulated - all data", "in diseased patients", "only CAD patients"), multiple = FALSE),
              selectInput(inputId = "miOnt" , label = "Choose ontology" , choices = c("Biological Process", "Cellular Component", "Molecular Function"), multiple = FALSE),
             # selectInput(inputId = "UDMI" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE),
             selectInput(inputId = "miAnnot" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
             
             numericInput(inputId = "mTpCO",value = 0.05 ,label = "p-value cutoff", min = 0.01 , max = 0.05,step = 0.005),
             numericInput(inputId = "mTfCO",value = 1 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
             
             selectInput(inputId = "miDB" , label = "Choose database to map" , choices = c("miRTarBase", "miRDB"), multiple = FALSE),
              
              checkboxInput(inputId = "miEXP", label = "Choose miRNA Targets retrieved experimentally (miRTarBase)",value = TRUE),
              checkboxInput(inputId = "miSign", label = "Choose only significant miRNAs to map",value = TRUE),
              
              
              br(),
              br(),
              actionButton(inputId = "clickMIF2",label = "Start Analysis", style = "color: #FFF; background-color: #000000; border-color: #FFF"),
              shinyjs::hidden(p(id = "textclickMIF2", "you have to check at least one checkbox!"))
            )
          ),
          #uiOutput("UIMI")
          tagList( 
            tabBox(
              tabPanel(
                id="t5","DE miRNA", 
                span("Significant differentially expressed miRNA",style = "color: black; font-size: 16px" ),
                br(),
                #downloadButton("dlDTMI", "Download table"),
                br(),
                selectInput(inputId = "dtmi_uda" , label = "Choose up- or downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE),
                dataTableOutput("DTMI")%>% withSpinner(color= "#000000"),
                tags$head(tags$style("#DTMI table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
                
                tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
                
              ),
              tabPanel(
                id="t6","Volcano Plot", 
                
                plotOutput("VPMI")%>% withSpinner(color= "#000000"),
                #downloadButton("dlVPMI", "Download plot in PNG format")
                downloadButton("dl_VP_miRNA1", "Download plot in PNG format"),
                downloadButton("dl_VP_miRNA2", "Download plot in PDF format")
              ),
              tabPanel(
                id="t7","HeatMap", 
                
                plotOutput("HMMI")%>% withSpinner(color= "#000000"),
                #downloadButton("dlHMMI", "Download plot in PNG format")
                downloadButton("dl_HM_miRNA1", "Download plot in PNG format"),
                downloadButton("dl_HM_miRNA2", "Download plot in PDF format")
              ),
              tabPanel(
                id="t8","MA plot", 
                
                plotOutput("MAMI")%>% withSpinner(color= "#000000"),
                #downloadButton("dlHMMI", "Download plot in PNG format")
                downloadButton("dl_MA_miRNA1", "Download plot in PNG format"),
                downloadButton("dl_MA_miRNA2", "Download plot in PDF format")
              ),
              tabPanel(
                id="t9","PCA plot", 
                
                plotOutput("PCAMI")%>% withSpinner(color= "#000000"),
                #downloadButton("dlHMMI", "Download plot in PNG format")
                downloadButton("dl_PCA_miRNA1", "Download plot in PNG format"),
                downloadButton("dl_PCA_miRNA2", "Download plot in PDF format")
              ),
              tabPanel(
                id="t10","Expression plot", 
                
                plotOutput("ExpMI")%>% withSpinner(color= "#000000"),
                #downloadButton("dlHMMI", "Download plot in PNG format")
                downloadButton("dl_EXP_miRNA1", "Download plot in PNG format"),
                downloadButton("dl_EXP_miRNA2", "Download plot in PDF format")
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
              selectInput(inputId = "UDMI" , label = "Choose upregulated / downregulated" , choices = c("upregulated", "downregulated"), multiple = FALSE),
              conditionalPanel("input.miOnt == 'Biological Process'",
                               prettySwitch(
                                 "gos_of_interest_MI",
                                 label = "Plot only the GO Terms of interest (containing the words platelet, hemostasis and coagulation):",
                                 value = FALSE,
                                 status = "default",
                                 slim = FALSE,
                                 fill = TRUE,
                                 bigger = FALSE,
                                 inline = FALSE,
                                 width = NULL
                               )            
              ),
              tabPanel(
                id= "tE","miRNA Targets",
                #downloadButton("dlMIDT", "Download table"),
                br(),br(),
                dataTableOutput("MIDT")%>% withSpinner(color= "#000000"),
                tags$head(tags$style("#DTMIT table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
                
                tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
              ),
              tabPanel(
                id="t10","miRNA Target Enrichment", 
                #selectInput(inputId = "MSIG" , label = "Choose significant or all mRNA" , choices = c("significant","all"), multiple = FALSE),
                #actionButton(inputId = "clickSIG",label = "Get mRNA ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
                span("The  GO terms" ,style = "color: black; font-size: 18px"),br(),
                #downloadButton("dlDTMIT", "Download table"),
                br(),br(),
                dataTableOutput("DTMIT")%>% withSpinner(color= "#000000"),
                tags$head(tags$style("#DTMIT table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
                
                tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                               border-top: 1px solid #000000;}",media="screen", type="text/css"))
                
              ),
              tabPanel(
                id="t11","miRNA Target Network Graph", 
                conditionalPanel(
                condition = "input.miEXP == true || input.miSign == true",
                span("the results show miRNAs (darkblue nodes) and their mapping genes, while green nodes are upregulated genes and red nodes are downregulated genes.",
                     style = "color: black; font-size: 18px"),
                br(),
                forceNetworkOutput("MInet", height = "700px")%>% withSpinner(color= "#000000")
                
                ),
                conditionalPanel(
                  condition = "input.miEXP == false && input.miSign == false",
                  span("Unfortunately, by not choosing any of the checkboxes, the resulting network size is too big & therefore will not be displayed",
                       style = "color: black; font-size: 18px"),
                  br(),
                  
                )
              ),
              
              tabPanel(
                id="t12","miRNA Target Enrichment Graph", 
                span("miRNA Target Enrichment Graph",style = "color: black; font-size: 24px"),br(),br(),
                plotOutput("MIGOnet", height = "850px")%>% withSpinner(color= "#000000"),
                downloadButton("dl_graph_miRNA1", "Download plot in PNG format"),
                downloadButton("dl_graph_miRNA2", "Download plot in PDF format")
                
              ),
              tabPanel(
                title = "percentage of target genes in GO categories", solidHeader = TRUE,# width = 12,
                plotlyOutput("MIpG",height = "800px")%>% withSpinner(color= "#000000"),
                #downloadButton("dlPNGGO", "Download plot in PNG format")
              ),
              tabPanel(
                title = "enrichment FDR per GO categories of categories", solidHeader = TRUE,# width = 12,
                #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
                #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
                #br(),
                plotlyOutput("MIpFE",height = "800px")%>% withSpinner(color= "#000000"),
                #downloadButton("dlPNGGO", "Download plot in PNG format")
              ),#MIGOp
              tabPanel(
                title = "Top 10 most significant GO categories", solidHeader = TRUE,# width = 12,
                #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
                #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
                #br(),
                plotlyOutput("MIGOp",height = "800px")%>% withSpinner(color= "#000000"),
                #downloadButton("dlPNGGO", "Download plot in PNG format")
              ),
              tabPanel(
                title = "GO categories with the 10 highest DE genes", solidHeader = TRUE,# width = 12,
                #numericInput(inputId = "GOFE",value = 1.0 ,label = "Choose fold enrichment threshhold", min = 0 , max = 10, step = 0.5),
                #actionButton(inputId = "clickGOFE",label = "get Plot ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000"),
                #br(),
                plotlyOutput("MIGOC",height = "800px")%>% withSpinner(color= "#000000"),
                #downloadButton("dlPNGGO", "Download plot in PNG format")
              ),
              tabPanel(
                title = "Scatterplot- p.adjusted vs. identity percentage", solidHeader = TRUE,# width = 12,
                #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
                #br(),
                #withMathJax(),
                #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
                plotlyOutput("MIGO_perc_term",height = "800px")%>% withSpinner(color= "#000000"),
                #downloadButton("dlPNGGO", "Download plot in PNG format")
                #downloadButton("dlPNGGO", "Download plot in PNG format")
              ),
              tabPanel(
                title = "Scatterplot- p.adjusted vs. gene logFC", solidHeader = TRUE,# width = 12,
                #span("The z-score, computed by the circ_dat() Method from the GOplot package, is a score meant to give an orientation to understand if a category is more likely to be increased (z-score = positive) or decreased (z-score = negative)",style = "color: black; font-size: 15px"),
                #br(),
                #withMathJax(),
                #helpText("$$ zscore= \\frac{(up-down)}{\\sqrt{count}}$$"),
                plotlyOutput("MIGO_padj_logFC",height = "800px")%>% withSpinner(color= "#000000"),
                #downloadButton("dlPNGGO", "Download plot in PNG format")
              )
            )
          )
        )
        
        
      ),

# ciRNA Analysis ----------------------------------------------------------

      tabItem(
        tabName = "cA",
        fluidRow(
          tabBox(
            title = "Filter the results",
            # The id lets us use input$tabset1 on the server to find the current tab
            
            id = "VPfilter3", #height = "200px", 
            # tabPanel(id = "t1", "circRNA Analysis - Guide",
            #          span("circRNA analysis - Guide", style = "color: black; font-size: 18px; font-weight: bold"),br(),br(),
            #          "This page shows the differential expression results of circular RNA. The results are based on transcriptome analysis of reticulated and mature platelets from 19 CAD patients.",br(),br(),
            #          "The results of the circRNA differential expression analysis can be browsed by simply selecting the adjusted p-value threshold and the Log2 fold-change threshold. This analysis is represented by a volcano plot and the significantly differentially expressed circRNAs are listed in a data table. The results are presented in a downloadable data table and plots.",br(),br(),
            #          "The circRNA interaction analysis can also be browsed by selecting the adjusted p-value threshold and Log2 fold-change threshold of the interacting gene expression. A network of interactions between significantly upregulated genes (green), significantly downregulated genes (red) and circRNAs (dark blue) can be displayed. Other analyses include a correlation analysis of genes and circRNAs, biotypes of genes whose expression is significantly correlated with circRNA expression, and a bar graph showing the number of circRNAs per host gene."
            #          ),
            tabPanel(id="t2","Get differentially expressed circRNA", 
                     #selectInput(inputId = "ciDS" , label = "Choose from patient datasets" , choices = c("MP vs. RP","MP vs. RP in CAD patients", "MP vs. RP in MI patients" , "CAD vs MI"), multiple = FALSE),
                     #pickerInput(inputId = "ciDS" , label = "Choose from patient datasets" , choices = c("MP vs. RP in CAD patients","MP vs. RP","MP vs. RP in MI patients","CAD vs MI"), multiple = FALSE,choicesOpt = list(
                     #  disabled = c("MP vs. RP in CAD patients","MP vs. RP","MP vs. RP in MI patients","CAD vs MI") %in% c("MP vs. RP","MP vs. RP in MI patients","CAD vs MI")
                     #)),
                     span("Filter differentially expressed circRNA", style = "color: black; font-size: 15px; font-weight: bold"),br(),
                     numericInput(inputId = "pCOcirc",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "fCOcirc",value =  0.5 , label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
                     
                     span("CircRNA interaction analysis: Filter differentially expressed genes", style = "color: black; font-size: 15px; font-weight: bold"), br(),
                     numericInput(inputId = "pCOcirc2",value = 0.05 , label = "p-value cutoff (gene expression)", min = 0.01 , max = 0.05, step = 0.005),
                     numericInput(inputId = "fCOcirc2",value =  0.5 , label = "log2(Foldchange) cutoff (gene expression)", min = 0 , max = 4.0, step = 0.5),
                     selectInput(inputId = "ciAnnot" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
                     br(),
                     actionButton(inputId = "clickCirc",label = "get Table ",icon = icon('arrow'), style = "color: #FFF; background-color:#000000; border-color: #000000")#,
                     
            #)#,
            #tabPanel(id="t2","ciRNA Interaction Analysis", 
                     #pickerInput(inputId = "ciDS2" , label = "Choose from patient datasets" , choices = c("MP vs. RP in CAD patients","MP vs. RP","MP vs. RP in MI patients","CAD vs MI"), multiple = FALSE,choicesOpt = list(
                     #  disabled = c("MP vs. RP in CAD patients","MP vs. RP","MP vs. RP in MI patients","CAD vs MI") %in% c("MP vs. RP","MP vs. RP in MI patients","CAD vs MI")
                     #)),                     
                     
                     ##actionButton(inputId = "clickCirc2",label = "get Table ",icon = icon('arrow'), style = "color: #FFF; background-color:#000000; border-color: #000000")#,
                     
            ),
            
          ),
          #fluidRow(
          box(
            span("circRNA analysis - Guide", style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,br(),br(),
            "This page shows the differential expression results of circular RNA. The results are based on transcriptome analysis of reticulated and mature platelets.",br(),br(),
            "The results of the circRNA differential expression analysis can be browsed by simply selecting the adjusted p-value threshold and the Log2 fold-change threshold. This analysis is represented by a volcano plot and the significantly differentially expressed circRNAs are listed in a data table. The results are presented in a downloadable data table and plots.",br(),br(),
            "The circRNA interaction analysis can also be browsed by selecting the adjusted p-value threshold and Log2 fold-change threshold of the interacting gene expression. A network of interactions between significantly upregulated genes (green), significantly downregulated genes (red) and circRNAs (dark blue) can be displayed. Other analyses include a correlation analysis of genes and circRNAs, biotypes of genes whose expression is significantly correlated with circRNA expression, and a bar graph showing the number of circRNAs per host gene."
            
            #,
            # tabPanel(
            #   id="t8","MA plot", 
            #   "under construction"
            #   #plotOutput("CIMA")%>% withSpinner(color= "#000000")#,
            #   #downloadButton("dlHMMI", "Download plot in PNG format")
            # )
            )
        ),
        fluidRow(
          #uiOutput("Crna")
          tagList(
            tabBox(
              width = 12, 
              height = "1000px",
            tabPanel(
              id="t6","Volcano Plot", 
              
              plotOutput("CIVP",height = "800px")%>% withSpinner(color= "#000000"),
              #downloadButton("dlVPMI", "Download plot in PNG format")
              downloadButton("dl_VP_ciRNA1", "Download plot in PNG format"),
              downloadButton("dl_VP_ciRNA2", "Download plot in PDF format")
              
            ),
            tabPanel(
              id="t5","List of differentially expressed circRNA", "",
              selectInput(inputId = "circUD" , label = "Choose up or downregulated ciRNAs" , choices = c("upregulated","downregulated"), multiple = FALSE),
              br(),
              actionButton(inputId = "clickCUD",label = "Get ciRNAs ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
              br(),
              #downloadButton("dlCTO", "Download Table"),
              br(),
              dataTableOutput("CTO"),
              tags$head(tags$style("#CTO table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
              
              tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                                    border-top: 1px solid #000000;}",media="screen", type="text/css"))
              
            ),
              #tabPanel(id = "CIVP1", "Network", "",
              #         br(),
                       #downloadButton("dlCIVP", "Download plot in PNG format"),
              #         br(),
              #         visNetworkOutput("ciNetwork", width = "100%", height = "800px") %>% withSpinner(color= "#000000")
                       
              #),
              
              tabPanel(id = "CR1", "List of differentially expressed genes & circRNAs", "",
                       #selectInput(inputId = "circUD" , label = "Choose up or downregulated ciRNAs" , choices = c("upregulated","downregulated"), multiple = FALSE),
                       #actionButton(inputId = "clickCUD",label = "Get ciRNAs ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
                       #br(),
                       #downloadButton("dlCTO", "Download Table"),
                       br(),
                       dataTableOutput("CircCor")%>% withSpinner(color= "#000000"),
                       tags$head(tags$style("#CircCor table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),

                       tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
                                    border-top: 1px solid #000000;}",media="screen", type="text/css"))

              ),
               tabPanel(id = "CR2", "Biotypes of correlating genes", "",
                        plotlyOutput("CircBP",height = "800px")%>% withSpinner(color= "#000000")
               ),
              tabPanel(id = "CR3", "circRNA-gene Network", #"under construction",
                       forceNetworkOutput("CInet", height = "700px")%>% withSpinner(color= "#000000")
              ),
              # tabPanel(id = "CR4", "circRNA correlation", #"under construction",
              #          radioGroupButtons(
              #            inputId = "circ_cor",
              #            label = "choose correlation statistic",
              #            choices = c("mscore", 
              #                        "partial correlation", "correlation"),
              #            checkIcon = list(
              #              yes = icon("ok",
              #                         lib = "glyphicon"))
              #          ),
              #          numericInput(inputId = "circ_cutoff",value =  0.0 , label = "statistic cutoff ((the absolute) displayed values are ≥ than the cutoff)", min = -1.0 , max = 1.0, step = 0.1),
              #     
              #          plotlyOutput("Circ_mscore",height = "800px")%>% withSpinner(color= "#000000")
              #          
              # ),
              tabPanel(id = "CR5", "number of circRNA per genes", "",
                       numericInput(inputId = "circ_cutoff_gene",value =  0.0 , label = "minimal amount of circRNA", min = 0 , max = 200, step = 5),
                       plotlyOutput("Circ_gene",height = "800px")%>% withSpinner(color= "#000000")
              )#,
              # tabPanel(id = "CR6", "GO enrichment", "",
              #          #selectInput(inputId = "circUD" , label = "Choose up or downregulated ciRNAs" , choices = c("upregulated","downregulated"), multiple = FALSE),
              #          #actionButton(inputId = "clickCUD",label = "Get ciRNAs ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
              #          #br(),
              #          #downloadButton("dlCTO", "Download Table"),
              #          br(),
              #          selectInput(inputId = "circGO_ont" , label = "Choose ontology" , choices = c("Biological Process", "Cellular Component", "Molecular Function"), multiple = FALSE),
              #          selectInput(inputId = "circGO_UDA" , label = "Choose Up/Downregulated" , choices = c("upregulated","downregulated"), multiple = FALSE, selected = "upregulated"),
              #          
              #          dataTableOutput("CircGO")%>% withSpinner(color= "#000000"),
              #          tags$head(tags$style("#CircCor table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
              #          
              #          tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
              #                       border-top: 1px solid #000000;}",media="screen", type="text/css"))
              #          
              # ),
              # tabPanel(id = "CR7", "GO enrichment plots", #"under construction",
              #          selectInput(inputId = "circGO_ont2" , label = "Choose ontology" , choices =c("Biological Process", "Cellular Component", "Molecular Function"), multiple = FALSE),
              #          selectInput(inputId = "circGO_UDA2" , label = "Choose Up/Downregulated" , choices = c("upregulated","downregulated"), multiple = FALSE, selected = "upregulated"),
              #          
              #          radioGroupButtons(
              #            inputId = "circ_GO_pl",
              #            label = "Choose GO plot",
              #            choices = c("percentage in GO category", 
              #                        "enrichment FDR", "Top 10 most significant GO categories", "GO categories with 10 highest DE genes"),
              #            checkIcon = list(
              #              yes = icon("ok",
              #                         lib = "glyphicon"))
              #          ),
              #          #numericInput(inputId = "circ_cutoff",value =  0.0 , label = "statistic cutoff ((the absolute) displayed values are ≥ than the cutoff)", min = -1.0 , max = 1.0, step = 0.1),
              #          
              #          plotlyOutput("CircGO_pl",height = "800px")%>% withSpinner(color= "#000000")
              #          
              # ),
            )
          )
          
        )
        #))
      ),
      # tabItem(
      #   tabName = "KPM",
      #   #width =12,
      #   fluidRow(
      #     tabBox(
      #       # width =12, 
      #       tabPanel(id="t2","Mature vs. Reticulated", 
      #                selectInput(inputId = "KTPM1" , label = "Choose TPM cutoff" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
      #                selectInput(inputId = "KD1" , label = "Choose on which on dataset of patients the KPM-Algorithm will be performed on their significantly differentially expressed datasets" , choices = c("healthy","diseased", "diseased CCS", "diseased MI"), multiple = FALSE),
      #                selectInput(inputId = "KUD1" , label = "Choose the upregulated/downregulated GOs" , choices = c("upregulated", "downregulated"), multiple = FALSE),
      #                selectInput(inputId = "KO1" , label = "Choose Ontology" , choices = c("biological process", "cellular component", "molecular function"), multiple = FALSE),
      #                selectInput(inputId = "KES1" , label = "Choose KPM network gene annotation" , choices = c("ENSEMBL", "gene symbol"), multiple = FALSE),
      #                
      #                actionButton(inputId = "clickK1",label = "get Results ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")#,
      #       ),     
      #       tabPanel(id="t3","CCS vs ACS", 
      #                selectInput(inputId = "KTPM2" , label = "Choose TPM cutoff" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
      #                selectInput(inputId = "KD2" , label = "Choose on which on dataset of patients the KPM-Algorithm will be performed on their significantly differentially expressed datasets" , choices = c("diseased", "diseased - mature platelets", "diseased - reticulated platelets"), multiple = FALSE),
      #                selectInput(inputId = "KUD2" , label = "Choose the upregulated/downregulated GOs" , choices = c("upregulated", "downregulated"), multiple = FALSE),
      #                selectInput(inputId = "KO2" , label = "Choose Ontology" , choices = c("biological process", "cellular component", "molecular function"), multiple = FALSE),
      #                selectInput(inputId = "KES2" , label = "Choose KPM network gene annotation" , choices = c("ENSEMBL", "gene symbol"), multiple = FALSE),
      #                
      #                actionButton(inputId = "clickK2",label = "get Results ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")#,
      #       )
      #     ),
      #     box(
      #       title = span("De novo functional enrichment analysis using KeyPathwayMiner - Guide",style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
      #       "The de novo functional enrichment analysis is performed here by processing the dataset of each patient group's differentially expressed genes with the KeyPathwayMiner algorithm in order to identify new networks between genes that are then functionally categorized in Gene Ontologies by ShinyGO.",br(),br(),
      #       'The results are first categorized in two filters: The first being the differentiation of mature platelets versus reticulated platelets and the second being the comparative analysis between two diseased groups of patients - patients with Chronic coronary syndrome and patients with Acute coronary syndrome. The results are based on the transcriptome analysis of reticulated and mature platelets from 4 healthy patients and 10 diseased patients, of which 6 were diagnosed with Chronic coronary syndrome (CCS) and 4 were diagnosed with Acute coronary syndrome (ACS).',br(),br(),
      #       'The browsability of these results is ensured by being able to filter by the TPM (transcript per million) read count normalization threshold of the differentially expressed genes, the desired patient dataset and the up- or downregulated genes of this dataset, that correspond to the different sets in the before mentioned categories, the desired Ontology and the desired gene ID format.', br(), br(),
      #       'The representation of these results is then shown in the results by data tables that show the mapping of the genes to the functional categories, different network plots that also visualize those mappings and the relationships between the functional categories, and a Gene Ontology hierarchy plot to visualize the most important functional categories and the relationship between each other.' 
      #     )
      #   ),
      #   fluidRow(
      #     uiOutput("KPMO")
      #   )
      # ),
      
      
      tabItem(
        tabName = "DAS",
        #width =12,
        fluidRow(
          #tabBox(
          # width =12, 
          box(id="t2",span("Configure your search",style = "color: black; font-size: 16px; font-weight: bold"), 
              #selectInput(inputId = "DAS1" , label = "Choose differential Alternative Splicing Dataset" , choices = c("Reticulated Platelets vs Mature  Platelets (healthy)","Mature Platelets vs Reticulated Platelets (stable CCS)", "healthy vs stable CCS (Mature Platelets)", "healthy vs stable CCS (Reticulated Platelets)"), multiple = FALSE),
              #pickerInput(inputId = "DAS1" , label = "Choose differential Alternative Splicing Dataset" , choices = c("Mature Platelets vs Reticulated Platelets (stable CCS)", "Reticulated Platelets vs Mature  Platelets (healthy)","healthy vs stable CCS (Mature Platelets)", "healthy vs stable CCS (Reticulated Platelets)"), multiple = FALSE,choicesOpt = list(
              #  disabled = c("Mature Platelets vs Reticulated Platelets (stable CCS)", "Reticulated Platelets vs Mature  Platelets (healthy)","healthy vs stable CCS (Mature Platelets)", "healthy vs stable CCS (Reticulated Platelets)") %in% c("Reticulated Platelets vs Mature  Platelets (healthy)","healthy vs stable CCS (Mature Platelets)", "healthy vs stable CCS (Reticulated Platelets)")
              #)),
              numericInput(inputId = "DASe",value = 0.2 ,label = paste0('Choose DELTAPSI = E(','\U0394','PSI) per LSV junction cutoff'), min = 0.2 , max = 1, step = 0.05),
              numericInput(inputId = "DASp",value = 0.95 ,label = paste0('Choose confidence cutoff'), min = 0.95 , max = 1 , step = 0.05),
              selectInput(inputId = "DASga" , label = "Choose gene annotation" , choices = c("ENSEMBL ID", "HGNC symbol"), multiple = FALSE),
              
              
              actionButton(inputId = "clickDAS1",label = "get Results ",icon = icon('arrow'), style = "color: #fff; background-color: #000000; border-color: #000000")#,
          ),     
          #),
          box(
            title = span("Differential Alternative Splicing - Guide",style = "color: white; font-size: 18px; font-weight: bold"), background = "black" , solidHeader = TRUE,
            "The differential alternative splicing analysis results presented on this page have been obtained by using MAJIQ and Voila on the different datasets, the most important of which can be browsed here. The dataset used throughout Platlas is based on transcriptome analysis of reticulated and mature platelets" ,br(),br(),
            "This analysis shows which local splicing variations per gene are differentially spliced between mature and reticulated platelets in CAD patients. Delta E(PSI) is used to filter the degree of differentiation, with a default value of 0.2. This means that only local splicing variations per gene that are 20% more (less) spliced in the corresponding data sets are selected. The filter P(|deltaPSI|>= 0.20) acts as a p-value threshold, where the default value is 0.95, but the higher the number, the lower the actual p-value." ,br(),br(),
            "The results are presented here in different plots showing the differences in E(PSI) between RPs and MPs for each splicing event, a list of significantly alternatively spliced genes and bar plots showing the number of significant local splicing variations per gene. "
          )
        ),
        fluidRow(
          #uiOutput("DASO")
          tabBox(
            width =12,
            height = "1500px",
            title = "Get the results",
            tabPanel(
              id= "tT", "DELTAPSI = E(PSI) per junction",
              # uiOutput("DTO")
              #selectInput(inputId = "UpDown" , label = "Choose up or downregulated genes" , choices = c("upregulated","downregulated"), multiple = FALSE),
              #actionButton(inputId = "click4",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
              br(),
              br(),
              plotlyOutput("DASBX", height = "800px"),
              #downloadButton("downloadDASBX", "Download table"),
            ),
            # tabPanel(
            #   id= "tH", 
            #   #uiOutput("HMO")
            #   
            #   title = "Differential Alternative Splicing Activity", background = "red" , solidHeader = TRUE,
            #   br(),
            #   downloadButton("dlDASVP", "Download plot in PNG format"),
            #   br(),
            #   plotOutput("DASVP", height = "800px"),
            #   #forceNetworkOutput("KPMnet", height = "850px")
            #   
            #   
            #   
            # ),
            tabPanel(
              id= "tH", 
              #uiOutput("HMO")
              
              title = "Significantly alternatively spliced Genes", background = "red" , solidHeader = TRUE,
              #textOutput("DASueber"),
              #selectInput(inputId = "DASUD", label = "Choose significantly alternatively spliced LSV list ", choices = c("Condition 1", "Condition 2"), multiple = FALSE),
              #actionButton(inputId = "clickDASDT",label = "Get Genes ",icon = icon('arrow'), style = "color: #FFF; background-color: #000000; border-color: #000000"),
              br(),
              br(),
             # downloadButton("dlDASDT", "Download Table"),
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
            # tabPanel(
            #   id= "tN",
            #   #uiOutput("HMO")
            # 
            #   title = "LSV in differentially expressed genes", background = "red" , solidHeader = TRUE,
            #   span("Choose the filters, which define differentially expressed genes and get a list of the LSVs in those Genes, that also shows if they are upregulated/downregulated",style = "color: black; font-size: 18px"),
            #   #selectInput(inputId = "DASTPM" , label = "Read count normalization: Choose TPM (Transcripts per million) threshholds" , choices = c("TPM > 0.1","TPM > 0.3", "TPM > 1", "TPM > 2"), multiple = FALSE),
            #   numericInput(inputId = "DASpc",value = 0.05 , label = "p-value cutoff", min = 0.01 , max = 0.05, step = 0.005),
            #   numericInput(inputId = "DASfc",value =  1.0 , label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
            #   actionButton(inputId = "DD_click",label = "Start Analysis", style = "color: #fff; background-color: #000000; border-color: #000000"),
            # 
            #   br(),
            #   br(),
            #   # fluidRow({width =12
            #   #column(6,plotlyOutput("DASDIF",height = "400px")),
            #   #column(6,plotlyOutput("DASDIFTOE",height = "400px")),
            #   #}),
            # 
            # 
            #   #br(),
            #   #br(),
            #   #downloadButton("dlDASDIFtab", "Download Table"),
            #   br(),
            #   dataTableOutput(outputId = "DASDIFtab"),
            #   tags$head(tags$style("#DASDIFtab table {background-color: #DCDCDC; color : #000000;  }", media="screen", type="text/css")),
            # 
            #   tags$head(tags$style(".table>tbody>tr>td, .table>tbody>tr>th, .table>tfoot>tr>td, .table>tfoot>tr>th, .table>thead>tr>td, .table>thead>tr>th {
            #                  border-top: 1px solid #000000;}",media="screen", type="text/css"))
            # 
            #   #forceNetworkOutput("KPMGOnet", height = "850px")
            #   #  br(),
            #   # downloadButton("dlKPMN", "Download plot in PNG format")
            # 
            # ),
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
                     numericInput(inputId = "MGIfCO1",value = 0.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
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
                     numericInput(inputId = "DSIfCO1",value = 0.5 ,label = "log2(Foldchange) cutoff", min = 0 , max = 4.0, step = 0.5),
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
