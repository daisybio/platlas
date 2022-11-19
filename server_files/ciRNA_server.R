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

output$Crna <- renderUI({
  CRUI()
})




# CircRNA -----------------------------------------------------------------

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

