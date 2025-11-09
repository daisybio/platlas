source("www/functions.R")
source("www/create_plots.R")
library(shiny)
server <- function(input,output,session){
  
  source(file.path("server_files", "difex_server.R"), local = TRUE)$value
  source(file.path("server_files", "fe_server.R"), local = TRUE)$value
  source(file.path("server_files", "miRNA_server.R"), local = TRUE)$value
  source(file.path("server_files", "ciRNA_server.R"), local = TRUE)$value
  source(file.path("server_files", "das_server.R"), local = TRUE)$value
  source(file.path("server_files", "download_server.R"), local = TRUE)$value
  
 # --------------------homepage
  study_doi <- "10.1093/eurheartj/ehaf694"                 # <-- set your DOI
  study_doi_url <- paste0("https://doi.org/", study_doi)
  study_pdf_rel <- "EHJ_paper.pdf"          # file at: www/papers/platlas_study.pdf
  
  #graphics.off()
  #dev.off()
  #output$slickr <- renderSlickR({
  #  imgs <- list.files("www/home/", pattern=".png", full.names = TRUE)
  #  slickR(imgs)
  #})
  output$start <- renderImage({
    list(src = "www/platlas-black.jpeg",
         contentType = 'image/jpeg',
         width = 650,
         height = 400,
         alt = "This is alternate text")
      }, deleteFile = FALSE)
  output$summary <- renderImage({
    list(src = "www/Platlas_graphical_summary2.jpeg",
         contentType = 'image/jpeg',
         width = 1200,
         height = 400,
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  observeEvent(input$home_preview_pdf, {
    pdf_full <- file.path("www", study_pdf_rel)
    if (file.exists(pdf_full)) {
      showModal(modalDialog(
        title = "Study preview (PDF)",
        size  = "l",
        easyClose = TRUE,
        footer = tagList(
          modalButton("Close"),
          tags$a(href = study_pdf_rel, class = "btn btn-primary",
                 download = "Platlas_study.pdf", "Download PDF")
        ),
        tags$iframe(src = study_pdf_rel,
                    style = "width:100%;height:70vh;border:0;")
      ))
    } else {
      showNotification("Local PDF not found in www/.", type = "error")
    }
  })
  # disable tabs Exposure, Covariate, and Construct on page load
  #shinyjs::disable(selector = '.navbar-nav a[data-value="DIFEX_t2"')
  #shinyjs::disable(selector = '.navbar-nav a[data-value="FEX_t2"')
  
}

#shinyApp(ui = ui, server = server)