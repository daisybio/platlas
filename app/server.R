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
  # disable tabs Exposure, Covariate, and Construct on page load
  #shinyjs::disable(selector = '.navbar-nav a[data-value="DIFEX_t2"')
  #shinyjs::disable(selector = '.navbar-nav a[data-value="FEX_t2"')
  
}

#shinyApp(ui = ui, server = server)