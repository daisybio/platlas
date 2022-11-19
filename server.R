source("www/functions.R")

server <- function(input,output,session){
  
  source(file.path("server_files", "difex_server.R"), local = TRUE)$value
  source(file.path("server_files", "fe_server.R"), local = TRUE)$value
  source(file.path("server_files", "miRNA_server.R"), local = TRUE)$value
  source(file.path("server_files", "ciRNA_server.R"), local = TRUE)$value
  source(file.path("server_files", "kpm_server.R"), local = TRUE)$value
  source(file.path("server_files", "das_server.R"), local = TRUE)$value
  source(file.path("server_files", "int_server.R"), local = TRUE)$value
  source(file.path("server_files", "download_server.R"), local = TRUE)$value
  
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