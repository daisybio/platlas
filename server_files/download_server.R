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

