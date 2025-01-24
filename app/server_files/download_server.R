#___________download_BUTTONS___________

# output$downloadAll <- downloadHandler(
#   filename = function() {
#     paste("Differentially_Expressed_Genes", ".csv", sep = "")
#   },
#   content = function(file) {
#     write.csv2(saved_ud_tab(), file, row.names = FALSE)
#   },
#   contentType = "text/csv"
# )
# output$downloadGO <- downloadHandler(
#   filename = function() {
#     paste("Functional_Enrichment_Results", ".csv", sep = "")
#   },
#   content = function(file) {
#     write.csv2(saved_GO(), file, row.names = FALSE)
#   },
#   contentType = "text/csv"
# )
# output$dlGOC <- downloadHandler(
#   filename = function() {
#     paste("circular_functional_enrichment_plot", ".png", sep = "")
#   },
#   content = function(file) {
#     png(file)
#     sav_cir()
#     dev.off()
#   },
#   contentType = 'image/png'
# )

# DIFEX DOWNLOAD ----------------------------------------------------------

output$dl_VP_DIFEX1 <- downloadHandler(
  filename = function() {
    paste("VolcanoPlot", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(savedVP_plot())
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_VP_DIFEX2 <- downloadHandler(
  filename = function() {
    paste("VolcanoPlot", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(savedVP_plot())
    dev.off()
  },
  contentType = 'application/pdf'
)

output$dl_heatmap_DIFEX1 <- downloadHandler(
  
  filename = function() {
    paste("Heatmap", ".png", sep = "")
  },
  content = function(file) {
    png(file)
   # pheatmap(saved_heatmap_plot())
    grid::grid.newpage()
    grid::grid.draw(saved_heatmap_plot())
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_heatmap_DIFEX2 <- downloadHandler(
  
  filename = function() {
    paste("Heatmap", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    grid::grid.newpage()
    grid::grid.draw(saved_heatmap_plot())
    dev.off()
  },
  contentType = 'application/pdf'
)
output$dl_MA_DIFEX1 <- downloadHandler(
  filename = function() {
    paste("MA_plot", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_MA_plot())
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_MA_DIFEX2 <- downloadHandler(
  filename = function() {
    paste("MA_plot", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_MA_plot())
    dev.off()
  },
  contentType = 'application/pdf'
)
output$dl_BT_DIFEX1 <- downloadHandler(
  filename = function() {
    paste("Biotypes", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_BT_plot())
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_BT_DIFEX2 <- downloadHandler(
  filename = function() {
    paste("Biotypes", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_BT_plot())
    dev.off()
  },
  contentType = 'application/pdf'
)
output$dl_EXP_DIFEX1 <- downloadHandler(
  filename = function() {
    paste("Expression_analysis", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_EXP_plot())
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_EXP_DIFEX2 <- downloadHandler(
  filename = function() {
    paste("Expression_analysis", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_EXP_plot())
    dev.off()
  },
  contentType = 'application/pdf'
)
output$dl_PCA_DIFEX1 <- downloadHandler(
  filename = function() {
    paste("PCA", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_pcabi())
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_PCA_DIFEX2 <- downloadHandler(
  filename = function() {
    paste("PCA", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_pcabi())
    dev.off()
  },
  contentType = 'application/pdf'
)


# FE DOWNLOAD -------------------------------------------------------------

#png
#pdf


#dl_graph_FE1
output$dl_graph_FE1 <- downloadHandler(
  filename = function() {
    paste("GO_graph", ".png", sep = "")
  },
  content = function(file) {
   png(file)
   #grid::grid.newpage()
   grid::grid.draw(goplot(saved_GO_graph()))
   dev.off()
  },
  contentType = 'image/png'
)
output$dl_graph_FE2 <- downloadHandler(
  filename = function() {
    paste("GO_graph", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    #grid::grid.newpage()
    grid::grid.draw(goplot(saved_GO_graph()))
    dev.off()
  },
  contentType = 'application/pdf'
)

#dl_graph_GSEA
output$dl_graph_GSEA1 <- downloadHandler(
  filename = function() {
    paste("GSEA_graph", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    #grid::grid.newpage()
    grid::grid.draw(my_gseaplot2(saved_GSEA_plot(), geneSetID = saved_GSEA_plot_index() , title = saved_GSEA_plot()$Description[saved_GSEA_plot_index()], pvalue_table = TRUE))
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_graph_GSEA2 <- downloadHandler(
  filename = function() {
    paste("GSEA_graph", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    #grid::grid.newpage()
    grid::grid.draw(my_gseaplot2(saved_GSEA_plot(), geneSetID = saved_GSEA_plot_index() , title = saved_GSEA_plot()$Description[saved_GSEA_plot_index()], pvalue_table = TRUE))
    dev.off()
  },
  contentType = 'application/pdf'
)

output$dl_graph_circGOplot1 <- downloadHandler(
  filename = function() {
    paste("circular_graph", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    #grid::grid.newpage()
    grid::grid.draw(GOChord(saved_chord(), space = 0.02, limit = c(3, 0), gene.order = 'logFC', gene.space = 0.25, gene.size = 3, ribbon.col = brewer.pal(n = 5, name = "Set2"), border.size = 0.01))
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_graph_circGOplot2 <- downloadHandler(
  filename = function() {
    paste("circular_graph", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    #grid::grid.newpage()
    grid::grid.draw(GOChord(saved_chord(), space = 0.02, limit = c(3, 0), gene.order = 'logFC', gene.space = 0.25, gene.size = 3, ribbon.col = brewer.pal(n = 5, name = "Set2"), border.size = 0.01))
    dev.off()
  },
  contentType = 'application/pdf'
)

output$dl_graph_GSEA_hm_1 <- downloadHandler(
  
  filename = function() {
    paste("GSEAHeatmap", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    # pheatmap(saved_heatmap_plot())
    grid::grid.newpage()
    grid::grid.draw(saved_GSEA_heatmap())
    dev.off()
  },
  contentType = 'image/png'
)
output$dl_graph_GSEA_hm_2 <- downloadHandler(
  
  filename = function() {
    paste("GSEAHeatmap", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    grid::grid.newpage()
    grid::grid.draw(saved_GSEA_heatmap())
    dev.off()
  },
  contentType = 'application/pdf'
)



# MIRNA DOWNLOAD ----------------------------------------------------------

#png
output$dl_VP_miRNA1 <- downloadHandler(
  filename = function() {
    paste("miRNA_VolcanoPlot", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_VPMI())
    dev.off()
  },
  contentType = 'image/png'
)

output$dl_HM_miRNA1 <- downloadHandler(
  filename = function() {
    paste("miRNA_Heatmap", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    # pheatmap(saved_heatmap_plot())
    grid::grid.newpage()
    grid::grid.draw(saved_HMMI())
    dev.off()
  },
  contentType = 'image/png'
)

output$dl_MA_miRNA1 <- downloadHandler(
  filename = function() {
    paste("miRNA_MAPlot", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_MAMI())
    dev.off()
  },
  contentType = 'image/png'
)

output$dl_PCA_miRNA1 <- downloadHandler(
  filename = function() {
    paste("miRNA_PCAPlot", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_PCAMI())
    dev.off()
  },
  contentType = 'image/png'
)

output$dl_EXP_miRNA1 <- downloadHandler(
  filename = function() {
    paste("miRNA_EXPPlot", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_ExpMI())
    dev.off()
  },
  contentType = 'image/png'
)

output$dl_graph_miRNA1 <- downloadHandler(
  filename = function() {
    paste("miRNA_GO_graph", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    #grid::grid.newpage()
    grid::grid.draw(goplot(saved_MIGOnet()))
    dev.off()
  },
  contentType = 'image/png'
)

#pdf
output$dl_VP_miRNA2 <- downloadHandler(
  filename = function() {
    paste("miRNA_VolcanoPlot", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_VPMI())
    dev.off()
  },
  contentType = 'application/pdf'
)

output$dl_HM_miRNA2 <- downloadHandler(
  filename = function() {
    paste("miRNA_Heatmap", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    # pheatmap(saved_heatmap_plot())
    grid::grid.newpage()
    grid::grid.draw(saved_HMMI())
    dev.off()
  },
  contentType = 'application/pdf'
)

output$dl_MA_miRNA2 <- downloadHandler(
  filename = function() {
    paste("miRNA_MAPlot", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_MAMI())
    dev.off()
  },
  contentType = 'application/pdf'
)

output$dl_PCA_miRNA2 <- downloadHandler(
  filename = function() {
    paste("miRNA_PCAPlot", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_PCAMI())
    dev.off()
  },
  contentType = 'application/pdf'
)

output$dl_EXP_miRNA2 <- downloadHandler(
  filename = function() {
    paste("miRNA_EXPPlot", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_ExpMI())
    dev.off()
  },
  contentType = 'application/pdf'
)

output$dl_graph_miRNA2 <- downloadHandler(
  filename = function() {
    paste("miRNA_GO_graph", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    #grid::grid.newpage()
    grid::grid.draw(goplot(saved_MIGOnet()))
    dev.off()
  },
  contentType = 'application/pdf'
)




# CIRNA DOWNLOAD ----------------------------------------------------------

#png
output$dl_VP_ciRNA1 <- downloadHandler(
  filename = function() {
    paste("ciRNA_VolcanoPlot", ".png", sep = "")
  },
  content = function(file) {
    png(file)
    plot(saved_CIVP())
    dev.off()
  },
  contentType = 'image/png'
)

#pdf

output$dl_VP_ciRNA2 <- downloadHandler(
  filename = function() {
    paste("ciRNA_VolcanoPlot", ".pdf", sep = "")
  },
  content = function(file) {
    pdf(file)
    plot(saved_CIVP())
    dev.off()
  },
  contentType = 'application/pdf'
)






# output$dlVPMI <- downloadHandler(
#   filename = function() {
#     paste("miRNA_Volcanoplot", ".png", sep = "")
#   },
#   content = function(file) {
#     png(file)
#     saved_VPMI()
#     dev.off()
#   },
#   contentType = 'image/png'
# )
# output$dlHMMI <- downloadHandler(
#   filename = function() {
#     paste("miRNA_Heatmap", ".png", sep = "")
#   },
#   content = function(file) {
#     png(file)
#     saved_HMMI()
#     dev.off()
#   },
#   contentType = 'image/png'
# )
# output$dlMIDT <- downloadHandler(
#   filename = function() {
#     paste("miRNA_Targets", ".csv", sep = "")
#   },
#   content = function(file) {
#     write.csv2(getMIPaths(), file, row.names = FALSE)
#   },
#   contentType = "text/csv"
# )
# output$dlDTMIT <- downloadHandler(
#   filename = function() {
#     paste("miRNA_Target_Enrichment", ".csv", sep = "")
#   },
#   content = function(file) {
#     write.csv2(saved_DTMIT(), file, row.names = FALSE)
#   },
#   contentType = "text/csv"
# )
# output$dlCIVP <- downloadHandler(
#   filename = function() {
#     paste("VolcanoPlot_circRNA", ".png", sep = "")
#   },
#   content = function(file) {
#     png(file)
#     saved_CIVP()
#     dev.off()
#   },
#   contentType = 'image/png'
# )
# output$dlCTO <- downloadHandler(
#   filename = function() {
#     paste("circRNAs", ".csv", sep = "")
#   },
#   content = function(file) {
#     write.csv2(cRNA_ud_dt(), file, row.names = FALSE)
#   },
#   contentType = "text/csv"
# )
# output$downloadKPMD <- downloadHandler(
#   filename = function() {
#     paste("Functional_Enrichment_KPM", ".csv", sep = "")
#   },
#   content = function(file) {
#     write.csv2(Kpm(), file, row.names = FALSE)
#   },
#   contentType = "text/csv"
# )
# output$dlDASVP <- downloadHandler(
#   filename = function() {
#     paste("VolcanoPlot_DAS", ".png", sep = "")
#   },
#   content = function(file) {
#     png(file)
#     saved_DASVP()
#     dev.off()
#   },
#   contentType = 'image/png'
# )
# output$dlDASDT <- downloadHandler(
#   filename = function() {
#     paste("DiffAlternative_Splicing", ".csv", sep = "")
#   },
#   content = function(file) {
#     write.csv2(DASDataTable(), file, row.names = FALSE)
#   },
#   contentType = "text/csv"
# )
# output$dlDASDIFtab <- downloadHandler(
#   filename = function() {
#     paste("DAS_and_DE", ".csv", sep = "")
#   },
#   content = function(file) {
#     write.csv2(DD_tab(), file, row.names = FALSE)
#   },
#   contentType = "text/csv"
# )
# 
