
#database <- read.csv2("C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/miRTarBase_MTI.CSV",header = TRUE, na.strings = "_")
#C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/
#miRNAfile <- unzip("www/miRTarBase_HUMAN.zip", files ="miRTarBase_HUMAN.CSV" )
#human <- read.csv2(miRNAfile,header = TRUE, na.strings = "_")
#strong_experimental_evidence <- human[!grepl("(Weak)",human[,8]),]
#weak_experimental_evidence <- human[grepl("(Weak)",human[,8]),]
#C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/
#genDF <- read.csv2("www/all_diseased_TPM-0-1_disease-phenotype-all.CSV",header = TRUE, na.strings = "_")
#C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/
library(data.table)
library(gridExtra)
library(grid)
library(readxl)
goterms <- read.table("www/goterms.txt", header = FALSE, sep = "\t", quote = "")
colnames(goterms) = c("id","term1")
path<-"www/all_"
names2links <- function(tab){
  hi <- data.frame()
  i <- 1
  for(g in tab[,1]){#
    link <- HTML(paste('<a href="https://www.ensembl.org/Homo_sapiens/Gene/Sequence?g=',g, '">',g, "</a>", sep = ""))
    #l <- paste("https://www.ensembl.org/Homo_sapiens/Gene/Sequence?g=",g,sep = "")
    #link <- a(href = l, g)
    tab[i,1] <- link
    i <- i + 1
  }
  # colnames(tab)[1] = "Gene_ID"
  # colnames(tab)[2] = "Log2(FC)"
  # colnames(tab)[3] = "adjusted p-value"
  # colnames(tab)[4] = "biotypes"
  return(tab)
}
fe_shiny_paths <- function(filsi,t,ds,o,uda){
  filt_ordner <- switch(filsi,"C" = "CCS/", "M" = "MPRP/")
  off_UDA <- switch(uda, "upregulated" = "up_","downregulated" = "down_","all" = "all_")
  norm_UDA <- switch(off_UDA, "up_" = "UP", "down_" = "DOWN", "all_" = "ALL")
  tpm <- switch(t, "TPM > 0.1" = "TPM-0-1","TPM > 0.3" = "TPM-0-3", "TPM > 1" = "TPM-1", "TPM > 2" = "TPM-2")
  off_zwischen <- switch(ds,"healthy" = "healthy","diseased"="diseased", "disease - stable CAD"="diseased" , "disease - MI"="diseased","disease - mature platelets"="diseased", "disease - reticulated platelets"="diseased")
  filt_zwischen <- switch(filsi, "M" = "phenotype-disease", "C" = "disease-phenotype")
  filt_zwischen_nom <- switch(filsi, "M" = "MP-RP", "C" = "CAD-MI")
  off_dataset <- switch(ds,"healthy" = "healthy","disease - stable CAD" = "stable-CAD" , "disease - MI" = "MI","diseased" = "all","disease - mature platelets" = "MP", "disease - reticulated platelets" = "RP")
  nom_dataset <- switch(off_dataset,"healthy" = "healthy", "stable-CAD"= "CAD", "all" = "diseased", "MI" = "MI", "MP" = "MP", "RP" = "RP" )
  nom_o <- switch(o, "biological process" = "BP","cellular component" = "CC", "molecular function" = "MF")
  ent <- ".txt"
  cir <- ".tsv"
  edg <- "_e.txt"
  im  <- "_p.png"
  entab_file = ""
  circ_file = ""
  edg_file = ""
  im_file = ""
  lis = list()
  if(filsi == "M" && off_zwischen == "healthy"){
    entab_file = paste(c("FE_results/",filt_ordner,"endtab/" ,off_UDA, off_zwischen, "_", tpm, ent), collapse = "")
    circ_file = paste(c("FE_results/",filt_ordner, "circ/","circ_",off_UDA, off_zwischen, "_", tpm, cir), collapse = "")
    edg_file = paste(c("FE_results/",filt_ordner,"edges/" ,nom_o, "_",filt_zwischen_nom, "_", tpm,"_", nom_dataset,"_" ,norm_UDA, edg), collapse = "")
    im_file = paste(c("FE_results/",filt_ordner,"plot/",nom_o, "_",filt_zwischen_nom, "_", tpm,"_", nom_dataset,"_" ,norm_UDA, im), collapse = "")
  }else{
    entab_file = paste(c("FE_results/",filt_ordner,"endtab/", off_UDA, off_zwischen, "_", tpm,"_" ,filt_zwischen,"-", off_dataset ,ent), collapse = "")
    circ_file = paste(c("FE_results/",filt_ordner,"circ/" ,"circ_",off_UDA, off_zwischen, "_", tpm,"_" ,filt_zwischen,"-", off_dataset, cir), collapse = "")
    edg_file = paste(c("FE_results/",filt_ordner,"edges/" ,nom_o, "_",filt_zwischen_nom, "_", tpm,"_", nom_dataset,"_" ,norm_UDA, edg), collapse = "")
    im_file = paste(c("FE_results/",filt_ordner,"plot/" ,nom_o, "_",filt_zwischen_nom, "_", tpm,"_", nom_dataset,"_" ,norm_UDA, im), collapse = "")
    
  }
  print(im_file)
  lis = list(entab_file,circ_file,edg_file,im_file)
  return(lis)
}

getGeneText <- function(t){
  et <- c()
  if(grepl(" ",t, fixed = TRUE)== TRUE){
    h <- strsplit(t, split = " ")#[[1]]
    et <- h
  }else if(grepl("\t",t, fixed = TRUE)== TRUE){
    h <- strsplit(t, split = "\t")#[[1]]
    et <- h
  }else if(grepl("\n",t, fixed = TRUE)== TRUE){
    h <- strsplit(t, split = "\n")#[[1]]
    et <- h
  }
  return(et)
}

getIDif <- function(tabin, p, f, ud,num){

  if (ud == "upregulated") {
    #print(f)
    endt <- tabin[(tabin[,num] >= as.numeric(f)) & (tabin[,(num + 4)] <= as.numeric(p)), ]

  }else{
    endt <- tabin[tabin[,num] <= as.numeric(-f) & tabin[,(num + 4)] <= as.numeric(p), ]
  }
  return(endt)
}

getMTAB <- function(pa){
  pas <- switch(pa,"Mature Platelets vs. Reticulated - all data" =  "rp_vs_mp_sign_diff_mir_counts.CSV", "in diseased patients" = "disease_diff_mir_counts.CSV", "only CAD patients"= "cad_rp_vs_mp_diff_mir_counts.CSV")
  enpa <- paste("www/", pas, sep = "")
  enpat <- read.csv2(enpa)
  return(enpat)
}


intsect <- function(eins, zwei){
  eins <- as.data.frame(eins)
  zwei <- as.vector(zwei)
  res <- eins[eins[,1] %in% zwei,]
  return(res)
}
names2GOlinks <- function(tab){
  hi <- data.frame()
  i <- 1
  for(g in tab[,1]){#
    link <- HTML(paste('<a href="http://amigo.geneontology.org/amigo/medial_search?q=',g, '">',g, "</a>", sep = ""))
    #l <- paste("https://www.ensembl.org/Homo_sapiens/Gene/Sequence?g=",g,sep = "")
    #link <- a(href = l, g)
    tab[i,1] <- link
    i <- i + 1
  }
  # colnames(tab)[1] = "GO-entry"
  #  colnames(tab)[6] = "categorization"
  
  
  return(tab)
}

# getGOs2tab <- function(tab){
#   tab$Functional_Category = tab$Functional_Category %>% make_clean_names()
#   goterms_clean <- goterms
#   goterms_clean$term1 <- goterms_clean$term1 %>% make_clean_names()
#   merged <- left_join(tab,goterms_clean, by = c("Functional.Category" = "term1"), keep = FALSE)
#   me_tab <- data.frame(merged$id, merged$Functional_Category, merged$Enrichment_FDR, merged$Genes_in_list, merged$Total_genes,merged$Genes,merged$V6,merged$V7)
#   colnames(me_tab)= c("id","Functional_Category","Enrichment_FDR", "Genes_in_list", "Total_genes","Genes" ,"V6", "V7")
#   return(me_tab)
# }

names2milinks <- function(tab){
  hi <- data.frame()
  i <- 1
  for(g in tab[,1]){#
    link <- HTML(paste('<a href="https://mirbase.org/results/?query=',g, '">',g, "</a>", sep = ""))
    tab[i,1] <- link
    i <- i + 1
  }
  colnames(tab)[1] = "miRNA"
  # colnames(tab)[6] = "categorization"
  return(tab)
}

miRNA_path <- function(ont, expe, netorenrich,ud,db,data,signif,datei){
  anfang <- "www/miRNA/"
  o <- switch(ont,"Biological Process"="BP_", "Cellular Component" = "CC_", "Molecular Function" = "MF_"  )
  e <- ""
  n <- ""
  s <- ""
  if(expe == TRUE){
    e <- "E_"
  }
  n <- switch(netorenrich, "N" = "N_", "G"= "G_")
  u <- switch(ud, "upregulated"="up", "downregulated"= "down")
  d <- switch(db,"miRTarBase"= "tar", "miRDB" = "rdb" )
  da <- switch(data, "Mature Platelets vs. Reticulated - all data" = "rp_mp","in diseased patients" = "dis","only CAD patients" = "cad"  )
  if(signif == TRUE){
    s <- "_sig"
  }
  letzt <- switch(datei, "1" ="_Liste_mit_Ontologies.txt", "2" = "_Nodes.txt", "3" = "_Edges.txt", "0" = ".tsv" )
  #c <- paste(b,letzt, sep = "")
  
  finalpath <- paste0(anfang,o,e,n,u,"_",d,"_",da,s,letzt)
  
}

circRNA_path <- function(inp){
  cp <- ""
  inp <- as.character(inp)
  cp <- switch(inp,"MP vs. RP"= "MP_RP.CSV","MP vs. RP in CAD patients"="CCS.csv","MP vs. RP in MI patients"="MPRP_MI.CSV", "CAD vs MI"= "CAD_MI.CSV"  )
  Wcp <- "www/circRNA/circRNA_"
  endCP <- paste(Wcp,cp,sep = "")
  print(Wcp)
  df <- read.csv(endCP)
  #df <- df[!is.na(df$padj),]
  return(df)
}
getImirnagene <- function(genlist, sigmirna, m2g){
  endt <- NULL
  endt <- m2g[(m2g[,1] %in% sigmirna) & (m2g[,3] %in% genlist),]
  colnames(endt) = c("miRNA ID", "hgnc symbol", "ensembl ID")
  return(endt)
}


paths2frames <- function(oneortwo,P,IT,IH,CT){
  end <- switch(CT,".CSV" = ".CSV",".tsv"= ".tsv" )
  tpm <- switch(as.character(IT),"TPM > 0.1"="TPM-0-1","TPM > 0.3"="TPM-0-3","TPM > 1"= "TPM-1", "TPM > 2"= "TPM-2", "TPM > 0.2"="TPM-0-2")
  hod <- ""
  finalpath <- ""
  if(oneortwo == "1"){
    if(as.character(IH)== "healthy"){
      pathWHod<- paste(P,as.character(IH), sep= "")
      newpath<- paste(pathWHod,tpm,sep = "_")
      finalpath <- paste(newpath,end,sep= "")
    }else{
      pathWHod<- paste(P,"diseased", sep= "")
      newpath<- paste(pathWHod,tpm,sep = "_")
      tail <- switch(as.character(IH),"diseased"="_phenotype-disease-all","disease - CAD"= "_phenotype-disease-stable-CAD","disease - ACS"= "_phenotype-disease-MI"  )
      tailplus <- paste(newpath, tail, sep = "")
      finalpath <- paste(tailplus,end ,sep= "")
    }
  }else{
    pathWHod<- paste(P,"diseased", sep= "")
    newpath<- paste(pathWHod,tpm,sep = "_")
    tail <- switch(as.character(IH),"diseased"="_disease-phenotype-all","disease - mature platelets"= "_disease-phenotype-MP", "disease - reticulated platelets" = "_disease-phenotype-RP" )
    tailplus <- paste(newpath, tail, sep = "")
    finalpath <- paste(tailplus,end,sep= "")
  }
  
  
  
  return(finalpath)
}

getData <- function(hum,sig,gDT, ud){
  Hum_sig<- hum[grep(paste(sig[,1],collapse="|"),hum[,2]),]
  colnames(sig)[1] = "miRNA"
  jointdataset <- merge(Hum_sig, sig, by = 'miRNA')
  entab <- jointdataset[ ,c(1,4,5,11)]
  neu_entab <- data.frame()
  if(ud == "u"){
    neu_entab <- subset(entab, entab[,4] > 0)
  }else{
    neu_entab <- subset(entab, entab[,4] < 0)
  }
  gen <- gDT[,11]
  gen =  gen[!(gen %in% "")]
  nentab <- neu_entab[neu_entab$Target.Gene %in% gen, ]
  geneList <- nentab[,4]
  names(geneList) <- as.character(nentab[,3])
  geneList <- sort(geneList, decreasing = TRUE)
  return(geneList)
}

KPMpaths2frames <- function(oneortwo,O,IT,IH,UD,CT,fil,ES){
  end <- switch(CT,".txt" = ".txt", ".png" = ".png", ".sif" = ".sif" )
  tpm <- switch(as.character(IT), "TPM > 0.1" = "TPM-0-1", "TPM > 0.3" = "TPM-0-3", "TPM > 1" = "TPM-1","TPM > 2" = "TPM-2"  )
  hod <- switch(as.character(IH),"healthy"="healthy","diseased"= "diseased", "diseased CCS"= "CCS","diseased MI" = "MI", "diseased - mature platelets" = "MP", "diseased - reticulated platelets" = "RP" )
  filter <- switch(oneortwo, "1" = "MP_vs_RP", "2" = "CCS_vs_MI" )
  ont <- switch(O,"biological process"="www/KPM/BP_KPM","cellular component"="www/KPM/CC_KPM", "molecular function"="www/KPM/MF_KPM","no"= "www/KPM/KPM" )
  ud <- switch(UD,"upregulated" = "UP", "downregulated"= "DOWN" )
  filnam <- switch(fil,"1"="Liste_mit_Ontologies", "2"= "Hierarchie_Baum", "3" ="Network", "4"="Nodes", "5"="Edges")
  if(fil == "6" && ES == "ENSEMBL"){
    filnam <- "Nodes_E"
  }else if(fil == "6" && ES == "gene symbol"){
    filnam <- "Nodes_S"
  }
  OUD <- paste(ont, ud, sep = "_")
  FHOD <- paste(filter, hod, sep = "_")
  TF <- paste(tpm,filnam, sep = "_")
  OF <- paste(OUD, FHOD, sep = "_")
  TFOF <- paste(OF, TF, sep = "_")
  finalpath <- paste(TFOF,end, sep = "")
  return(finalpath)
}

#mat2<- data.frame()
#C:/Users/Leonora/Desktop/Studium/6.Semester/Bachelor-Thema__implementierung_eines_Platelet_Atlas/R-programme-shiny/www/
path1 <-"www/heat_"
pathDAS <- "www/"
pickGeneAnnot <- function(geneID,tab){ #right gene id annotation
  endt <- NULL
  if(geneID == "e"){
    endt <- tab[,-c(2,3)]
  }else{
    endt <- tab[,-c(1,3)]
    
  }
  
  endt <- endt[!(startsWith(endt[,1]," ")),]
  endt<- endt[!(duplicated(endt)),]

  row.names(endt) = endt[,1]
  endt = endt[,-1]
  return(endt)
}
metapath <- "www/simulated_sample_names_all.tsv"
metatab <- as.data.frame(fread(metapath, header = TRUE, sep = "\t", na.strings = "NA", quote = "" ))
row.names(metatab) = metatab$sample_name
metatab = metatab[,-1]
#colnames(metatab) = c("medical_condition","platelet_condition")
metatab$CCS<- NULL
metatab$CCS[metatab$condition == "CCS"] <- 1
metatab$CCS[metatab$condition != "CCS"] <- 0
metatab$ACS<- NULL
metatab$ACS[metatab$condition == "ACS"] <- 1
metatab$ACS[metatab$condition != "ACS"] <- 0
metatab$MP<- NULL
metatab$MP[metatab$platelet == "MP"] <- 1
metatab$MP[metatab$platelet != "MP"] <- 0
metatab$RP<- NULL
metatab$RP[metatab$platelet == "RP"] <- 1
metatab$RP[metatab$platelet != "RP"] <- 0

metatab$healthy[metatab$condition == "healthy"] <- 1
metatab$healthy[metatab$condition != "healthy"] <- 0
#metatab$platelet_condition [metatab$medical_condition == "healthy"] <- 0

getDASdata <- function(ds){
  p <- switch(ds,"Reticulated Platelets vs Mature  Platelets (healthy)"= "healthy-RP_healthy-MP.deltapsi.tsv","Mature Platelets vs Reticulated Platelets (stable CCS)"="DAS_neu.xlsx" ,"healthy vs stable CCS (Mature Platelets)"= "healthy-MP_stable-MP.deltapsi.tsv","healthy vs stable CCS (Reticulated Platelets)" = "healthy-RP_stable-RP.deltapsi.tsv")
  left <- switch(ds,"Reticulated Platelets vs Mature  Platelets (healthy)"= "Reticulated Platelets in healthy patients","Mature Platelets vs Reticulated Platelets (stable CCS)"="Mature Platelets in patients with stable CCS" ,"healthy vs stable CCS (Mature Platelets)"= "healthy patients (Mature Platelets)","healthy vs stable CCS (Reticulated Platelets)" = "healthy patients (Reticulated Platelets)" )
  right <- switch(ds,"Reticulated Platelets vs Mature  Platelets (healthy)"= "Mature Platelets in healthy patients" ,"Mature Platelets vs Reticulated Platelets (stable CCS)"= "Reticulated Platelets in patients with stable CCS","healthy vs stable CCS (Mature Platelets)"="patients with stable CCS (Mature Platelets)" ,"healthy vs stable CCS (Reticulated Platelets)" = "patients with stable CCS (Reticulated Platelets)")
  endpat<- paste(pathDAS,p, sep = "")
  DF <- read_excel(endpat)#, sep = "\t" , na.strings = "NA", header = FALSE
  #View(DF)
  res <- list(DF,left,right) # zugriff res[[1]]
  return(res)
}

getValidEasy <- function(atab,e,p){ ##signifikant alternative splicing
  #
  out <- as.data.frame(atab)#data.frame(matrix(NA,ncol = 6,nrow = 0))
  #View(out)
  out <- subset(out,(abs(as.numeric(out[,6]))>=e & as.numeric(out[,3]) >=p))
  #View(valideasy)
  return(out)
}
miUD <- function(tab,f,p){
  U <- tab[tab[,3] >= f & tab[,6] >= p,]
  D <- tab[tab[,3] <= (-f) & tab[,6] >= p,]
  return(list(U,D))
}
getDIFmi <- function(name,tpm,logfc,pval){
  t <- ""
  if(tpm == "TPM > 0.1"){
    t <- "TPM-0-1"
  }else if(tpm == "TPM > 0.3"){
    t <- "TPM-0-3"
  }else if(tpm == "TPM > 1"){
    t <- "TPM-1"
  }else if(tpm == "TPM > 2"){
    t <- "TPM-2"
  }
  d1<- NULL
  if(name == "Mature Platelets vs. Reticulated - all data"){
    na1 <- c("www/all_diseased_",t,"_phenotype-disease-all.CSV")
    na2 <- c("www/all_healthy_",t,".CSV")
    na1 <- paste0(na1, collapse = "")
    na2 <- paste0(na2, collapse = "")
    print(na1)
    dx <- read.csv2(na1)
    dy <- read.csv2(na2)
    d1 <- rbind(dx,dy)
    
  }else if(name == "in diseased patients"){
    na <- c("www/all_diseased_",t,"_phenotype-disease-all.CSV")
    na <- paste0(na, collapse = "")
    d1 <- read.csv2(na)
    
  }else if(name == "only CAD patients"){
    na <- c("www/all_diseased_",t,"_phenotype-disease-stable-CAD.CSV")
    na <- paste0(na, collapse = "")
    d1 <- read.csv2(na)
  }
  UundD <-miUD(d1,logfc,pval)
  U <- UundD[[1]]
  D <- UundD[[2]]
  u <- U[,c(1,11)]
  d <- D[,c(1,11)]
  u <- u[!duplicated(u),]
  d <- d[!duplicated(d),]
  return(list(u,d))
}
umformen <- function(t){
  um <- t
  um <- um %>% mutate(V2 = as.numeric(V2)*10)
  um <- um %>% mutate(V3 = abs(as.numeric(V3) - 1))
  return(um)
}

f <- list(
  family = "Arial",
  size = 14,
  color = "#00000"
)

f_foldenrich <- function(a, num){
  b <- data.frame()
  b <- a[a[,3] > num, ]
  return(b)
}

getScatter <-  function(non,sig){
  zeilen <- NROW(non)
  endtab <- data.frame(matrix(data = NA,nrow = 0,ncol = 3))
  for(z in 1:zeilen){
    s <- non[z,1]
    en <- non[z,2]
    if(s %in% sig[,1] & en %in% sig[,2]){
      endtab[nrow(endtab) + 1,] = list(s,en,"significantly alternatively spliced junctions")
    }else{
      endtab[nrow(endtab) + 1,] = list(s,en,"all junctions")
    }
  }
  return(endtab)
}

getEorS <- function(eors,kptab){
  #endtab <- NULL
  colnames(kptab)[c(7,8)] = c("ensembl_gene_id","hgnc_symbol")
  if(eors == "e"){
    kptab <- kptab[,-c(6,8)]
  }else if(eors == "s"){
    kptab <- kptab[,-c(6,7)]
  }
  return(kptab)
}



#new gseaplot 
my_gseaplot2 <- function(x, geneSetID, title = "", color="green", base_size = 11,
                      rel_heights=c(1.5, .5, 1), subplots = 1:3,
                      pvalue_table = FALSE, ES_geom="line") {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  
  geneList <- position <- NULL ## to satisfy codetool
  
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))
  
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                          size=1)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                           size=1, data = subset(gsdata, position == 1))
  }
  
  p.res <- p + es_layer +
    theme(legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))
  
  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  if (length(geneSetID) == 1) {
    ## geneList <- gsdata$geneList
    ## j <- which.min(abs(geneList))
    ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
    ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]
    
    ## v <- sort(c(v1, v2))
    ## inv <- findInterval(geneList, v)
    
    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,
           xmax=~xmax,
           ymin=~ymin,
           ymax=~ymax,
           fill=~I(col)),
      data=d,
      alpha=.9,
      inherit.aes=FALSE)
  }
  
  ## p2 <- p2 +
  ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
  ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
  ## theme(legend.position="none") +
  ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                            color="grey")
  p.pos <- p.pos + ylab("Ranked List Metric") +
    xlab("Rank in Ordered Dataset") +
    theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("NES", "pvalue", "p.adjust")]
    rownames(pd) <- ""
    x <- pd
    # pd <- round(pd, 4)
    for (i in seq_len(ncol(pd))) {
      pd[, i] <- format(pd[, i], digits = 4)
    }
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }
  
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())
  
  if (length(subplots) == 1)
    return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                    l=.2, unit="cm")))
  
  
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  
  # aplot::plot_list(gglist = plotlist, ncol=1, heights=rel_heights)
 aplot::gglist(gglist = plotlist, ncol=1, heights=rel_heights)
}


gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

gseaScores <- getFromNamespace("gseaScores", "DOSE")


tableGrob2 <- function(d, p = NULL) {
  # has_package("gridExtra")
  d <- d[order(rownames(d)),]
  #check_installed('gridExtra', 'for `tableGrob2()`.')
  tp <- gridExtra::tableGrob(d)
  if (is.null(p)) {
    return(tp)
  }
  
  # Fix bug: The 'group' order of lines and dots/path is different
  p_data <- ggplot_build(p)$data[[1]]
  # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  p_data <- p_data[order(p_data[["group"]]), ]
  pcol <- unique(p_data[["colour"]])
  ## This is fine too
  ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
  j <- which(tp$layout$name == "rowhead-fg")
  
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
  }
  return(tp)
}

