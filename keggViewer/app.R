library(shiny)
library(mapview)
library(png)
library(raster)
library(leaflet)
library(biomaRt)
library(XML)
library(DT)

CONTRASTS <- c("ERCvsWT0", "WT24vsERC", "WT24vsWT0")
PATHWAYS <- list(
 "Carbon metabolism"="dme01200",
 "Glycolysis / Gluconeogenesis"="dme00010",
 "Citrate cycle (TCA cycle)"="dme00020",
 "Pentose phosphate pathway"="dme00030",
 "Pentose and glucuronate interconversions"="dme00040",
 "Fructose and mannose metabolism"="dme00051",
 "Galactose metabolism"="dme00052",
 "Pyruvate metabolism"="dme00620",
 "Oxidative phosphorylation"="dme00190",
 "Protein export"="dme03060",
 "Protein processing in endoplasmic reticulum"="dme04141",
 "SNARE interactions in vesicular transport"="dme04130",
 "ABC transporters"="dme02010",
 "MAPK signaling pathway - fly"="dme04013",
 "Wnt signaling pathway"="dme04310",
 "Notch signaling pathway"="dme04330",
 "Hedgehog signaling pathway - fly"="dme04341",
 "TGF-beta signaling pathway"="dme04350",
 "Hippo signaling pathway - fly"="dme04391",
 "Hippo signaling pathway -multiple species"="dme04392",
 "FoxO signaling pathway"="dme04068",
 "Phosphatidylinositol signaling system"="dme04070",
 "mTOR signaling pathway"="dme04150",
 "Endocytosis"="dme04144",
 "AGE-RAGE signaling pathway in diabetic complications"="dme04933"
)

ui <- fluidPage(
  fluidRow(
    column(4, selectInput("pathway", "Pathway:", choices=PATHWAYS, selected=PATHWAYS[2])),
    column(4, selectInput("contrast", "Contrast:", choices=CONTRASTS, selected=CONTRASTS[1])),
    column(4, sliderInput("plotSize", "Plot size:", min=400, max=1000, value=600, step=100))
  ),
  uiOutput("mapUI"),
  h3("Genes involved in this pathway:"),
  DT::dataTableOutput("DeTable")
)

server <- function(input, output, session) {
 
  output$mapUI <- renderUI({
    leafletOutput("map", height=input$plotSize)
  })
   
  ##
  ## read DE results
  ##
  withProgress({
    f <- "de.RData"
    if(file.exists(f)) {
      # do nothing, it will be loaded later
    } else {
      files <- c("/fsimb/groups/imb-bioinfocf/projects/jgu/jgu_berger_2016_01_kaiser_RNA-Seq/with_UMIs/results/DE_DESeq2/ERCvsWT0.csv",
                 "/fsimb/groups/imb-bioinfocf/projects/jgu/jgu_berger_2016_01_kaiser_RNA-Seq/with_UMIs/results/DE_DESeq2/WT24vsERC.csv",
                 "/fsimb/groups/imb-bioinfocf/projects/jgu/jgu_berger_2016_01_kaiser_RNA-Seq/with_UMIs/results/DE_DESeq2/WT24vsWT0.csv")
      
      de <- lapply(files, read.csv, check.names=FALSE)
      names(de) <- c("ERCvsWT0", "WT24vsERC", "WT24vsWT0")
    }
  }, value=0, message="Loading differential expression data")
  
  ##
  ## convert gene_id to CGID
  ##
  withProgress({
    f <- "de.RData"
    if(file.exists(f)) {
      load(f)
    } else {
      de <- lapply(de, function(x) {
          ann <- getBM(attributes=c("ensembl_gene_id", "entrezgene", "flybasecgid_gene"),
                       filters="ensembl_gene_id", values=x$gene_id,
                       mart=useMart("ensembl", dataset="dmelanogaster_gene_ensembl"))
          x$entrez <- ann$entrezgene[match(x$gene_id, ann$ensembl_gene_id)]
          x$cgid   <- ann$flybasecgid_gene[match(x$gene_id, ann$ensembl_gene_id)]
          x
      })
      save(de, file=f)
    }
  }, value=0, message="Annotating genes with Biomart (may take a while)")
 
  ##
  ## download and parse the kegg pathway id
  ##
  downloadKeggXml <- function(pathway) {
    
    # get from kegg's SOAP
    kegg <- xmlToList(xmlParse(paste0("http://rest.kegg.jp/get/", pathway, "/kgml")))
    
    # parse
    kegg <- lapply(kegg, function(x) {
      if(class(x) == "list") {
        if(!is.null(x$.attrs) & !is.null(x$graphics)) {
          if(x$.attrs["type"] == "gene") {
             data.frame(id=x$.attrs["id"],
                        cgid=sub("^dme:Dmel_", "", unlist(strsplit(x$.attrs["name"], " "))),
                        x=as.integer(x$graphics["x"]),
                        y=as.integer(x$graphics["y"]),
                        w=as.integer(x$graphics["width"]),
                        h=as.integer(x$graphics["height"]),
                        row.names=NULL)
          }
        }
      }
    })
    
    do.call(rbind, c(kegg[!sapply(kegg, is.null)], make.row.names=FALSE))
  }
  
  ##
  ## download kegg image
  ##
  downloadKeggImage <- function(pathway){
    
    download.file(paste0("http://rest.kegg.jp/get/", pathway, "/image"), paste0("/tmp/", pathway, ".png"))
    png <- readPNG(paste0("/tmp/", pathway, ".png"))
      
    blue  <- raster(png[, , 1])
    green <- raster(png[, , 2])
    red   <- raster(png[, , 3])
    
    # create the leaflet map from the rastered image
    img <- brick(red, green, blue)
    m <- viewRGB(img, maxpixels=ncol(png) * nrow(png))   # Red-Green-Blue plot of a multi-layered Raster object
    list(m=m, png=png)
  }
  
  ## 
  ## make the kegg map
  ##
  doPlot <- function(ohs, img, de) {
    
    # add DE info and collapse rows (genes) pointing to the same kegg box
    ohs$gene <- de[match(ohs$cgid, de$cgid), "gene_name"]
    ohs$fc   <- round(de[match(ohs$cgid, de$cgid), grepl("log2 fold change", colnames(de))], digits=2)
    ohs$fdr  <- de[match(ohs$cgid, de$cgid), grepl("BH adjusted p-values", colnames(de))]
    ohs$fdr  <- ifelse(is.na(ohs$fdr), 1, ohs$fdr)
    ohs <- by(ohs, ohs$id, function(x){
      out <- x[1, , drop=F]
      x$fdr2 <- format(x$fdr, scientific=TRUE, digits=3)
      out$label <- paste(apply(x, 1, function(x) paste0(x["gene"], " (", x["fc"], ", ", x["fdr2"], ")")), collapse=", ")
      if(sum(x$fdr < .01) > 0) {  # put red/green labels if any gene of the cluster is significant
        i <- which.max(abs(x$fc[x$fdr < .01]))
        out$dir <- ifelse(x$fc[i] < -1, "red", ifelse(x$fc[i] > 1, "green", "grey"))
      } else {
        out$dir <- "grey"
      }
      out
    })
    ohs <- do.call(rbind, ohs)
    
    # calculate the gene coordinates in the new space
    ohs.scaled <- ohs
    ohs.scaled$x <- (ohs$x / ncol(img$png)) * img$m@object[[1]]@extent@xmax
    ohs.scaled$y <- ((nrow(img$png) - ohs$y) / nrow(img$png)) * img$m@object[[1]]@extent@ymax
    ohs.scaled$w <- (ohs$w / ncol(img$png)) * img$m@object[[1]]@extent@xmax
    ohs.scaled$h <- (ohs$h / nrow(img$png)) * img$m@object[[1]]@extent@ymax
    ohs.scaled$xini <- ohs.scaled$x - ohs.scaled$w / 2
    ohs.scaled$xend <- ohs.scaled$x + ohs.scaled$w / 2
    ohs.scaled$yini <- ohs.scaled$y - ohs.scaled$h / 2
    ohs.scaled$yend <- ohs.scaled$y + ohs.scaled$h / 2
    
    # draw the object
    img$m@map %>%
    addRectangles(lng1=ohs.scaled$xini, lng2=ohs.scaled$xend, lat1=ohs.scaled$yini, lat2=ohs.scaled$yend, color=ohs$dir, popup=ohs$label, opacity=1, stroke=TRUE, fillOpacity=0)
  }
  
  ##
  ## Download pathway XML description
  ##
  pathwayXML <- reactive({
    # load pathway description file from kegg
    f <- paste0(input$pathway, ".xml.RData")
    if(file.exists(f)) {
      load(f)
    } else {
      withProgress(ohs <- downloadKeggXml(input$pathway), value=0, message="Downloading pathway description from Kegg")
      save(ohs, file=f)
    }
    
    ohs
  })
  
  ##
  ## Download pathway IMG file
  ##
  pathwayIMG <- reactive({
    # download and raster pathway image from kegg
    f <- paste0(input$pathway, ".img.RData")
    if(file.exists(f)) {
      load(f)
    } else {
      withProgress(img <- downloadKeggImage(input$pathway), value=0, message="Downloading pathway picture from Kegg")
      save(img, file=f)
    }
    
    img
  })   
  
  ##
  ## select contrast
  ##
  DE <- reactive({
    de[[input$contrast]]
  })
  
  ##
  ## render output plot
  ##
  output$map <- renderLeaflet({
    input$plotSize
    doPlot(ohs=pathwayXML(), img=pathwayIMG(), de=DE())
  })
  
  ##
  ## render output table
  ##
  output$DeTable <- renderDataTable({
    
    # select genes in the contrast which are involved in this pathway
    de=DE()
    ohs=pathwayXML()
    DT::datatable(de[de$cgid %in% ohs$cgid, ])
  })
}

shinyApp(ui, server)