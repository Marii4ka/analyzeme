library(shiny)
library(jsonlite)
library(digest)
library(SNPRelate)
library(gwascat)
library(ggbio)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

options(shiny.maxRequestSize = -1)

shinyServer(function(input, output) {
  temp.file <- digest(date())
  
  load.file.fn <- reactive({
    data <- input$inFile
    path <- as.character(data$datapath)
    print(typeof(path))
    if(length(path) == 0) {
      stop()
    }
    print(paste("zcat < ", path, " > ./temp/", temp.file, sep=""))
    system(paste("unzip -p ", path, " > ./temp/", temp.file, sep=""))
    d <- read.table(paste("./temp/", temp.file, sep=""), sep="\t", header=FALSE, colClasses=c("character", "character", "numeric", "character"), col.names=c("rsid", "chrom", "position", "genotype"))
    
  })

  output$snpsByChr <- renderPlot({
    withProgress(message="Graphing data", value=0.1, {
        #     d <- tryCatch({
        #       load.file.fn()
        #     }, error = function(err) {
        d <- read.table("./static/example.txt", sep="\t", header=FALSE, colClasses=c("character", "character", "numeric", "character"), col.names=c("rsid", "chrom", "position", "genotype"))
        #     })
        incProgress(0.1)
        d$chrom = ordered(d$chrom, levels=c(seq(1, 22), "X", "Y", "MT"))
        incProgress(0.1)
        setProgress(0.75)
        ggplot(d, aes(x=chrom, fill=factor(chrom))) + ggtitle("Number of mutations per chromosome") + geom_bar() + theme(legend.position = "none")
        
      })
    })
  output$riskTable <- renderDataTable({
    withProgress(message="Finding SNP associations", value=0.1, {
      d <- tryCatch({
        load.file.fn()
      }, error = function(err) {
        d <- read.table("./static/example.txt", sep="\t", header=FALSE, colClasses=c("character", "character", "numeric", "character"), col.names=c("rsid", "chrom", "position", "genotype"))
      })
      incProgress(0.1)
      d$chrom = ordered(d$chrom, levels=c(seq(1, 22), "X", "Y", "MT"))
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
      levels(d$chrom) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")
      my.snps <- with(d, GRanges(seqnames=chrom, 
                                 IRanges(start=position, width=1), 
                                 rsid=rsid, genotype=genotype))
      incProgress(0.1)
      data(gwrngs19)
      gwrngs.emd <- as.data.frame(elementMetadata(gwrngs19))
      incProgress(0.1)
      dm <- merge(d, gwrngs.emd, by.x="rsid", by.y="SNPs")
      risk.alleles <- gsub("[^\\-]*-([ATCG?])", "\\1", dm$Strongest.SNP.Risk.Allele)
      i.have.risk <- mapply(function(risk, mine) {
        risk %in% unlist(strsplit(mine, ""))
      }, risk.alleles, dm$genotype)
      incProgress(0.1)
      dm$i.have.risk <- i.have.risk
      my.risk <- dm[dm$i.have.risk, ]
      rel.cols <- c(colnames(d), "Disease.Trait", "Risk.Allele.Frequency",
                    "p.Value", "i.have.risk", "X95..CI..text.", "Initial.Sample.Size")
      my.risk.ord <- my.risk[order(my.risk$Risk.Allele.Frequency), rel.cols]
      risk.table <- my.risk.ord[,c("genotype", "Disease.Trait", "p.Value", "Risk.Allele.Frequency", "Initial.Sample.Size")]
      setProgress(1)
      colnames(risk.table) <- c("Genotype", "Disease Trait", "P Value","Risk Allele Frequency", "Population Sampled")
      risk.table
    })
  }, options = list(pageLength = 10))
  
  
})