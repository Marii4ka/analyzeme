library(shiny)
library(jsonlite)
library(digest)
library(SNPRelate)
library(gwascat)
library(ggbio)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

options(shiny.maxRequestSize = -1)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
cols <- c("#006BA4","#FF800E","#A2C8EC","#898989","#ABABAB","#595959","#5F9ED1","#CFCFCF","#FFBC79","#C85200","#0063A4","#FF400E","#A228EC","#fd325a","#ABAB4B","#595459","#579ED1","#C6CFCF","#FF6C79","#C55200","#FF1C79","#C81200","#006314","#F1400E","#1228EC","#7dff18","#138503")
cols.alpha <- add.alpha(cols,0.7)
#palette(cols.alpha)

shinyServer(function(input, output) {
  temp.file <- digest(date())
  
  load.file.fn <- reactive({
    data <- input$inFile
    path <- as.character(data$datapath)
    if(length(path) == 0) {
      stop()
    }
    system(paste("unzip -p ", path, " > ./temp/", temp.file, sep=""))
    d <- read.table(paste("./temp/", temp.file, sep=""), sep="\t", header=FALSE, colClasses=c("character", "character", "numeric", "character"), col.names=c("rsid", "chrom", "position", "genotype"))
    
  })
  
  unzip.file.fn <- reactive({
    data <- input$inFile
    path <- as.character(data$datapath)
    if(length(path) == 0) {
      stop()
    }
    system(paste("unzip -p ", path, " > ./temp/", temp.file, sep=""))
    d <- paste("./temp/", path, sep="")
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
  
  output$allPCA.image <- renderUI({
    images <- c("http://www.i2symbol.com/images/abc-123/o/white_smiling_face_u263A_icon_256x256.png")
    #images <- c("./static/pca.all.png")
    tags$img(src= images)
    })
    
  output$allPCA.plot <- renderPlot({
    withProgress(message="Performing principle component analysis", value=0, {
        d <- tryCatch({
          unzip.file.fn()
        }, error = function(err) {
          d <- "./static/example.txt"
        })
        temp.path <- paste("./temp/", temp.file, sep="")
        system(paste("time ./bin/plink --23file ", d, " --snps-only no-DI --chr 22 --make-bed --out ", temp.path, sep=""))
        incProgress(0.1)
        system(paste("time ./bin/plink --bfile ./static/all.chr.22.clean --bmerge ",temp.path, ".bed ",temp.path, ".bim ",temp.path, ".fam --make-bed --out ", temp.path, ".merge", sep=""))
        bed.fn <- paste(temp.path,".merge.bed", sep="")
        fam.fn <- paste(temp.path,".merge.fam", sep="")
        bim.fn <- paste(temp.path,".merge.bim", sep="")
        snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, paste(temp.path, ".merge.gds", sep=""))
        incProgress(0.1)
        genofile <- snpgdsOpen(paste(temp.path, ".merge.gds", sep=""))
        incProgress(0.1)
        pca <- snpgdsPCA(genofile, num.thread=2) #long
        incProgress(0.4)
        pc.percent <- pca$varprop*10
        
        pops<-read.csv("./static/popids.csv", header=T, stringsAsFactors = FALSE)
        selectedRows <- (pops$Sample %in% pca$sample.id) 
        pops.present <- pops[selectedRows,]
        pops.present <- rbind(pops.present, c("You", "You"))
        pop<-factor(c(pops.present$Population))#[match(pops.present$Sample,pca$sample.id)]
        pop[2505] <- "You"
        tab <- data.frame(sample.id = pca$sample.id,
                          EV1 = pca$eigenvect[,1],    # the first eigenvector
                          EV2 = pca$eigenvect[,2],    # the second eigenvector
                          pop = pop,
                          stringsAsFactors = FALSE)
        
        ### all plot --------
        setIncrement(1)
        plot(tab$EV2, tab$EV1, xlab="PC 2", ylab="PC 1", main="Population clustering and you", pch=19, col=as.integer(tab$pop))
        legend("topleft", legend=unique(tab$pop), pch=19, border=NA, col=as.integer(unique(tab$pop)), cex=0.7)
        arrows(tab$EV2[2505]-1e-2, tab$EV1[2505]+1e-2, tab$EV2[2505]-1.1e-4, tab$EV1[2505]-1.1e-4, col=add.alpha("#EE0000",0.7),lwd=1.5)
        text(tab$EV2[2505]-1.1e-2, tab$EV1[2505]+1.1e-2,"You")
    })
  })
  
  
})