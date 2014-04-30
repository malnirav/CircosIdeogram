#### Make ideogram from datasets ####
library(ggbio)
library(gridExtra)
library(hiAnnotator)
library(BSgenome.Hsapiens.UCSC.hg18)
library(biovizBase)

freeze <- "hg18"

#### get data/tracks ####
## get point data aka sites ##
## the data is processed from PMID: 24489369 & 21561906 ##
load("sites.RData")
sites.gr <- makeGRanges(sites, soloStart=TRUE, chromCol="Chr", freeze='hg18')
sites.gr$Alias <- ifelse(grepl("HIV",sites.gr$setName), "HIV", "MLV")
sites.gr$y <- runif(length(sites.gr), 0, 3)

## get genes ##
genes <- getUCSCtable("refFlat","RefSeq Genes")
genes <- unique(makeGRanges(genes, freeze='hg18'))
genes <- keepSeqlevels(genes,seqlevels(sites.gr))

## get genome data ##
hg18Ideo <- getIdeogram("hg18",cytoband=TRUE)
hg18Ideo <- keepSeqlevels(hg18Ideo, seqlevels(sites.gr))

## get gc perc ##
gcperc <- getUCSCtable("gc5Base","GC Percent")
gcperc$strand <- "*"
gcperc <- makeGRanges(unique(gcperc), chromCol="chrom", freeze='hg18')
gcperc <- keepSeqlevels(gcperc, seqlevels(sites.gr))
gcperc$ratio <- gcperc$sumData/gcperc$validCount

#### format names/data ####
seqlevelsStyle(hg18Ideo) <- "NCBI"
seqlevelsStyle(sites.gr) <- "NCBI"
seqlevelsStyle(genes) <- "NCBI"
seqlevelsStyle(gcperc) <- "NCBI"

#### create reduced tracks ####
## reduced gene track ##
strand(genes) <- "*"
genes.reduced <- reduce(genes, min.gapwidth = 100000)
mcols(genes.reduced)$counts <- countOverlaps(genes.reduced, genes, 
                                             ignore.strand=TRUE)
mcols(genes.reduced)$logCount <- log(mcols(genes.reduced)$counts)

## reduced GC track ##
gc.reduced <- do.call(c,tile(reduce(hg18Ideo),width=1000000))
seqinfo(gc.reduced) <- seqinfo(hg18Ideo)
gc.reduced <- getFeatureCounts(gc.reduced, gcperc, colnam="gc", 
                               widths=c(0), weightsColname="ratio")
gc.reduced$counts <- rescale(gc.reduced$gc.0bp, to=c(0,1))

#### circle plot ####
hg18Ideo$y <- 1
cytobandColor <- getOption("biovizBase")$cytobandColor
p <- ggplot() + layout_circle(reduce(hg18Ideo), geom = "ideo", fill = "gray70",
                              radius = 34, trackWidth = 2.5) 
p <- p + layout_circle(hg18Ideo, geom = "rect", aes(fill=gieStain, y=y), color=NA,
                       radius = 34, trackWidth = 2.4, group.selfish=FALSE, 
                       stat="identity", alpha=.8) +
  scale_fill_manual(values=cytobandColor)
p <- p + layout_circle(data=reduce(hg18Ideo), geom = "text", aes(label = seqnames), 
                       vjust = 0, radius = 37, trackWidth = 1)
p <- p + layout_circle(data=gc.reduced, aes(y=counts), geom = "bar", 
                       fill = "#4daf4a", color = NA, radius = 29, 
                       trackWidth = 5.5, size=0.2)
p <- p + layout_circle(data=genes.reduced, aes(y=logCount), geom = "bar", 
                       fill = "#e41a1c", color = NA, radius = 24, 
                       trackWidth = 5.5, size=0.2)
p <- p + layout_circle(data=subset(sites.gr,Alias=="HIV"), aes(y=y), 
                        color="#ff7f00", geom = "point", radius = 17, 
                        trackWidth = 7.2, size=0.3)
p <- p + layout_circle(data=subset(sites.gr,Alias!="HIV"), aes(y=y), 
                         color="#984ea3", geom = "point", radius = 10, 
                         trackWidth = 7.2, size=0.3)

dataset <- c("GC"="#4daf4a", "Genes"="#e41a1c", "HIV"="#ff7f00", "MLV"="#984ea3")
p2 <- qplot(data=data.frame(Track=factor(names(dataset),
                                         levels=names(dataset))),1,Track,
            fill=Track, geom="tile") + scale_fill_manual(name="",values=dataset)

g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
mylegend <- g_legend(p2)

# Ciros_plot #
grid.arrange(arrangeGrob(p + theme(legend.position="none")),
                   mylegend, ncol=2, widths=c(10, 1))
