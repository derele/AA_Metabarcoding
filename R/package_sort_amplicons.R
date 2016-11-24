## library(devtools)
devtools::install_github("derele/MultiAmplicon")

library(MultiAmplicon)
library(dada2)
library(phyloseq)
library(pheatmap)
library(ggplot2)
library(reshape)
library(DECIPHER)
library(parallel)
library(phangorn)

a.files <- list.files(path="/SAN/Metabarcoding/Hyena/second/fastq_raw/",
                      pattern=".fastq.gz", full.names=TRUE)

Ffq.file <- a.files[grepl("R1", a.files)]
Rfq.file <- a.files[grepl("R2", a.files)]

## asses the Quality
## plotQualityProfile(Rfq.file[[20]])

samples <- gsub("/SAN/Metabarcoding/Hyena/second/fastq_raw/(.*?)\\.fastq\\.gz", "\\1", Ffq.file)

filt_path <- "/SAN/Metabarcoding/Hyena/second/DaDafilt"
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))

## ## Filter # only run when new filtered data is needed
## for(i in seq_along(Ffq.file)) {
##     fastqPairedFilter(c(Ffq.file[i], Rfq.file[i]), c(filtFs[i], filtRs[i]),
##                       truncLen=c(200,200), 
##                       maxN=0, maxEE=2, truncQ=2, 
##                       compress=TRUE, verbose=TRUE)
## }

## obtain the primer data
file.ptable <- "/SAN/Metabarcoding/Hyena/second/sort_amplicons_input.csv"
ptable <- read.csv(file.ptable, sep=" ", header=TRUE)

primerF <- as.character(ptable[,3])
primerR <- as.character(ptable[,4])

names(primerF) <- as.character(ptable[,1])
names(primerR) <- as.character(ptable[,2])

rownames(ptable) <- apply(ptable[,3:4], 1, paste, collapse=":")
ptable$pRnames <- apply(ptable[,1:2], 1, paste, collapse=":")

filtFs.e <- filtFs[file.exists(filtFs)&file.exists(filtRs)]
filtRs.e <- filtRs[file.exists(filtFs)&file.exists(filtRs)]

files <- PairedReadFileSet(filtFs.e, filtRs.e)
primers <- PrimerPairsSet(primerF, primerR)

MA.Hy <- MultiAmplicon(primers, files)
MA.Hy <- sortAmplicons(MA.Hy)

## name it when putting the below in the package:
derepEmpty <- function(x){
    d <- derepFastq(x[file.info(x)$size>21])
    names(d) <- x[file.info(x)$size>21]
    return(d)
}

## different amplicons (in rows) dereplicated 
MultiDerepF  <- apply(MA.Hy@FstratifiedFiles, 1, derepEmpty)
names(MultiDerepF) <- names(MA.Hy@PrimerPairsSet)

## check empty are removed
MultiDerepF <- MultiDerepF[unlist(lapply(MultiDerepF, length)>0)]

## same for reverse reads
MultiDerepR  <- apply(MA.Hy@RstratifiedFiles, 1, derepEmpty)
names(MultiDerepR) <- names(MA.Hy@PrimerPairsSet)

## check empty are removed
MultiDerepR <- MultiDerepR[unlist(lapply(MultiDerepR, length)>0)]

## list of dada sample inference objects
MultiDadaF  <- lapply(MultiDerepF, function(x){
    dada(x, err=NULL, selfConsist=TRUE)
})

## list of dada sample inference objects
MultiDadaR  <- lapply(MultiDerepR, function(x){
    dada(x, err=NULL, selfConsist=TRUE)
})


mergers <- lapply(seq_along(MultiDadaF), function (i){
    mergePairs(MultiDadaF[[i]], MultiDerepF[[i]],
               MultiDadaR[[i]], MultiDerepR[[i]],
               justConcatenate=TRUE, verbose=TRUE)
})
names(mergers) <- names(MultiDadaF)

ST <- lapply(mergers,  makeSequenceTable)

lapply(ST, function (x){
    dim(x)
})

STnoC <- mclapply(ST, removeBimeraDenovo, verbose=TRUE,
                  mc.cores=20)

## tough fix fot the filenames at this point
STnoC <- lapply(STnoC, function(x){
   n <- gsub(".*?(.{4}_S\\d+)_L001_.*",  "\\1" ,
             rownames(x),
             perl=TRUE)
   rownames(x) <- n
   return(x)
})

lapply(STnoC, dim)

lapply(seq_along(ST), function (i){
    sum(STnoC[[i]]/sum(ST[[i]]))
})


lapply(ST, nrow)
lapply(STnoC, nrow)
lapply(STnoC, ncol)

sumSample <- lapply(STnoC, rowSums)
dadaMapped <- melt(sumSample)
dadaMapped$sample <- unlist(lapply(sumSample, names))

Dmap <- t(cast(dadaMapped, sample~L1))
Dmap[is.na(Dmap)] <- 0

pdf("figures/primers_MA_sorted.pdf",
    width=15, height=15, onefile=FALSE)
pheatmap(log10(MA.Hy@rawCounts+.1))
dev.off()

pdf("figures/primers_dadaMap.pdf",
    width=15, height=15, onefile=FALSE)
pheatmap(log10(Dmap+.1))
dev.off()

plot.exclude.clusters <- function(OTUs, k=2) {
    distOTUs <- hclust(dist(t(log10(OTUs+1))))
    clusCatOTU <- factor(unlist(cutree(distOTUs, k = k)))
    ## plot
    plot(distOTUs)
    ## add cluster highlighting
    rect.hclust(distOTUs, k=k)
    return(clusCatOTU)
}

pdf("figures/dada_sample_exclusion_test.pdf", width=15, height=15)
cluster.table <- plot.exclude.clusters(Dmap, 2)
dev.off()

Samples.to.exclude <- names(cluster.table[cluster.table==1])

########## Remove failed samples ###########################
library(dplyr)

tSTnoC <- lapply(STnoC, function(x) as.data.frame(t(x)))
all.otu.counts <- bind_rows(tSTnoC)
all.otu.counts[is.na(all.otu.counts)] <- 0

pdf("figures/full_otu_dadaMap.pdf", width=15, height=15)
pheatmap(log10(all.otu.counts+.1))
dev.off()


## phylogenetic trees 
align.seqtab <- function (seqtab){
    seqs <- getSequences(seqtab)
    names(seqs) <- seqs 
    alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
}

alignments <- mclapply(STnoC, align.seqtab, mc.cores=20)

get.tree.from.alignment <- function (alignment){
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm) # Note, tip order != sequence order
    fit <- pml(treeNJ, data=phang.align)
    ## negative edges length changed to 0!
    fitGTR <- update(fit, k=4, inv=0.2)
    fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                        rearrangement = "stochastic", control = pml.control(trace = 0))
    return(fitGTR)
}

trees <- mclapply(alignments, get.tree.from.alignment, mc.cores=20)

num.seq <- lapply(trees, function (x) length(x$tree$tip.label))
sum.edge.length <- lapply(trees, function (x) sum(x$tree$edge.length))
number.edges <- lapply(trees, function (x) length(x$tree$edge.length))
l.per.edge <- lapply(trees, function (x) sum(x$tree$edge.length)/length(x$tree$edge.length))

primer.overview <- as.data.frame(cbind(
    num.seq=unlist(num.seq),
    sum.edge.length=unlist(sum.edge.length),
    l.per.edge=unlist(l.per.edge)))


primer.overview <- merge(ptable, primer.overview, by=0)
rownames(primer.overview) <- primer.overview$Row.names
primer.overview$Row.names <- NULL

## Or load at this point for complete data up to here
## load("/home/ele/Sunday.Rdata")


##### TAXONOMY assignment
assign.full.tax <- function(seqtab){
    seqs <- getSequences(seqtab)
    taxa <- assignTaxonomy(seqs, "/SAN/db/RDP/Silva_123/SILVA_123_dada2.fasta")
    taxaS <- addSpecies(taxa, "/SAN/db/RDP/silva_species_assignment_v123.fa.gz", verbose=TRUE)
    return(taxaS)
}

set.seed(100) # Initialize random number generator for reproducibility
tax.l <- mclapply(STnoC, assign.full.tax, mc.cores=20)

## save(tax.l, file="/SAN/Metabarcoding/Hyena/second/taxa.Rdata")

lapply(tax.l, function (x) {
    rownames(x) <- NULL
    head(x)
})

## how many are assigned to Genus level

tax.tab <- sapply(c("Phylum", "Class", "Order", "Family", "Genus"), function (y){
    sapply(tax.l, function (x, level=y) {
        ## not just reporting the level above and not NA
        assigned <- !grepl("_\\w{2}", x[, level]) & !is.na(x[, level])
        length(assigned[assigned])/length(assigned)
    })
})

rownames(tax.tab) <- names(ps.l)
tax.tab[order(tax.tab[, "Genus"]), ]

primer.overview <- merge(primer.overview, tax.tab, by=0)
rownames(primer.overview) <- primer.overview$Row.names
primer.overview$Row.names <- NULL

######### FROM HERE completely Hyena specific ###################

Hyena.Cat <- read.csv("/home/ele/Dropbox/Animal_number_variable_codes_EH.csv")
Hyena.Cat <- as.data.frame(apply(Hyena.Cat, 2, function (x) gsub(" ", "", x)))

Hyena.Cat$V1 <- ifelse(Hyena.Cat$V1==1, "Male", "Female")
Hyena.Cat$V2 <- ifelse(Hyena.Cat$V2==2, "Cub", "Adult")
Hyena.Cat$V3 <- ifelse(Hyena.Cat$V3==1, "high", "low")
Hyena.Cat$V3[is.na(Hyena.Cat$V3)] <- "low" ## watchout not TRUE!!!

names(Hyena.Cat) <- c("animal", "ID", "sex", "age", "rank", "pack")
rownames(Hyena.Cat) <- Hyena.Cat$ID
Hyena.Cat$ID <- Hyena.Cat$animal <- NULL


get.sample.data <- function(seqtab.nochim){
    samples.out <- rownames(seqtab.nochim)
    samples.out <- samples.out[!samples.out%in%Samples.to.exclude]
    subject.ids <- gsub(".*?(.{4})_S\\d+.*?$", "\\1", samples.out)
    subject.df <- as.data.frame(cbind(Subject=samples.out, Subject.ids=subject.ids))
    subject.df <- subject.df[subject.df$Subject.ids%in%rownames(Hyena.Cat), ]
    subject.df <- merge(subject.df, Hyena.Cat,
                        by.x="Subject.ids", by.y=0)
    subject.df <- subject.df[!duplicated(subject.df$Subject), ]
    subject.df <- subject.df[!is.na(subject.df$Subject), ]
    rownames(subject.df) <- subject.df$Subject
    seqtab.phylo <- seqtab.nochim[!duplicated(rownames(seqtab.nochim)), ]
    seqtab.phylo <- seqtab.phylo[rownames(seqtab.phylo)%in%subject.df$Subject, ]
    ## order the two identically
    subject.df <- subject.df[rownames(seqtab.phylo), ]
    subject.df$rep <- ifelse(duplicated(subject.df$Subject.id), "r2", "r1")
    list(subject.df, seqtab.phylo)
}

sub.df.l <- lapply(STnoC, get.sample.data)

## Construct phyloseq object (straightforward from dada2 outputs)

ps.l <- lapply(seq_along(sub.df.l), function (i) {
    ps <- phyloseq(otu_table(sub.df.l[[i]][[2]], taxa_are_rows=FALSE),
                   sample_data(sub.df.l[[i]][[1]]),
                   phy_tree(trees[[i]]$tree),
                   tax_table(tax.l[[i]]))
    return(ps)
})

names(ps.l) <- names(sub.df.l)

Shannon.l <- lapply(ps.l, estimate_richness, measures="Shannon")

Shannon.sample.l <- lapply(seq_along(ps.l), function(i){
    merge(sample_data(ps.l[[i]]), Shannon.l[[i]], by=0)
})


sex.Shannon <- lapply(Shannon.sample.l, function(x){
    t.test(Shannon ~ sex, data=x)$p.value
})


age.Shannon <- lapply(Shannon.sample.l, function(x){
    t.test(Shannon ~ age, data=x)$p.value
})


rank.Shannon <- lapply(Shannon.sample.l, function(x){
    t.test(Shannon ~ rank, data=x)$p.value
})


rep.Shannon <- lapply(Shannon.sample.l, function(x){
    t.test(Shannon ~ rep, data=x)$p.value
})

Shannon.df <- cbind(sex=unlist(sex.Shannon),
                    age=unlist(age.Shannon),
                    rank=unlist(rank.Shannon),
                    rep=unlist(rep.Shannon))

rownames(Shannon.df) <- names(ps.l)

annotation <- merge(primer.overview, Shannon.df, by=0)
annotation$is.16.S <- ifelse(grepl("ADM|ACM|Klin", annotation$pRname), "16S", "18S")

pdf("figures/shannon.pdf", onefile=FALSE)
pheatmap(Shannon.df, cluster_cols=FALSE, cluster_rows=TRUE,
         color = colorRampPalette(c("firebrick3", "navy"))(4),
         breaks = 10^(0:-4),
         show_rownames=FALSE,
         legend = FALSE,
         display_numbers=TRUE,
         number_format = "%.1e",
         annotation_row = annotation[, c(6, 7, 13)])
dev.off()


primer.plot <- lapply(ps.l, function(x){
    R18S <- x
    R18S <- filter_taxa(R18S, function(x) mean(x)>0.2, TRUE)
    UF <- UniFrac(R18S, weighted=TRUE)
    UF[is.na(UF)] <- 0
    UFO <- ordinate(R18S, method="PCoA", distance=UF)
    UFOplot <- plot_ordination(R18S, UFO,
                               "samples", color="rank",
                               label="Subject.ids", axes=1:2) +
        ##    geom_point(size=5) +
        geom_path() +
        scale_colour_hue(guide = FALSE) +
        theme_bw()
    return(UFOplot)
})


library(gridExtra)

pdf("figures/UFO.pdf")
do.call("grid.arrange", c(primer.plot, ncol=3))
dev.off()

ps.bak1 <- ps.l[rownames(annotation)[annotation$is.16.S%in%"16S"]][[1]]

pdf("figures/Bacteria_richnesSex.pdf", width=14, height=8)
plot_richness(ps.bak1, x="sex") + theme_bw()  + geom_boxplot()
dev.off()

pdf("figures/Bacteria_richnesAge.pdf", width=14, height=8)
plot_richness(ps.bak1, x="age") + theme_bw()  + geom_boxplot()
dev.off()

pdf("figures/Bacteria_richnesRank.pdf", width=14, height=8)
plot_richness(ps.bak1, x="rank") + theme_bw()  + geom_boxplot()
dev.off()

pdf("figures/Bacteria_richnesPack.pdf", width=14, height=8)
plot_richness(ps.bak1, x="pack") + theme_bw()  + geom_boxplot()
dev.off()


############################### BARCHARTS
top60 <- names(sort(taxa_sums(ps.l[[1]]), decreasing=TRUE))[1:60]
ps.top60 <- transform_sample_counts(ps.l[[1]], function(OTU) OTU/sum(OTU))
ps.top60 <- prune_taxa(top60, ps.top60)

## plot_bar(ps.top20, x="rep", fill="Family") + facet_wrap(~V2, scales="free_x")

## plot_bar(ps.top20, x="V1", fill="Family") + facet_wrap(~V2, scales="free_x")

## plot_bar(ps.top20, x="V2", fill="Family") + facet_wrap(~V2, scales="free_x")


###########################################

enterotype <- subset_taxa(ps, Genus!= "-1")

### setting up distance methods to be tested

dist_methods <- unlist(distanceMethodList)

## 1:3 require tree
## dist_methods <- dist_methods[-(1:3)]

dist_methods <- dist_methods[-which(dist_methods=="ANY")]
plist <- vector("list", length(dist_methods))
names(plist) <- dist_methods

iDist <- distance(ps.l[[1]], method=dist_methods[[1]])
iMDS  <- ordinate(ps.l[[1]], "MDS", distance=iDist)



for( i in dist_methods[c(4, 5, 6, 11, 16, 22)] ){ # other give error
    ## Calculate distance matrix
    iDist <- distance(enterotype, method=i)
    ## Calculate ordination
    iMDS  <- ordinate(enterotype, "MDS", distance=iDist)
    ## Make plot
    ## Don't carry over previous plot (if error, p will be blank)
    p <- NULL
    ## Create plot, store as temp variable, p
    p <- plot_ordination(enterotype, iMDS, color="age", shape="rank")
    ## Add title to each plot
    p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
                                        # Save the graphic to file.
    plist[[i]] <- p
}






