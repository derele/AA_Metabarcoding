## library(devtools)
## devtools::install_github("derele/MultiAmplicon")

library(MultiAmplicon)
library(phyloseq)
library(dada2)
library(pheatmap)
## library(ggplot2)
library(reshape)
library(DECIPHER)
library(parallel)
library(phangorn)
## library(plyr)
library(dplyr)
## library(gridExtra)
## library(nlme)
## library(structSSI)

a.files <- list.files(path="/SAN/Metabarcoding/AA_combi/all_fastq_raw",
                      pattern=".fastq.gz", full.names=TRUE)

Ffq.file <- a.files[grepl("R1", a.files)]
Rfq.file <- a.files[grepl("R2", a.files)]

## asses the Quality
## plotQualityProfile(Rfq.file[[20]])

samples <- gsub("/SAN/Metabarcoding/AA_combi/all_fastq_raw/(S\\d+).*?\\.fastq\\.gz", "\\1", Ffq.file)

filt_path <- "/SAN/Metabarcoding/AA_combi/DaDaFiltN"
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples

filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## ## Filter # only run when new filtered data is needed
## for(i in seq_along(Ffq.file)) {
##     fastqPairedFilter(c(Ffq.file[i], Rfq.file[i]), c(filtFs[i], filtRs[i]),
##                       truncLen=c(170,170), 
##                       maxN=0, maxEE=2, truncQ=2, 
##                       compress=TRUE, verbose=TRUE)
## }

## obtain the primer data
file.ptable <- "/SAN/Metabarcoding/AA_combi/sort_amplicons_input.csv"
ptable <- read.csv(file.ptable, sep=",", header=TRUE)

primerF <- as.character(ptable[,3])
primerR <- as.character(ptable[,4])

names(primerF) <- as.character(ptable[,1])
names(primerR) <- as.character(ptable[,2])

rownames(ptable) <- apply(ptable[,3:4], 1, paste, collapse=":")
ptable$pRnames <- apply(ptable[,1:2], 1, paste, collapse=":")

files <- PairedReadFileSet(filtFs, filtRs)

primers <- PrimerPairsSet(primerF, primerR)

MA <- MultiAmplicon(primers, files)

MA1 <- sortAmplicons(MA)

pdf("figures/primers_MA_sorted.pdf", 
    width=25, height=15, onefile=FALSE)
cluster <- plot_Amplicon_numbers(rawCounts(MA1))
dev.off()

MA2 <- derepMulti(MA1)

MA3 <- dadaMulti(MA2, err=NULL, selfConsist=TRUE,
                 multithread=TRUE)

MA4 <- MultiAmplicon:::mergeMulti(MA3, justConcatenate=TRUE)
                 
MA5 <- MultiAmplicon:::sequenceTableMulti(MA4)

MA6 <- MultiAmplicon:::noChimeMulti(MA5, mc.cores=20)


########## From here: work with MA package finished #### 

names(MA6@sequenceTableNoChime) <- rownames(MA6)
STnoC <- MA6@sequenceTableNoChime

all.dada.seq <- DNAStringSet(unlist(lapply(STnoC, colnames)))

## writeFasta(all.dada.seq, "/SAN/Metabarcoding/AA_combi/all_dada.fasta")

## the level percent of non-Bimera sequences
lapply(seq_along(MA6@sequenceTable), function (i){
    sum(STnoC[[i]]/sum(MA6@sequenceTable[[i]]))
})

sumSample <- lapply(STnoC, rowSums)
dadaMapped <- melt(sumSample)
dadaMapped$sample <- unlist(lapply(sumSample, names))

Dmap <- t(cast(dadaMapped, sample~L1))
Dmap[is.na(Dmap)] <- 0

pdf("figures/primers_dadaMap.pdf",
    width=15, height=15, onefile=FALSE)
pheatmap(log10(Dmap+1))
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

pdf("figures/dada_sample_exclusion_test.pdf", width=24, height=12)
cluster.table.sample <- plot.exclude.clusters(Dmap, 2)
dev.off()

Samples.to.exclude <- names(cluster.table.sample[cluster.table.sample==2])


pdf("figures/dada_primer_exclusion_test.pdf", width=24, height=12)
cluster.table.primer <- plot.exclude.clusters(t(Dmap), 2)
dev.off()

Primers.to.exclude <- names(cluster.table.primer[cluster.table.primer==2])



########## Remove failed samples ###########################    
tSTnoC <- lapply(STnoC, function(x) as.data.frame(t(x)))
all.otu.counts <- bind_rows(tSTnoC)
all.otu.counts[is.na(all.otu.counts)] <- 0


## this figure is disturbing... see what happened!
## is it just the "data crunching" above?
pdf("figures/full_otu_dadaMap.pdf", width=15, height=15)
pheatmap(log10(all.otu.counts+1))
dev.off()

## follow this 
## https://f1000research.com/articles/5-1492/v2

## phylogenetic trees 
align.seqtab <- function (seqtab){
    seqs <- getSequences(seqtab)
    names(seqs) <- seqs 
    alignment <- AlignSeqs(RNAStringSet(DNAStringSet(seqs)), anchor=NA)
}

alignments <- mclapply(STnoC, align.seqtab, mc.cores=20)



get.tree.from.alignment <- function (alignment){
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm) # Note, tip order != sequence order
    fit <- pml(treeNJ, data=phang.align)
    ## negative edges length changed to 0!
    fitGTR <- update(fit, k=4, inv=0.2)
    fitGTR <- optim.pml(fitGTR, model="GTR",
                        optInv=TRUE, optGamma=TRUE,
                        rearrangement = "stochastic",
                        control = pml.control(trace = 0))
    return(fitGTR)
}


tree.l <- mclapply(alignments, get.tree.from.alignment, mc.cores=20)

##### TAXONOMY assignment
assign.full.tax <- function(seqtab){
    seqs <- getSequences(seqtab)
    taxa <- assignTaxonomy(seqs, "/SAN/db/RDP/Silva_123/SILVA_123_dada2_exp.fasta")
    return(taxa)
    ## taxaS <- addSpecies(taxa, "/SAN/db/RDP/silva_species_assignment_v123.fa.gz", verbose=TRUE)
    ## return(taxaS)
}

set.seed(100) # Initialize random number generator for reproducibility
tax.l <- mclapply(STnoC, assign.full.tax, mc.cores=20)

lapply(tax.l, function (x) {
    rownames(x) <- NULL
    head(x)
})


## save(MA6, file="/SAN/Metabarcoding/allMA.Rdata")

## save(Samples.to.exclude, file="/SAN/Metabarcoding/exclude.Rdata")
## save(STnoC, file="/SAN/Metabarcoding/table.Rdata")
## save(tax.l, file="/SAN/Metabarcoding/taxa.Rdata")
save(alignments, file="/SAN/Metabarcoding/align.Rdata")
## save(tree.l, file="/SAN/Metabarcoding/trees.Rdata")

