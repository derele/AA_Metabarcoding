library(ggplot2)
library(pheatmap)
library(data.table)
library(RSvgDevice)

#################### Analyse the OTU data alone #####################

get.otu.files.combi <- function (path){
    OTU.files <- list.files(path=path,
                            pattern="*.fastq.otu_table.txt", full.names=TRUE)
    OTU.all <- lapply(OTU.files, read.table, header=TRUE)
    regex <- paste0(path, "/(.*?)\\.fastq\\.otu_table\\.txt")
    names(OTU.all) <- gsub(regex,
                           "\\1", OTU.files)
    OTU.all <- lapply(names(OTU.all), function (name){
        OTUID <- paste(OTU.all[[name]]$OTUID, name, sep="|")
        cbind(OTUID, OTU.all[[name]][, -1])
    })
    all.samples <- unique(unlist(lapply(OTU.all, colnames)))
    ## fill all columns for all samples
    OTU.l <- lapply(OTU.all, function (table){
        new.table <- as.data.frame(matrix(table$OTUID, ncol=1))
        new.table[, all.samples] <- 0
        new.table[,colnames(table)] <- table
        return(new.table[, all.samples])
    })
    ## Make it an all numeric data
    OTU.all <- do.call(rbind, OTU.l)
    rownames(OTU.all) <- OTU.all$OTUID
    OTU.all$OTUID <- NULL
    return(OTU.all)
}

OTUs.uc  <- get.otu.files.combi("/SAN/Metabarcoding/Hyena/second/sorted_amps/usearch/")
OTUs.Luis.uc  <- get.otu.files.combi("/SAN/Metabarcoding/Hyena/first_good/sorted_amplicons/usearch/")

sum(OTUs.uc)
sum(OTUs.Luis.uc)

pdf("figures/Hyena_OTU_heat.pdf", width=14, height=14)
pheatmap(log10(OTUs.uc+1),
         show_rownames=FALSE,
         show_colnames=FALSE,
         annotation_col=data.frame(row.names=colnames(OTUs.uc),
                                   is.control=as.numeric(
                                       grepl("H2O|Argave|Wolf|Paramix",
                                             colnames(OTUs.uc)))),
         annotation_legend = FALSE)
dev.off()

png("figures/Hyena_OTU_heat.png", res=300, width = 1480, height = 1480)
pheatmap(log10(OTUs.uc+1),
         show_rownames=FALSE,
         show_colnames=FALSE,
         treeheight_row=0,
         treeheight_col=0,
         annotation_col=data.frame(row.names=colnames(OTUs.uc),
                                   is.control=as.numeric(
                                       grepl("H2O|Argave|Wolf|Paramix",
                                             colnames(OTUs.uc)))),
         annotation_legend = FALSE)
dev.off()


## a crude background reduction by setting  all counts below an outlier detection to zero
## http://stats.stackexchange.com/questions/56402/detecting-outliers-in-count-data
bg.correct <- function (OTU){
    out.z <- function(x){
        trans <- log10(as.numeric(x))
        ## a trick to not assess the distribution of zeros, ones 
        ## assumed here to be true negatives
        NN <- which(trans>0)
        rob.z <- (trans-median(trans[NN]))/mad(trans[NN])
        z.outl <- which(!rob.z>quantile(rob.z[NN], 0.05, na.rm=TRUE))
        return(z.outl)
    }
    for(i in 1:nrow(OTU)){
        OTU[i, out.z(OTU[i,])] <- 0
    }
    return(OTU)
}

OTUs <- bg.correct(OTUs.uc)
## 90,427 removed 06/09/2016
sum(OTUs.uc)-sum(OTUs)
## 139,211 removed 25/10/2016

OTUs.Luis <- bg.correct(OTUs.Luis.uc)
sum(OTUs.Luis.uc)-sum(OTUs.Luis)
## 77,099 removed 25/10/2016

png("figures/Hyena_OTU_heat_BCcor.png", res=300, width = 1480, height = 1480)
pheatmap(log10(OTUs+1),
         show_rownames=FALSE,
         show_colnames=TRUE,
         treeheight_row=0,
         treeheight_col=0,
         annotation_col=data.frame(row.names=colnames(OTUs),
                                   is.control=as.numeric(
                                       grepl("H2O|Argave|Wolf|Paramix",
                                             colnames(OTUs)))),
         annotation_legend = FALSE)
dev.off()


png("figures/Hyena_OTU_heat_BCcor_Luis.png", res=300, width = 1480, height = 1480)
pheatmap(log10(OTUs.Luis+1),
         show_rownames=FALSE,
         show_colnames=TRUE,
         treeheight_row=0,
         treeheight_col=0,
         annotation_legend = FALSE)
dev.off()





######## Analyse just the READ counts per amplicon and sample #####
sum(OTUs)
amplicon <- gsub("OTU\\d+\\|", "", rownames(OTUs))
SUM.amp <- do.call(cbind, by(OTUs, amplicon, colSums))

devSVG("figures/Hyena_Primer_heat.svg", width=14, height=14)
pheatmap(log10(t(SUM.amp)+1),
         show_rownames=FALSE,
         show_colnames=FALSE,
         treeheight_row=0,
         treeheight_col=0,
         annotation_col=data.frame(row.names=rownames(SUM.amp),
                                   is.control=as.numeric(
                                       grepl("H2O|Argave|Wolf|Paramix",
                                             rownames(SUM.amp)))),
         annotation_legend = FALSE)
dev.off()


get.cluster.tree.df <- function(data, k){
   clusdat <- get.scaled.and.clustered(data)
   Cluster <- factor(unlist(cutree(clusdat, k = k)))
   as.data.frame(Cluster)
}


plot.exclude.clusters <- function(OTUs) {
    distOTUs <- hclust(dist(t(log10(OTUs+1))))
    clusCatOTU <- factor(unlist(cutree(distOTUs, k = 2)))
    ## plot
    plot(distOTUs)
    ## add cluster highlighting
    rect.hclust(distOTUs, k=2)
    return(clusCatOTU)
}

pdf("figures/sample_exclusion_test.pdf", width=15, height=15)
cluster.table <- plot.exclude.clusters(OTUs)
dev.off()


## for our second run samples exclude cluster 2
Samples.to.exclude <- names(cluster.table[cluster.table==2])
OTUs <- OTUs[, !colnames(OTUs)%in%Samples.to.exclude]



## For Luis's run nothing excluded
pdf("figures/sample_exclusion_test_Luis.pdf", width=15, height=15)
cluster.table.Luis <- plot.exclude.clusters(OTUs.Luis)
dev.off()

## Just sum the replicates

##    At some point use more scruteny (from
##    https://support.bioconductor.org/p/61234/): Thanks to the
##    underlying theory behind dispersion estimation, you can easily
##    test whether your technical replicates really do represent
##    technical replicates. Specifically, read counts in technical
##    replicates should follow a Poisson distribution, which is a
##    special case of the negative binomial with zero dispersion. So,
##    simply fit a model using edgeR or DESeq2 with a separate
##    coefficient for each group of technical replicates. Thus all the
##    experimental variation will be absorbed into the model
##    coefficients and the only thing left will be the technical
##    variability of the replicates. For true technical replicates,
##    the dispersion should be zero for all genes. So if you estimate
##    dispersions using this model, and plotBCV/plotDispEsts shows the
##    dispersion very near to zero, then you can be confident that you
##    really have technical replicates. If the dispersion is nonzero,
##    then there is some additional source of unaccounted-for
##    variation. I have used this method on a pilot dataset with
##    several technical replicates for each condition. edgeR said the
##    dispersion was something like 10^-3 or less for all genes except
##    for the very low-expressed genes. -Ryan

##    idea: use this to exclude OTUs (genes) that are not credible
##    i.e. show more than technical variation. 
   
sum.columns <- function(x){
    reps <- as.factor(substring(colnames(x), 1, 4))
    y <- do.call(rbind, by(t(x), reps, colSums))
    t(y)
}
   
## mean between replicates
SOTUs <- sum.columns(OTUs)
   
## Luis does not have replicats
SOTUs.Luis <- OTUs.Luis
