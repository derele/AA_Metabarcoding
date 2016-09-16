library(ggplot2)
library(pheatmap)
library(reshape2)
library(parallel)
library(data.table)
library(RSvgDevice)

#################### Analyse the OTU data alone #####################


OTU.files <- list.files(path="/SAN/Metabarcoding/Hyena/second/sorted_amps/usearch/",
                        pattern="*.fastq.otu_table.txt", full.names=TRUE)

OTU.all <- mclapply(OTU.files, read.table, header=TRUE)

names(OTU.all) <- gsub("/SAN/Metabarcoding/Hyena/second/sorted_amps/usearch//(.*?)\\.fastq\\.otu_table\\.txt",
                       "\\1", OTU.files)


OTU.all <- lapply(names(OTU.all), function (name){
    OTUID <- paste(OTU.all[[name]]$OTUID, name, sep="|")
    cbind(OTUID, OTU.all[[name]][, -1])
})


## fill all columns for all samples
all.samples <- unique(unlist(lapply(OTU.all, colnames)))

OTU.all <- lapply(OTU.all, function (table){
    new.table <- as.data.frame(matrix(table$OTUID, ncol=1))
    new.table[, all.samples] <- 0
    new.table[,colnames(table)] <- table
    return(new.table[, all.samples])
})


OTUs <- do.call(rbind, OTU.all)

## Make it an all numeric data
rownames(OTUs) <- OTUs$OTUID
OTUs$OTUID <- NULL

sum(OTUs)

devSVG("figures/Hyena_OTU_heat.svg", width=14, height=14)
pheatmap(log10(OTUs+1),
         show_rownames=FALSE,
         show_colnames=FALSE,
         treeheight_row=0,
         treeheight_col=0,
         annotation_col=data.frame(row.names=colnames(OTUs),
                                   is.control=as.numeric(
                                       grepl("H2O|Argave|Wolf|Paramix",
                                             colnames(OTUs)))),
         annotation_legend = FALSE)
dev.off()

png("figures/Hyena_OTU_heat.png", res=300, width = 1480, height = 1480)
pheatmap(log10(OTUs+1),
         show_rownames=FALSE,
         show_colnames=FALSE,
         treeheight_row=0,
         treeheight_col=0,
         annotation_col=data.frame(row.names=colnames(OTUs),
                                   is.control=as.numeric(
                                       grepl("H2O|Argave|Wolf|Paramix",
                                             colnames(OTUs)))),
         annotation_legend = FALSE)
dev.off()


## a crude background reduction by setting  all counts below an outlier detection to zero
## http://stats.stackexchange.com/questions/56402/detecting-outliers-in-count-data
out.z <- function(x){
    trans <- log10(as.numeric(x))
    ## a trick to not assess the distribution of zeros, ones and twos
    ## assumed here to be true negatives
    NN <- which(trans>0)
    rob.z <- (trans-median(trans[NN]))/mad(trans[NN])
    z.outl <- which(!rob.z>quantile(rob.z[NN], 0.05, na.rm=TRUE))
}

for(i in 1:nrow(OTUs)){
    OTUs[i, out.z(OTUs[i,])] <- 0
}

sum(OTUs)

## 90,427 removed 06/09/2016

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

sum(OTUs)

amplicon <- gsub("OTU\\d+\\|", "", rownames(OTUs))

######## Analyse just the READ counts per amplicon and sample #####
SUM.amp <- do.call(cbind, by(OTUs, amplicon, colSums))

## nummber of reads overall
sum(SUM.amp)


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

TAX.raw <-
    read.csv("/SAN/Metabarcoding/Hyena/second/sorted_amps/usearch/ALL_outs_nt.taxtable",
             sep = ",")

TAX.raw$query <- gsub(".fastq.otus.fa", "", TAX.raw$query)

TAX.raw$amplicon <- gsub("OTU\\d+\\|", "", TAX.raw$query)

## unique betst bitscore for every query
T.l <- by(TAX.raw, TAX.raw$query,  function (x) {
     b.bit <- max(x$bitscore)
     all.best <- x[x$bitscore==b.bit, ]
     ### A little last common ancestor play here... BUT wait ... 
     ## a shortcut throwing out OTUs that don't agree at least on the
     ## class level and allowing only the best hit afterwards
     u.family <- unique(all.best$class)
     if(length(u.family)==1){
         return(all.best)
     }
})


TAX <- do.call(rbind, T.l)
tail(TAX[order(TAX$amplicon), ])
rownames(TAX) <- NULL



## Only consider Euks now
TAX <- TAX[TAX$superkingdom%in%"Eukaryota", ]

## remove some really weird stuff FIND later out where the errors
## are!! Database errors... 
table(TAX$phylum)

## should be fixed in database at some point... but now as a shortkut
## here
TAX <- TAX[!TAX$phylum%in%c("Cnidaria", "Porifera",
                            "Bacillariophyta",  ## maybe okay?
                            "Eustigmatophyceae" ## maybe okay?
                            ) ,]

## ## now use only best hit
TAX.best <- TAX[!duplicated(TAX$query), ]
## tail(TAX[order(TAX$amplicon), ])

## table(TAX$phylum)

### Summarizing by class
foo <- merge(TAX, OTUs, by.x="query", by.y=0)
foobar <- do.call(rbind, by(foo[, 18:ncol(foo)], foo$class, colSums))
foobar <- foobar[order(rowSums(foobar), decreasing=TRUE), ]
foobar <- foobar[!rownames(foobar)%in%c("", "undef"), ]
foobar <- foobar[, !grepl("H2O|Argave|Wolf", colnames(foobar))]

mean.columns <- function(x){
  reps <- as.factor(gsub("_S\\d+$", "", colnames(x)))
  y <- do.call(rbind, by(t(x), reps, colMeans))
  t(y)
}

## mean between replicates
baz <- mean.columns(foobar)
## removing stuff with very low support from only one replicate
baz[baz<1] <- 0

## remove lowly represented classes
baz <- baz[rowSums(baz)>20, ]


devSVG("figures/Hyena_class_sumNT_heat.svg", width=7, height=7)
pheatmap(log10(baz+1),
         show_rownames=TRUE,
         show_colnames=FALSE,
         treeheight_row=0,
         treeheight_col=0)
dev.off()


### preliminary single species
best.merge <- merge(TAX.best, OTUs, by.x="query", by.y=0)


## Ancylostoma
get.level <- function(regex, level="species"){
    subs <- best.merge[grepl(regex, best.merge[, level]), ]
    subs <- by(subs[,18:ncol(subs)], as.character(subs[, level]), colSums)
    do.call(rbind, subs)
}


teaser <- rbind(get.level("Ancylo"), 
                get.level("isospora"),
                get.level("Sarcocystis buffalonis|Sarcocystis .*canis$"))


write.table(teaser, "~/Dropbox/Hyena_2ndDataset_teaser_reps.csv", sep=",")

write.table(mean.columns(teaser), "~/Dropbox/Hyena_2ndDataset_teaser_merged.csv", sep=",")

write.table(rownames(mean.columns(teaser)), "~/Dropbox/Hyena_lables.txt")
