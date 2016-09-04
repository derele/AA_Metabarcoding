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

TAX.raw <- read.csv("/SAN/Metabarcoding/Hyena/second/sorted_amps/usearch/ALL_outs.taxtable",
                sep = ",")

TAX.raw$query <- gsub(".fastq.otus.fa", "", TAX.raw$query)

TAX.raw$amplicon <- gsub("OTU\\d+\\|", "", TAX.raw$query)

## unique betst bitscore for every query
T.l <- by(TAX.raw, TAX.raw$query,  function (x) {
     b.bit <- max(x$bitscore)
     all.best <- x[x$bitscore==b.bit, ]
     ### A little last common ancestor play here... BUT wait ... 
     ## a shortcut throwing out OTUs that don't agree at least on the
     ## family level and allowing only the best hit afterwards
     u.family <- unique(all.best$family)
     if(length(u.family)==1){
         return(all.best)
     }
})


TAX <- do.call(rbind, T.l)
head(TAX[order(TAX$amplicon), ])
rownames(TAX) <- NULL

## Only consider Euks now
TAX <- TAX[TAX$superkingdom%in%"Eukaryota", ]

## remove some really weird stuff FIND later out where the errors
## are!! Database errors... 
table(TAX$phylum)

TAX <- TAX[!TAX$phylum%in%c("Cnidaria", "Porifera",
                            "Bacillariophyta",  ## maybe okay?
                            "Eustigmatophyceae" ## maybe okay?
                            ) ,]

table(TAX$phylum)

### Summarizing by class
foo <- merge(TAX, OTUs, by.x="query", by.y=0)
foobar <- do.call(rbind, by(foo[, 18:ncol(foo)], foo$class, colSums))
foobar <- foobar[order(rowSums(foobar), decreasing=TRUE), ]
foobar <- foobar[!rownames(foobar)%in%c("", "undef"), ]
foobar <- foobar[, !grepl("H2O|Argave|Wolf", colnames(foobar))]

foobar <- foobar[rowSums(foobar)>116, ]

mean.columns <- function(x){
  reps <- as.factor(gsub("_S\\d+$", "", colnames(x)))
  y <- do.call(rbind, by(t(x), reps, colMeans))
  t(y)
}

baz <- mean.columns(foobar)
baz[baz<1] <- 0

devSVG("figures/Hyena_class_sum_heat.svg", width=7, height=7)
pheatmap(log10(baz+1),
         show_rownames=TRUE,
         show_colnames=FALSE,
         treeheight_row=0,
         treeheight_col=0)
dev.off()


