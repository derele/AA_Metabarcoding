library(ggplot2)
## library(DESeq)
library(pheatmap)
library(reshape2)

#################### Analyse the OTU data alone #####################
OTUs <- read.csv("ALL_otus_samples_replaced.csv", sep = ";")
OTUs$X <- NULL # remove an empty variable  

## Make it an all numeric data
rownames(OTUs) <- OTUs$OTU 
OTUs$OTU <- NULL

## in the first Row we get a operational taxonomic unit (OTU) with a
## arbitrary number and the name of the amplicon from the multi
## amplicon PCR approach
amplicon <- gsub("OTU_\\d+_", "", rownames(OTUs))


######## Analyse just the READ counts per amplicon and sample #####
SUM.amp <- do.call(cbind, by(OTUs, amplicon, colSums))

## nummber of reads overall
sum(SUM.amp)

## display how many reads per Amplicon were obtained
pheatmap(SUM.amp)

## pdf("overall_heat.pdf", width=14, height=14)
pheatmap(log10(t(SUM.amp+1)))
## dev.off()

## There are some samples that cluster with water controls and seem to
## be failed regarding their overall read counts we can select from
## the hierarchical distance clustering to get these
d <- dist(log10(SUM.amp+1))
c <- hclust(d)

## pdf("without_nomralization_cluster.pdf", width=18, height=14)
plot(c)
## dev.off()

## FAILED samples
clusters <- cut(as.dendrogram(c), h=7)
str(clusters)
str(clusters[[2]][1])

failed <- clusters[[2]][1]
failed <- rownames(SUM.amp)[unlist(failed)]

## A matrix is easier to work with when you want to multiply over all
## rows using R's vectorization
D.mat <- as.matrix(t(OTUs))
per.sample <- rowSums(D.mat)
per.sample[order(per.sample)]

## a simple normalization by multiplying with the fraction of the top
## (what about median??) count per sample
norm.factor <- max(per.sample)/per.sample
D.mat.snorm <- D.mat*norm.factor

SUM.amp.norm <- do.call(cbind, by(t(D.mat.snorm), amplicon, colSums))

## display how many reads per Amplicon were obtained
pheatmap(SUM.amp.norm)

## pdf("overall_norm_heat.pdf", width=14, height=14)
pheatmap(log10(SUM.amp.norm+1))
## dev.off()

d.snorm <- dist(log10(SUM.amp.norm+1))
c.snorm <- hclust(d.snorm)
plot(c.snorm)

## What about the failed samples after normalization:
clusters <- cut(as.dendrogram(c.snorm), h=15)
str(clusters)
str(clusters[[2]][2])

failed2 <- clusters[[2]][2]
failed2 <- rownames(SUM.amp.norm)[unlist(failed2)]

all(failed2%in%failed)
all(failed%in%failed2)

cbind(per.sample[order(per.sample)],
      names(per.sample[order(per.sample)])%in%failed, 
      names(per.sample[order(per.sample)])%in%failed2)

#################### Now back to the FULL DATA: OTU counts in samples ###################
d.snorm <- dist(D.mat.snorm[!rownames(D.mat.snorm)%in%failed2, ])
c.snorm <- hclust(d.snorm)
plot(c.snorm)

c.snorm <- hclust(log10(d.snorm), method="average")
plot(c.snorm)

## pdf("simple_nomralization_cluster.pdf", width=14, height=14)
plot(c.snorm)
## dev.off()

## Let's look how replicates behave #####################################################
OTUs.snorm <- as.data.frame(t(D.mat.snorm)) 

reshape.two.cols <- function (x){
    D <- reshape(x, direction="long",
                 varying=1:ncol(x),
                 times=as.character(colnames(x)),
                 timevar="library",
                 v.names=c("OTUc"),
                 ids = row.names(x))

    rownames(D) <- NULL
    D$sample <- gsub("_S\\d+", "", D$library)

    lib.per.sample <- unique(D[,c("library", "sample")])

    lib.per.sample <- do.call(rbind,
                              by(lib.per.sample, lib.per.sample$sample,
                                 function (x) cbind(x, rep=1:nrow(x))))

    R <- merge(D, lib.per.sample)

    ## select only replicate 1 and 2 for replicated samples
    R <- R[R$rep<3, ]
    DR <- reshape(R, direction="wide",
                  v.names = c("OTUc", "library"), timevar= "rep",
                  idvar = c("id", "sample"))

    ## remove non-repeated samples
    DR <- DR[!is.na(DR$OTUc.1) & !is.na(DR$OTUc.2), ]

    DR$F1 <- DR$library.1%in%failed | DR$library.2%in%failed
    DR$F2 <- DR$library.1%in%failed2 | DR$library.2%in%failed2
    return(DR)
}


rep.OTUs <- reshape.two.cols(OTUs)

rep.OTUs.snorm <- reshape.two.cols(OTUs.snorm)

ggplot(rep.OTUs, aes(OTUc.1, OTUc.2, color=F1, shape=F2)) + geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~sample)

## adding regression equation and R2 using code from github
## http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph
library(devtools)
source_gist("524eade46135f6348140")


pdf("replicates_raw.pdf", width=30, height = 30)
ggplot(rep.OTUs, aes(log10(OTUc.1+1), log10(OTUc.2+1), color=F1, shape=F2)) +
    geom_point() +
    stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    geom_smooth(method="lm",se=FALSE) +
    facet_wrap(~sample, scale="free")
dev.off()

pdf("replicates_snorm.pdf", width=30, height = 30)
ggplot(rep.OTUs.snorm, aes(log10(OTUc.1+1), log10(OTUc.2+1),
                           color=F1, shape=F2)) +
    geom_point() +
    stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    geom_smooth(method="lm",se=FALSE) +
    facet_wrap(~sample, scale="free")
dev.off()

################ YOUR TASK: FIND A BETTER NORMALIZATION OR PERFORME RAREFRACTION ######################
## You can try: a) normalize to a lower number and use floor() to make
## all low counts zero the R package phyloseq has a normalization
## routine b) this is based on the normalization in the package
## DESeq(2) which you can also try to use directly alternatively c)
## you can do performe rarefraction (there should be R packages, or
## you can try non-R software with exported data)

#######################################################################################
## read in the BLAST taxonomy assignment
BLASTs <- read.csv("ALL_blast_results_119.csv", sep = ";")

## do we have a common row of ids to merge 
length(intersect(rownames(OTUs), BLASTs$qu))

## yes but in the BLAST some are double
length(unique(BLASTs$qu))
length(BLASTs$qu)

## The Blast similarity report, however, contains duplicated entries.
## Part of the Problem is completely duplicated lines (bug in an
## upstream script?). Another part is multiple HSPs for a blast hit.
head(BLASTs[duplicated(BLASTs$qu), ])

## because blast results are sorted by evalue (lower to higher) we can
## in a very crude effort just take the first occurence of a query
BLASTs <- BLASTs[!duplicated(BLASTs$qu),]
head(BLASTs)

OTUs$amplicon <- amplicon

## now merge the two talbles
DATA <- merge(OTUs, BLASTs, by.x = 0, by.y = "qu") ## all=TRUE)

table(DATA$amplicon)

## How many OTUs are found in only one sample
table(rowSums(DATA[,c(2:104)]>0)>1)

################ Summation by taxonomic assignment #########################
TAX.sum <- by(DATA, DATA$taxid, function (x) {
    colSums(x[,2:105])
})

TAX.sum <- do.call(rbind,  TAX.sum)

TAX.sum.mat <- t(TAX.sum)

d <- dist(TAX.sum.mat[!rownames(TAX.sum.mat)%in%failed2, ])
c <- hclust(d)
plot(c)

## a simple normalization by multiplying with the fraction of the top
## count per sample
TAX.sum.mat.snorm <- TAX.sum.mat*norm.factor 


d <- dist(TAX.sum.mat.snorm[!rownames(TAX.sum.mat.snorm)%in%failed2, ])
c <- hclust(d)

plot(c)

TAX.sum.snorm <- t(TAX.sum.mat.snorm)

TAX.sum.snorm.rep <- reshape.two.cols(as.data.frame(TAX.sum.snorm))

pdf("replicates_tax_sum_snorm.pdf", width=30, height = 30)
ggplot(TAX.sum.snorm.rep, aes(log10(OTUc.1+1), log10(OTUc.2+1),
                           color=F1, shape=F2)) +
    geom_point() +
    stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    geom_smooth(method="lm",se=FALSE) +
    facet_wrap(~sample, scale="free")
dev.off()


### Using a crude way of resolving what should be resolved by a more
### proper normalization #############################

table(TAX.sum.snorm.rep$OTUc.1>0)
table(TAX.sum.snorm.rep$OTUc.2>0)
table(TAX.sum.snorm.rep$OTUc.1>0 & TAX.sum.snorm.rep$OTUc.2>0)

TAX.sum.snorm.rep.trunc <- subset(TAX.sum.snorm.rep,
                                  OTUc.1 >0 & OTUc.2 >0)

pdf("replicates_tax_sum_snorm_trunc.pdf", width=30, height = 30)
ggplot(TAX.sum.snorm.rep.trunc, aes(log10(OTUc.1), log10(OTUc.2),
                           color=F1, shape=F2)) +
    geom_point() +
    stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
    geom_smooth(method="lm",se=FALSE) +
    facet_wrap(~sample, scale="free")
dev.off()

## add species name instead of taxid
taxid.species <- unique(DATA[,c("taxid", "species")])

TAX.sum.snorm <- merge(TAX.sum.snorm, taxid.species, by.x = 0, by.y = "taxid")

TAX.sum.snorm$Row.names <- NULL

rownames(TAX.sum.snorm) <- make.unique(as.character(TAX.sum.snorm$species))
TAX.sum.snorm$species <- NULL
