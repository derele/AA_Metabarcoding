library(phyloseq)

load(file="/SAN/Metabarcoding/exclude.Rdata") ## -> Samples.to.exclude
## load(MA6, file="/SAN/Metabarcoding/allMA.Rdata") ## -> MA6 just for
## control 
load(file="/SAN/Metabarcoding/table.Rdata") ## -> STnoC
load(file="/SAN/Metabarcoding/taxa.Rdata") ## tax.l
load(file="/SAN/Metabarcoding/trees.Rdata") ## tree.l

######### FROM HERE completely Hyena specific ###################
Hyena.Cat <- read.csv("/home/ele/Dropbox/Hyena_Hartmann_MS/Animal_number_variable_codes_EH.csv")
Hyena.Cat <- as.data.frame(apply(Hyena.Cat, 2, function (x) gsub(" ", "", x)))

Hyena.Cat$V1 <- ifelse(Hyena.Cat$V1==1, "Male", "Female")
Hyena.Cat$V2 <- ifelse(Hyena.Cat$V2==2, "Cub", "Adult")
Hyena.Cat$V3 <- ifelse(Hyena.Cat$V3==1, "high", "low")
## Hyena.Cat$V3[is.na(Hyena.Cat$V3)] <- "low" ## watchout not TRUE!!!
Hyena.Cat <- Hyena.Cat[!is.na(Hyena.Cat$V3), ]

names(Hyena.Cat) <- c("animal", "Hyena.ID", "sex", "age", "rank", "pack")
Hyena.Cat$ID <- Hyena.Cat$animal <- NULL

Exp.Cat <- read.table("/SAN/Metabarcoding/AA_combi/sample_table.csv", header=FALSE)
names(Exp.Cat) <- c("Hyena.ID", "Sample.ID")

## might be needed to add this to samples to exclude... "I436"

subject.df <- merge(Exp.Cat, Hyena.Cat, by="Hyena.ID")
subject.df <- subject.df[!subject.df$Sample.ID%in%Samples.to.exclude, ]
subject.df$rep <- ifelse(duplicated(subject.df$Sample.ID), "r2", "r1")
rownames(subject.df) <- subject.df$Sample.ID

all.samples <- unique(unlist(lapply(STnoC, rownames)))

STnoC.filled <- lapply(STnoC, function (foo){
    missing.samples <- all.samples[!all.samples%in%rownames(foo)]
    if(length(missing.samples)>0){
        bar <- matrix(0, nrow=length(missing.samples), ncol=ncol(foo))
        rownames(bar) <- missing.samples
        foobar <- rbind(foo, bar)
    } else {foobar <- foo}
    foobar[all.samples, ]
})

Sync.STnoC.Hy <- lapply(STnoC.filled, function (x) {
    x[rownames(subject.df), ]
})


## Construct phyloseq object (straightforward from dada2 outputs)
ps.l <- lapply(seq_along(Sync.STnoC.Hy), function (i) {
    ps <- phyloseq(otu_table(Sync.STnoC.Hy[[i]], taxa_are_rows=FALSE),
                   sample_data(subject.df),
    ## phy_tree(tree.l[[i]]$tree),
                   tax_table(tax.l[[i]]))
    return(ps)
})

save(ps.l, file="/SAN/Metabarcoding/phlyoSeq_list.Rdata")

all.tax.counts <- Reduce(cbind, Sync.STnoC.Hy)
all.tax <- Reduce(rbind, tax.l)

PSA <- phyloseq(otu_table(all.tax.counts, taxa_are_rows=FALSE),
                sample_data(subject.df),
                tax_table(all.tax))

## correct number of taxa
sum(unlist(lapply(ps.l, function (x) nrow(tax_table(x)))))

save(PSA, file="/SAN/Metabarcoding/phlyoSeq_cat.Rdata")

## not happy with phyloseq merge_samples so coding this manually
sum.columns <- function(ps){
    SD <- sample_data(ps)$Hyena.ID
    OT <- otu_table(ps)
    do.call(rbind, by(OT, SD, colSums))
}

## mean between replicates
SumRSVs <- sum.columns(PSA)

nrow()

sum.subject.df <- subject.df[!duplicated(subject.df$Hyena.ID), ]
sum.subject.df$Sample.ID <- NULL
sum.subject.df$rep <- NULL
rownames(sum.subject.df) <- sum.subject.df$Hyena.ID

PS <- phyloseq(otu_table(SumRSVs, taxa_are_rows=FALSE),
               sample_data(sum.subject.df),
               tax_table(all.tax))

save(PS, file="/SAN/Metabarcoding/phlyoSeq_Hy.Rdata")
