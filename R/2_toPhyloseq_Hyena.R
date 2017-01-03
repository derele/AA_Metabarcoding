
load(file="/SAN/Metabarcoding/exclude.Rdata") ## -> Samples.to.exclude
## load(MA6, file="/SAN/Metabarcoding/allMA.Rdata") ## -> MA6 just for
## control 
load(file="/SAN/Metabarcoding/table.Rdata") ## -> STnoC
load(file="/SAN/Metabarcoding/taxa.Rdata") ## tax.l
## load(file="/SAN/Metabarcoding/trees.Rdata") ## tree.l

######### FROM HERE completely Hyena specific ###################
Hyena.Cat <- read.csv("/home/ele/Dropbox/Animal_number_variable_codes_EH.csv")
Hyena.Cat <- as.data.frame(apply(Hyena.Cat, 2, function (x) gsub(" ", "", x)))

Hyena.Cat$V1 <- ifelse(Hyena.Cat$V1==1, "Male", "Female")
Hyena.Cat$V2 <- ifelse(Hyena.Cat$V2==2, "Cub", "Adult")
Hyena.Cat$V3 <- ifelse(Hyena.Cat$V3==1, "high", "low")
## Hyena.Cat$V3[is.na(Hyena.Cat$V3)] <- "low" ## watchout not TRUE!!!
Hyena.Cat <- Hyena.Cat[!is.na(Hyena.Cat$V3), ]

names(Hyena.Cat) <- c("animal", "ID", "sex", "age", "rank", "pack")
rownames(Hyena.Cat) <- Hyena.Cat$ID
Hyena.Cat$ID <- Hyena.Cat$animal <- NULL

Exp.Cat <- read.table("/SAN/Metabarcoding/AA_combi/sample_table.csv", header=FALSE)

get.sample.data <- function(seqtab.nochim){
    samples.out <- rownames(seqtab.nochim)
    samples.out <- samples.out[!samples.out%in%Samples.to.exclude]
    subject.df <- Exp.Cat[Exp.Cat$V1%in%rownames(Hyena.Cat), ]
    subject.df <- merge(subject.df, Hyena.Cat,
                        by.x="V1", by.y=0)
    rownames(subject.df) <- subject.df$V2
    seqtab.phylo <- seqtab.nochim[rownames(seqtab.nochim)%in%rownames(subject.df), ]
    ## order the two identically
    subject.df <- subject.df[rownames(seqtab.phylo), ]
    subject.df$rep <- ifelse(duplicated(subject.df$V1), "r2", "r1")
    list(subject.df, seqtab.phylo)
}

sub.df.l <- lapply(STnoC, get.sample.data)

## Construct phyloseq object (straightforward from dada2 outputs)

ps.l <- lapply(seq_along(sub.df.l), function (i) {
    ps <- phyloseq(otu_table(sub.df.l[[i]][[2]], taxa_are_rows=FALSE),
                   sample_data(sub.df.l[[i]][[1]]),
                   ##                   phy_tree(tree.l[[i]]$tree),
                   tax_table(tax.l[[i]]))
    return(ps)
})

names(ps.l) <- names(STnoC)
