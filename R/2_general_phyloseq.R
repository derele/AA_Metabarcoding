library(phyloseq)
library(reshape)

load(file="/SAN/Metabarcoding/table.Rdata") ## -> STnoC
load(file="/SAN/Metabarcoding/taxa.Rdata") ## tax.l

ps.l <- lapply(seq_along(STnoC), function (i) {
    ps <- phyloseq(otu_table(STnoC[[i]], taxa_are_rows=FALSE),
                   ## phy_tree(tree.l[[i]]$tree),
                   tax_table(tax.l[[i]]))
    return(ps)
})

names(ps.l) <- names(STnoC)

ps.l <- ps.l[!grepl("Klin|ADM", names(ps.l))]


num.taxa <-sapply(c("Genus", "Family", "Order", "Class", "Phylum"), function (rank){
    lapply(ps.l, function (x) 
        length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads <- unlist(lapply(ps.l, function (x) 
    sum(otu_table(x))))

PrimTax <- as.data.frame(cbind(num.taxa, num.reads))

PrimTax <- as.data.frame(apply(PrimTax, 2, unlist))

ggplot(PrimTax, aes(Genus, num.reads)) +
    geom_point() #+
##    scale_y_log10()

ggplot(PrimTax, aes(Order, num.reads)) +
    geom_point() +
    scale_y_log10()

ggplot(PrimTax, aes(Phylum, num.reads)) +
    geom_point()


unique.taxa.cumsum <- function(ps.list, rank){
    taxonNa <- lapply(ps.list, function (x) 
        get_taxa_unique(x, taxonomic.rank = rank))
    names(taxonNa) <- names(ps.list)
    ## taxonNa <- taxonNa[order(unlist(lapply(taxonNa, length)),
    ##                         decreasing=TRUE)]
    unique.taxonNa.cumsum <- sapply(0:length(taxonNa),
                                    function (i) length(unique(unlist(taxonNa[0:i]))))
    return(unique.taxonNa.cumsum)
}

ps.numord.l <- ps.l[order(PrimTax$num.reads, decreasing=TRUE)]


cum.tax <- sapply(c("Genus", "Family", "Class", "Order", "Phylum"), function (rank){
    unique.taxa.cumsum(ps.numord.l, rank)
})

colnames(cum.tax) <- paste0("cum.", colnames(cum.tax))

cum.tax <- as.data.frame(cbind(1:nrow(cum.tax), cum.tax))


cum.plot <- ggplot(cum.tax) +
    geom_step(aes(V1, cum.Genus), color="red") +
    geom_text(aes(20, 1100, label="Genera"), color="red") +
    geom_step(aes(V1, cum.Family), color="purple") +
    geom_text(aes(20, 610, label="Families"), color="purple") +
    geom_step(aes(V1, cum.Order), color="blue") +
    geom_text(aes(20, 310, label="Orders"), color="blue") +
    geom_step(aes(V1, cum.Class), color="seagreen") +
    geom_text(aes(20, 135, label="Classes"), color="seagreen") +
    geom_step(aes(V1, cum.Phylum), color="green") +
    geom_text(aes(20, 63, label="Phyla"), color="green") +
    scale_y_log10("Cummulative count of taxa") +
    scale_x_continuous("Number of primers considerd (starting with the one with highest read count)") + 
    annotation_logticks(sides="l") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

pdf("figures/Cummulative_Taxa.pdf")
cum.plot
dev.off()

before33 <- lapply(ps.numord.l[1:32], function (x) {
    get_taxa_unique(x, "Phylum")
})

before33 <- unique(unlist(before33))
at33 <- get_taxa_unique(ps.numord.l[[33]], "Phylum")
at33[!at33%in%before33]

before28 <- lapply(ps.numord.l[1:28], function (x) {
    get_taxa_unique(x, "Phylum")
})

before28 <- unique(unlist(before28))
at28 <- get_taxa_unique(ps.numord.l[[28]], "Phylum")
at28[!at28%in%before28]


### putting it all together (throwing away the primer data)
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


all.tax.counts <- Reduce(cbind,
                         STnoC.filled[!grepl("Klin|ADM", names(STnoC.filled))])
all.tax <- Reduce(rbind, tax.l[!grepl("Klin|ADM", names(STnoC.filled))])

PSA <- phyloseq(otu_table(all.tax.counts, taxa_are_rows=FALSE),
                tax_table(all.tax))

## Just to find out what phyla of low throughput primers are 
ordered.throughput <- num.reads[order(num.reads, decreasing=TRUE)]

cum.throu <- cumsum(ordered.throughput)

one.sample <- colSums(otu_table(PSA))

one.sample <- matrix(one.sample, ncol=length(one.sample))

colnames(one.sample) <- colnames(otu_table(PSA))

OneSample <- phyloseq(otu_table(one.sample, taxa_are_rows=FALSE),
                      tax_table(PSA))

set.seed(124)
rare.at.seq <- sapply(cum.throu, function (x) {
    rarefy_even_depth(OneSample, sample.size = x, 
    verbose = FALSE, replace = TRUE)
})


cum.tax$rare.Genus <- c(0, length(unique(get_taxa_unique(PSA.rare, "Genus"))))

rare.tax <- sapply(c("Genus", "Family", "Class", "Order", "Phylum"), function (rank){
    sapply(rare.at.seq, function (x) {
        length(unique(get_taxa_unique(x, rank)))
    })
})

rare.tax <- rbind(c(0, 0, 0, 0, 0), rare.tax)

rownames(rare.tax) <- NULL

colnames(rare.tax) <- paste0("rare.", colnames(rare.tax))

Tax.eval <- cbind(cum.tax, rare.tax)

Eval.plot <- cum.plot + 
    geom_step(data=Tax.eval, aes(V1, rare.Genus),
              color="red", linetype=2) + 
    geom_step(data=Tax.eval, aes(V1, rare.Family),
              color="purple", linetype=2) +
    geom_step(data=Tax.eval, aes(V1, rare.Order),
              color="blue", linetype=2) +
    geom_step(data=Tax.eval, aes(V1, rare.Class),
              color="seagreen", linetype=2) +
    geom_step(data=Tax.eval, aes(V1, rare.Phylum),
              color="green", linetype=2)

pdf("figures/Eval_Taxa.pdf")
Eval.plot
dev.off()
