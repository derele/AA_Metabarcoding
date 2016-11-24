library(phyloseq); packageVersion("phyloseq")

source("R/Hyena_dada2_bacteria.R")
samples.out <- rownames(seqtab.nochim)

subject.ids <- substring(samples.out, 1, 4)

subject.df <- as.data.frame(cbind(Subject=samples.out, Subject.ids=subject.ids))

subject.df <- subject.df[subject.df$Subject.ids%in%Hyena.Cat$ID.Hyena, ]

subject.df <- merge(subject.df, Hyena.Cat,
                    by.x="Subject.ids", by.y="ID.Hyena")

subject.df <- subject.df[!duplicated(subject.df$Subject), ]
subject.df <- subject.df[!is.na(subject.df$Subject), ]
rownames(subject.df) <- subject.df$Subject

seqtab.phylo <- seqtab.nochim[!duplicated(rownames(seqtab.nochim)), ]
seqtab.phylo <- seqtab.phylo[rownames(seqtab.phylo)%in%subject.df$Subject, ]

## order the two identically
subject.df <- subject.df[rownames(seqtab.phylo), ]
subject.df$rep <- ifelse(duplicated(subject.df$Subject.id), "r2", "r1")

subject.df$V3 <- NULL

## Construct phyloseq object (straightforward from dada2 outputs)

ps <- phyloseq(otu_table(seqtab.phylo, taxa_are_rows=FALSE),
               sample_data(subject.df), 
               tax_table(taxa))


pdf("figures/Bacteria_richnesV1.pdf", width=14, height=8)
plot_richness(ps, x="V1") + theme_bw() # + geom_jitter()
dev.off()

pdf("figures/Bacteria_richnesV2.pdf", width=14, height=8)
plot_richness(ps, x="V2") + theme_bw()#  + geom_jitter()
dev.off()

## pdf("figures/Bacteria_richnesV3.pdf", width=14, height=8)
## plot_richness(ps, x="V3") + theme_bw() + geom_jitter()
## dev.off()

pdf("figures/Bacteria_richnesV4.pdf", width=14, height=8)
plot_richness(ps, x="V4") + theme_bw()# + geom_jitter()
dev.off()

## ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")

## plot_ordination(ps, ord.nmds.bray, color="V4", title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

## plot_bar(ps.top20, x="rep", fill="Family") + facet_wrap(~V2, scales="free_x")

## plot_bar(ps.top20, x="V1", fill="Family") + facet_wrap(~V2, scales="free_x")

## plot_bar(ps.top20, x="V2", fill="Family") + facet_wrap(~V2, scales="free_x")

enterotype <- subset_taxa(ps, Genus != "-1")

dist_methods <- unlist(distanceMethodList)

## 1:3 require tree
dist_methods <- dist_methods[-(1:3)]

dist_methods <- dist_methods[-which(dist_methods=="ANY")]

plist <- vector("list", length(dist_methods))

names(plist) <- dist_methods


for( i in dist_methods[c(1:3, 8, 13, 19)] ){ # other give error
    ## Calculate distance matrix
    iDist <- distance(enterotype, method=i)
    ## Calculate ordination
    iMDS  <- ordinate(enterotype, "MDS", distance=iDist)
    ## Make plot
    ## Don't carry over previous plot (if error, p will be blank)
    p <- NULL
    ## Create plot, store as temp variable, p
    p <- plot_ordination(enterotype, iMDS, color="SeqTech", shape="Enterotype")
    ## Add title to each plot
    p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
                                        # Save the graphic to file.
    plist[[i]] <- p
}



library(plyr)

df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

p <- ggplot(df, aes(Axis.1, Axis.2, color=V4, shape=V1)) + theme_bw()
p <- p + geom_point(size=3, alpha=0.5)
p <- p + facet_wrap(~distance, scales="free")
p <- p + ggtitle("MDS on various distance metrics for bacterial metabarcoding dataset")
p


