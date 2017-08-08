library(grid)
library(gridExtra)

############ Log count based weighted unifrac distances ################
ps.l.log <- lapply(ps.l.clean, function (PS) {
    transform_sample_counts(PS, function(x) log(1 + x))
})

ps.l.log <- lapply(ps.l.log, function (PS) {
    p <- subset_taxa(PS, colSums(otu_table(PS))>0)
    subset_samples(p, rowSums(otu_table(p))>0)
})

ord.l.log <- mclapply(ps.l.log, function (x){
    ordinate(x, method = "MDS", distance = "wunifrac")
}, mc.cores=20)

ord.l.log <- ord.l.log[unlist(lapply(ord.l.log, length))==5]

evals.l.log <- lapply(ord.l.log, function (x) x[["values"]][["Eigenvalues"]])

ord.plot.l <- lapply(1:length(ord.l.log), function (i){
    plot_ordination(ps.l.log[[i]], ord.l.log[[i]], color = "species", label = "V1") +
        labs(col = "species", title=names(ord.l.log)[[i]]) +
        coord_fixed(sqrt(evals.l.log[[i]][2] / evals.l.log[[i]][1])) 
})

pdf("figures/unifac_log_MDS.pdf", width=46, height=46)
do.call(grid.arrange, ord.plot.l)
dev.off()

### bray
ord.l.log <- mclapply(ps.l.log, function (x){
    ordinate(x, method = "MDS", distance = "bray")
}, mc.cores=20)

ord.l.log <- ord.l.log[unlist(lapply(ord.l.log, length))==5]

evals.l.log <- lapply(ord.l.log, function (x) x[["values"]][["Eigenvalues"]])

ord.plot.l <- lapply(1:length(ord.l.log), function (i){
    plot_ordination(ps.l.log[[i]], ord.l.log[[i]], color = "species", label = "V1") +
        labs(col = "species", title=names(ord.l.log)[[i]]) +
        coord_fixed(sqrt(evals.l.log[[i]][2] / evals.l.log[[i]][1])) 
})

pdf("figures/bray_log_MDS.pdf", width=46, height=46)
do.call(grid.arrange, ord.plot.l)
dev.off()

############ Rank based weighted unifrac distances ################

ps.l.rank <- lapply(ps.l.clean, function (PS) {
    transform_sample_counts(PS, function(x) rank(x))
})

ps.l.rank <- lapply(ps.l.rank, function (PS) {
    p <- subset_taxa(PS, colSums(otu_table(PS))>0)
    subset_samples(p, rowSums(otu_table(p))>0)
})

ord.l.rank <- mclapply(ps.l.rank, function (x){
    ordinate(x, method = "MDS", distance = "wunifrac")
}, mc.cores=20)

ord.l.rank <- ord.l.rank[unlist(lapply(ord.l.rank, length))==5]

evals.l.rank <- lapply(ord.l.rank, function (x) x[["values"]][["Eigenvalues"]])

ord.plot.l <- lapply(1:length(ord.l.rank), function (i){
    plot_ordination(ps.l.rank[[i]], ord.l.rank[[i]], color = "species", label = "V1") +
        labs(col = "species", title=names(ord.l.rank)[[i]]) +
        coord_fixed(sqrt(evals.l.rank[[i]][2] / evals.l.rank[[i]][1])) 
})

pdf("figures/unifac_rank_MDS.pdf", width=46, height=46)
do.call(grid.arrange, ord.plot.l)
dev.off()


########### Combined data ############################

PSA.log <- transform_sample_counts(PSA.clean, function(x) log(1 + x))

PSA.log <- subset_taxa(PSA.log, colSums(otu_table(PSA.log))>0)
PSA.log <- subset_samples(PSA.log, rowSums(otu_table(PSA.log))>0)

ord.log <- ordinate(PSA.log, method = "MDS", distance = "bray")

evals.log <- ord.log[["values"]][["Eigenvalues"]]

ord.plot.log <- 
    plot_ordination(PSA.log, ord.log, color = "species", label = "V1") +
        labs(col = "species", title="combined from all primers") +
        coord_fixed(sqrt(evals.log[2] / evals.log[1]))

pdf("figures/bray_log_MDS_combi.pdf", width=8, height=8)
print(ord.plot.log)
dev.off()

### bray
PSA.rank <- transform_sample_counts(PSA.clean, function(x) log(1 + x))

PSA.rank <- subset_taxa(PSA.rank, colSums(otu_table(PSA.rank))>0)
PSA.rank <- subset_samples(PSA.rank, rowSums(otu_table(PSA.rank))>0)

ord.rank <- ordinate(PSA.rank, method = "MDS", distance = "bray")

evals.rank <- ord.rank[["values"]][["Eigenvalues"]]

ord.plot.rank <- 
    plot_ordination(PSA.rank, ord.rank, color = "species", label = "V1") +
        labs(col = "species", title="combined from all primers") +
        coord_fixed(sqrt(evals.rank[2] / evals.rank[1]))

pdf("figures/bray_rank_MDS_combi.pdf", width=8, height=8)
print(ord.plot.rank)
dev.off()


########## Some diversity estimates ## not very useful
rich.plot <- lapply(ps.l.clean, function (ps){
     plot_richness(ps, x="species", measures=c("Observed", "Chao1"))
})

pdf("figures/richness_all.pdf", width=46, height=46)
do.call(grid.arrange, rich.plot)
dev.off()



