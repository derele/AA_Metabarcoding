library(nlme)
library(ggplot2)
library(exactRankTests)

################# From here Phyloseq ######################
## load(file="/SAN/Metabarcoding/phlyoSeq_list.Rdata") ## -> ps.l

## load(file="/SAN/Metabarcoding/phlyoSeq_cat.Rdata") ## -> PSA

load(file="/SAN/Metabarcoding/phlyoSeq_Hy.Rdata") ## -> PS

PS.genus <- tax_glom(PS, "Genus", NArm = TRUE)

PS.genus.na <- tax_glom(PS, "Genus", NArm = FALSE)

Rich_plot <- plot_richness(PS.genus, x="rank") + theme_bw()  
Rich_plot$layers <- NULL
Rich_plot + geom_boxplot()

## What phyla are detected

table(tax_table(PS)[, "Phylum"])

euk_phyla <- names(table(tax_table(subset_taxa(PS, Kingdom %in% "Eukaryota"))[, "Phylum"]))

undef <- c("", "uncultured","undef")
clear.prey.phyla <- c("Vertebrata", "Chordata")
plant.phyla <- c("Chlorophyta", "Chlorophyta_ph", "Streptophyta",
                 "Phragmoplastophyta", "Ochrophyta", "Klebsormidiophyceae",
                 "Eustigmatophyceae", "Bacillariophyta")
clear.parasite.phyla <- c("Nematoda", "Apicomplexa",
                          "Platyhelminthes", "Microsporidia")

## how many taxa in clear parasite phyla 
subset_taxa(PS, Phylum %in% clear.parasite.phyla)
subset_taxa(PS.genus, Phylum %in% clear.parasite.phyla)

## how many taxa in clear prey phyla 
subset_taxa(PS, Phylum %in% clear.prey.phyla)
subset_taxa(PS.genus, Phylum %in% clear.prey.phyla)

## most very unlikely to be correctly annotated at the genus level
unname(tax_table(subset_taxa(PS.genus, Phylum %in% clear.prey.phyla))[, "Genus"])

test.Chao <- function(ps){
    df <- merge(sample_data(ps),
                estimate_richness(ps, measures="Chao1"),
                by=0)
    test <- wilcox.exact(Chao1 ~ rank, data = df, exact = TRUE)
    test
}

test.Chao(PS)
test.Chao(PS.genus)
test.Chao(PS.genus.na)

test.Chao(subset_taxa(PS,
                         !Phylum %in% c(clear.prey.phyla, undef) &
                         Kingdom %in% "Eukaryota"
                         ))

test.Chao(subset_taxa(PS.genus.na,
                         !Phylum %in% c(clear.prey.phyla, undef) &
                         Kingdom %in% "Eukaryota"
                         ))

test.Chao(subset_taxa(PS.genus,
                         !Phylum %in% c(clear.prey.phyla, undef) &
                         Kingdom %in% "Eukaryota"
                         ))


test.Chao(subset_taxa(PS, Phylum %in% clear.parasite.phyla))

test.Chao(subset_taxa(PS.genus.na, Phylum %in% clear.parasite.phyla))

test.Chao(subset_taxa(PS.genus,Phylum %in% clear.parasite.phyla))


test.Chao(subset_taxa(PS, Phylum %in% clear.prey.phyla))

test.Chao(subset_taxa(PS, Phylum %in% plant.phyla))


subset_taxaE<- function (physeq, subset, ...) {
    if (is.null(tax_table(physeq))) {
        cat("Nothing subset. No taxonomyTable in physeq.\n")
        return(physeq)
    }
    else {
        oldMA <- as(tax_table(physeq), "matrix")
        oldDF <- data.frame(oldMA)
        newDF <- subset(oldDF, subset, ...)
        newMA <- as(newDF, "matrix")
        if (inherits(physeq, "taxonomyTable")) {
            return(tax_table(newMA))
        }
        else {
            tax_table(physeq) <- tax_table(newMA)
            return(physeq)
        }
    }
}


Chao.phyla <- lapply(euk_phyla, function(x){
    ps <- subset_taxaE(PS.genus, tax_table(PS.genus)[, "Phylum"] %in% x)
    ps
    tryCatch(test.Chao(ps), error = function (e) NULL)
})

names(Chao.phyla) <- euk_phyla

p.Chao <- lapply(Chao.phyla, function (x) {
    x[c("statistic", "p.value")]
})

names(p.Chao) <- euk_phyla

p.Chao <- p.Chao[!unlist(lapply(p.Chao, is.null))]

p.Chao.table <- data.frame(do.call(rbind, p.Chao))

## p.Chao.table$FDR <- p.adjust(p.Chao.table$p.value, method="BH")

p.Chao.table[p.Chao.table$p.value<0.05, ]

## working with worm counts
WC <- read.csv("/home/ele/Dropbox/Hyena_Hartmann_MS/Hyena_WormCounts_fixed.csv")

colnames(WC)[2:ncol(WC)] <- paste0("count_", colnames(WC)[2:ncol(WC)])


## For selected genera

get.ps.counts <- function (ps, level){
    Tax.tab <- otu_table(ps)
    tax.cols <- unname(tax_table(ps)[, level])
    colnames(Tax.tab) <- tax.cols
    All.data <- cbind(sample_data(ps), Tax.tab)
    list(merge(WC, All.data, by.x="ID.Hyena", by.y="Hyena.ID"), tax.cols)
}
    


############# Untargeted #################
get.top.cors <- function(ps, level){
    psc<- get.ps.counts(ps, level)
    foo <- psc[[1]]
    tax.cols <- psc[[2]]
    count.col <- grep("count_", colnames(foo), value=TRUE)
    foo.cor <- data.frame(t(cor(foo[, count.col], foo[, tax.cols],
                                use="pairwise.complete.obs")))
    top.cors <- lapply(count.col, function (x){
        head(foo.cor[order(foo.cor[,x], decreasing=TRUE), x, drop=FALSE], n=20)
    })
    names(top.cors) <- count.col
    return(top.cors)
}

## runs for very long...
## get.top.cors(PSA.genus, "Genus")

get.top.cors(PS.genus, "Genus")

get.top.cors(subset_taxa(PS.genus, Phylum %in%
                                   c("Apicomplexa", "Platyhelminthes",
                                     "Nematoda", "Microsporidia")),
             "Genus")



##count_Mite Hemialges  1.0000000

##                             count_Cystoisospora2
## Thoreauomyces                          0.8296620
## Rurikoplites                           0.8296620
## Coelomomyces                           0.8250951
## Turicibacter                           0.8005079
## Crenosoma                              0.7881375
## Haemonchus                             0.7784763
## Zoniolaimus                            0.7495413
## Isospora                               0.7398843

##                  count_Dipilydium
## Adenocephalus         0.891580404
## Amidostomum           0.750319752
## Cryptosporidium       0.690625689
## Leidyana              0.542524472
## Gregarina             0.372994050
## Dipylidium            0.323629188

##                 count_Spirurida
## Strongyloides        0.73909046

##                  count_U_Mesocestoide
## Caryospora                0.926864282
## Eimeria                   0.925575434
## Eimeria.1                 0.925575434

##                   count_Taenia
## Strongyloides       0.69137862


##                   count_Taenia.45.
## Anoplocephala           1.00000000
## Apharyngostrigea        1.00000000
## Austrodiplostomum       1.00000000
## Neoctenotaenia          1.00000000

##                  count_Cystoisospora2
## Crenosoma                  0.78813753
## Haemonchus                 0.77847627
## Zoniolaimus                0.74954128
## Isospora                   0.73988429

##                count_Cystoisospora_small
## Geneiorhynchus                0.86393186
## Toxoplasma                    0.86392121
## Besnoitia                     0.86290245
## Besnoitia.1                   0.86290245
## Gregarina                     0.71671086


##                     count_Spirometra
## Spirometra               0.679405146
## Diphyllobothrium         0.591959397

##                   count_Ancylostoma2
## Anoplocephala             0.98986803
## Apharyngostrigea          0.98986803
## Austrodiplostomum         0.98986803
## Neoctenotaenia            0.98986803
## Moniezia                  0.98985348
## Paraschneideria           0.98390109
## Cylicocyclus              0.98205917
## Petrovinema               0.90560644

##                   count_Ancylostoma_small
## Teladorsagia                    0.7734495
## Enterocytozoon                  0.7734495
## Ancylostoma                     0.4450958
## Bothridium                      0.2963130



ggplot(Tax.tab, aes(Ancylostoma+1, Ancylostoma_small+Ancylostoma2+1,
                    label=V1)) +
    geom_point() +
    ##    geom_text()+
    scale_x_log10("Sequence Abundance") +
    scale_y_log10("Flotation Counts") +
    ggtitle("Ancylostoma") 

ggplot(Tax.tab, aes(Spirometra.x+1, Spirometra.y+1)) +
    geom_point() +
    scale_x_log10("Sequence Abundance") +
    scale_y_log10("Flotation Counts") +
    ggtitle("Spirometra")

### for all genera



## inspect the reported phyla
lapply(ps.l, function(ps){
    table(tax_table(ps)[, "Phylum"], exclude = NULL)
})

## remove otus without phylum level reported
Phyl.ps.l <- lapply(ps.l, function (ps){
    subset_taxa(ps, !is.na(Phylum) &
                    !Phylum %in% c("", "undef", "uncharacterized"))
})


## Are there phyla that are comprised of mostly low-prevalence
## features? Compute the total and average prevalences of the features
## in each phylum.
prevdf.l <- lapply(Phyl.ps.l, function (ps){
   prev <- apply(X = otu_table(ps),
                 MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
   prevTax <- data.frame(Prevalence = prev, 
                         TotalAbundance = taxa_sums(ps),
                         tax_table(ps))
   return(prevTax)
})


lapply(prevdf.l, function(prevTax){
    allAbu<- sum(prevTax$TotalAbundance)
    plyr::ddply(prevTax, "Phylum", function(df1){
        cbind(meanPrevalence=mean(df1$Prevalence),
              maxPrevalence=max(df1$Prevalence),
              TotAbund=sum(df1$TotalAbundance),
              PercAbund=round((sum(df1$TotalAbundance)/allAbu)*100, 3)
              )
    })
})


## invest a bit more time here to look through those

## exclude everything that apperars in less than prevalenceTreshold
## samples (here 2) samples
prevalenceThreshold <- 2

## 
TSps.l <- lapply(seq_along(prevdf.l), function(i){
    keep <- rownames(prevdf.l[[i]])[(prevdf.l[[i]]$Prevalence >= prevalenceThreshold)]
    prune_taxa(keep, Phyl.ps.l[[i]])
})


Gen.l <- mclapply(TSps.l, function (ps) {
    tax_glom(ps, "Genus", NArm = TRUE)},
    mc.cores=20)

Dist.l <- mclapply(TSps.l, function (ps) {
    tip_glom(ps, h=0.1)},
    mc.cores=20)

get.num.taxa <- function(x, level="Genus"){
    length(get_taxa_unique(x, taxonomic.rank = level))
}

cbind(raw = lapply(ps.l, get.num.taxa),
      wPhylum = lapply(Phyl.ps.l, get.num.taxa),
      wPrev = lapply(TSps.l, get.num.taxa),
      GenAgg = lapply(Gen.l, get.num.taxa),
      DistAgg = lapply(Dist.l, get.num.taxa)
      )

## the latter are so few that we can  have a direct look:
lapply(Gen.l, function(x) get_taxa_unique(x, taxonomic.rank = "Genus"))

lapply(Dist.l, function(x) get_taxa_unique(x, taxonomic.rank = "Genus"))

get.all.shannon <- function(ps){
    Shan <- estimate_richness(ps, measures="Shannon")
    x <- merge(sample_data(ps), Shan, by=0)
    s <- t.test(Shannon ~ sex, data=x)$p.value
    a <- t.test(Shannon ~ age, data=x)$p.value
    r <- t.test(Shannon ~ rank, data=x)$p.value
    n <- t.test(Shannon ~ rep, data=x)$p.value
    Shannon.df <- cbind(sex=s, age=a,
                        rank=r, rep=n)
    return(Shannon.df)
}


do.call(rbind, lapply(ps.l, get.all.shannon))
do.call(rbind, lapply(Phyl.ps.l, get.all.shannon))
do.call(rbind, lapply(TSps.l, get.all.shannon))
do.call(rbind, lapply(Gen.l, get.all.shannon))
do.call(rbind, lapply(Dist.l, get.all.shannon))


Shannon.df <- do.call(rbind, lapply(Phyl.ps.l, get.all.shannon))
rownames(Shannon.df) <- names(ps.l)

annotation <- merge(primer.overview, Shannon.df, by=0)
rownames(annotation) <- annotation$Row.names
annotation$Row.names <- NULL
annotation$is.16.S <- ifelse(grepl("ADM|ACM|Klin", annotation$pRname), "16S", "18S")

pdf("figures/shannon.pdf", onefile=FALSE)
pheatmap(Shannon.df, cluster_cols=FALSE, cluster_rows=TRUE,
         color = colorRampPalette(c("firebrick3", "navy"))(4),
         breaks = 10^(0:-4),
         show_rownames=FALSE,
         legend = FALSE,
         display_numbers=TRUE,
         number_format = "%.1e",
         annotation_row = annotation[, c(6, 7, 13)])
dev.off()


## one example of one bacterial amplicon
ps.bak1 <- ps.l[[3]]

pdf("figures/Bacteria_richnesSex.pdf", width=14, height=8)
plot_richness(ps.bak1, x="sex") + theme_bw()  + geom_boxplot()
dev.off()

pdf("figures/Bacteria_richnesAge.pdf", width=14, height=8)
plot_richness(ps.bak1, x="age") + theme_bw()  + geom_boxplot()
dev.off()

pdf("figures/Bacteria_richnesRank.pdf", width=14, height=8)
plot_richness(ps.bak1, x="rank") + theme_bw()  + geom_boxplot()
dev.off()

pdf("figures/Bacteria_richnesPack.pdf", width=14, height=8)
plot_richness(ps.bak1, x="pack") + theme_bw()  + geom_boxplot()
dev.off()


all.div.plots <- function(ps, what, annotation){
 pl <- plot_richness(ps, x=what, measures="Shannon") +
     theme_bw()  +
     geom_violin() +
     ggtitle(annotation$pRnames)
 if(annotation$is.16.S%in%"18S"){
     p <- pl + theme(strip.background = element_rect(fill="orange"))
 } else{
     p <- pl + theme(strip.background = element_rect(fill="blue"))
 }
 if(as.numeric(as.character(annotation[, what]))<0.01){
     p <- p + theme(panel.background = element_rect(fill = "palevioletred2"))
 } else {
     p
 }
}

pdf("figures/richness_Rank.pdf", width=60, height=20)
do.call("grid.arrange",
        c(lapply(seq_along(ps.l), function (i) {
            all.div.plots(ps.l[[i]], "rank", annotation[i,])
        }),
        ncol=10))
dev.off()


pdf("figures/richness_Age.pdf", width=60, height=20)
do.call("grid.arrange",
        c(lapply(seq_along(ps.l), function (i) {
            all.div.plots(ps.l[[i]], "age", annotation[i,])
        }),
        ncol=10))
dev.off()



get.female.adult.subset <- function(ps){
    dat <- sample_data(ps)
    keep <- dat$sex%in%"Female" & dat$age%in%"Adult"
    subset_samples(ps, keep)
}

FAps.l <- lapply(ps.l, get.female.adult.subset)
FAPhyl.ps.l <- lapply(Phyl.ps.l, get.female.adult.subset)
FATSps.l <- lapply(TSps.l, get.female.adult.subset)
FAGen.l <- lapply(Gen.l, get.female.adult.subset)
FADist.l <- lapply(Dist.l, get.female.adult.subset)


pdf("figures/FArichness_Rank.pdf", width=60, height=20)
do.call("grid.arrange",
        c(lapply(seq_along(FAPhyl.ps.l), function (i) {
            all.div.plots(FAPhyl.ps.l[[i]], "rank", annotation[i,])
        }),
        ncol=10))
dev.off()


get.shannon.lme <- function(ps){
    Shan <- estimate_richness(ps, measures="Shannon")
    x <- merge(sample_data(ps), Shan, by=0)
    test <- lme(fixed = Shannon ~ rank + pack, data = x,
                random = ~ 1 | Subject.ids)
    s.test <- summary(test)
    r.effect <- s.test$tTable["ranklow", "Value"]
    r.pval <- s.test$tTable["ranklow", "p-value"]
    a.test <- anova(test)
    c.pval <- a.test$"p-value"[[3]]
    Shannon.df <- cbind(r.low.effect=r.effect, r.pval=r.pval,
                        c.pval=c.pval)
    return(Shannon.df)
}


do.call(rbind, lapply(FAps.l, get.shannon.lme))
do.call(rbind, lapply(FAPhyl.ps.l, get.shannon.lme))
do.call(rbind, lapply(FATSps.l, get.shannon.lme))
do.call(rbind, lapply(FAGen.l, get.shannon.lme))
do.call(rbind, lapply(FADist.l, get.shannon.lme))

Shannon.df <- do.call(rbind, lapply(FAps.l, get.shannon.lme))
rownames(Shannon.df) <- names(ps.l)

annotation <- merge(primer.overview, Shannon.df, by=0)
rownames(annotation) <- annotation$Row.names
annotation$Row.names <- NULL
annotation$is.16.S <- ifelse(grepl("ADM|ACM|Klin", annotation$pRname), "16S", "18S")

lme.div.plots <- function(ps, annotation){
 pl <- plot_richness(ps, x="rank", measures="Shannon") +
     theme_bw()  +
     geom_boxplot() +
     ggtitle(annotation$pRnames)
 if(annotation$is.16.S%in%"18S"){
     p <- pl + theme(strip.background = element_rect(fill="orange"))
 } else{
     p <- pl + theme(strip.background = element_rect(fill="blue"))
 }
 if(as.numeric(as.character(annotation$r.pval))<0.05){
     p <- p + theme(panel.background = element_rect(fill = "palevioletred2"))
 } else {
     p
 }
}

pdf("figures/LMErichness_Rank.pdf", width=60, height=20)
do.call("grid.arrange",
        c(lapply(seq_along(FAps.l), function (i) {
            lme.div.plots(FAps.l[[i]], annotation[i,])
        }),
        ncol=10))
dev.off()

lapply(FAps.l,  function(ps){
    tab <- table(tax_table(ps)[,2])
    tab[order(tab)]
})


Para.ps.l <- lapply(FAps.l, function(ps){
    tryCatch(
        subset_taxa(ps, Phylum %in% c("Apicomplexa", "Platyhelminthes",
                                      "Nematoda")),
        error = function(e) NULL)
})

do.call(rbind, lapply(Para.ps.l, function(ps){
    tryCatch(get.shannon.lme(ps),
             error = function (e) NA)
}))


pdf("figures/Parasite_richnesRank43.pdf", width=14, height=8)
plot_richness(Para.ps.l[[43]], x="rank") + theme_bw()  + geom_boxplot()
dev.off()

pdf("figures/Parasite_richnesRank45.pdf", width=14, height=8)
plot_richness(Para.ps.l[[45]], x="rank") + theme_bw()  + geom_boxplot()
dev.off()

NonPrey.ps.l <- lapply(FAps.l, function(ps){
    tryCatch(
        subset_taxa(ps, !Phylum %in% c("Chordata", "Vertebrata", 
        "Chlorophyta_ph", "Phragmoplastophyta",
        "Streptophyta", "Chlorophyta", "Arthropoda",
        "", "undef")),
        error = function(e) NULL)
})

do.call(rbind, lapply(NonPrey.ps.l, function(ps){
    tryCatch(get.shannon.lme(ps),
             error = function (e) NA)
}))

Prey.ps.l <- lapply(FAps.l, function(ps){
    tryCatch(
        subset_taxa(ps, Phylum %in% c("Chordata", "Vertebrata")),
        error = function(e) NULL)
})

do.call(rbind, lapply(Prey.ps.l, function(ps){
    tryCatch(get.shannon.lme(ps),
             error = function (e) NA)
}))

sapply(get_taxa_unique(FAps.l[[43]], "Phylum"), function (p) {
    p <- "Apicomplexa"
    p.sub <- subset_taxa(FAps.l[[43]], Phylum%in%"Apicomplexa")
    sum(otu_table(p.sub))
})


################### ORDINATIONS ##################

## log transformaitons on sample counts 
## lets go direclty for the subsetted data
log.ps.l <- lapply(FAps.l, transform_sample_counts, function(x) log(1 + x))
log.Phyl.ps.l <- lapply(FAPhyl.ps.l, transform_sample_counts, function(x) log(1 + x))
log.TSps.l <- lapply(FATSps.l, transform_sample_counts, function(x) log(1 + x))
log.Gen.l <- lapply(FAGen.l, transform_sample_counts, function(x) log(1 + x))
log.Dist.l <- lapply(FADist.l, transform_sample_counts, function(x) log(1 + x))

## some of the amplicons give errors
ord.works <- c(1:3, 5:10, 12:18, 20:30, 32, 34:36, 40, 41, 43:48)

## for the moment only interested in 16s
amp.16S <- which(annotation$is.16.S%in%"16S")
# for  4 it does not work somehow...
amp.16S <- c(3, 5,13)


get.ord.plot <- function(ps, meth = "MDS", dist = "bray"){
    or <- ordinate(ps, method = meth, distance = dist)
    evals <- or$eig
    plot_ordination(ps, or, color = "rank",
                    shape = "pack") +
        theme_bw()
}

logBC.l <- lapply(log.ps.l[amp.16S], get.ord.plot)
logSSBC.l <- lapply(log.Phyl.ps.l[amp.16S], get.ord.plot)
logTSBC.l <- lapply(log.TSps.l[amp.16S], get.ord.plot)
logGenBC.l <- lapply(log.Gen.l[amp.16S], get.ord.plot)
logDistBC.l <- lapply(log.Dist.l[amp.16S], get.ord.plot)

pdf("figures/ordinate1.pdf", width=15, height=7)
do.call("grid.arrange", c(logBC.l, ncol=2))
dev.off()

pdf("figures/ordinate2.pdf", width=15, height=7)
do.call("grid.arrange", c(logSSBC.l, ncol=2))
dev.off()

pdf("figures/ordinate3.pdf", width=15, height=7)
do.call("grid.arrange", c(logTSBC.l, ncol=2))
dev.off()

pdf("figures/ordinate4.pdf", width=15, height=7)
do.call("grid.arrange", c(logGenBC.l, ncol=2))
dev.off()

pdf("figures/ordinate5.pdf", width=15, height=7)
do.call("grid.arrange", c(logDistBC.l, ncol=2))
dev.off()

logWUF.l <- lapply(log.ps.l[amp.16S], get.ord.plot,
                   meth="PCoA" , dist="wunifrac")

logSSWUF.l <- lapply(log.Phyl.ps.l[amp.16S], get.ord.plot,
                     meth="PCoA" , dist="wunifrac")

logTSWUF.l <- lapply(log.TSps.l[amp.16S], get.ord.plot,
                     meth="PCoA" , dist="wunifrac")

logGenWUF.l <- lapply(log.Gen.l[amp.16S], get.ord.plot,
                      meth="PCoA" , dist="wunifrac")

logDistWUF.l <- lapply(log.Dist.l[amp.16S], get.ord.plot,
                       meth="PCoA" , dist="wunifrac")

pdf("figures/ordinate6.pdf", width=15, height=7)
do.call("grid.arrange", c(logWUF.l, ncol=2))
dev.off()

pdf("figures/ordinate7.pdf", width=15, height=7)
do.call("grid.arrange", c(logSSWUF.l, ncol=2))
dev.off()

pdf("figures/ordinate8.pdf", width=15, height=7)
do.call("grid.arrange", c(logTSWUF.l, ncol=2))
dev.off()

pdf("figures/ordinate9.pdf", width=15, height=7)
do.call("grid.arrange", c(logGenWUF.l, ncol=2))
dev.off()

pdf("figures/ordinate10.pdf", width=15, height=7)
do.call("grid.arrange", c(logDistWUF.l, ncol=2))
dev.off()


## ## Check ccpna
## ordX <- lapply(log.ps.l, function(ps){
##     ps.ccpna <- ordinate(ps, "CCA", formula = ps ~ rank + pack)
##     ps.scores <- vegan::scores(ps.ccpna)
##     sites <- data.frame(ps.scores$sites)
##     sites$SampleID <- rownames(sites)
##     sites <- merge(sites, sample_data(ps), by=0, all.x=TRUE)
##     species <- data.frame(ps.scores$species)
##     species$otu_id <- seq_along(colnames(otu_table(ps)))
##     species <- left_join(species, tax)
##     evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
##     ggplot() +
##         geom_point(data = row_scores,
##                    aes(x = li.Axis1, y = li.Axis2), shape = 2) +
##         geom_point(data = col_scores,
##                    aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
##                    size = .3, alpha = 0.6) +
##         scale_color_brewer(palette = "Set2") +
##         facet_grid(~ rank) +
##         guides(col = guide_legend(override.aes = list(size = 3))) +
##         labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
##              y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
##         coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
##         theme(panel.border = element_rect(color = "#787878",
##                                           fill = alpha("white", 0)))
## })

## ## Supervised learning
library(caret)

predict.pls <- function(pslog, method){
    dataMatrix <- data.frame(rank = sample_data(pslog)$rank, otu_table(pslog))
    ## take 20 mice at random to be the training set, and the remaining 24
    ## the test set
    trainingHyena <- sample(unique(sample_data(pslog)$Subject), size = 20)
    inTrain <- which(sample_data(pslog)$Subject %in% trainingHyena)
    training <- dataMatrix[inTrain,]
    testing <- dataMatrix[-inTrain,]
    plsFit <- train(rank ~ ., data = training,
                    method = method, preProc = "center")
    plsClasses <- predict(plsFit, newdata = testing)
    table(plsClasses, testing$rank)
}

PLS <- lapply(log.Gen.l[c(1, 16:18)], predict.pls, "pls")

PLS

## [[1]]
          
## plsClasses high low
##       high   16  12
##       low     5   8

## [[2]]
          
## plsClasses high low
##       high    7   4
##       low    17  13

## [[3]]
          
## plsClasses high low
##       high   14   9
##       low     8  10

## [[4]]
          
## plsClasses high low
##       high   14  12
##       low     8   7



RF <- lapply(log.TSps.l[c(1, 16:18)], predict.pls, "rf")

RF

## [[1]]
          
## plsClasses high low
##       high   13   4
##       low    13  11

## [[2]]
          
## plsClasses high low
##       high    0   0
##       low    26  15

## [[3]]
          
## plsClasses high low
##       high   11  10
##       low    12   8

## [[4]]
          
## plsClasses high low
##       high    6   4
##       low    17  14



################## hierarchical testing ###############
get.hier.test <- function(ps){
    ## ## Playing with DEseq2
    ## ps_dds <-  phyloseq_to_deseq2(ps, ~ rank + pack)
    ## ts <- counts(ps_dds)
    ## geoMeans <- apply(ts, 1, function(row) {
    ##     if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
    ## })
    ## ps_dds <- estimateSizeFactors(ps_dds, geoMeans=geoMeans)
    ## ps_dds <- estimateDispersions(ps_dds)
    ## abund <- getVarianceStabilizedData(ps_dds)
    abund <- t(otu_table(ps))
    short_names <- make.names(substr(rownames(abund), 1, 5),
                              unique=TRUE)
    rownames(abund) <- short_names
    el <- phy_tree(ps)$edge
    el0 <- el
    el0 <- el0[nrow(el):1, ]
    el_names <- c(short_names, seq_len(phy_tree(ps)$Nnode))
    el[, 1] <- el_names[el0[, 1]]
    el[, 2] <- el_names[as.numeric(el0[, 2])]
    unadj_p <- treePValues(el, abund, sample_data(ps)$rank)
    ## Warning: essentially perfect fit: summary may be unreliable
    hfdr_res <- hFDR.adjust(unadj_p, el, 0.99)
    tax <- data.frame(tax_table(ps)[, c("Family", "Genus")])
    tax$seq <- short_names
    hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
    lj <- left_join(tax, hfdr_res@p.vals)
    head(arrange(lj, adjp), n=50)
}

## for a full test
##
## XYZ <- mclapply(log.ps.l, tryCatch(get.hier.test, error = function
##                 (e) NA), mc.cores=20)
##
## which might not make sense as we are interested in things with a
## clear annotation here anyways

tests.gen.coll <- mclapply(log.Gen.l, tryCatch(get.hier.test, error = function (e) NA),
                           mc.cores=20)

lapply(tests.gen.coll, function (x) x[x$adjp<0.05& !is.na(x$adjp),])

## ## Not so
## [[17]]
##                      Family                        Genus      seq       unadjp        adjp
## 1           Eggerthellaceae                      Slackia CAGTG.24 0.0001681410 0.000336282
## 2       Erysipelotrichaceae  Erysipelotrichaceae_UCG-001  CAGTA.6 0.0005669904 0.001133981
## 3     Peptostreptococcaceae                   Romboutsia CAGTG.30 0.0067334885 0.013466977
## 4 Bacteroidales_S24-7_group Bacteroidales_S24-7_group_ge CAGTG.38 0.0135579717 0.027115943
## 5       Erysipelotrichaceae  Erysipelotrichaceae_UCG-004  CAGTA.7 0.0159865623 0.031973125

## Slackia:
## Couple of papers, among them: Characterization of Slackia exigua
## isolated from human wound infections, including abscesses of
## intestinal origin

## Erysipelotrichaceae: 
## Wikipedia Erysipelotrichia are a class of bacteria of the phylum
## Firmicutes. Species of this class are known to be common in the gut
## microbiome, as they have been isolated from swine manure [1] and
## increase in composition of the mouse gut microbiome for mice
## switched to diets high in fat.[2]

## Romboutsia
##  bacterium isolated from the right human colon by colonoscopy in a
## 63-year-old French man with severe anaemia with melaena.

## Bacteroidales_S24-7:
## representatives constitute a substantial component of the murine
## gut microbiota, as well as being present within the human gut



table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]], Genus%in%"Slackia")))>0,
      sample_data(log.Gen.l[[4]])$rank)

##       high low
## FALSE    9  20
## TRUE    24   8

table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]], Genus%in%"Slackia"))),
      sample_data(log.Gen.l[[4]])$rank)

##                  high low
## 0                   9  20
## 1.09861228866811    1   0
## 1.6094379124341     2   1
## 2.07944154167984    1   0
## 2.19722457733622    3   1
## 2.30258509299405    3   2
## 2.39789527279837    1   2
## 2.56494935746154    0   1
## 2.63905732961526    2   1
## 2.77258872223978    3   0
## 2.99573227355399    1   0
## 3.29583686600433    1   0
## 3.3322045101752     1   0
## 3.40119738166216    2   0
## 3.43398720448515    1   0
## 3.52636052461616    1   0
## 3.76120011569356    1   0

table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]],
                                    Genus%in%"Erysipelotrichaceae_UCG-001")))>0,
      sample_data(log.Gen.l[[4]])$rank)

##       high low
## FALSE   15  24
## TRUE    18   4


table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]],
                                    Genus%in%"Erysipelotrichaceae_UCG-001"))),
      sample_data(log.Gen.l[[4]])$rank)

##                  high low
## 0                  15  24
## 1.6094379124341     2   0
## 1.79175946922805    2   1
## 1.94591014905531    4   1
## 2.07944154167984    1   2
## 2.19722457733622    1   0
## 2.30258509299405    3   0
## 2.484906649788      2   0
## 2.56494935746154    2   0
## 2.77258872223978    1   0


table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]],
                                    Genus%in%"Romboutsia")))>0,
      sample_data(log.Gen.l[[4]])$rank)

##       high low
## FALSE   31  19
## TRUE     2   9

table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]],
                                    Genus%in%"Romboutsia"))),
      sample_data(log.Gen.l[[4]])$rank)

##                  high low
## 0                  31  19
## 2.07944154167984    0   1
## 2.19722457733622    0   1
## 2.70805020110221    1   0
## 2.83321334405622    0   1
## 2.99573227355399    1   0
## 3.09104245335832    0   1
## 3.13549421592915    0   1
## 3.40119738166216    0   1
## 3.66356164612965    0   1
## 3.87120101090789    0   1
## 3.97029191355212    0   1



table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]],
                                    Genus%in%"Bacteroidales_S24-7_group_ge")))>0,
      sample_data(log.Gen.l[[4]])$rank)

##       high low
## FALSE   32  22
## TRUE     1   6


table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]],
                                    Genus%in%"Bacteroidales_S24-7_group_ge"))),
      sample_data(log.Gen.l[[4]])$rank)

##                   high low
## 0                   32  22
## 0.693147180559945    1   1
## 1.79175946922805     0   2
## 2.30258509299405     0   1
## 2.70805020110221     0   1
## 3.36729582998647     0   1



table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]],
                                    Genus%in%"Erysipelotrichaceae_UCG-004")))>0,
      sample_data(log.Gen.l[[4]])$rank)

##       high low
## FALSE   29  18
## TRUE     4  10

table(rowSums(otu_table(subset_taxa(log.Gen.l[[17]],
                                    Genus%in%"Erysipelotrichaceae_UCG-004"))),
      sample_data(log.Gen.l[[4]])$rank)

##                  high low
## 0                  29  18
## 1.38629436111989    0   2
## 1.6094379124341     2   0
## 1.79175946922805    0   1
## 2.19722457733622    0   1
## 2.484906649788      1   1
## 2.56494935746154    0   1
## 2.63905732961526    1   0
## 3.2188758248682     0   1
## 3.29583686600433    0   1
## 4.4188406077966     0   1
## 4.49980967033027    0   1



tests.dist.coll <- mclapply(log.Dist.l, tryCatch(get.hier.test, error = function (e) NA),
                            mc.cores=20)

lapply(tests.dist.coll, function (x) x[x$adjp<0.05& !is.na(x$adjp),])

the.test <- lapply(tests.gen.coll,
                   function (x) x[x$adjp<0.05& !is.na(x$adjp),])
## Linguatula




