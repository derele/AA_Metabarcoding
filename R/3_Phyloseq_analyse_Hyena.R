library(nlme)
library(ggplot2)
library(exactRankTests)
library(phyloseq)
library(RSvgDevice)

################# From here Phyloseq ######################
load(file="/SAN/Metabarcoding/phlyoSeq_list.Rdata") ## -> ps.l

## Briefly, oberall number of reads for different amplicons
sum.seq <- lapply(ps.l, function (x){
    sum(otu_table(x))
})

sum.seq <- unlist(sum.seq)

sum(sum.seq)
mean(sum.seq, rm.na=TRUE)

range(sum.seq)


## oberall number of reads
num.RSVs <- lapply(ps.l, function (x){
    ncol(otu_table(x))
})

num.RSVs <- unlist(num.RSVs)

sum(num.RSVs)
mean(num.RSVs, rm.na=TRUE)

range(num.RSVs)

## Read numbersbetween samples before normalization and merging of technical replicates
load(file="/SAN/Metabarcoding/phlyoSeq_cat.Rdata") ## -> PSA

sample.seq <- rowSums(otu_table(PSA))
range(sample.seq)


## After merging of technical replicates and normalization 
load(file="/SAN/Metabarcoding/phlyoSeq_Hy.Rdata") ## -> PS

keep.PS <- colSums(otu_table(PS))>0

PS <- subset_taxa(PS, keep.PS)
PS <- subset_taxa(PS, Kingdom %in% c("Bacteria", "Eukaryota"))

PS.genus <- tax_glom(PS, "Genus", NArm = TRUE)

Genus_table <- cbind(tax_table(PS.genus), t(otu_table(PS.genus)))
rownames(Genus_table) <- NULL

write.csv(Genus_table,
          "/home/ele/Genus_abundance_table_export.csv", quote=FALSE,
          row.names=FALSE)


PS.phylum <- tax_glom(PS, "Phylum", NArm = TRUE)

PS.class <- tax_glom(PS, "Class", NArm = TRUE)

PS.order <- tax_glom(PS, "Order", NArm = TRUE)

PS.family <- tax_glom(PS, "Family", NArm = TRUE)

## Genera and RSVs Bacteria 
subset_taxa(PS.genus, Kingdom %in% "Bacteria")
subset_taxa(PS, Kingdom %in% "Bacteria")

## Genera and RSVs Eukaryota 
subset_taxa(PS, Kingdom %in% "Eukaryota")
subset_taxa(PS.genus, Kingdom %in% "Eukaryota")

## Bacteria read numbes 
sum(otu_table(subset_taxa(PS, Kingdom %in% "Bacteria")))
sum(otu_table(subset_taxa(PS.genus, Kingdom %in% "Bacteria")))

## Eukaryote read numbers 
sum(otu_table(subset_taxa(PS, Kingdom %in% "Eukaryota")))
sum(otu_table(subset_taxa(PS.genus, Kingdom %in% "Eukaryota")))

## ## Abundance vs. Intensity
## Genus.Abu <- unname(colSums(otu_table(PS.genus)))
## Genus.Prev <- unname(apply(otu_table(PS.genus), 2,
##                            function (x) sum(x>0)))

## Prev.Abu.Genus <- data.frame(cbind(
##     Abu=Genus.Abu,
##     Prev=Genus.Prev,
##     tax_table(PS.genus)))

## unicellular_phy<- c("Apicomplexa",  "Dinoflagellata", "Microsporidia")

## multicellular_phy <- c("Arthropoda",  "Chordata",
##                        "Platyhelminthes",  "Vertebrata")

## Prev.Abu.Genus <- Prev.Abu.Genus[Prev.Abu.Genus$Phylum%in%c(unicellular_phy,
##                                                             multicellular_phy), ]# |
## ##                                 Prev.Abu.Genus$Kingdom%in%"Bacteria", ]
                                                            
## Prev.Abu.Genus$Prev <- as.numeric(as.character(Prev.Abu.Genus$Prev))
## Prev.Abu.Genus$Abu <- as.numeric(as.character(Prev.Abu.Genus$Abu))

## ggplot(Prev.Abu.Genus, aes(Abu, Prev, color=Phylum)) +
##     geom_point() +
## ##     scale_y_log10() +
##     scale_x_log10() +
##     geom_smooth(method="lm", se=FALSE)


## RSV.abu <- unname(colSums(otu_table(PS)))


## What phyla are detected
table(tax_table(PS)[, "Phylum"])

tax.frame <- data.frame(tax_table(PS.genus))
rownames(tax.frame) <- NULL

## SOME FIXES FOR TAXONOMY
table(tax.frame$Phylum)
tax.frame$Phylum[tax.frame$Phylum%in%"Vertebrata"] <-
    "Chordata"

tax.frame$Phylum[tax.frame$Phylum%in%"Chlorophyta_ph"] <-
    "Chlorophyta"

tax.frame$Phylum[tax.frame$Phylum%in%"Euglenida"] <-
    "Euglenozoa"           

tax.frame$Phylum[tax.frame$Phylum%in%"Spirochaetes"] <-
    "Spirochaetae"         

tax.frame$counts <- unname(colSums(otu_table(PS.genus)))

spurious.phyla <- names(table(tax.frame$Phylum)[table(tax.frame$Phylum)<2])

## removes spurious phyla for which only one 
tax.frame <- tax.frame[!tax.frame$Phylum%in%spurious.phyla, ]
tax.frame <- tax.frame[!tax.frame$Phylum%in%"undef", ]

tax.frame$Phylum <- droplevels(tax.frame$Phylum)

Phy.roles <- read.csv("./Euk_phyla.csv")
names(Phy.roles) <- c("Phylum", "role")

tax.frame <- merge(tax.frame, Phy.roles, all.x=TRUE, by="Phylum")

## here what we have in total RSVs
as.data.frame(table(tax_table(PS)[, "Kingdom"]))

RSVcounts <- as.data.frame(table(tax_table(PS)[, "Phylum"]))
names(RSVcounts) <- c("Phylum", "numRSVs")

RSVcounts$Phylum[RSVcounts$Phylum%in%"Vertebrata"] <-
    "Chordata"

RSVcounts$Phylum[RSVcounts$Phylum%in%"Chlorophyta_ph"] <-
    "Chlorophyta"

RSVcounts$Phylum[RSVcounts$Phylum%in%"Euglenida"] <-
    "Euglenozoa"           

RSVcounts$Phylum[RSVcounts$Phylum%in%"Spirochaetes"] <-
    "Spirochaetae"         

RSVcounts[RSVcounts$Phylum%in%euk_phyla& !RSVcounts$Phylum%in%spurious.phyla ,]

## the rest are not assigned to phylum
sum(RSVcounts[RSVcounts$Phylum%in%bak_phyla,"numRSVs"])
sum(RSVcounts[RSVcounts$Phylum%in%euk_phyla,"numRSVs"])


sum(otu_table(subset_taxa(PS, Kingdom %in% "Bacteria")))-
sum(otu_table(subset_taxa(PS.genus, Kingdom %in% "Bacteria")))


sum(otu_table(subset_taxa(PS, Kingdom %in% "Eukaryota")))-
sum(otu_table(subset_taxa(PS.genus, Kingdom %in% "Eukaryota")))

tax.frame.bac <- tax.frame[tax.frame$Kingdom%in%"Bacteria", ]
tax.frame.euk <- tax.frame[tax.frame$Kingdom%in%"Eukaryota", ]


devSVG("figures/figures_Hyena/Figure1a_taxon_dist.svg", width=12, height=6)
ggplot(tax.frame, aes(x=Phylum, color=role, y=..count..)) +
    geom_bar() +
    coord_polar() +
    scale_y_continuous("Genera annotated per phylum") +
    facet_grid(.~Kingdom, scales="free_x", drop=TRUE) +
    theme_bw()
dev.off()

devSVG("figures/figures_Hyena/Figure1b_taxon_weighted.svg", width=12, height=6)
ggplot(tax.frame, aes(x=Phylum, y=..count.., weight=counts, color=role)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Seqeuncing read counts per phylum") +
    facet_grid(.~Kingdom, scales="free_x", drop=TRUE)+
    theme_bw()
dev.off()

showP <- cbind(as.character(unique(sort(tax.frame$Phylum[tax.frame$Kingdom%in%"Eukaryota"]))))

write.csv(showP, "Euk_phyla.csv", row.names=FALSE, quote=FALSE)

library(xtable)

table1_Tax <- cbind(
    Kingdom =tapply(tax.frame$Kingdom, tax.frame$Phylum,
                    function (x) as.character(unique(x))),
    genera = table(tax.frame$Phylum),
    reads = tapply(tax.frame$counts, tax.frame$Phylum, sum))


euk_phyla <- names(table(tax_table(subset_taxa(PS, Kingdom %in% "Eukaryota"))[, "Phylum"]))

bac_phyla <- names(table(tax_table(subset_taxa(PS, Kingdom %in% "Bacteria"))[, "Phylum"]))

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

## working with worm counts ####################################################
######## WC - Worm Counts ######################################################
################################################################################
WC <- read.csv("/home/ele/Documents/Hyena_Hartmann_MS/Hyena_WormCounts_fixed.csv",
               as.is=TRUE)

colnames(WC)[2:ncol(WC)] <- paste0("count_", colnames(WC)[2:ncol(WC)])

WC[, 2:ncol(WC)] <- sapply(WC[, 2:ncol(WC)],
                           function(x) as.numeric(as.character(x)))

WC$count_Cystoisospora<- round(WC$count_Cystoisospora_small +
                               WC$count_Cystoisospora2)

WC$count_Ancylostoma<- round(WC$count_Ancylostoma_small+ WC$
                             count_Ancylostoma2)

WC$count_Spirometra <- round(WC$count_Spirometra)

## extract only what is usable
WC <- WC[, c("ID.Hyena", "count_Cystoisospora", "count_Cystoisospora_small",
             "count_Cystoisospora2",
             "count_Ancylostoma", "count_Ancylostoma_small", "count_Ancylostoma2",
             "count_Spirometra")]

## use only the samples we have seqdata for

WC <- WC[WC$ID.Hyena%in%sample_data(PS)$Hyena.ID,]

apply(WC, 2, function (x) summary.factor(as.numeric(as.character(x))>0))

## For selected genera and whole Phyla

Genus.tab <- otu_table(PS.genus)
Genus.cols <- make.names(make.unique(unname(tax_table(PS.genus)[, "Genus"])))
colnames(Genus.tab) <- Genus.cols

Family.tab <- otu_table(PS.family)
Family.cols <- make.names(make.unique(unname(tax_table(PS.family)[, "Family"])))
colnames(Family.tab) <- Family.cols

Order.tab <- otu_table(PS.order)
Order.cols <- make.names(make.unique(unname(tax_table(PS.order)[, "Order"])))
colnames(Order.tab) <- Order.cols

Class.tab <- otu_table(PS.class)
Class.cols <- make.names(make.unique(unname(tax_table(PS.class)[, "Class"])))
colnames(Class.tab) <- Class.cols

Phylum.tab <- otu_table(PS.phylum)
Phylum.cols <- make.names(make.unique(unname(tax_table(PS.phylum)[, "Phylum"])))
colnames(Phylum.tab) <- Phylum.cols

tax.tab <- cbind(Phylum.tab, Class.tab, Order.tab, Family.tab, Genus.tab)
tax.cols <- make.names(make.unique(colnames(tax.tab)))

All.data <- data.frame(cbind(sample_data(PS.genus), tax.tab))
All.data <- merge(WC, All.data, by.x="ID.Hyena", by.y="Hyena.ID")

## The raw taxonomy information to compare against
tax.raw <- data.frame(tax_table(PS))

get.lower.in.hihger.tax <- function(level, higher, lower.level){
    ss <- tax.raw[tax.raw[, level]%in%higher, ]
    names(table(as.character(ss[, lower.level])))
}

get.lower.in.hihger.tax("Order", "Rhabditida", "Genus")
    
############# Untargeted #################
count.col <- grep("count_", colnames(All.data), value=TRUE)

spear.cor <- data.frame(t(cor(All.data[, count.col], All.data[, tax.cols],
                              method="spearman",
                              use="pairwise.complete.obs")))

top.cors.spear <- lapply(count.col, function (x){
        head(spear.cor[order(spear.cor[,x], decreasing=TRUE), x, drop=FALSE], n=10)
})

names(top.cors.spear) <- count.col

## Ancylostoma
cor.test(All.data$Ancylostoma, All.data$count_Ancylostoma,
         use="pairwise.complete.obs", method="spearman")

cor.test(All.data$Rhabditida, All.data$count_Ancylostoma,
    use="pairwise.complete.obs", method="spearman")

cor.test(All.data$Ancylostoma + All.data$Haemonchus, ## + All.data$Ostertagia,
         All.data$count_Ancylostoma,
         use="pairwise.complete.obs", method="spearman")

cor.test(All.data$Ancylostoma +  All.data$Ostertagia, 
         All.data$count_Ancylostoma,
         use="pairwise.complete.obs", method="spearman")

devtools::source_gist("524eade46135f6348140",
                      filename = "ggplot_smooth_func.R")

pdf("figures/figures_Hyena/Figure2a_Ancylostoma.pdf")
ggplot(All.data, aes(Rhabditida+1 , count_Ancylostoma+1,
                     label=ID.Hyena)) +
    geom_point() +
    stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE) +
    stat_smooth(method="lm", se=FALSE) +
    scale_x_log10("Sequence abundance") +
    scale_y_log10("FEC") +
    ggtitle("Ancylostoma (Rhabditida)") +
    theme_bw()
dev.off()

lm(log10(Rhabditida+1)~log10(count_Ancylostoma+1), data=All.data)

glm(Rhabditida~count_Ancylostoma, data=All.data, family="poisson")

Rhab.genera <- get.lower.in.hihger.tax("Order", "Rhabditida", "Genus")

Rhab.counts <- colSums(All.data[, Rhab.genera])

round(sort(Rhab.counts/sum(Rhab.counts)*100, decreasing=TRUE), 2)


## Spirometra 
cor.test(All.data$Diphyllobothriidea, All.data$count_Spirometra,
         use="pairwise.complete.obs", method="spearman")

cor.test(All.data$Spirometra, All.data$count_Spirometra,
         use="pairwise.complete.obs", method="spearman")

cor.test(All.data$Diphyllobothrium, All.data$count_Spirometra,
         use="pairwise.complete.obs", method="spearman")

pdf("figures/figures_Hyena/Figure2b_Spirometra.pdf")
ggplot(All.data, aes(Diphyllobothriidea+1, count_Spirometra+1, label=ID.Hyena)) +
    geom_point() +
    stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE) +
    stat_smooth(method="lm", se=FALSE) +
    scale_x_log10("Sequence Abundance") +
    scale_y_log10("FEC") +
    ggtitle("Spirometra (Diphyllobothriidea)")+
    theme_bw()
dev.off()

Dip.genera <- get.lower.in.hihger.tax("Order", "Diphyllobothriidea", "Genus")

Dip.counts <- colSums(All.data[, Dip.genera])

round(sort(Dip.counts/sum(Dip.counts)*100, decreasing=TRUE), 2)

range(All.data[All.data$count_Spirometra==0 & All.data$Diphyllobothriidea>0,
               "Diphyllobothriidea"], na.rm=TRUE)

## Best combinations for "_small"

All.data$Cocci.comb <-
    All.data$Eimeriidae +
    All.data$Besnoitia +
    All.data$Toxoplasma

cor.test((All.data$Coccidia),
         All.data$count_Cystoisospora_small,
         use="pairwise.complete.obs", method="spearman")

cor.test((All.data$Cocci.comb),
         All.data$count_Cystoisospora_small,
         use="pairwise.complete.obs", method="spearman")

### Best combinations for "2"
cor.test(All.data$Eimeria.1,
         All.data$count_Cystoisospora2,
         use="pairwise.complete.obs", method="spearman")

pdf("figures/figures_Hyena/Figure2c_Coccidia_small.pdf")
ggplot(All.data, aes(Cocci.comb +1, 
                     count_Cystoisospora_small+1, label=ID.Hyena)) +
    geom_point() +
    stat_smooth_func(geom="text", method="lm", hjust=0, parse=TRUE) +
    stat_smooth(method="lm", se=FALSE) +
    scale_x_log10("Sequence Abundance") +
    scale_y_log10("FEC") +
    ggtitle("Coccidia (small oocyst types)") +
    theme_bw()
dev.off()

summary(lm(log10(count_Cystoisospora_small+1)~log10(Cocci.comb +1),
           data=All.data))

pdf("figures/figures_Hyena/Figure2d_Coccidia_larger.pdf")
ggplot(All.data, aes(Eimeria.1+1, count_Cystoisospora2+1,
                     label=ID.Hyena)) +
    geom_point()+
    stat_smooth_func(geom="text", method="lm", hjust=-1.5, parse=TRUE) +
    stat_smooth(method="lm", se=FALSE) +
    scale_x_log10("Sequence Abundance") +
    scale_y_log10("FEC") +
    ggtitle("Coccidia (large oocyst types)") +
    theme_bw()
dev.off()

################ Richness, diversity and evenness ######################
load(file="/SAN/Metabarcoding/phlyoSeq_Hy_rare.Rdata") # -> PS.rare,
PS.r.genus <- tax_glom(PS.rare, "Genus", NArm = TRUE)

test.Chao.rank <- function(ps){
    mes <- c("Observed", "Chao1", "Shannon")
    est <- estimate_richness(ps, measures=mes)
    df <- merge(sample_data(ps), est, by=0)
    df$pilou <- df$Shannon/log(df$Observed)

    test.ric <- wilcox.exact(Observed ~ rank, data = df, exact = TRUE)
    med.ric <- tapply(df$Observed, df$rank, median)
    test.dif <- wilcox.exact(Chao1 ~ rank, data = df, exact = TRUE)
    med.div <- tapply(df$Chao1, df$rank, median)
    test.eve <- wilcox.exact(pilou ~ rank, data = df, exact = TRUE)
    med.eve <- tapply(df$pilou, df$rank, median)
    list(Obs.ric = list(W=test.ric$statistic, p.value=test.ric$p.value,
                        med.high=med.ric["high"], med.low=med.ric["low"]),
         Cha.div = list(W=test.dif$statistic, p.value=test.dif$p.value,
                        med.high=med.div["high"], med.low=med.div["low"]),
         Pil.eve = list(W=test.eve$statistic, p.value=test.eve$p.value,
                        med.high=med.eve["high"], med.low=med.eve["low"]))
}


test.Chao.age <- function(ps){
    mes <- c("Observed", "Chao1", "Shannon")
    est <- estimate_richness(ps, measures=mes)
    df <- merge(sample_data(ps), est, by=0)
    df$pilou <- df$Shannon/log(df$Observed)

    test.ric <- wilcox.exact(Observed ~ age, data = df, exact = TRUE)
    test.dif <- wilcox.exact(Chao1 ~ age, data = df, exact = TRUE)
    test.eve <- wilcox.exact(pilou ~ age, data = df, exact = TRUE)

    test.ric <- wilcox.exact(Observed ~ age, data = df, exact = TRUE)
    med.ric <- tapply(df$Observed, df$age, median)
    test.dif <- wilcox.exact(Chao1 ~ age, data = df, exact = TRUE)
    med.div <- tapply(df$Chao1, df$age, median)
    test.eve <- wilcox.exact(pilou ~ age, data = df, exact = TRUE)
    med.eve <- tapply(df$pilou, df$age, median)
    list(Obs.ric = list(W=test.ric$statistic, p.value=test.ric$p.value,
                        med.Adult=med.ric["Adult"], med.Cub=med.ric["Cub"]),
         Cha.div = list(W=test.dif$statistic, p.value=test.dif$p.value,
                        med.Adult=med.div["Adult"], med.Cub=med.div["Cub"]),
         Pil.eve = list(W=test.eve$statistic, p.value=test.eve$p.value,
                        med.Adult=med.eve["Adult"], med.Cub=med.eve["Cub"]))
}

do.call(rbind, test.Chao.rank(PS.rare))
do.call(rbind, test.Chao.rank(PS.r.genus))

PS.rare.E <- subset_taxa(PS.rare,
                         !Phylum %in% c(clear.prey.phyla, undef) &
                         Kingdom %in% "Eukaryota")

PS.r.g.E <- subset_taxa(PS.r.genus,
                        !Phylum %in% c(clear.prey.phyla, undef) &
                        Kingdom %in% "Eukaryota")


do.call(rbind, test.Chao.rank(PS.rare.E))

do.call(rbind, test.Chao.rank(PS.r.g.E))



PS.rare.est <- estimate_richness(PS.rare.E, measures=c("Observed", "Chao1", "Shannon"))
PS.rare.df <- merge(sample_data(PS.rare.E), PS.rare.est, by=0)
PS.rare.df$pilou <- PS.rare.df$Shannon/log(PS.rare.df$Observed)
PS.rare.mdf = reshape2::melt(PS.rare.df,
                             measure.vars = c("Observed", "Chao1", "pilou"))

pdf("figures/figures_Hyena/Figure4a_euk_rank_diversity_RSVs.pdf")
ggplot(PS.rare.mdf, aes(rank, value)) +
    geom_boxplot(na.rm = TRUE) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0))+
    ylab("Measure")+    
    facet_wrap(~variable, nrow = 1, scales="free_y")+
    ggtitle("Richness, diversity and evenness for rarified RSV counts")
dev.off()


do.call(rbind, test.Chao.rank(
                   subset_taxa(PS.r.genus,
                               !Phylum %in% c(clear.prey.phyla, undef) &
                               Kingdom %in% "Eukaryota"
                               )))

PS.r.g.E <- subset_taxa(PS.r.genus,
                         !Phylum %in% c(clear.prey.phyla, undef) &
                         Kingdom %in% "Eukaryota")

do.call(rbind, test.Chao.rank(PS.r.g.E))

PS.r.g.est <- estimate_richness(PS.r.g.E, measures=c("Observed", "Chao1", "Shannon"))
PS.r.g.df <- merge(sample_data(PS.r.g.E), PS.r.g.est, by=0)
PS.r.g.df$pilou <- PS.r.g.df$Shannon/log(PS.r.g.df$Observed)
PS.r.g.mdf = reshape2::melt(PS.r.g.df,
                             measure.vars = c("Observed", "Chao1", "pilou"))

pdf("figures/figures_Hyena/Figure4b_euk_rank_diversity_genera.pdf")
ggplot(PS.r.g.mdf, aes(rank, value)) +
    geom_boxplot(na.rm = TRUE) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0))+
    ylab("Measure")+    
    facet_wrap(~variable, nrow = 1, scales="free_y")+
    ggtitle("Richness, diversity and evenness for rarified genus counts")
dev.off()

do.call(rbind, test.Chao.rank(
                   subset_taxa(PS.rare,
                               Kingdom %in% "Bacteria"
                               )))

do.call(rbind, test.Chao.rank(subset_taxa(
                   PS.r.genus,
                   Kingdom %in% "Bacteria"
               )))

do.call(rbind, test.Chao.rank(subset_taxa(
                   PS.rare, Phylum %in% clear.parasite.phyla)))

do.call(rbind,
        test.Chao.rank(subset_taxa(
            PS.r.genus, Phylum %in% clear.parasite.phyla)))

do.call(rbind,
        test.Chao.rank(subset_taxa(
            PS.r.genus, Phylum %in% clear.prey.phyla)))

do.call(rbind,
        test.Chao.rank(subset_taxa(
            PS.r.genus, Phylum %in% plant.phyla)))

## a post-hoc test for the phyla contributing to diversity differences

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

test.phyla <- euk_phyla[!euk_phyla%in%spurious.phyla]

Chao.phyla <- lapply(test.phyla, function(x){
    ps <- subset_taxaE(PS.r.genus, tax_table(PS.r.genus)[, "Phylum"] %in% x)
    tryCatch(do.call(rbind,test.Chao.rank(ps)), error = function (e) NULL)
})

names(Chao.phyla) <- test.phyla

obs.order <- order(unlist(lapply(Chao.phyla, function (x) x["Obs.ric","p.value"])))

names(Chao.phyla)[obs.order]

chao.order <- order(unlist(lapply(Chao.phyla, function (x) x["Cha.div","p.value"])))

names(Chao.phyla)[chao.order]


########################## Bacteria for age

## test.Chao.age(subset_taxa(PS,
##                          Kingdom %in% "Bacteria"
##                          ))

## reporting for genera
do.call(rbind, test.Chao.age(subset_taxa(PS.rare,
                                         Kingdom %in% "Bacteria"
                                         )))

do.call(rbind, test.Chao.age(subset_taxa(PS.r.genus,
                                         Kingdom %in% "Bacteria"
                                         )))

do.call(rbind, test.Chao.age(subset_taxa(PS.r.genus,
                                         Kingdom %in% "Eukaryota"
                                         )))


PS.r.g.B <- subset_taxa(PS.r.genus,
                         Kingdom %in% "Bacteria")

do.call(rbind, test.Chao.rank(PS.r.g.B))

PS.r.gB.est <- estimate_richness(PS.r.g.B, measures=c("Observed", "Chao1", "Shannon"))
PS.r.gB.df <- merge(sample_data(PS.r.g.B), PS.r.gB.est, by=0)
PS.r.gB.df$pilou <- PS.r.gB.df$Shannon/log(PS.r.gB.df$Observed)
PS.r.gB.mdf = reshape2::melt(PS.r.gB.df,
                             measure.vars = c("Observed", "Chao1", "pilou"))

pdf("figures/figures_Hyena/Figure3a_bak_age_diversity_Generas.pdf")
ggplot(PS.r.gB.mdf, aes(age, value)) +
    geom_boxplot(na.rm = TRUE) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0))+
    ylab("Measure")+    
    facet_wrap(~variable, nrow = 1, scales="free_y")+
    ggtitle("Richness, diversity and evenness for rarified genus counts")
dev.off()




bak_phyla <- names(table(tax_table(subset_taxa(PS.genus,
                                               Kingdom %in% "Bacteria"))[, "Phylum"]))

Chao.phyla.bak <- lapply(bak_phyla, function(x){
    ps <- subset_taxaE(PS.genus, tax_table(PS.genus)[, "Phylum"] %in% x)
    tryCatch(do.call(rbind, test.Chao.age(ps)), error = function (e) NULL)
})

names(Chao.phyla.bak) <- bak_phyla

obs.order.bak <- order(unlist(lapply(Chao.phyla.bak, function (x) x["Obs.ric","p.value"])))

names(Chao.phyla.bak)[obs.order.bak]

chao.order.bak <- order(unlist(lapply(Chao.phyla.bak, function (x) x["Cha.div","p.value"])))

names(Chao.phyla.bak)[chao.order.bak]

Chao.phyla.bak[c("Tenericutes", "Actinobacteria","Bacteroidetes", "Firmicutes")]
## Wow: Ten less in Cubs, Act more in Cubs, Bact less in Cubs


## Bacteria for Sex ... not reported as not able to disentangle from
## age

### Testing individual taxa for correlation with diversity.

## Eukaryotes
EDiv <- PS.rare.est
names(EDiv) <- paste0("Euk_", names(EDiv))

All.data <- merge(All.data, EDiv, by.x= "ID.Hyena", by.y=0) 

## Bacteria
BDiv <- PS.r.gB.est
names(BDiv) <- paste0("Bac_", names(BDiv))

All.data <- merge(All.data, BDiv, by.x= "ID.Hyena", by.y=0) 


div.col <- grep("Euk_|Bac_", colnames(All.data), value=TRUE)

div.spear.cor <- data.frame(t(cor(All.data[, div.col], All.data[, c(tax.cols, count.col)],
                                  method="spearman",
                                  use="pairwise.complete.obs")))


All.female <- All.data[All.data$sex%in%"Female", ]

cor.test(All.female$count_Spirometra, All.female$Euk_Observed, method="spearman")

cor.test(All.female$count_Spirometra, All.female$Bac_Observed, method="spearman")

cor.test(All.female$Spirometra, All.female$Bac_Observed, method="spearman")

cor.test(All.female$Diphyllobothriidae, All.female$Bac_Observed, method="spearman")


cor.test(All.female$count_Ancylostoma, All.female$Euk_Observed, method="spearman")

cor.test(All.female$count_Ancylostoma, All.female$Bac_Observed, method="spearman")

cor.test(All.female$Ancylostoma, All.female$Euk_Observed, method="spearman")

cor.test(All.female$Ancylostoma, All.female$Bac_Observed, method="spearman")

cor.test(All.female$Rhabditida, All.female$Bac_Observed, method="spearman")

cor.test(All.female$Spirometra, All.female$Bac_Observed, method="spearman")

cor.test(All.female$Ancylostoma, All.female$Bac_Observed, method="spearman")

cor.test(All.female$Rhabditida, All.female$Bac_Observed, method="spearman")


## This is a result of Juveniles having lower Dipylidium and lower
## Bacterial diversity:
cor.test(All.female$Dipylidiidae, All.female$Bac_Observed, method="spearman")

cor.test(All.female$Dipylidiidae, All.female$Euk_Observed, method="spearman")


## This might be reason to worry???!!! 
cor.test(All.female$Vertebrata, All.female$Euk_Observed, method="spearman")

cor.test(All.female$Vertebrata, All.female$Bac_Observed, method="spearman")



################### ORDINATIONS ##################
## log transformaitons on sample counts 
## lets go direclty for the subsetted data

## BACTERIA #### FOR AGE:
logBAC <- transform_sample_counts(
    subset_taxa(PS.genus, Kingdom%in%"Bacteria"),
    function(x) log10(1+x))

out.bc.log <- ordinate(logBAC, method = "NMDS", distance = "bray")

pdf("figures/figures_Hyena/Figure3b_Ord_bac.pdf")
plot_ordination(logBAC, out.bc.log, color="age", shape="age") +
    labs(col = "Binned Age") +
    scale_shape_discrete(solid=FALSE) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    ggtitle("Ordination on log 1+x transformed abundance of bacterial genera")
dev.off()

## ## Supervised learning
library(caret)

dataMatrixBac <- data.frame(age = sample_data(logBAC)$age,
                            otu_table(logBAC))
plsFitBac <- train(age ~ ., data = dataMatrixBac,
                   method = "pls", preProc= "center",
                   trControl = trainControl(method = "LOOCV"))

pls_biplot_Bac <- list("loadings" = loadings(plsFitBac$finalModel),
                       "scores" = scores(plsFitBac$finalModel))
class(pls_biplot_Bac$scores) <- "matrix"

pls_biplot_Bac$scores <- data.frame(sample_data(logBAC),
                                    pls_biplot_Bac$scores)

tax <- data.frame(tax_table(logBAC), stringsAsFactors = FALSE)
main_phyla <- sort(c("Tenericutes", "Actinobacteria",
                     "Bacteroidetes", "Firmicutes", "Proteobacteria"))

tax$Phylum[!(tax$Phylum %in% main_phyla)] <- "Z_Other"
tax$Phylum <- factor(tax$Phylum, levels = c(main_phyla, "Z_Other"))
class(pls_biplot_Bac$loadings) <- "matrix"
pls_biplot_Bac$loadings <- data.frame(tax, pls_biplot_Bac$loadings)

pdf("figures/figures_Hyena/Figure3c_pls_bi_bac.pdf")
ggplot() +
    geom_point(data = pls_biplot_Bac$loadings,
               aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Phylum),
               size = 2) +
    geom_point(data = pls_biplot_Bac$scores,
               aes(x = Comp.1, y = Comp.2, shape = age), size = 3) +
    scale_color_brewer(palette = "Set1") +
    scale_shape_discrete(solid=FALSE) +
    labs(x = "Axis1", y = "Axis2", col = "Phylum") +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    theme_bw()
##    theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
dev.off()



scores.vector <- (pls_biplot_Bac$loadings$Comp.1*-1) *
    pls_biplot_Bac$loadings$Comp.2 

pls.ordered <- pls_biplot_Bac$loadings[order(scores.vector), ]

q25 <- head(pls.ordered, n = (nrow(pls.ordered)*0.25))
q75 <- tail(pls.ordered, n = (nrow(pls.ordered)*0.25))

tar.Phylum <- "Tenericutes"

fisher.test(pls.ordered$Genus %in%
            pls.ordered$Genus[pls.ordered$Phylum%in%tar.Phylum],
            pls.ordered$Genus %in%
            q75$Genus)

tar.Phylum <- "Bacteroidetes"

fisher.test(pls.ordered$Genus %in%
            pls.ordered$Genus[pls.ordered$Phylum%in%tar.Phylum],
            pls.ordered$Genus %in%
            q75$Genus, alternative="greater")


tar.Phylum <- "Actinobacteria"

fisher.test(pls.ordered$Genus %in%
            pls.ordered$Genus[pls.ordered$Phylum%in%tar.Phylum],
            pls.ordered$Genus %in%
            q25$Genus)


## get unnormalized data to do the DESeq stuff

load(file="/SAN/Metabarcoding/phlyoSeq_Hy_raw.Rdata") ## -> PS.raw

PS.raw.genus <- tax_glom(PS.raw, "Genus")

## differences by phyla
library(DESeq2)
Bac.diagdds <- phyloseq_to_deseq2(subset_taxa(PS.raw.genus, Kingdom%in%"Bacteria"), ~ age)
Bac.diagdds <- DESeq(Bac.diagdds, test="LRT", fitType="parametric", reduced= ~ 1)

Bac.res <- results(Bac.diagdds, cooksCutoff = FALSE)
alpha <- 0.05
Bac.sigtab <- Bac.res[which(Bac.res$padj < alpha), ]
Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS.genus)[rownames(Bac.sigtab), ], "matrix"))
rownames(Bac.sigtab) <- NULL
Bac.sigtab

library(scales)

get.sigtab.plot <- function (sigtab){
    sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels=main_phyla)
    ## Genus order
    x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
    x <- sort(x, TRUE)
    sigtab$Genus <- factor(as.character(sigtab$Genus), levels=names(x))
    ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum, label=scientific(padj))) +
        geom_point(size=6) +
        geom_text(vjust=1.6, color="black")+
        theme_bw()+
        scale_color_brewer(palette = "Set1") +
        theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust=1))
}

plot.bacterial.diff <- get.sigtab.plot(Bac.sigtab)

devSVG("figures/figures_Hyena/Figure3d_single_diff.svg")
plot.bacterial.diff
dev.off()

pdf("figures/figures_Hyena/Figure3d_single_diff.pdf")
plot.bacterial.diff
dev.off()

######################### Eukaryotes ############################

logEuk <- transform_sample_counts(
    subset_taxa(PS.genus, Kingdom%in%"Eukaryota"&
                    !Phylum%in%clear.prey.phyla),
    function(x) log10(1+x))


dataMatrixEuk <- data.frame(age = sample_data(logEuk)$age,
                            rank = sample_data(logEuk)$rank,
                            otu_table(logEuk))

plsFitEuk <- train(rank ~ ., data = dataMatrixEuk,
                   method = "pls", preProc= "center",
                   trControl = trainControl(method = "LOOCV"))

pls_biplot_Euk <- list("loadings" = loadings(plsFitEuk$finalModel),
                       "scores" = scores(plsFitEuk$finalModel))
class(pls_biplot_Euk$scores) <- "matrix"

pls_biplot_Euk$scores <- data.frame(sample_data(logEuk),
                                    pls_biplot_Euk$scores)

tax <- data.frame(tax_table(logEuk), stringsAsFactors = FALSE)

main_phyla <- sort(c("Basidiomycota", "Ascomycota", "Platyhelminthes",
                     "Arthropoda", "Apicomplexa", "Ochrophyta",
                     "Blastocladiomycota", "Nematoda"))

tax$Phylum[!(tax$Phylum %in% main_phyla)] <- "Z_Other"

tax$Phylum <- factor(tax$Phylum, levels = c(main_phyla, "Z_Other"))

class(pls_biplot_Euk$loadings) <- "matrix"
pls_biplot_Euk$loadings <- merge(tax, pls_biplot_Euk$loadings, by=0)

colnames(pls_biplot_Euk$loadings) <-
    make.names(colnames(pls_biplot_Euk$loadings))

pdf("figures/figures_Hyena/Figure4c_pls_bi_euk.pdf", width=5, height=10)
ggplot() +
    geom_jitter(data = pls_biplot_Euk$loadings,
                aes(x="1", y = 25 * Comp.1 , col = Phylum), size = 2) +
    geom_jitter(data = pls_biplot_Euk$scores,
               aes(x="1",y = Comp.1, shape=rank), size=4) +
    scale_color_brewer(palette = "Set1") +
    labs(y = "Sinlge Axis 1", col = "Phylum", x = "") +
    scale_shape_discrete(solid=FALSE) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    ##    theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
    theme_bw()
dev.off()

pdf("/home/ele/Figure4c_labeled.pdf", width=18, height=4)
ggplot() +
    geom_jitter(data = pls_biplot_Euk$loadings,
                aes(y="1", x = 25 * Comp.1 , col = Phylum), size = 2) +
    geom_label(data = pls_biplot_Euk$scores,
               aes(y="1",x = Comp.1, label=Hyena.ID), size=2, position="jitter") +
    scale_color_brewer(palette = "Set1") +
    labs(y = "Random scatter", col = "Phylum", x = "") +
    scale_shape_discrete(solid=FALSE) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    ##    theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
    theme_bw()
dev.off()

write.csv(data[order(data$Comp.1), ], "/home/ele/PLS_differences_Figure4c.csv",
          row.names = FALSE)


pls.ordered <- pls_biplot_Euk$loadings[order(pls_biplot_Euk$loadings$Comp.1), ]

q25 <- head(pls.ordered, n = (nrow(pls.ordered)*0.25))
q75 <- tail(pls.ordered, n = (nrow(pls.ordered)*0.25))


tar.Phylum <- "Basidiomycota"

fisher.test(pls.ordered$Genus %in%
            pls.ordered$Genus[pls.ordered$Phylum%in%tar.Phylum],
            pls.ordered$Genus %in%
            q75$Genus)


tar.Phylum <- "Apicomplexa"

fisher.test(pls.ordered$Genus %in%
            pls.ordered$Genus[pls.ordered$Phylum%in%tar.Phylum],
            pls.ordered$Genus %in%
            q25$Genus)

fisher.test(pls.ordered$Genus %in%
            pls.ordered$Genus[pls.ordered$Phylum%in%tar.Phylum],
            pls.ordered$Genus %in%
            c(q25$Genus, q75$Genus))


## Eukaryote single genus differences
PS.female.euk <- subset_samples(
    subset_taxa(PS.raw.genus, Kingdom%in%"Eukaryota"))## , age%in%"Adult")

Euk.diagdds <- phyloseq_to_deseq2(PS.female.euk, ~ rank)
Euk.diagdds <- DESeq(Euk.diagdds, test="LRT", fitType="parametric", reduced= ~ 1)

Euk.res <- results(Euk.diagdds, cooksCutoff = FALSE)
alpha <- 0.1
Euk.sigtab <- Euk.res[which(Euk.res$padj < alpha), ]
Euk.sigtab <- cbind(as(Euk.sigtab, "data.frame"), as(tax_table(PS.raw.genus)[rownames(Euk.sigtab), ], "matrix"))
rownames(Euk.sigtab) <- NULL
Euk.sigtab

plot.euk.diff <- get.sigtab.plot(Euk.sigtab)

devSVG("figures/figures_Hyena/Figure4d_single_diff.svg")
plot.euk.diff
dev.off()


pdf("figures/figures_Hyena/Figure4d_single_diff.pdf")
plot.euk.diff
dev.off()


NSplsFitEuk.age <- train(age ~ ., data = dataMatrixEuk,
                         method = "pls", preProc= "center",
                         trControl = trainControl(method = "LOOCV"))
