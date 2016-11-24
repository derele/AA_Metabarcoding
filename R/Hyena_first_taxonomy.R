source("R/Hyena_first_rough.R")

read.taxtable <- function (path) {
    TAX.raw <- read.csv(path, sep = ",")
    TAX.raw$query <- gsub(".fastq.otus.fa", "", TAX.raw$query)
    TAX.raw$amplicon <- gsub("OTU\\d+\\|", "", TAX.raw$query)
    return(TAX.raw)
}

TAX.raw <-
    read.taxtable("/SAN/Metabarcoding/Hyena/second/sorted_amps/usearch/ALL_outs_nt.taxtable")

TAX.raw.Luis <-
    read.taxtable("/SAN/Metabarcoding/Hyena/first_good/sorted_amplicons/usearch/ALL_outs_NT.taxtable")


get.TAX.from.raw <- function(tax.raw, level="class"){
    ## unique betst bitscore for every query only if this is resolved at a
    ## certain level 
    T.l <- by(tax.raw, tax.raw$query, function (x, level=level) {
        b.bit <- max(x$bitscore)
        all.best <- x[x$bitscore==b.bit, ]
        ## A little last common ancestor play here... BUT wait ... 
        ## a shortcut throwing out OTUs that don't agree at least on the
        ## class level and allowing only the best hit afterwards
        u.family <- unique(all.best[, "class"])
        if(length(u.family)==1){
            ## get the best entry only
            return(all.best[1,])
        }
    })
    TAX <- do.call(rbind, T.l)
    rownames(TAX) <- NULL
    return(TAX)
}

TAX <- get.TAX.from.raw(TAX.raw, level="class")

TAX.Luis <- get.TAX.from.raw(TAX.raw.Luis, level="class")

## Only consider Euks now
TAX.Euk <- TAX[TAX$superkingdom%in%"Eukaryota", ]
TAX.Euk.Luis <- TAX.Luis[TAX.Luis$superkingdom%in%"Eukaryota", ]

## remove some really weird stuff FIND later out where the errors
## are!! Database errors...  should be fixed in database at some
## point... but now as a shortkut here
excl.phyla <- c("Cnidaria", "Porifera",
                "Bacillariophyta",  ## maybe okay?
                "Eustigmatophyceae") ## maybe okay?

TAX <- TAX[!TAX$phylum%in%excl.phyla ,]
TAX.Luis <- TAX.Luis[!TAX.Luis$phylum%in%excl.phyla ,]

### Summarizing by some taxon 
sumAmpByTax <- function(tax, OTUs, level,
                        exclude.sample.regex="H2O|Argave|Wolf"){
    mer <- merge(tax, OTUs, by.x="query", by.y=0)    
    mergedT.l<- by(mer[, (ncol(TAX)+1):ncol(mer)], mer[, level], colSums)
    mergedT <- do.call(rbind, mergedT.l)
    mergedT <- mergedT[order(rowSums(mergedT), decreasing=TRUE), ]
    ## remove empty and undef values
    mergedT <- mergedT[!rownames(mergedT)%in%c("", "undef"), ]
    ## remove some sample if wanted
    mergedT <- mergedT[, !grepl(exclude.sample.regex, colnames(mergedT))]
    return(mergedT)
}

ClassAbu.rep <- sumAmpByTax(TAX.Euk, OTUs, "class")
SpecAbu.rep  <- sumAmpByTax(TAX.Euk, OTUs, "species")

ClassAbu <- sumAmpByTax(TAX.Euk, SOTUs, "class")
SpecAbu  <- sumAmpByTax(TAX.Euk, SOTUs, "species")
GenusAbu  <- sumAmpByTax(TAX.Euk, SOTUs, "genus")

ClassAbu.Luis <- sumAmpByTax(TAX.Luis, SOTUs.Luis, "class")
GenusAbu.Luis  <- sumAmpByTax(TAX.Luis, SOTUs.Luis, "genus")
SpecAbu.Luis  <- sumAmpByTax(TAX.Luis, SOTUs.Luis, "species")

write.csv(SpecAbu.Luis, "~/Dropbox/Hyena_Species_abundance.csv")
write.csv(GenuAbu.Luis, "~/Dropbox/Hyena_Genus_abundance.csv")
write.csv(ClassAbu.Luis, "~/Dropbox/Hyena_Class_abundance.csv")

devSVG("figures/Hyena_class_sumNT_heat.svg", width=12, height=12)
pheatmap(log10(ClassAbu+1),
         show_rownames=TRUE,
         show_colnames=FALSE,
         treeheight_row=0,
         treeheight_col=0)
dev.off()

devSVG("figures/Hyena_species_sumNT_heat.svg", width=12, height=12)
pheatmap(log10(SpecAbu+1),
         show_rownames=FALSE,
         show_colnames=TRUE,
         treeheight_row=0,
         treeheight_col=0)
dev.off()

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

write.table(baz, "~/Dropbox/Hyena_2ndDataset_Classes_table_20161012.csv", sep=",")


write.csv(OTUs, "~/Dropbox/Hyena_2ndDataset_CompleteOTU_table_20161012.csv")

write.table(foo, "~/Dropbox/Hyena_2ndDataset_TaxaBlast_merge_20161012.csv", sep=",", row.names = FALSE)

