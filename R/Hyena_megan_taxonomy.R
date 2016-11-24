
foo <- 



library(stringr)


MeDat <- readLines("/SAN/Metabarcoding/Hyena/second/sorted_amps/usearch/ALL_outs_vs_NT.megantax")
foo <- lapply(MeDat, function (line){
    all.col <- strsplit(line, ";")[[1]]
    all.col <- lapply(all.col, function (x) gsub(" *", "", x))
    all.col <- all.col[unlist(lapply(all.col, function (x) nchar(x)>0))]
    regex <- paste("d__(?<domain>.*)",
                   "k__(?<kingdom>.*)",
                   "p__(.*)|c__(.*)|o__(.*)|f__(.*)|g__(.*)|s__(.*)", sep="|")
    gregexpr(regex, all.col, perl=TRUE)
    ##     col.nice <- gsub(,
    ##                     "\\1\\2\\3\\4\\5\\6\\7\\8",
    ##                     all.col)


    ## addding the NAs for missing
    col.nice <- col.nice[1:17]
})
bar <- as.data.frame(do.call(rbind,foo))
names(bar) <- c("OTU", "domain", "s.domain", "kingdom", "s.kingdom",
                "phylum", "s.phylum", "class", "s.class",
                "order", "s.order", "family", "s.family",
                "genus", "s.genus", "species", "s.species")
