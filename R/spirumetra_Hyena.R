Spiru.ps.l <- lapply(FAps.l, function(ps){
    tryCatch(
        subset_taxa(ps, Order %in% "Spirurida"),
        error = function(e) NULL)


Taenia.ps.l <- lapply(FAps.l, function(ps){
    tryCatch(
        subset_taxa(ps, Genus %in% "Taenia"),
        error = function(e) NULL)
})

Spirom.ps.l <- lapply(FAps.l, function(ps){
    tryCatch(
        subset_taxa(ps, Genus %in% "Spirometra"),
        error = function(e) NULL)
})

lapply(Spirom.ps.l, function (x){
    tryCatch(
        table(rowSums(otu_table(x))>0),
        error = function(e) NULL)
})

rowSums(otu_table(Spirom.ps.l[[47]])>0)
    
    
