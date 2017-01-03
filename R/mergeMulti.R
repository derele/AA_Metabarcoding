mergeMulti <- function(MA, ...){
    MA@mergers <- lapply(seq_along(MultiDadaF), function (i){
        mergePairs(MA@dadaF[[i]], MA@derepF[[i]],
                   MA@dadaR[[i]], MA@derepR[[i]],
                   ...)
    })
    return(MA)
}
