library(Biostrings)
library(ShortRead)
## library(dada2)


Ffq.file <- "/SAN/Metabarcoding/Hyena/second/fastq_raw/I247_S121_L001_R1_001.fastq.gz"

Rfq.file <- "/SAN/Metabarcoding/Hyena/second/fastq_raw/I247_S121_L001_R2_001.fastq.gz"

file.ptable <- "/SAN/Metabarcoding/Hyena/second/sort_amplicons_input.csv"

ptable <- read.csv(file.ptable, sep=" ", header=TRUE)

sort.into.amplicons <- function(Ffls, Rfls, primerF, primerR,
                                nameF, nameR, n=1e6, ...){
    nPF <- length(primerF) 
    nPR <- length(primerR)
    if (nPF!=nPR) {
        stop (paste("different number of forward and reverse primers:",
                    length(primerF), length(primerR),"provide primer pairs!"))
    }    
    data <- data.frame(primerF, primerR, nameF, nameR,
                       matrix(0, nrow=nPF, ncol=length(Ffls)))
    names(data)[5:ncol(data)] <- basename(Ffls)
    for(i in seq_along(Ffls)) {
        tmpbaseF <- paste0(tempfile(), basename(Ffls[[i]]))
        tmpbaseR <- paste0(tempfile(), basename(Rfls[[i]]))
        out <- get.primer.map(primerF, primerR, nameF, nameR)
        Fprimer <- out$pF
        Rprimer <- out$pR
        f <- FastqStreamer(Ffls[[i]], n = n)
        r <- FastqStreamer(Rfls[[i]], n = n)
        ## request forward and reverse file simultaneously
        while(length(suppressWarnings(Ffq <- yield(f))) &&
              length(suppressWarnings(Rfq <- yield(r)))){
                  fM <- lapply(Fprimer, function(x){
                      as.vector(isMatchingStartingAt(x, sread(Ffq),
                                                     fixed=FALSE))
                  })
                  rM <- lapply(Rprimer, function(x){
                      as.vector(isMatchingStartingAt(x, sread(Rfq),
                                                     fixed=FALSE))
                  })
                  matches <- apply(out$map, 1, function(x){
                      map.primerF <- as.numeric(x["map.pF"])
                      map.primerR <- as.numeric(x["map.pR"])
                      select <- fM[[map.primerF]] & rM[[map.primerR]]
                      tmppathF <- paste0(tmpbaseF, x[["nF"]], ":", x[["nR"]], ".fastq")
                      tmppathR <- paste0(tmpbaseR, x[["nF"]], ":", x[["nR"]], ".fastq")
                      writeFastq(Ffq[select], file=tmppathF, mode="a")
                      writeFastq(Rfq[select], file=tmppathR, mode="a")
                      number.matches <- length(select[select==TRUE])
                      return(number.matches)
                  })
                  ## need to add over the while loop because of fastq streaming 
                  data[, basename(Ffls[[i]])] <- data[, basename(Ffls[[i]])] + matches
              }
    }
    return(data)
}

foo <- sort.into.amplicons(Ffq.file, Rfq.file,
                           as.character(ptable$TS.SequenceF),
                           as.character(ptable$TS.SequenceR),
                           as.character(ptable$corrected.NameF),
                           as.character(ptable$corrected.NameR))

get.primer.map <- function(pF, pR, nF, nR){
    ## test for equal length of primer F and R    
    l.pF <- nchar(pF)
    l.pR <- nchar(pR)
    u.pF <- sort(unique(pF))
    map.pF <-  as.numeric(factor(pF))
    Fprimer <- DNAStringSet(u.pF)
    u.pR <- sort(unique(pR))
    map.pR <-  as.numeric(factor(pR))
    Rprimer <- DNAStringSet(u.pR)
    map <- cbind(map.pF, map.pR, l.pF, l.pR, nF=nF, nR=nR)
    return(list(map=map, pF=Fprimer, pR=Rprimer))
}



primer.to.regex <- function(primer){
    if (!is.vector(primer) || !is.character(primer)){
        stop("please provide primer sequence as character vector")
    }
    translator <- list(A = "A",
                       C = "C",
                       G = "G",
                       T = "T",
                       W = "(A|T)",
                       S = "(C|G)",
                       K = "(T|G)",
                       M = "(A|C)",
                       Y = "(C|T)",
                       R = "(A|G)",
                       V = "(A|C|G)",
                       D = "(A|T|G)",
                       B = "(T|G|C)",
                       H = "(A|T|C)",
                       N = "(A|T|C|G)")
    primer.l<- strsplit(primer, "")[[1]]
    if (any(!primer.l%in%names(translator))){
        stop(paste0("found non IUPAC symbol in primer sequence: ", primer))
    }
    regex.l <- lapply(primer.l, function (x) translator[[x]])
    paste(regex.l, collapse="")
}

