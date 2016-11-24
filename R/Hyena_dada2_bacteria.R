library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")


path <- "/SAN/Metabarcoding/Hyena/second/sorted_amps/"

ampF <- c("ADM330_Klin0785_CR_R1.fastq", "ACM_008_Klin0341_CR_R1.fastq")
ampR <- c("ADM330_Klin0785_CR_R2.fastq", "ACM_008_Klin0341_CR_R2.fastq")

files <- c(file.path(path, ampF),
           file.path(path, ampR))

fastqs <- files[grepl(".fastq$", files)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # reverse

## plotQualityProfile(fnFs[[1]])
## plotQualityProfile(fnRs[[1]])
## uhhh ... got to cut at 200 latest


filt.path <- file.path(path, "filteredDada")
if(!file_test("-d", filt.path)) dir.create(filt.path)
filtFs <- file.path(filt.path, paste0(ampF, "_F_filt.fastq.gz"))
filtRs <- file.path(filt.path, paste0(ampR, "_R_filt.fastq.gz"))


for(i in seq_along(fnFs)) {
    fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                      trimLeft=c(10, 10), truncLen=c(200,200),
                      maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                      compress=TRUE, verbose=TRUE)
}


derepMultiSampleFastq <- function(fls,
                                  headerCaptureRegex="^(?:.*sample=)?(\\w*+).*",
                                  n=1e6, ...){
    samplefiles <- character()
    for(file in fls) {
        tmpbase <- paste0(tempfile(), basename(file))
        f <- FastqStreamer(file, n = n)
        while( length(suppressWarnings(fq <- yield(f))) ){
            samples <- sub(headerCaptureRegex, "\\1", id(fq), perl=TRUE)
            ## get the uniqes samples in sequential order
            uni.samples <- sort(unique(samples))
            ## factor to map unique values to sequential integers
            map <- as.numeric(factor(samples))
            for (i in seq_along(uni.samples)){                
                tmppath <- paste0(tmpbase, uni.samples[[i]], ".fastq")
                ## lazy append for further stream   
                writeFastq(fq[map%in%i], file=tmppath, mode="a")
                samplefiles <- append(samplefiles, tmppath)
            }
        }
    }
    derepFastq(samplefiles, ...)
}

derepFs <- derepMultiSampleFastq(filtFs,
                                 "RECORD:(.*?)_L001.*", verbose=TRUE) 

derepRs <- derepMultiSampleFastq(filtRs,
                                 "RECORD:(.*?)_L001.*", verbose=TRUE)


sample.names <- grep("R1", list.files(tempdir()), value=TRUE)
                            
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=NULL, selfConsist = TRUE)
dadaRs <- dada(derepRs, err=NULL, selfConsist = TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE,
                      returnRejects=TRUE, justConcatenate=TRUE)

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

taxa <- assignTaxonomy(seqtab.nochim, "/SAN/db/RDP/rdp_train_set_14.fa.gz")
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

