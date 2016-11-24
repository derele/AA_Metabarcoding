source("R/Hyena_first_rough.R")

library(dada2); packageVersion("dada2")
## library(ShortRead); packageVersion("ShortRead")

read.ESeq <- function (file){
  content <- readLines(file)
  header.lines <- grep( "^>", content)
  start.lines <- header.lines+1
  end.lines <- c(header.lines[-1]-1, length(content))
  sq <- sapply(1:length(start.lines), function (i) {
    list(content[start.lines[i]:end.lines[i]])})
  names(sq) <- substr(content[header.lines],2, nchar(content[header.lines]))
  sq <- unlist(lapply(sq, paste, collapse=""))
  sq <- sq[nchar(sq)>0]
  class(sq) <- "ESeq"
  return(sq)
}

OTU.seq.file <- "/SAN/Metabarcoding/Hyena/second/sorted_amps/usearch/ALL_outs.fa"
OTU.seq.plain <- read.ESeq(OTU.seq.file)

OTU.seq.bak <- OTU.seq.plain[grep("ADM|ACM|Klin", names(OTU.seq.plain), value=TRUE)]

set.seed(100) # Initialize random number generator for reproducibility
taxa <- assignTaxonomy(OTU.seq.bak, "/SAN/db/RDP/silva_nr_v123_train_set.fa.gz")

unname(taxa)



library(plyr)

df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

p <- ggplot(df, aes(Axis.1, Axis.2, color=rank, shape=age)) + theme_bw()
p <- p + geom_point(size=3, alpha=0.5)
p <- p + facet_wrap(~distance, scales="free")
p <- p + ggtitle("MDS on various distance metrics for bacterial metabarcoding dataset")
p

