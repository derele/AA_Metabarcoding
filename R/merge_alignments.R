library(Biostrings)

if(!exists("aln", mode="matrix")){
    source("R/4_entropy.R")
}


load(file="/SAN/Metabarcoding/align.Rdata") ## -> alignments

## reformat alignment from Silva  
RNAref <- RNAStringSet(DNAStringSet(apply(aln, 1, paste, collapse="")))

dbConn <- dbConnect(SQLite(), ":memory:")

Seqs2DB(RNAref, "RNAStringSet", dbConn, "Ref", tblName="Ref")

for(i in 1:length(alignments){
    Seqs2DB(alignaments[[i]], "RNAAStringSet", dbConn, "AA1", tblName="AA1")
}

Seqs2DB(alignaments[[1]], "RNAAStringSet", dbConn, "amp1", tblName="amp1")

## read it as DNA as it has Ts
moreDB <- readDNAStringSet("/SAN/Metabarcoding/Hyena/second/dada_nr_hits.fasta")
## align as RNA

Seqs2DB(RNAStringSet(moreDB), "XStringSet", dbConn, "moreDB", tblName="moreDB")



## AlignDB(dbConn, tblName=c("moreDB", "moreDB"), add2tbl="moreDB")


Seqs2DB(moreDB, "RNAStringSet", dbConn, "oldDB", tblName="oldDB")


AlignDB(dbConn, tblName=c("Ref", "oldDB"), add2tbl="AA")


dbDisconnect(dbConn)

