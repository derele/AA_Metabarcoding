# AA_Metabarcoding

## How to process

## stratify the data into amplicons and samples
## all samples in one file with sample identifier:

## sort into amplicons


## 
parallel path/to/AA_Metabarcoding/scripts/pipeLine.sh {} ::: /path/to/*.fastq

gawk '{if(/^>/){print $0"|"FILENAME} else{print $0}}' *.fastq.otus.fa > ALL_outs.fa

blastn -query ALL_outs.fa -evalue 1e-20 -perc_identity 97 -db /SAN/db/blastdb/silvaEukarya/SilvaEuk.fast -outfmt 11 > ALL_outs_vs_SilvaEuk.asn1
