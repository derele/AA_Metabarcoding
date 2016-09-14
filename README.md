# AA_Metabarcoding


[![DOI](https://zenodo.org/badge/67146407.svg)](https://zenodo.org/badge/latestdoi/67146407)

An overview of a metabarcoding pipeline for multiple marker data from
the Fluidigm Access Array. The present version is intended for a
overview of the preliminary processing of data. Code can be reviewed
in the /scripts and /R folders. The processing is not reproducible at
the present version as raw data files cannot be accessed (yet). 

## Preprocessing

### stratify the data into amplicons and samples
#### all samples in one file with sample identifier:
```shell
scripts/combine_and_lable_samples.pl *R1.fastq.gz > ALL_R1.fastq
scripts/combine_and_lable_samples.pl *R2.fastq.gz > ALL_R2.fastq
```

#### sort into amplicons
`scripts/sort_amplicons.pl`

#### run the usearch pipeline 
`parallel scripts/pipeLine.sh {} ::: /path/to/sorted/*.fastq`

## a little fix befor import of otu tables  table for R
`sed -i  's/#OTU ID/OTUID/ig' *.fastq.otu_table.txt`

#### concatenate all otus
`gawk '{if(/^>/){print $0"|"FILENAME} else{print $0}}' *.fastq.otus.fa > ALL_outs.fa`

#### blast all otus against a silva subset of nt
`blastn -query ALL_outs.fa -evalue 1e-20 -perc_identity 97 -db /SAN/db/blastdb/nt/nt -gilist /SAN/db/blastdb/silvaEukarya/SilvaEuk.gi  -outfmt 11 > ALL_outs_vs_SilvaEuk.asn1`

#### get the taxonomy based on these blasts
`scripts/blast2alltax_outfmt11.pl ALL_outs_vs_SilvaEuk.asn1 > ALL_outs.taxtable`

