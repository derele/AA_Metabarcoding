#!/bin/bash

## filter and tirm in one step ## harsh
pfix=$(basename $1)

usearch -fastq_filter $1 -fastq_trunclen 250 -fastq_maxee 0.5 -fastq_truncqual 15 -fastaout $pfix.filtered.fasta 

## get a "sample identifier" to start right after the >
sed -i 's/>RECORD:/>/ig' $pfix.filtered.fasta
sed -i 's/-/_/ig' $pfix.filtered.fasta
sed -i 's/_L001_/;/ig' $pfix.filtered.fasta

## dereplicate
usearch -derep_fulllength $pfix.filtered.fasta -fastaout $pfix.uniques.fasta -sizeout 

## cluster
usearch -cluster_otus $pfix.uniques.fasta -otus $pfix.otus.fa -uparseout $pfix.out.up -relabel OTU -minsize 2

usearch -usearch_global $pfix.filtered.fasta -db $pfix.otus.fa -strand plus -id 0.97 -otutabout $pfix.otu_table.txt -biomout $pfix.otu_table.json

## the OTU sequences have tax=xxx; annotations, these will be included
## as an extra column in the tabbed file or as taxonomy metadata in
## the BIOM file. These annotations can be generated using the
## -fastaout option of the utax command, e.g.:

## usearch -utax otus.fa -db 16s.udb -strand both -fastaout otus_tax.fa
