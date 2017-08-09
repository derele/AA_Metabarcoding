# download the SILVA taxa mapping
# ## wget http://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_123.txt

# cut -f3 tax_slv_ssu_123.txt | sort | uniq -c
# what a mess

# reordered by taxonomic level
# 	8275 genus
# 	  26 subfamily
# 	1412 family
# 	  10 superfamily
# 	  21 suborder
# 	1138 order
# 	  16 superorder
# 	   6 infraclass
# 	  58 subclass
# 	 541 class
# 	   1 superclass
# 	   1 infraphylum
# 	  33 subphylum
# 	 239 phylum
# 	   5 superphylum
# 	   1 infrakingdom
# 	   4 subkingdom
# 	  15 kingdom
# 	   2 superkingdom
# 	   4 major_clade
# 	   3 domain
#      1       # out of all the labels, they leave Escherichia blank??

# start from scratch with the silva.nr_v123.align headers
# grep '>' silva.nr_v123.align | cut -f1,3 | cut -f2 -d'>' > silva.nr_v123.full

#!/usr/bin/R

map.in <- read.table("/SAN/db/RDP/Silva_123/tax_slv_ssu_123.txt",header=F,sep="\t",stringsAsFactors=F)
map.in <- map.in[,c(1,3)]
colnames(map.in) <- c("taxlabel","taxlevel")

#fix Escherichia nonsense
map.in$taxlevel[which(map.in$taxlabel=="Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia;")] = "genus"

taxlevels <- c("root","domain","major_clade","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus")
taxabb <- c("ro","do","mc","pk","ki","bk","ik","pp","ph","bp","ip","pc","cl","bc","ic","po","or","bo","pf","fa","bf","ge")
tax.mat <- matrix(data="",nrow=nrow(map.in),ncol=length(taxlevels))
tax.mat[,1] <- "root"
colnames(tax.mat) <- taxlevels

#change this if you want just the default 6
## outlevels <- taxlevels
## default 6
outlevels <- c("domain","phylum","class","order","family","genus")

#put all the taxa into a matrix

for (i in 1:nrow(map.in)) {
	taxname <- unlist(strsplit(as.character(map.in[i,1]), split=';'))
	print(taxname);
	while ( length(taxname) > 0) {
		#regex to look for exact match
		tax.exp <- paste(paste(taxname,collapse=";"),";",sep="")
		tax.match <- match(tax.exp,map.in$taxlabel)
		tax.mat[i,map.in[tax.match,2]] <- tail(taxname,1)
		taxname <- head(taxname,-1)
	}
}


for (i in 1:nrow(tax.mat)) {
	#this fills in the empty gaps by using the closest higher taxonomic level appended with an abbreviation for the current taxonomic level
	#if you don't want this behavior, cut it out
	for (j in 1:ncol(tax.mat)) {
		if(tax.mat[i,j] < 0) { tax.mat[i,j] <- paste(tmptax,taxabb[j],sep="_")}
		else { tmptax <- tax.mat[i,j]}
	}
	#this maps the new name to the input taxonomic levels
	map.in[i,"taxout"] <- paste(paste(tax.mat[i,outlevels],collapse=";"),";",sep="")
}

# replace spaces with underscores
map.in$taxout <- gsub(" ","_",map.in$taxout)

# bring in the old taxonomic levels from SILVA and remap them using the new levels

tax.in <- read.table("/SAN/db/RDP/Silva_123/silva.nr_v123.full",header=F,stringsAsFactors=F,sep="\t")
colnames(tax.in) <- c("taxid","taxlabel")

tax.in$taxlabel <- gsub("[[:space:]]+;", ";", tax.in$taxlabel) #fix extra space in "Eukaryota;Opisthokonta;Nucletmycea;Fungi;Dikarya;Ascomycota;Pezizomycotina;Dothideomycetes;Pleosporales;Phaeosphaeriaceae;Parastagonospora ;"

tax.in$id <- 1:nrow(tax.in)

tax.write <- merge(tax.in,map.in,all.x=T,sort=F)
tax.write <- tax.write[order(tax.write$id),]

## write.table(tax.write[,c("taxid","taxout")],
##             file="/SAN/db/RDP/Silva_123/silva.nr_v123.full.tax",
##             sep="\t", row.names=F, quote=F, col.names=F)

sids <- strsplit(tax.write$taxid, "\\.")
tax.write$seqacc <- unlist(lapply(sids, "[", 1))

library(Biostrings)

FAS <- Biostrings::readRNAStringSet("/SAN/db/RDP/Silva_123/SILVA_123_SSURef_Nr99_tax_silva.fasta")

nnames <- strsplit(names(FAS), " ")
nnames <- unlist(lapply(nnames, "[", 1))
nnames <- strsplit(nnames, "\\.")
nnames <- unlist(lapply(nnames, "[", 1))

FAS.taxed <- FAS[nnames%in%tax.write$seqacc]
nnames.taxed <- nnames[nnames%in%tax.write$seqacc]

FAS.taxed <- FAS.taxed[match(tax.write$seqacc, nnames.taxed)]
names(FAS.taxed) <- tax.write$taxout
FAS.taxed <- DNAStringSet(FAS.taxed)

## This would be a very first quite "incomplete" database.
## Biostrings::writeXStringSet(FAS.taxed, "/SAN/db/RDP/Silva_123/SILVA_123_dada2.fasta", format="fasta")


## now need to run blast and the blast2alltax scirpt



FUZZ.tax <- read.csv("/SAN/Metabarcoding/AA_combi/all_dada_vs_nt.taxtable",
                     as.is=TRUE)
## first see that we have them more then once
FUZZ.tax <- FUZZ.tax[duplicated(FUZZ.tax$subject), ]
## then exclude duplicates
FUZZ.tax <- FUZZ.tax[!duplicated(FUZZ.tax$subject), ]
## then exclude undef at genus level
FUZZ.tax <- FUZZ.tax[!grepl("undef", FUZZ.tax$genus), ]

write.table(FUZZ.tax$subject,
            "/SAN/Metabarcoding/AA_combi/selected_dada_nr_hits.acc",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

## now run blastcmd to get the fasta file

## blastdbcmd -db /SAN/db/blastdb/nt/nt -entry_batch selected_dada_nr_hits.acc -outfmt '%a|%s' > dada_nr_hits.fasta

## plus a little bit of formatting magic
## tr '|' '\n' < dada_nr_hits.fasta > tmp; mv tmp dada_nr_hits.fasta
## awk '{if(NR%2){print ">"$0} else {print}}' dada_nr_hits.fasta > tmp; mv tmp dada_nr_hits.fasta

FUZZ <- Biostrings::readDNAStringSet("/SAN/Metabarcoding/AA_combi/dada_nr_hits.fasta")

names(FUZZ) <- gsub("\\.\\d+$", "", names(FUZZ))

newname <- sapply(names(FUZZ), function(acc){
    focus <- FUZZ.tax[FUZZ.tax$subject%in%acc,]
    ## names corresponding to 
    ## outlevels <- c("domain","phylum","class","order","family","genus")
    paste(focus[,c("superkingdom", "phylum", "class", 
                   "order", "family", "genus")], collapse=";")
})

names(FUZZ) <- newname

## now get rid of whole genomes and other really long sequences for
## efficiency of RDP
FUZZ <- FUZZ[(width(FUZZ)<5000 & width(FUZZ)>1400)]

## properly named
FUZZ <- FUZZ[(nchar(names(FUZZ))>10)]
FUZZ <- FUZZ[!grepl("character\\(0\\)", names(FUZZ))]

allFUZZ <- DNAStringSet(c(FAS.taxed, FUZZ))

## Biostrings::writeXStringSet(FUZZ,
##                            "/SAN/Metabarcoding/AA_combi/dada2_WHexp.fasta",
##                            format="fasta")


## Biostrings::writeXStringSet(allFUZZ,
##                             "/SAN/db/RDP/Silva_123/SILVA_123_dada2_WHexp.fasta",
##                             format="fasta")
