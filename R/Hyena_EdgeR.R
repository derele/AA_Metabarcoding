library(edgeR)
library(anytime)

source("R/Hyena_first_rough.R")
source("R/Hyena_first_taxonomy.R")

Hyena.Cat <- read.csv("/home/ele/Dropbox/Animal_number_variable_codes_EH.csv")
Hyena.Cat <- as.data.frame(apply(Hyena.Cat, 2, function (x) gsub(" ", "", x)))

Hyena.Cat$V1 <- ifelse(Hyena.Cat$V1==1, "Male", "Female")
Hyena.Cat$V2 <- ifelse(Hyena.Cat$V2==2, "Cub", "Adult")
Hyena.Cat$V3 <- ifelse(Hyena.Cat$V3==1, "high", "low")
Hyena.Cat$V3[is.na(Hyena.Cat$V3)] <- "low" ## watchout not TRUE!!!

names(Hyena.Cat) <- c("animal", "ID", "sex", "age", "rank", "pack")
rownames(Hyena.Cat) <- Hyena.Cat$ID
Hyena.Cat$ID <- Hyena.Cat$animal <- NULL

## compare with Luis' Hyenas:
Lu.Cat <- read.csv("/home/ele/Dropbox/Hyena_data_Luis.csv",header=TRUE)

Lu.Cat$Collection.Date <- ifelse(nchar(Lu.Cat$Collection.Date)==5,
                                 paste0("200", Lu.Cat$Collection.Date),
                                 paste0("20", Lu.Cat$Collection.Date))
Lu.Cat$CollectionDate <- anydate(Lu.Cat$Collection.Date)

Lu.Cat <- merge(Lu.Cat, Hyena.Cat)

## the OTUs
clean.SOTUs <- SOTUs[, colnames(SOTUs)%in%rownames(Hyena.Cat)]
## get them into the same order
clean.SOTUs <- clean.SOTUs[, rownames(Hyena.Cat)]

Hyena.design <- model.matrix(~Hyena.Cat$sex +
                                 Hyena.Cat$age +
                                 Hyena.Cat$rank +
                                 Hyena.Cat$pack)

colnames(Hyena.design) <- gsub("design.frame\\$", "", colnames(Hyena.design))

get.my.models <- function(otus, design, cutoff=50){
    cutoff <- 50
    keep <- rowSums(otus)>cutoff
    DGEList <- DGEList(otus[keep,])
    DGEList <- calcNormFactors(DGEList )
    DGEList <- estimateDisp(DGEList, design = design, robust = TRUE)
    fit <- glmFit(DGEList, design)
    contrasts.V1 <- makeContrasts(V11vs2 = V12,
                                  V22vs3 = V23,
                                  V4IvsM = V4M,
                                  V4IvsP = V4P,
                                  V4MvsP = V4M-V4P, 
                                  levels=design)
    fitLRT.list <- lapply(colnames(contrasts.V1),
                          function (x) glmLRT(fit,
                                              contrast = contrasts.V1[,x]))
    names(fitLRT.list) <- colnames(contrasts.V1)

    top.list <- lapply(fitLRT.list, function (x){
        topTags(x, 1000000)[[1]]
    })
    names(top.list) <- colnames(contrasts.V1)
    OTU.list <- lapply(top.list, function(x) {
        rownames(x[x$FDR<0.05,]) ####### ## CHANGE FDR!!!!!!!!!!
    })
    names(OTU.list) <- colnames(contrasts.V1)
    return(list(OTU.list, top.list, fitLRT.list))
}

Smodels <- get.my.models(SOTUs, Hyena.design)
SpeciesSmodels <- get.my.models(SpecAbu, Hyena.design)
ClassSmodels <- get.my.models(ClassAbu, Hyena.design)
GenusSmodels <- get.my.models(GenusAbu, Hyena.design)

########## Not a good plot

pdf("figures/mdsV4Color.pdf")
plotMDS(DGEList,
        labels=paste(design.frame$V1, design.frame$V2, design.frame$V3, design.frame$V4, sep=""), 
        col=ifelse(design.frame$V4%in%"I", "red",
            ifelse(design.frame$V4%in%"P", "blue", "green")))
dev.off()

pdf("figures/mdsV1Color.pdf")
plotMDS(DGEList,
        labels=paste(design.frame$V1, design.frame$V2, design.frame$V3, design.frame$V4, sep=""), 
        col=ifelse(design.frame$V1%in%"1", "red", "blue"))
dev.off()

pdf("figures/mdsV2Color.pdf")
plotMDS(DGEList,
        labels=paste(design.frame$V1, design.frame$V2, design.frame$V3, design.frame$V4, sep=""), 
        col=ifelse(design.frame$V2%in%"2", "red", "blue"))
dev.off(which=)

pdf("figures/mdsV3Color.pdf")
plotMDS(DGEList,
        labels=paste(design.frame$V1, design.frame$V2, design.frame$V3, design.frame$V4, sep=""), 
        col=ifelse(design.frame$V3%in%"2", "red", "blue"))
dev.off()

##################################

## Look at
TAX[TAX$query%in%OTU.list[["V11vs2"]],]

## only bacteria
###

TAX[TAX$query%in%OTU.list[["V22vs3"]],]
V2 <- c("Ascogregarina", "Dipylidium")
##

TAX[TAX$query%in%OTU.list[["V31vs2"]],]

V3 <- c("Linguatula", "Cochliomyia", "Cylicostephanus", "Protostrongylus")

#### #### 
TAX[TAX$query%in%OTU.list[["V4IvsP"]],]

V4IvsP <- c("Cochliomyia", "Lucilia")

## "Linguatula" was differentially abundant when analysing the
## technical pseudoreplicates

##
TAX[TAX$query%in%OTU.list[["V4MvsP"]],]

V4IvsP <- c("Cochliomyia", "Lucilia")

annotation.columns <- design.frame[,2:5]

special.pheatmap <- function(OTUs, OTU.list, test="V11vs2", level="species"){
    data <- log(SOTUs[OTU.list[[test]], ]+0.1)
    taxdata <- TAX.best[TAX.best$query%in%rownames(data), ]
    rownames(taxdata) <- taxdata$query
    labels.row <- as.character(taxdata[rownames(data), level])
    pheatmap(data, 
             annotation_col=annotation.columns,
             labels_row=labels.row)
}

pdf("figures/testedSpecAbu_pheatmatV4.pdf", height=14, width=14)
pheatmap(log(SpecAbu[unique(c(OTU.list[["V4MvsP"]], OTU.list[["V4IvsP"]])), ]+0.1),
         annotation_col=annotation.columns)
dev.off()

pdf("figures/tested_pheatmatV2.pdf", height=14, width=14)
special.pheatmap("V22vs3", "species")
dev.off()

pdf("figures/tested_pheatmatV3.pdf", height=14, width=14)
pheatmap(log(SOTUs[OTU.list[["V31vs2"]], ]+0.1),
         annotation_col=annotation.columns)
dev.off()

pdf("figures/tested_pheatmatV4IP.pdf", height=14, width=14)
pheatmap(log(SOTUs[OTU.list[["V4IvsP"]], ]+0.1),
         annotation_col=annotation.columns)
dev.off()

pdf("figures/tested_pheatmatV4MP.pdf", height=14, width=14)
pheatmap(log(SOTUs[OTU.list[["V4MvsP"]], ]+0.1),
         annotation_col=annotation.columns)
dev.off()


## get the same order
SpecAbu <- SpecAbu[, rownames(design.frame)]

## the same for the species. 
cutoff <- 23
keep <- rowSums(SpecAbu)>cutoff

DGEList <- DGEList(SpecAbu[keep,])
DGEList <- calcNormFactors(DGEList )

DGEList <- estimateDisp(DGEList, design = Hyena.design, robust = TRUE)

fit <- glmFit(DGEList, Hyena.design)

contrasts.V1 <- makeContrasts(V11vs2 = V12,
                              V22vs3 = V23,
                              V4IvsM = V4M,
                              V4IvsP = V4P,
                              V4MvsP = V4M-V4P, 
                              levels=Hyena.design)

fitLRT.list <- lapply(colnames(contrasts.V1),
                      function (x) glmLRT(fit,
                                          contrast = contrasts.V1[,x]))

names(fitLRT.list) <- colnames(contrasts.V1)

top.list <- lapply(fitLRT.list, function (x){
    topTags(x, 1000000)[[1]]
})


names(top.list) <- colnames(contrasts.V1)

OTU.list <- lapply(top.list, function(x) {
    rownames(x[x$FDR<0.01,]) ####### ## CHANGE FDR!!!!!!!!!!
})

names(OTU.list) <- colnames(contrasts.V1)


########## Not a good plot

pdf("figures/mdsV4Color.pdf")
plotMDS(DGEList, labels=paste(design.frame$V1, design.frame$V2, design.frame$V3, design.frame$V4, sep=""), 
        col=ifelse(design.frame$V4%in%"I", "red",
            ifelse(design.frame$V4%in%"P", "blue", "green")))
dev.off()

pdf("figures/mdsV1Color.pdf")
plotMDS(DGEList, labels=paste(design.frame$V1, design.frame$V2, design.frame$V3, design.frame$V4, sep=""), 
        col=ifelse(design.frame$V1%in%"1", "red", "blue"))
dev.off()

pdf("figures/mdsV2Color.pdf")
plotMDS(DGEList, labels=paste(design.frame$V1, design.frame$V2, design.frame$V3, design.frame$V4, sep=""), 
        col=ifelse(design.frame$V2%in%"2", "red", "blue"))
dev.off()

pdf("figures/mdsV3Color.pdf")
plotMDS(DGEList, labels=paste(design.frame$V1, design.frame$V2, design.frame$V3, design.frame$V4, sep=""), 
        col=ifelse(design.frame$V3%in%"2", "red", "blue"))
dev.off()


pdf("figures/Dipylidium_pheatmatV2.pdf", height=14, width=14)
pheatmap(log(SpecAbu[grep("Dip", rownames(SpecAbu)), ]+1),
                     annotation_col=annotation.columns)
dev.off()
