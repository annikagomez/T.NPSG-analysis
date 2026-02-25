library(DESeq2)

##Set WD to source file location
setwd("INSERTFILEPATH")

#read in raw counts of t. erythraeum transcripts
rawCounts <- read.csv("metaT_files/tery_raw_gene_counts.csv",header=TRUE)
head(rawCounts)
geneID <- rawCounts$gene
sampleIndex <- grepl("SS0\\d+",colnames(rawCounts))
rawCounts <- round(as.matrix(rawCounts[,sampleIndex]))
rownames(rawCounts) <- geneID
head(rawCounts)

#metadata file indicates which samples are from the slick and which are from colonies
metaData <- read.csv("metaT_files/slick_metadata.csv",header=TRUE)
head(metaData)
rownames(metaData) <- metaData$id
metaData$id <- factor(metaData$id)
head(metaData)

rawCounts <- rawCounts[,unique(rownames(metaData))]
all(colnames(rawCounts) == rownames(metaData))
metaData$source <- factor(metaData$source, levels=c("slick", "nonslick"))

dds <- DESeqDataSetFromMatrix(countData=rawCounts, 
                              colData=metaData, 
                              design=~source)
dds <- DESeq(dds)

deseq2Results <- results(dds,contrast=c("source","slick",'nonslick'))
summary(deseq2Results)


plotMA(deseq2Results)


library(scales) # needed for oob parameter

# Coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# Examine this data frame
head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .05, "Significant", NA)

write.csv(deseq2ResDF, "metaT_files/t_e_deseq_resdf.csv")


