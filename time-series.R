counts <- read.delim("/mnt/projects/veronika/results/anduril/execute/_deseq_ERvsNTd15_counts_array1/array/_index")
samples <- read.delim("/mnt/projects/veronika/data/samples.csv")

# prepare sample table
sampleTable <- merge(samples, countFiles, by.x="Alias", by.y="Key")
names(sampleTable)[ncol(sampleTable)] <- "fileName"
sampleTable <- sampleTable[,c(c("Alias", "fileName"), names(sampleTable)[!names(sampleTable) %in% c("Alias", "fileName")])]

sampleTable$day <- NA
sampleTable$day[grepl("3_", sampleTable$Alias)]
sampleTable$day[grepl("3_", sampleTable$Alias)] <- "day3"
sampleTable$day[grepl("8_", sampleTable$Alias)] <- "day8"
sampleTable$day[grepl("15_", sampleTable$Alias)] <- "day15"
sampleTable$day <- factor(sampleTable$day, levels=c("day3", "day8", "day15"))

sampleTable$treatment <- NA
sampleTable$treatment[grepl("Kd", sampleTable$Alias)]
sampleTable$treatment[grepl("Kd", sampleTable$Alias)] <- "knockdown"
sampleTable$treatment[grepl("Nt", sampleTable$Alias)] <- "non-targeting"
sampleTable$treatment <- factor(sampleTable$treatment, levels=c("non-targeting", "knockdown"))

# DESeq
library(DESeq2)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "/", design = ~ treatment + day + treatment:day)
dds <- DESeq(dds, test="LRT", reduced = ~ treatment + day)
res <- results(dds)

#plotdata <- data.frame(gene=character(0), sample=character(0), count=numeric(0), day=character(0), treatment=character(0))
#g <- "ENSG00000002586"
plot.top.n <- 30
pdf("/mnt/projects/veronika/results/time-series.pdf")
for (g in res@rownames[order(res$padj)][1:plot.top.n]) {
  data <- plotCounts(dds, g, intgroup=c("day","treatment"), returnData=TRUE)
  p <- ggplot(data, aes(x=day, y=count, color=treatment, group=treatment)) +
    theme_bw() +
    geom_point() + 
    stat_smooth(se=FALSE,method="loess") + 
    scale_y_log10() +
    ggtitle(sprintf("%s FC=%.1f padj=%.2g", g, res[g,"log2FoldChange"], res[g,"padj"]))
  print(p)
}
dev.off()

# pairwise comparison
sampleTable.3to8 <- sampleTable[sampleTable$day %in% c("day3", "day8"),] ; sampleTable.3to8$day <- factor(sampleTable.3to8$day, levels=c("day3", "day8"))
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable.3to8, directory = "/", design = ~ treatment + day + treatment:day)
dds <- DESeq(dds, test="LRT", reduced = ~ treatment + day)
res <- results(dds)
res <- res[order(res$padj),]

# significant up
up <- as.data.frame(res[!is.na(res$padj) & res$padj < 0.01 & res$log2FoldChange > 0,])
dn <- as.data.frame(res[!is.na(res$padj) & res$padj < 0.01 & res$log2FoldChange < 0,])

#g <- "ENSG00000111145"
plot.top.n <- 30
pdf("/mnt/projects/veronika/results/time-series.day3to8.pdf")
for (g in res@rownames[order(res$padj)][1:plot.top.n]) {
  data <- plotCounts(dds, g, intgroup=c("day","treatment"), returnData=TRUE)
  p <- ggplot(data, aes(x=day, y=count, color=treatment, group=treatment)) +
    theme_bw() +
    geom_point() + 
    stat_smooth(se=FALSE,method="loess") + 
    scale_y_log10() +
    ggtitle(sprintf("%s FC=%.1f padj=%.2g", g, res[g,"log2FoldChange"], res[g,"padj"]))
  print(p)
}
dev.off()
