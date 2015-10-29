d3 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseq_ERvsNTd3/results.csv")
d8 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseq_ERvsNTd8/results.csv")
d15 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseq_ERvsNTd15/results.csv")

pdf("/mnt/projects/veronika/results/expression-pairwise-correlations-exp2.pdf", width=12, height=12)
par(mfcol=c(3,3))

# day 3
plot(d3$KdD3_1, d3$KdD3_2, log="xy", cex=0.1, xlab="KdD3 (replicate 1)", ylab="KdD3 (replicate 2)", main="Exp #2 D3 E/R KD replicate 1 vs. 2") ; abline(0, 1)
plot(d3$NtD3_1, d3$NtD3_2, log="xy", cex=0.1, xlab="NtD3 (replicate 1)", ylab="NtD3 (replicate 2)", main="Exp #2 D3 NT control replicate 1 vs. 2") ; abline(0, 1)
plot(d3$KdD3_1, d3$NtD3_1, log="xy", cex=0.1, xlab="KdD3 (replicate 1)", ylab="NtD3 (replicate 1)", main="Exp #2 D3 E/R KD vs. NT control") ; abline(0, 1)

# day 8
plot(d8$KdD8_1, d8$KdD8_2, log="xy", cex=0.1, xlab="KdD8 (replicate 1)", ylab="KdD8 (replicate 2)", main="Exp #2 D8 E/R KD replicate 1 vs. 2") ; abline(0, 1)
plot(d8$NtD8_1, d8$NtD8_2, log="xy", cex=0.1, xlab="NtD8 (replicate 1)", ylab="NtD8 (replicate 2)", main="Exp #2 D8 NT control replicate 1 vs. 2") ; abline(0, 1)
plot(d8$KdD8_1, d8$NtD8_1, log="xy", cex=0.1, xlab="KdD8 (replicate 1)", ylab="NtD8 (replicate 1)", main="Exp #2 D8 E/R KD vs. NT control") ; abline(0, 1)

# day 15
plot(d15$KdD15_1, d15$KdD15_2, log="xy", cex=0.1, xlab="KdD15 (replicate 1)", ylab="KdD15 (replicate 2)", main="Exp #2 D15 E/R KD replicate 1 vs. 2") ; abline(0, 1)
plot(d15$NtD15_1, d15$NtD15_2, log="xy", cex=0.1, xlab="NtD15 (replicate 1)", ylab="NtD15 (replicate 2)", main="Exp #2 D15 NT control replicate 1 vs. 2") ; abline(0, 1)
plot(d15$KdD15_1, d15$NtD15_1, log="xy", cex=0.1, xlab="KdD15 (replicate 1)", ylab="NtD15 (replicate 1)", main="Exp #2 D15 E/R KD vs. NT control") ; abline(0, 1)

dev.off()

pdf("/mnt/projects/veronika/results/expression-pairwise-correlations-exp1-vs-exp2.pdf", width=12, height=4)
par(mfcol=c(1,3))

plot(d3$KdD3_3, d3$KdD3_4, log="xy", cex=0.1, xlab="KdD3 (exp #1, rep #1)", ylab="KdD3 (exp #1, rep #2)", main="Exp #1 D3 E/R KD Rep #1 vs. #2", ylim=c(1,60000)) ; abline(0, 1)
plot(d3$KdD3_1, d3$KdD3_2, log="xy", cex=0.1, xlab="KdD3 (exp #2, rep #1)", ylab="KdD3 (exp #2, rep #2)", main="Exp #2 D3 E/R KD Rep #1 vs. #2") ; abline(0, 1)
plot(d3$KdD3_1, d3$KdD3_3, log="xy", cex=0.1, xlab="KdD3 (exp #2, rep #1)", ylab="KdD3 (exp #1, rep #1)", main="D3 E/R KD Exp #1 vs. Exp #2") ; abline(0, 1)

dev.off()
