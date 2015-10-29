fuka.d13 <- read.delim("/mnt/projects/veronika/data/fuka/telamlKD.esetnsF.onlyG_13.annot.xls")
names(fuka.d13) <- paste0(names(fuka.d13), ".fuka.d13")
fuka.d20plus <- read.delim("/mnt/projects/chrisi/results/fuka/matAnn.telamlKD.REHandAT2.esetnsF.REH.AT2.balanced.annot.tsv")
names(fuka.d20plus) <- paste0(names(fuka.d20plus), ".fuka.d20plus")

d3 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd3/table.csv")
d3 <- merge(d3, fuka.d13[,c("syms.fuka.d13", "logFC.fuka.d13")], by.x="Gene", by.y="syms.fuka.d13", all.x=T)
d3 <- merge(d3, fuka.d20plus[,c("syms.fuka.d20plus", "logFC.fuka.d20plus")], by.x="Gene", by.y="syms.fuka.d20plus", all.x=T)
d8 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd8/table.csv")
d8 <- merge(d8, fuka.d13[,c("syms.fuka.d13", "logFC.fuka.d13")], by.x="Gene", by.y="syms.fuka.d13", all.x=T)
d8 <- merge(d8, fuka.d20plus[,c("syms.fuka.d20plus", "logFC.fuka.d20plus")], by.x="Gene", by.y="syms.fuka.d20plus", all.x=T)
d15 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd15/table.csv")
d15 <- merge(d15, fuka.d13[,c("syms.fuka.d13", "logFC.fuka.d13")], by.x="Gene", by.y="syms.fuka.d13", all.x=T)
d15 <- merge(d15, fuka.d20plus[,c("syms.fuka.d20plus", "logFC.fuka.d20plus")], by.x="Gene", by.y="syms.fuka.d20plus", all.x=T)

qCutoffPoint <- 1e-3
qCutoffLabel <- 1e-40

# FUKA DAY 13 ------------------------------------------------------------------------------------

pdf("/mnt/projects/veronika/results/correlation-fuka-day13.pdf", height=12, width=12)
par(mfrow=c(2,2))

# day 3
data <- d3[!is.na(d3$qValue) & d3$qValue < qCutoffPoint & !is.na(d3$logFC.fuka.d13),]
fit <- lm(logFC.fuka.d13 ~ fc, data)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(data$fc, data$logFC.fuka.d13, cex=0.2, xlab="FC Veronika day 3", ylab="FC Fuka day 13", main=sprintf("E/R KD Veronika D3 vs. Fuka D13\nq<%.2g n=%d p=%.2g R2=%.2g", qCutoffPoint, nrow(data), p, R))
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- data[(!is.na(data$qValue) & data$qValue < qCutoffLabel) | abs(data$logFC.fuka.d13) >= 1,]
text(sig$fc, sig$logFC.fuka.d13-0.05, sig$Gene, cex=0.5, col="red")

# day 8
data <- d8[!is.na(d8$qValue) & d8$qValue < qCutoffPoint & !is.na(d3$logFC.fuka.d13),]
fit <- lm(logFC.fuka.d13 ~ fc, data)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(data$fc, data$logFC.fuka.d13, cex=0.2, xlab="FC Veronika day 8", ylab="FC Fuka day 13", main=sprintf("E/R KD Veronika D8 vs. Fuka D13\nq<%.2g n=%d p=%.2g R2=%.2g", qCutoffPoint, nrow(data), p, R))
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- data[(!is.na(data$qValue) & data$qValue < qCutoffLabel) | abs(data$logFC.fuka.d13) >= 1,]
text(sig$fc, sig$logFC.fuka.d13-0.05, sig$Gene, cex=0.5, col="red")

# day 15
data <- d15[!is.na(d15$qValue) & d15$qValue < qCutoffPoint & !is.na(d3$logFC.fuka.d13),]
fit <- lm(logFC.fuka.d13 ~ fc, data)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(data$fc, data$logFC.fuka.d13, cex=0.2, xlab="FC Veronika day 15", ylab="FC Fuka day 13", main=sprintf("E/R KD Veronika D15 vs. Fuka D13\nq<%.2g n=%d p=%.2g R2=%.2g", qCutoffPoint, nrow(data), p, R))
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- data[(!is.na(data$qValue) & data$qValue < qCutoffLabel) | abs(data$logFC.fuka.d13) >= 1,]
text(sig$fc, sig$logFC.fuka.d13-0.05, sig$Gene, cex=0.5, col="red")

dev.off()

# FUKA DAY 20, 25 und 30 (REH & AT2) -----------------------------------------------------------------------

pdf("/mnt/projects/veronika/results/correlation-fuka-day20plus.pdf", height=12, width=12)
par(mfrow=c(2,2))

# day 3
data <- d3[!is.na(d3$qValue) & d3$qValue < qCutoffPoint & !is.na(d3$logFC.fuka.d20plus),]
fit <- lm(logFC.fuka.d20plus ~ fc, data)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(data$fc, data$logFC.fuka.d20plus, cex=0.2, xlab="FC Veronika day 3", ylab="FC Fuka day 20+", main=sprintf("E/R KD Veronika D3 vs. Fuka D20+\nq<%.2g n=%d p=%.2g R2=%.2g", qCutoffPoint, nrow(data), p, R))
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- data[(!is.na(data$qValue) & data$qValue < qCutoffLabel) | abs(data$logFC.fuka.d20plus) >= 1,]
text(sig$fc, sig$logFC.fuka.d20plus-0.05, sig$Gene, cex=0.5, col="red")

# day 8
data <- d8[!is.na(d8$qValue) & d8$qValue < qCutoffPoint & !is.na(d3$logFC.fuka.d20plus),]
fit <- lm(logFC.fuka.d20plus ~ fc, data)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(data$fc, data$logFC.fuka.d20plus, cex=0.2, xlab="FC Veronika day 8", ylab="FC Fuka day 20+", main=sprintf("E/R KD Veronika D8 vs. Fuka D20+\nq<%.2g n=%d p=%.2g R2=%.2g", qCutoffPoint, nrow(data), p, R))
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- data[(!is.na(data$qValue) & data$qValue < qCutoffLabel) | abs(data$logFC.fuka.d20plus) >= 1,]
text(sig$fc, sig$logFC.fuka.d20plus-0.05, sig$Gene, cex=0.5, col="red")

# day 15
data <- d15[!is.na(d15$qValue) & d15$qValue < qCutoffPoint & !is.na(d3$logFC.fuka.d20plus),]
fit <- lm(logFC.fuka.d20plus ~ fc, data)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(data$fc, data$logFC.fuka.d20plus, cex=0.2, xlab="FC Veronika day 15", ylab="FC Fuka day 20+", main=sprintf("E/R KD Veronika D15 vs. Fuka D20+\nq<%.2g n=%d p=%.2g R2=%.2g", qCutoffPoint, nrow(data), p, R))
abline(fit, col="red")
abline(v=0, lty=3)
abline(h=0, lty=3)
sig <- data[(!is.na(data$qValue) & data$qValue < qCutoffLabel) | abs(data$logFC.fuka.d20plus) >= 1,]
text(sig$fc, sig$logFC.fuka.d20plus-0.05, sig$Gene, cex=0.5, col="red")

dev.off()
