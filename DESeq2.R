#Author: Kiran Bharat Gaikwad
#Date: 21/06/2023
#Purpose: Estimate differentially expressed ion channels from rna-seq datasets of patients with breast cancer (DESeq2)

#Data import and preprocessing
file = "htseq_all_sample_count.tsv",sep="\t", row.names = FALSE)

#Differential expression analysis
library("DESeq2")
file1 = read.csv("htseq_all_sample_count.tsv", sep='\t', header = TRUE, row.names = 1, check.names =TRUE)
dim(file1)
head(file1)
countdata1 = as.matrix(file1)
dim(countdata1)
(colnames(countdata1))
condition1 = factor(c(rep("control",84),rep("ktcn",104)))
condition1
coldata1 = data.frame(row.names = colnames(countdata1), condition1)
coldata1
ddsFull = DESeqDataSetFromMatrix(countData=countdata1, colData= coldata1, design =~ condition1)
ddsFull
dds = DESeq(ddsFull)
dds
res = results(dds)
res
summary(res)
res_ordered = res[order(res$pvalue),]
res_d = as.data.frame(res)
res_d
res_d$diffexpressed <- "NO"


#Upregulated and downregulated identification
#if log2FoldChange > 1 and padj < 0.05, set as "UP"
res_d$diffexpressed[res_d$log2FoldChange > 1.0 & res_d$padj < 0.05] <- "UP"
#if log2FoldChange < -1 and padj < 0.05, set as "DOWN"
res_d$diffexpressed[res_d$log2FoldChange < -1.0 & res_d$padj < 0.05] <- "DOWN"
res_d <- cbind(rownames(res_d), data.frame(res_d, row.names=NULL))
res_d
colnames(res_d)[1] <- "hgnc_symbol"
res_d

a = res_d[which(res_d$diffexpressed == 'UP'),]
b = res_d[which(res_d$diffexpressed == 'DOWN'),]


#Export results to files 
write.table(a, file='PRJNA799648_up.txt', sep='\t', row.names = FALSE, col.names=TRUE)
write.table(b, file='PRJNA799648_down.txt', sep='\t', row.names = FALSE, col.names=TRUE)
write.table(res_d, file='all_genes.tsv', sep='\t', row.names = FALSE, col.names = TRUE)


#Processing for plotting
#BiocManager::install("ggplot2")
library(ggplot2)
res_d = read.csv('all_genes.tsv', sep='\t', header = TRUE, check.names =TRUE)
res_d
plot1 = ggplot(res_d, aes(x = log2FoldChange, y = -log10(pvalue)))+
  geom_point(aes(colour=diffexpressed), size=2, alpha=1)+
  scale_color_manual(values = c("green","grey","red"))+
  labs(title = "Keratoconus - Differentially expressed ion channels")+
  geom_vline(xintercept=c(-1,1), lty=4, col="black", lwd=0.8)+
  geom_hline(yintercept=1.301, lty=4, col="black", lwd=0.8)

plot1 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black"), plot.title = element_text(hjust = 0.5))
