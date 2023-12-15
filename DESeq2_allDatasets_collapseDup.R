setwd('D:\\IOB\\Projects\\keratoconus_072322\\rebuttal_all_genes_epithelial_samples\\replicates\\all_datasets_all_genes')

library("DESeq2")
library(dplyr)
file1 = read.csv("keratoconus_all_epithelial_051222.txt", sep='\t', header = TRUE, row.names = 1, check.names =TRUE)
dim(file1)
head(file1)
countdata1 = as.matrix(file1)
dim(countdata1)
(colnames(countdata1))
condition1 = factor(c(rep("control",84),rep("keratoconus",104)))
condition1
samples = c(c(rep("KR_38",5),rep("KR_50",3),rep("KR_51",3),rep("KR_53",6),rep("KR_56",4),rep("KR_57",4),rep("KR_61",3),rep("KR_64",3),rep("KR_08",1),rep("KR_13",1),rep("KR_14",1),rep("KR_15",1),rep("KR_18",1),rep("KR_19",4),rep("KR_21",4),rep("KR_23",1),rep("KR_24",3),rep("KR_25",3),rep("KR_26",3),rep("KR_27",1),rep("KR_28",1),rep("KR_33",1),rep("KR_35",1),rep("KR_48",1),rep("KR_52",1),rep("con_43",1),rep("con_48",1),rep("con_49",1),rep("con_51",1),rep("con_55",1),rep("con_S2",1),rep("con_S7",1),rep("con_S8",1),rep("con_S9",1),rep("con_S10",1),rep("GSM4587111",1),rep("GSM4587112",1),rep("GSM4587113",1),rep("GSM4587114",1),rep("GSM4587115",1),rep("GSM4587116",1),rep("GSM4587117",1),rep("CON_T_1",1),rep("CON_T_2",1),rep("CON_T_3",1),rep("CON_T_4",1),rep("CON_T_5",1),rep("CON_T_6",1),rep("CON_T_7",1),rep("KC_35",5),rep("KC_37",2),rep("KC_39",3),rep("KC_40",2),rep("KC_43",3),rep("KC_45",4),rep("KC_47",3),rep("KC_52",3),rep("KC_15",4),rep("KC_16",4),rep("KC_17",3),rep("KC_18",3),rep("KC_19",3),rep("KC_20",4),rep("KC_21",3),rep("KC_22",3),rep("KC_23",3),rep("KC_25",1),rep("KC_26",3),rep("KC_28",3),rep("KC_30",1),rep("KC_31",1),rep("KC_32",1),rep("KC_34",2),rep("KC_55",1),rep("GSE112155_KC_17",1),rep("GSE112155_KC_18",1),rep("GSE112155_KC_19",1),rep("GSE112155_KC_20",1),rep("GSE112155_KC_21",1),rep("kc_s1",1),rep("kc_s3",1),rep("kc_s4",1),rep("kc_s5",1),rep("kc_s6",1),rep("GSM4587118",1),rep("GSM4587119",1),rep("GSM4587120",1),rep("GSM4587121",1),rep("GSM4587122",1),rep("GSM4587123",1),rep("GSM4587124",1),rep("GSM4587125",1),rep("GSM4587126",1),rep("GSM4587127",1),rep("GSM4587128",1),rep("GSM4587129",1),rep("GSM4587130",1),rep("GSM4587131",1),rep("GSM4587132",1),rep("GSM4587133",1),rep("GSM4587134",1),rep("GSM4587135",1),rep("GSM4587136",1),rep("KC_T_1",1),rep("KC_T_2",1),rep("KC_T_3",1),rep("KC_T_4",1),rep("KC_T_5",1),rep("KC_T_6",1),rep("KC_T_7",1)))


run = colnames(countdata1)
coldata1 = data.frame(condition1, run, samples)

#coldata1$samples = row.names(coldata1)
rownames(coldata1) = 1:nrow(coldata1)
dim(coldata1)
head(coldata1)

ddsFull = DESeqDataSetFromMatrix(countData=countdata1, colData= coldata1, design =~ condition1)
ddsFull$samples
ddsFull$run
ddsCollapsed <- collapseReplicates( ddsFull,
                                    groupby = ddsFull$samples,
                                    run = ddsFull$run )
ddsFull$samples
ddsFull$run
confirm = colData(ddsCollapsed)

write.table(confirm, file='to_confirm_collapsing.txt', sep='\t', row.names = FALSE, col.names=TRUE)

original <- rowSums( counts(ddsFull)[ , ddsFull$sample == "KC_15" ] )
all( original == counts(ddsCollapsed)[ ,"KC_15" ] )

dds = DESeq(ddsCollapsed)
dds

res = results(dds)
res
summary(res)
res_ordered = res[order(res$padj),]
res_d = as.data.frame(res)
res_d
res_d$diffexpressed <- "NO"
#applying criteria
#if log2FoldChange > 1 and padj < 0.05, set as "UP"
res_d$diffexpressed[res_d$log2FoldChange > 1.0 & res_d$padj <= 0.05] <- "UP"
#if log2FoldChange < -1 and padj < 0.05, set as "DOWN"
res_d$diffexpressed[res_d$log2FoldChange < -1.0 & res_d$padj <= 0.05] <- "DOWN"
res_d <- cbind(rownames(res_d), data.frame(res_d, row.names=NULL))
res_d
colnames(res_d)[1] <- "ensembl_gene_id"
res_d

a = res_d[which(res_d$diffexpressed == 'UP'),]
nrow(a)
a
a$ensembl_gene_id
b = res_d[which(res_d$diffexpressed == 'DOWN'),]
b$ensembl_gene_id
b
nrow(b)
#setwd('D:\\IOB\\Projects\\keratoconus\\DESeq\\GSE77938')
write.table(a, file='replicatesCollapsed_allGenes_up.txt', sep='\t', row.names = FALSE, col.names=TRUE)
write.table(b, file='replicatesCollapsed_allGenes_down.txt', sep='\t', row.names = FALSE, col.names=TRUE)
