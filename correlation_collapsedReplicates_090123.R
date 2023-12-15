setwd("D:\\IOB\\Projects\\keratoconus_072322\\rebuttal_all_genes_epithelial_samples\\replicates\\all_datasets_all_genes\\ic\\collapsedRep\\correlation")
data1 = read.table("ic_expression_121323.txt", sep = '\t',header = TRUE)

#install.packages("WGCNA")
#1. prepare data
head(data1)
#rownames(data) <- data$hgnc_symbol
library("dplyr")
library("tibble")
data1 <- data1 %>% remove_rownames %>% column_to_rownames(var="hgnc_symbol")
head(data1)
dim(data1)

#2. QC - outlier detection ------------------------------------------------
library("WGCNA")
datExpr = data.matrix(data1)
dim(datExpr)
gsg <- goodSamplesGenes(t(datExpr))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)      
data1 <- data1[gsg$goodGenes == TRUE,]
outlier = data1[gsg$goodGenes == FALSE,]
write.table(outlier, file='outlier_ICs.tsv')
head(data1)
nrow(data1)

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset
library('DESeq2')
condition1 = factor(c(rep("control",84),rep("keratoconus",104)))
condition1
samples = c(c(rep("KR_38",5),rep("KR_50",3),rep("KR_51",3),rep("KR_53",6),rep("KR_56",4),rep("KR_57",4),rep("KR_61",3),rep("KR_64",3),rep("KR_08",1),rep("KR_13",1),rep("KR_14",1),rep("KR_15",1),rep("KR_18",1),rep("KR_19",4),rep("KR_21",4),rep("KR_23",1),rep("KR_24",3),rep("KR_25",3),rep("KR_26",3),rep("KR_27",1),rep("KR_28",1),rep("KR_33",1),rep("KR_35",1),rep("KR_48",1),rep("KR_52",1),rep("con_43",1),rep("con_48",1),rep("con_49",1),rep("con_51",1),rep("con_55",1),rep("con_S2",1),rep("con_S7",1),rep("con_S8",1),rep("con_S9",1),rep("con_S10",1),rep("GSM4587111",1),rep("GSM4587112",1),rep("GSM4587113",1),rep("GSM4587114",1),rep("GSM4587115",1),rep("GSM4587116",1),rep("GSM4587117",1),rep("CON_T_1",1),rep("CON_T_2",1),rep("CON_T_3",1),rep("CON_T_4",1),rep("CON_T_5",1),rep("CON_T_6",1),rep("CON_T_7",1),rep("KC_35",5),rep("KC_37",2),rep("KC_39",3),rep("KC_40",2),rep("KC_43",3),rep("KC_45",4),rep("KC_47",3),rep("KC_52",3),rep("KC_15",4),rep("KC_16",4),rep("KC_17",3),rep("KC_18",3),rep("KC_19",3),rep("KC_20",4),rep("KC_21",3),rep("KC_22",3),rep("KC_23",3),rep("KC_25",1),rep("KC_26",3),rep("KC_28",3),rep("KC_30",1),rep("KC_31",1),rep("KC_32",1),rep("KC_34",2),rep("KC_55",1),rep("GSE112155_KC_17",1),rep("GSE112155_KC_18",1),rep("GSE112155_KC_19",1),rep("GSE112155_KC_20",1),rep("GSE112155_KC_21",1),rep("kc_s1",1),rep("kc_s3",1),rep("kc_s4",1),rep("kc_s5",1),rep("kc_s6",1),rep("GSM4587118",1),rep("GSM4587119",1),rep("GSM4587120",1),rep("GSM4587121",1),rep("GSM4587122",1),rep("GSM4587123",1),rep("GSM4587124",1),rep("GSM4587125",1),rep("GSM4587126",1),rep("GSM4587127",1),rep("GSM4587128",1),rep("GSM4587129",1),rep("GSM4587130",1),rep("GSM4587131",1),rep("GSM4587132",1),rep("GSM4587133",1),rep("GSM4587134",1),rep("GSM4587135",1),rep("GSM4587136",1),rep("KC_T_1",1),rep("KC_T_2",1),rep("KC_T_3",1),rep("KC_T_4",1),rep("KC_T_5",1),rep("KC_T_6",1),rep("KC_T_7",1)))

run = colnames(countdata1)
coldata1 = data.frame(condition1, run, samples)
countdata1 = as.matrix(data1)

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
ddsFull$condition1 = relevel(ddsFull$condition1, ref = "control")

dds1 = DESeq(ddsCollapsed)

dds1
rownames(dds1)
dds_norm <- varianceStabilizingTransformation(dds1)
colnames(dds_norm)
norm.counts <- assay(dds_norm)
colnames(norm.counts)
#t_norm_counts <- t(norm.counts)
#nrow(t_norm_counts)
dim(norm.counts)
head(norm.counts)
control = norm.counts[,c(1:17,23:29, 86:110)]
dim(control)
kc = norm.counts[,c(18:22, 30:85)]
colnames(kc)
cor_data = t(control)
dim(control)

cor = cor(cor_data, method = "pearson")
write.table(cor, file='control_epithelial_cor_data.tsv', sep='\t')
library(corrplot)
library(psych)
#cor
#corrplot(cor, method='circle', type="upper", order="hclust", xlab='Ion channels', ylab='Ion channels', tl.col="black", tl.srt=90)
#dev.off()
cor_test_mat <- corr.test(cor_data)$p
 write.table(cor_test_mat, file='kc_epithelial_samples_correlation_IC_pValue_121323.tsv', sep='\t')
pdf("control_epithelial_corrplot.pdf", width = 200, height = 200)
corrplot(cor, method='circle', p.mat = cor_test_mat, type="lower", xlab='Ion channels', ylab='Ion channels', tl.col="black", tl.srt=90)
dev.off()
#corrplot(cor, order="hclust", xlab='Ion channels', ylab='Ion channels', tl.col="black", tl.srt=90)


