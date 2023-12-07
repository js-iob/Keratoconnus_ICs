setwd("D:\\IOB\\Projects\\FED_062723\\DESeq2")
data1 = read.table("DE_ic_expression_profile_081323.txt", sep = '\t',header = TRUE)

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
run = factor(c(rep("KR_38",5),rep("KR_50",3),rep("KR_51",3),rep("KR_53",6),rep("KR_56",4),rep("KR_57",4),rep("KR_61",3),rep("KR_64",3),rep("KR_08",1),rep("KR_13",1),rep("KR_14",1),rep("KR_15",1),rep("KR_18",1),rep("KR_19",4),rep("KR_21",4),rep("KR_23",1),rep("KR_24",3),rep("KR_25",3),rep("KR_26",3),rep("KR_27",1),rep("KR_28",1),rep("KR_33",1),rep("KR_35",1),rep("KR_48",1),rep("KR_52",1),rep('Con_43',1),rep('Con_48',1),rep('Con_49',1),rep('Con_51',1),rep('Con_55',1),rep('Con_S2',1),rep('Con_S7',1),rep('Con_S8',1),rep('Con_S9',1),rep('Con_S10',1),rep('Control_1',1),rep('Control_2',1),rep('Control_3',1),rep('Control_4',1),rep('Control_5',1),rep('Control_6',1),rep('Control_7',1),rep('CON_T_1',1),rep('CON_T_2',1),rep('CON_T_3',1),rep('CON_T_4',1),rep('CON_T_5',1),rep('CON_T_6',1),rep('CON_T_7',1),rep("KC_35",5),rep("KC_37",2),rep("KC_39",3),rep("KC_40",2),rep("KC_43",3),rep("KC_45",4),rep("KC_47",3),rep("KC_52",3),rep("KC_15",4),rep("KC_16",4),rep("KC_17",3),rep("KC_18",3),rep("KC_19",3),rep("KC_20",4),rep("KC_21",3),rep("KC_22",3),rep("KC_23",3),rep("KC_25",1),rep("KC_26",3),rep("KC_28",3),rep("KC_30",1),rep("KC_31",1),rep("KC_32",1),rep("KC_34",2),rep("KC_55",1),rep("GSM3058950",1),rep("GSM3058951",1),rep("GSM3058952",1),rep("GSM3058953",1),rep("GSM3058954",1),rep("GSM3058955",1),rep("GSM3058956",1),rep("GSM3058957",1),rep("GSM3058958",1),rep("GSM3058959",1),rep("Keratoconus_patient_7",1),rep("Keratoconus_patient_8",1),rep("Keratoconus_patient_9",1),rep("Keratoconus_patient_10",1),rep("Keratoconus_patient_11",1),rep("Keratoconus_patient_12",1),rep("Keratoconus_patient_13",1),rep("Keratoconus_patient_14",1),rep("Keratoconus_patient_15",1),rep("Keratoconus_patient_16",1),rep("Keratoconus_patient_17",1),rep("Keratoconus_patient_18",1),rep("Keratoconus_patient_19",1),rep("Keratoconus_patient_1",1),rep("Keratoconus_patient_2",1),rep("Keratoconus_patient_3",1),rep("Keratoconus_patient_4",1),rep("Keratoconus_patient_5",1),rep("Keratoconus_patient_6",1),rep("KC_T_1",1),rep("KC_T_2",1),rep("KC_T_3",1),rep("KC_T_4",1),rep("KC_T_5",1),rep("KC_T_6",1),rep("KC_T_7",1)))
run
coldata1 = data.frame(row.names =seq(1:188), colnames(countdata1), condition1, run)
countdata1 = as.matrix(data1)
coldata1$samples = row.names(coldata1)
rownames(coldata1) = 1:nrow(coldata1)
dim(coldata1)
head(coldata1)
ddsFull = DESeqDataSetFromMatrix(countData=countdata1, colData= coldata1, design =~ condition1)

ddsCollapsed <- collapseReplicates( ddsFull,
                                    groupby = ddsFull$run,
                                    run = ddsFull$samples )
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
control = norm.counts[,c(1:24,86:110)]
kc = norm.counts[,c(25:85)]
colnames(kc)
cor_data = t(kc)
dim(cor_data)

cor = cor(cor_data, method = "pearson")
write.table(cor, file='control_epithelial_cor_data.tsv', sep='\t')
library(corrplot)

#cor
#corrplot(cor, method='circle', type="upper", order="hclust", xlab='Ion channels', ylab='Ion channels', tl.col="black", tl.srt=90)
#dev.off()
cor_test_mat <- corr.test(cor_data)$p
write.table(cor_test_mat, file='control_epithelial_samples_correlation_IC_homogeneity_090223.tsv', sep='\t')
pdf("kc_epithelial_corrplot.pdf", width = 200, height = 200)
corrplot(cor, method='circle', p.mat = cor_test_mat, type="upper", xlab='Ion channels', ylab='Ion channels', tl.col="black", tl.srt=90)
dev.off()
#corrplot(cor, order="hclust", xlab='Ion channels', ylab='Ion channels', tl.col="black", tl.srt=90)


