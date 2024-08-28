library(GEOquery)
library(limma)
library(affy)
library(dplyr)
library(hugene10sttranscriptcluster.db)

#reading affy data (.cel files) in R
eset=ReadAffy(celfile.path = "./GSM1203269_RAW")
length(eset)

#reading annotation file in R
GPL6244 = getGEO("GPL6244", destdir = ".")
GPL6244 = Table(GPL6244)

#reannotating probes to latest ID
x <- hugene10sttranscriptclusterSYMBOL
ID <- as.character(GPL6244$ID)# Get the probe identifiers that are mapped to a gene symbol
Symbol <- as.list(x[ID])# Convert to a list
GPL6244 = cbind(Symbol,ID)
GPL6244 = as.data.frame(GPL6244)
View(GPL6244)

#backgroud correction/log transformation
eset=rma(eset,normalize = TRUE,background = TRUE)
exp_mat = exprs(eset) #extracting expression matrix
exp_mat_1 = as.data.frame(exp_mat)
exp_mat_1$ID_REF = rownames(exp_mat)
eset_merged = merge.data.frame(x = exp_mat_1,y = GPL6244,by.x = 'ID_REF',by.y = 'ID')

features = read.csv("comorbidity_features_list.csv",header = FALSE)
data_67_genes = eset_merged[eset_merged$Symbol %in% features$V1, ]
data_67_genes$Symbol = unlist(data_67_genes$Symbol)
data_67_genes = data_67_genes[,-1]
data_67_genes_mean = data_67_genes %>% group_by(Symbol) %>%
  summarise(across(everything(), list(mean = mean)))
data_67_genes_mean = as.data.frame(data_67_genes_mean)

genes_missing = features$V1[-which(features$V1 %in% data_67_genes_mean$Symbol)]
df_missing = matrix(nrow = length(genes_missing),ncol = length(data_67_genes_mean))
df_missing[,1] = genes_missing
df_missing = as.data.frame(df_missing)
df_missing[is.na(df_missing)] <- 0
colnames(df_missing) = colnames(data_67_genes_mean)

final_data = rbind(data_67_genes_mean,df_missing)
colnames(final_data)[1] = "genes"
final_data = as.data.frame(final_data)
final_data_1 = final_data[order(final_data$genes,decreasing = FALSE),]

write.csv(final_data_1,"input_data.csv",row.names = FALSE)
