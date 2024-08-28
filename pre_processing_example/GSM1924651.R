library(GEOquery)
library(affy)
library(dplyr)
library(illuminaHumanv4.db)

#############
#reading limma data in R
eset = lumiR("GSM1924651.csv")
length(eset)

#reading annotation file in R
GPL10558 = getGEO("GPL10558", destdir = ".")
GPL10558 = Table(GPL10558)

#reannotating probes to latest ID
x <- illuminaHumanv4SYMBOL
mapped_probes <- mappedkeys(x)
# Get the probe identifiers that are mapped to a gene symbol
ID <-as.character(GPL10558$ID)
# Convert to a list
Symbol <- as.list(x[ID])
GPL10558 = cbind(Symbol,ID)
GPL10558 = as.data.frame(GPL10558)
View(GPL10558)

#backgroud correction/normalization/log transformation
eset_B = lumiB(eset, method='forcePositive')
eset_BT = lumiT(eset_B, method='log2')
eset_BTN = lumiN(eset_BT, method='quantile')

exp_mat = exprs(eset_BTN)#extracting expression matrix
exp_mat_1 = as.data.frame(exp_mat)
exp_mat_1$ID_REF = rownames(exp_mat)
eset_merged = merge.data.frame(x = exp_mat_1,y = GPL10558,by.x = 'ID_REF',by.y = 'ID')

features = read.csv("comorbidity_features_list.csv",header = FALSE)
data_67_genes = eset_merged[eset_merged$Symbol %in% features$V1, ]
data_67_genes$Symbol = unlist(data_67_genes$Symbol)
data_67_genes = data_67_genes[,-1]
data_67_genes_mean = data_67_genes %>% group_by(Symbol) %>%  summarise(across(everything(), list(mean = mean)))
data_67_genes_mean = as.data.frame(data_67_genes_mean)

genes_missing = features$V1[-which(features$V1 %in% data_67_genes_mean$Symbol)]
df_missing = matrix(nrow = length(genes_missing),ncol = length(data_67_genes_mean))
df_missing[,1] = genes_missing
df_missing = as.data.frame(df_missing)
df_missing[is.na(df_missing)] <- 0

final_data = rbind(data_67_genes_mean,df_missing)
colnames(final_data)[1] = "genes"
final_data = as.data.frame(final_data)
final_data_1 = final_data[order(final_data$genes,decreasing = FALSE),]

write.csv(final_data,"input_data_GSM1924651.csv",row.names = FALSE)
