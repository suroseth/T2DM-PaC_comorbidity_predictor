library(dplyr)

#reading normalized sample
exp_mat = read.table("Control_sample_28.txt",sep = "\t",header = TRUE)
eset_merged = as.data.frame(exp_mat)

features = read.csv("comorbidity_features_list.csv",header = FALSE)
data_67_genes = eset_merged[eset_merged$X %in% features$V1, ]

data_67_genes_mean = data_67_genes %>% group_by(X) %>%
  summarise(across(everything(), list(mean = mean)))
data_67_genes_mean = as.data.frame(data_67_genes_mean)

genes_missing = features$V1[-which(features$V1 %in% data_67_genes_mean$X)]
df_missing = matrix(nrow = length(genes_missing),ncol = length(data_67_genes_mean ))
df_missing[,1] = genes_missing
df_missing = as.data.frame(df_missing)
df_missing[is.na(df_missing)] <- 0
colnames(df_missing) = colnames(data_67_genes_mean )

final_data = rbind(data_67_genes_mean ,df_missing)
colnames(final_data)[1] = "genes"
final_data = as.data.frame(final_data)
final_data_1 = final_data[order(final_data$genes,decreasing = FALSE),]
write.csv(final_data_1,"input.csv",row.names = FALSE)
