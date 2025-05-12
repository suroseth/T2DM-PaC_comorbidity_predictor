#ComBat
library(sva)

genes_67_features = read.csv("comorbidity_features.csv")

eset_gse15932 = read.csv("result_t2dpac_healthy_GSE15932_preprocessing.csv",row.names = 1)
eset_gse74629 = read.csv("result_t2dpac_healthy_GSE74629_preprocessing.csv",row.names = 1)

t2dpac_sample_id = c(colnames(eset_gse15932)[1:8],colnames(eset_gse74629)[1:14])
healthy_sample_id = c(colnames(eset_gse15932)[9:16],colnames(eset_gse74629)[15:28])

# Load your data (replace with your actual data)
data_platform1 <- as.data.frame(eset_gse15932)
data_platform2 <- as.data.frame(eset_gse74629)

x=intersect(rownames(eset_gse15932),rownames(eset_gse74629))
length(x)

xx1 = sapply(x,function(i){which(rownames(eset_gse15932) == i)})
data_platform1 = data_platform1[xx1,]
xx2 = sapply(x,function(i){which(rownames(eset_gse74629) == i)})
data_platform2 = data_platform2[xx2,]

identical(row.names(data_platform1),row.names(data_platform2))

combined_data = cbind(data_platform1,data_platform2)

batch <- factor(c(rep("Platform1", ncol(data_platform1)),
                  rep("Platform2", ncol(data_platform2))
))
#normalized_data_combat_param = read.csv("t2d_nor_param.csv",row.names = 1)
#normalized_data_combat_param$old_mu_avg7 = rowMeans(normalized_data_combat_param[,c(1,5,9,13,17,21)])
#normalized_data_combat_param$old_phi_avg7 = rowMeans(normalized_data_combat_param[,c(2,6,10,14,18,22)])
#normalized_data_combat_param$new_mu_avg7 = rowMeans(normalized_data_combat_param[,c(3,7,11,15,19,23)])
#normalized_data_combat_param$new_phi_avg7 = normalized_data_combat_param$new_phi_avg1
#1=184050, 2= 21321, 3=69528, 4=125158, 5=49641, 6=15932 7=avg

new_data_result_GSE74629 = adding_parameters(combined_data[,c(17:44)],
                                             normalized_data_combat_param,
                                             batch_number = 7,genes_67_features$Genes)
new_data_result_GSE15932 = adding_parameters(combined_data[,c(1:16)],
                                             normalized_data_combat_param,
                                             batch_number = 7,genes_67_features$Genes)
data_new = cbind(new_data_result_GSE15932,new_data_result_GSE74629)

#data_new_2 = data_new[which(rownames(data_new) %in% genes_67_features$Genes),]

#data_new_t2dpac = data_new_2[,c(1:8,17:30)]
data_new_t2dpac = data_new[,t2dpac_sample_id]
write.csv(data_new_t2dpac,"t2dpac_data_t2d_model_input_avg.csv")

#data_new_healthy = data_new_2[,c(9:16,31:44)]
data_new_healthy = data_new[,healthy_sample_id]
write.csv(data_new_healthy,"healthy_data_t2d_model_input_avg.csv")
