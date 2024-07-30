###########################################################################################
# To normalize T2DM-PaC comorbidity
# Input: given preprocessed sample files or your preprocessed gene expression data in CSV format
# OUTPUT: normalized gene expression values for 67 gene comorbidity features.  
###########################################################################################

library(sva)
input_data =  "input_file_name.csv"
if(length(input_data) >= 1)
{
  eset_gse15932 = read.csv("result_preprocessing_gse15932_3.csv",row.names = 1)
  eset_gse74629 = read.csv("result_preprocessing_GSE74629_3.csv",row.names = 1)
  input_data = read.csv("input_data.csv",row.names = 1)
  
  # Load your data (replace with your actual data)
  data_platform1 <- as.data.frame(eset_gse15932)
  data_platform2 <- as.data.frame(eset_gse74629)
  input_platform = as.data.frame(input_data)
  
  x=intersect(rownames(eset_gse15932),rownames(eset_gse74629))
  x1 = intersect(x,rownames(input_data))
  
  xx = sapply(x1,function(i){which(rownames(eset_gse15932) == i)})
  data_platform1 = data_platform1[xx,]
  
  xx2 = sapply(x1,function(i){which(rownames(eset_gse74629) == i)})
  data_platform2 = data_platform2[xx2,]

  xx3 = sapply(x1,function(i){which(rownames(input_data) == i)})
  data_platform3 = data_platform2[xx3,]
  
  identical(row.names(data_platform1),row.names(data_platform2))
  combined_data = cbind(data_platform1,data_platform2,input_data)
  
  #combined_data_x = merge.data.frame(data_platform1,data_platform2)
  batch <- factor(c(rep("Platform1", ncol(data_platform1)),
                    rep("Platform2", ncol(data_platform2)),
                    rep("Platform2", ncol(data_platform3))))
  normalized_data_combat <- ComBat(combined_data, 
                                   batch = batch, 
                                   mod = NULL, 
                                   par.prior = TRUE, 
                                   prior.plots = FALSE)
}else
  {
  eset_gse15932 = read.csv("result_preprocessing_gse15932_3.csv",row.names = 1)
  eset_gse74629 = read.csv("result_preprocessing_GSE74629_3.csv",row.names = 1)
  
  # Load your data (replace with your actual data)
  data_platform1 <- as.data.frame(eset_gse15932)
  data_platform2 <- as.data.frame(eset_gse74629)
  
  x=intersect(rownames(eset_gse15932),rownames(eset_gse74629))
  xx = sapply(x,function(i){which(rownames(eset_gse15932) == i)})
  data_platform1 = data_platform1[xx,]
  
  xx2 = sapply(x,function(i){which(rownames(eset_gse74629) == i)})
  data_platform2 = data_platform2[xx2,]
  
  identical(row.names(data_platform1),row.names(data_platform2))
  combined_data = cbind(data_platform1,data_platform2)
  
  #combined_data_x = merge.data.frame(data_platform1,data_platform2)
  batch <- factor(c(rep("Platform1", ncol(data_platform1)), rep("Platform2", ncol(data_platform2))))
  normalized_data_combat <- ComBat(combined_data, 
                                   batch = batch, 
                                   mod = NULL, 
                                   par.prior = TRUE, 
                                   prior.plots = FALSE)
  
}
features = read.csv("comorbidity_features_list.csv",header = FALSE)

data_67_genes = normalized_data_combat[row.names(normalized_data_combat) %in% features$V1, ]
if(length(data_67_genes) != length(features$V1))
{print("Input data does not have all the 67 gene featurres")  }
data_67_genes = cbind(row.names(data_67_genes),data_67_genes)
colnames(data_67_genes)[1] = "genes"
write.csv(data_67_genes,"normalized_data_combat_67features.csv",row.names = FALSE)
