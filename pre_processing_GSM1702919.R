library(GEOquery)

GPL10558=getGEO("GPL10558")
GPL10558_table = GPL10558@dataTable@table
##View(GPL10558_table)

GSM1702919=getGEO("GSM1702919")
GSM1702919 = Table(GSM1702919)
#View(GSM1702919)

x <- illuminaHumanv4SYMBOL
ID <-as.character(GPL10558_table$ID)
Symbol <- as.list(x[ID])
gpl10558 = cbind(Symbol,ID)
gpl10558 = as.data.frame(gpl10558)
#View(gpl10558)

eset_merged = merge.data.frame(x = GSM1702919,y = gpl10558,by.x = 'ID_REF',by.y = 'ID')

features = read.csv("comorbidity_features_list.csv",header = FALSE)
data_67_genes = eset_merged[eset_merged$Symbol %in% features$V1, ]

duplicated_id = as.data.frame(data_67_genes[which(duplicated(data_67_genes$Symbol)),])
unique_id = as.data.frame(data_67_genes[-which(duplicated(data_67_genes$Symbol)),])
duplicated_id$Symbol <- unlist(duplicated_id$Symbol)
unique_id$Symbol <- unlist(unique_id$Symbol)
for (i in 1:length(duplicated_id$Symbol)) 
{
  x = mean(data_67_genes[which(data_67_genes$Symbol == duplicated_id$Symbol[i]),2])
  unique_id[which(unique_id$Symbol==duplicated_id$Symbol[i]),2] = x
}
#View(unique_id)
unique_id_1 = unique_id[,c(3,2)]
colnames(unique_id_1) = c("genes","gsm1702919")

genes = features$V1[-which(features$V1 %in% unique_id$Symbol)]
gsm1702919 = rep(0,length(genes))

miss_genes = data.frame(genes,gsm1702919)
final_data = rbind(unique_id_1, miss_genes)

final_data_1 = final_data[order(final_data$genes,decreasing = FALSE),]
write.csv(final_data_1,"GSM1702919.csv",row.names = FALSE)
