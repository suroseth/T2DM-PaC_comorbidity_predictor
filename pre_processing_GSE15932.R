library(GEOquery)
library(limma)
library(umap)
library(WGCNA)
library(lumi)
library(affy)
library(DExMA)
library(DExMAdata)
library(dplyr)
library(tidyr)
library(stringr)
library(hgu133plus2.db)

#reading affy data (.cel files) in R
eset=ReadAffy(celfile.path = "./GSE15932_RAW")
length(eset)

#reading annotation file in R
gp570 = getGEO("GPL570", destdir = ".")
gp570 = Table(gp570)


#reannotating probes to latest ID
x <- hgu133plus2SYMBOL
# Get the probe identifiers that are mapped to a gene symbol
ID <- gp570$ID
# Convert to a list
Symbol <- as.list(x[ID])
gpl570 = cbind(Symbol,ID)
gpl570 = as.data.frame(gpl570)
View(gpl570)
#gpl570$Symbol[16] = "NA2"

#backgroud correction/normalization/log transformation
eset=rma(eset,normalize = TRUE,background = TRUE)
exp_mat = exprs(eset)#extracting expression matrix

#arranging gene ID accoring to probe ID in expression matrix and taking average of duplicate
aa=unlist(lapply(rownames(exp_mat),function(x){z=which(x==gpl570$ID);return(z)}));
gene_symbol=gpl570$Symbol[aa];
tempa=collapseRows(exp_mat, rowGroup=gene_symbol,rowID=rownames(exp_mat), method="Average")
avg_exp_mat=tempa$datETcollapsed
View(avg_exp_mat)

#designing case and control
gsms <- "1111111100000000"
sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
eset_f = avg_exp_mat[ ,sel]#extracting required samples
View(eset_f)

eset_gse15932  = as.matrix(eset_f)
colnames(eset_gse15932)
col_name = c( "GSM399666","GSM399667","GSM399669","GSM399676","GSM399679",
              "GSM399682","GSM399685","GSM399689","GSM399717","GSM399718",
              "GSM399719","GSM399761","GSM399762","GSM399763","GSM399765",
              "GSM399766") 
colnames(eset_gse15932) = col_name
View(eset_gse15932)
write.csv(eset_gse15932,"result_preprocessing_gse15932_3.csv")
