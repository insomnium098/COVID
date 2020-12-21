library(GEOquery)
library(Biobase)

gset <- getGEO("GSE162113", GSEMatrix =TRUE)

###GPL21290 contains data for Homo Sapiens
###GPL21493 contains data fort Mus musculos


data_human <- gset$`GSE162113-GPL21290_series_matrix.txt.gz`

data_mouse <- gset$`GSE162113-GPL21493_series_matrix.txt.gz`


###PhenoData

phenodata_human <- pData(data_human)
phenodata_mouse <- pData(data_mouse)



##Expression data
##Expression is obtained from the provided matrix at GEO


