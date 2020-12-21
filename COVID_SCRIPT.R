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
datos_raton_dia3 <- read.delim("Raw_data/GSE162113_day3_raw_counts.txt")

###El chido es el dia 7
datos_raton_dia7 <- read.delim("Raw_data/GSE162113_day7_raw_counts.txt")

###Los controles son los eGFP

#######FUNCIONES AUXILIARES QUE RECIBEN EL DF Y  EL TEJIDO DE ORIGEN, DEVUELVEN
###UN DF CON EL CONTROL DEL TEJIDO DE ORIGEN Y EL CASO CON EL TEJIDO DE ORIGEN

##El tejido de origen puede ser: "Heart","Kidney","Lung","Spleen"


###Para los controles
aux_control <- function(datos, tejido_origen){
  
  index_tejido <- grep(tejido_origen, colnames(datos))
  index_control <- grep("eGFP", colnames(datos))
  
  index_final <- intersect(index_control, index_tejido)
  return(datos[,index_final])
}

###Para los casos

aux_casos <- function(datos, tejido_origen){
  
  index_tejido <- grep(tejido_origen, colnames(datos))
  index_casos <- grep("hACE2", colnames(datos))
  
  index_final <- intersect(index_casos, index_tejido)
  return(datos[,index_final])
}


###Test de las funciones auxiliares
##El tejido de origen puede ser: "Heart","Kidney","Lung","Spleen"
##NOT RUN
#test <- aux_control(datos_raton_dia7,"Kidney")
#test2 <- aux_casos(datos_raton_dia7,"Spleen")
