library(GEOquery)
library(Biobase)
library(DESeq2)
library(dplyr)

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


####Diff exp between organs

#########Heart
#########
#########
#########
heart_controles <- aux_control(datos_raton_dia7,"Heart")
heart_casos <- aux_casos(datos_raton_dia7,"Heart")

heart_counts <- as.matrix(cbind(heart_controles, heart_casos))

condition <- factor(c( rep(c("control","enfermos"), 
                           c(ncol(heart_controles),
                             ncol(heart_casos)))))
coldata <- data.frame(row.names=colnames(heart_counts), condition)
dds <- DESeqDataSetFromMatrix(countData=heart_counts,
                              colData=coldata, 
                              design=~condition)
dds$condition <- relevel(dds$condition, ref="control")
dds <- DESeq(dds, parallel = TRUE)

####Obtener todos los resultados y quedarse solo con los significativos
res_heart <- as.data.frame(results(dds))
res_heart <- filter(res_heart, padj < 0.05)

####################END HEART
####################
####################


################Kidney
controles <- aux_control(datos_raton_dia7,"Kidney")
casos <- aux_casos(datos_raton_dia7,"Kidney")

counts <- as.matrix(cbind(controles, casos))

condition <- factor(c( rep(c("control","enfermos"), 
                           c(ncol(controles),
                             ncol(casos)))))
coldata <- data.frame(row.names=colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=coldata, 
                              design=~condition)
dds$condition <- relevel(dds$condition, ref="control")
dds <- DESeq(dds, parallel = TRUE)

####Obtener todos los resultados y quedarse solo con los significativos
res <- as.data.frame(results(dds))
res_kidney <- filter(res, padj < 0.05)


####################END Kidney
####################
####################


################Lung
controles <- aux_control(datos_raton_dia7,"Lung")
casos <- aux_casos(datos_raton_dia7,"Lung")

counts <- as.matrix(cbind(controles, casos))

condition <- factor(c( rep(c("control","enfermos"), 
                           c(ncol(controles),
                             ncol(casos)))))
coldata <- data.frame(row.names=colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=coldata, 
                              design=~condition)
dds$condition <- relevel(dds$condition, ref="control")
dds <- DESeq(dds, parallel = TRUE)

####Obtener todos los resultados y quedarse solo con los significativos
res <- as.data.frame(results(dds))
res_lung <- filter(res, padj < 0.05)


####################END Lung
####################
####################

################Spleen
controles <- aux_control(datos_raton_dia7,"Spleen")
casos <- aux_casos(datos_raton_dia7,"Spleen")

counts <- as.matrix(cbind(controles, casos))

condition <- factor(c( rep(c("control","enfermos"), 
                           c(ncol(controles),
                             ncol(casos)))))
coldata <- data.frame(row.names=colnames(counts), condition)
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=coldata, 
                              design=~condition)
dds$condition <- relevel(dds$condition, ref="control")
dds <- DESeq(dds, parallel = TRUE)

####Obtener todos los resultados y quedarse solo con los significativos
res <- as.data.frame(results(dds))
res_spleen <- filter(res, padj < 0.05)


####################END Spleen
####################
####################


####Guardar resultados de diff exp en csv

write.csv(res_heart,"res_heart.csv")
write.csv(res_kidney,"res_kidney.csv")
write.csv(res_lung,"res_lung.csv")
write.csv(res_spleen,"res_spleen.csv")


###Venn diagrams will be plotted using Venny 2.1 
###(https://bioinfogp.cnb.csic.es/tools/venny/)


#####Transcritos compartidos en todos los organos

t_compartidos <- Reduce(intersect, list(rownames(res_heart),
                       rownames(res_kidney),
                       rownames(res_lung),
                       rownames(res_spleen)))

#####Transcritos unicos de Heart
heart_unicos <- setdiff(rownames(res_heart),
                        c(rownames(res_kidney),
                          rownames(res_lung),
                          rownames(res_spleen)))
#####Transcritos unicos de Kidney
kidney_unicos <- setdiff(rownames(res_kidney),
                        c(rownames(res_heart),
                          rownames(res_lung),
                          rownames(res_spleen)))

#####Transcritos unicos de Lung
lung_unicos <- setdiff(rownames(res_lung),
                         c(rownames(res_heart),
                           rownames(res_kidney),
                           rownames(res_spleen)))

#####Transcritos unicos de Spleen
spleen_unicos <- setdiff(rownames(res_spleen),
                       c(rownames(res_heart),
                         rownames(res_kidney),
                         rownames(res_lung)))


####La red se hará con los datos de MouseNet V2
##https://www.inetbio.org/mousenet/



