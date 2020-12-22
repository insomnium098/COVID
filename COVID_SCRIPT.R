library(GEOquery)
library(Biobase)
library(DESeq2)
library(dplyr)
library(org.Mm.eg.db)

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

write.csv(res_heart,"Diff_Exp/res_heart.csv")
write.csv(res_kidney,"Diff_Exp/res_kidney.csv")
write.csv(res_lung,"Diff_Exp/res_lung.csv")
write.csv(res_spleen,"Diff_Exp/res_spleen.csv")


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

datos_network <- read.delim("Raw_data/MouseNetV2_symbol.txt",
                            header = FALSE)

###Red de Heart con todos los diff exp en res_heart

red_heart <- filter(datos_network,
                    V1 %in% rownames(res_heart))
red_heart <- filter(red_heart, V2 %in% rownames(res_heart))

write.csv(red_heart,"Networks_files/All_genes/network_heart.csv")

###Red de Kidney con todos los diff exp en res_kidney

red_kidney <- filter(datos_network,
                    V1 %in% rownames(res_kidney))
red_kidney <- filter(red_kidney, V2 %in% rownames(res_kidney))

write.csv(red_kidney,"Networks_files/All_genes/network_kidney.csv")

###Red de Lung con todos los diff exp en res_lung

red_lung <- filter(datos_network,
                     V1 %in% rownames(res_lung))
red_lung <- filter(red_lung, V2 %in% rownames(res_lung))

write.csv(red_lung,"Networks_files/All_genes/network_lung.csv")

###Red de Spleen con todos los diff exp en res_spleen

red_spleen <- filter(datos_network,
                   V1 %in% rownames(res_spleen))
red_spleen <- filter(red_spleen, V2 %in% rownames(res_spleen))

write.csv(red_spleen,"Networks_files/All_genes/network_spleen.csv")






#############
#############
#############KEGG

#####Se obtendran solo aquellos genes relacionados con la respuesta inmune
###Primero se obtendran los genes de las vias de señalizacion de kegg
kegg_mouse <- kegg.gsets(species = "mmu", id.type = "entrez", check.new=FALSE)
kegg_pathways <- kegg_mouse$kg.sets

##Hacemos un dataframe para que sea mas facil buscar las pathways
df_pathways <- as.data.frame(unlist(kegg_pathways))
colnames(df_pathways) <- "Entrez"
df_pathways$Pathway <- rownames(df_pathways)

###Definimos las pathways relacionadas con la respuesta inmune,
###estas son obtenidas de: 
#https://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=mmu

vias_inmune <- c("Hematopoietic cell lineage",
                 "Complement and coagulation cascades",
                 "Platelet activation",
                 "Toll-like receptor signaling pathway",
                 "NOD-like receptor signaling pathway",
                 "RIG-I-like receptor signaling pathway",
                 "Cytosolic DNA-sensing pathway",
                 "C-type lectin receptor signaling pathway",
                 "Natural killer cell mediated cytotoxicity",
                 "Antigen processing and presentation",
                 "T cell receptor signaling pathway",
                 "Th1 and Th2 cell differentiation",
                 "Th17 cell differentiation",
                 "IL-17 signaling pathway",
                 "B cell receptor signaling pathway",
                 "Fc epsilon RI signaling pathway",
                 "Fc gamma R-mediated phagocytosis",
                 "Leukocyte transendothelial migration",
                 "Intestinal immune network for IgA production",
                 "Chemokine signaling pathway")

####Filtramos y solo nos quedamos con las vias inmunes
for (i in 1:length(vias_inmune)){
  via <- vias_inmune[i]
  index_via <- grep(via, df_pathways$Pathway)
  if(!exists("index_via_final")){
    index_via_final <- index_via
  } else {
    index_via_final <- c(index_via_final, index_via)
  }
}

vias_chidas <- df_pathways[index_via_final,]


###Los genes en vias_chidas estan en ENTREZID, hay que convertirlas a 
##symbol

resultados_symbol = mapIds(org.Mm.eg.db,
                           keys=vias_chidas$Entrez, 
                           column="SYMBOL",#En esta opcion se selecciona que valor se quiere obtener de la conversion
                           keytype="ENTREZID",#En esta opcion se selecciona en que formato estan los datos de origen
                           multiVals="first")

vias_chidas$SYMBOL <- resultados_symbol

###Limpiar vias_chidas para que la pathway este bien

for (i in 1:length(vias_inmune)){
  via <- vias_inmune[i]
  index_via_limpia <- grep(via, vias_chidas$Pathway)
  vias_chidas[index_via_limpia,2] <- via
}

###Guardamos las vias inmunes y sus genes
write.csv(vias_chidas,"KEGG_IMMUNE_PATHWAYS_AND_GENES.csv")

#######
######
#####Finalmente filtramos las redes de los organos para quedarnos
#####Solamente con genes que esten involucrados en vias inmunes

red_heart_inmune <- filter(red_heart, V1 %in% vias_chidas$SYMBOL)
colnames(red_heart_inmune) <- c("Source","Target","Score")
####Hacemos merge con las vias chidas para saber a que via
###pertenece cada source node

##Funciona pero se hacen repetidos ya que los genes source pertenecen
## a varias vias, de momento no se utilizara
#test <- merge(red_heart_inmune, vias_chidas,
#              by.x = "Source",
#              by.y = "SYMBOL")

write.csv(red_heart_inmune,"Networks_files/Only_Immune_genes/red_heart_inmune.csv")


##
red_kidney_inmune <- filter(red_kidney, V1 %in% vias_chidas$SYMBOL)
colnames(red_kidney_inmune) <- c("Source","Target","Score")
write.csv(red_kidney_inmune,"Networks_files/Only_Immune_genes/red_kidney_inmune.csv")
##
red_lung_inmune <- filter(red_lung, V1 %in% vias_chidas$SYMBOL)
colnames(red_lung_inmune) <- c("Source","Target","Score")
write.csv(red_lung_inmune,"Networks_files/Only_Immune_genes/red_lung_inmune.csv")
##
red_spleen_inmune <- filter(red_spleen, V1 %in% vias_chidas$SYMBOL)
colnames(red_spleen_inmune) <- c("Source","Target","Score")
write.csv(red_spleen_inmune,"Networks_files/Only_Immune_genes/red_spleen_inmune.csv")



