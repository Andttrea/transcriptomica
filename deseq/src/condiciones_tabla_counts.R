library(dplyr)

# Cargamos la tabla de metadatos 
# NOTA: Usamos check.names = FALSE porque los nombres de las columnas tienen caracteres especiales 
metadata_table <- read.csv("/export/storage/users/andreavg/transcriptomica/deseq/data/metadata/GSE132040_MACA_Bulk_metadata.csv", check.names = FALSE)

# Creamos una variable que contenga el nombre de las muestras, el sexo y la edad. 
metadata <- metadata_table %>%
  # Seleccionamos solo las columnas que nos interesan, que son el nombre de la muestra, el sexo y la edad 
  select(srr_id = "raw file", edad = "characteristics: age", sexo = "characteristics: sex") %>%
  # Creamos el nombre de la condición concatenando el sexo y la edad
  mutate(condicion = paste0("m",sexo, "_", edad))

# Hacemos una función que va a asignarle a cada muestra su condición correspondiente en hisat2 y star
asignar_condicion <- function(path_archivo, metadata) {
  # Leemos la tabla de conteos con el nombre de las muestras
  tabla_counts <- read.table(path_archivo, header = TRUE, sep = "\t", check.names = FALSE)
  # Extraemos el nombre de las columnas , exceptuando la primera columna que es gene_id
  srr_names <- colnames(tabla_counts)[-1]
  # Asignamos la condición a cada muestra usando la función match para encontrar la posición de cada nombre de muestra en la tabla de metadatos 
  condition_names <- metadata$condicion[match(srr_names, metadata$srr_id)] 
  # Reemplazamos los nombres de las columnas en la tabla de conteos por las condiciones correspondientes  
  colnames(tabla_counts) <- c("gene_id", make.unique(condition_names, sep = "_rep"))
  
  # Retornamos la tabla de conteos con los nombres de las columnas actualizados a las condiciones correspondientes
  return(tabla_counts)
}

# Hacemos la función para asignar la condición a cada muestra en el objeto txi de salmon
asignar_condicion_txi <- function(path_archivo, metadata) {
  # Leemos el objeto txi generado por Salmon para los datos de paired-end o single-end
  txi <- readRDS(path_archivo)
  # Quitamos las versiones decimales de los rownames (ENSMUSG00001.5 -> ENSMUSG00001)
  # Lo hacemos en todas las matrices para mantener la consistencia
  rownames(txi$counts) <- sub("\\..*", "", rownames(txi$counts))
  rownames(txi$abundance) <- sub("\\..*", "", rownames(txi$abundance))
  rownames(txi$length) <- sub("\\..*", "", rownames(txi$length)) 
  
  # Obtenemos los nombres actuales (SRRs)
  srr_names <- colnames(txi$counts)
  # Asignamos la condición a cada muestra usando la función match para encontrar la posición de cada nombre de muestra en la tabla de metadatos
  condition_names <- metadata$condicion[match(srr_names, metadata$srr_id)]
  # Creamos nombres únicos para las réplicas
  new_names <- make.unique(condition_names, sep = "_rep")
  
  # Reemplazamos nombres en todas las matrices del objeto
  colnames(txi$counts)    <- new_names
  colnames(txi$abundance) <- new_names
  colnames(txi$length)    <- new_names
  # Retornamos el objeto txi con los nombres de las columnas actualizados a las condiciones correspondientes
  return(txi)
}

# Vamos a generar la tabla de counts para hisat2 y star, para paired-end y single-end, usando la función creada anteriormente
tabla_condicion_hisat_pe <- asignar_condicion("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_counts.tsv", metadata)
tabla_condicion_hisat_se <- asignar_condicion("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_counts.tsv", metadata)
tabla_condicion_star_pe <- asignar_condicion("/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_counts.tsv", metadata)
tabla_condicion_star_se <- asignar_condicion("/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_counts.tsv", metadata)

# Para salmon
txi_condicion_pe <- asignar_condicion_txi("/export/storage/users/andreavg/transcriptomica/deseq/results/salmon/paired_end/txi/txi_PE.rds", metadata)
txi_condicion_se <- asignar_condicion_txi("/export/storage/users/andreavg/transcriptomica/deseq/results/salmon/single_end/txi/txi_SE.rds", metadata)

# Guardamos las tablas de conteos con las condiciones asignadas en formato tsv 
write.table(tabla_condicion_hisat_pe, "/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_counts_condicion.tsv", sep="\t", quote=F, row.names=F)
write.table(tabla_condicion_hisat_se, "/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_counts_condicion.tsv", sep="\t", quote=F, row.names=F)
write.table(tabla_condicion_star_pe, "/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_counts_condicion.tsv", sep="\t", quote=F, row.names=F)
write.table(tabla_condicion_star_se, "/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_counts_condicion.tsv", sep="\t", quote=F, row.names=F)

# Guardamos los objetos txi con las condiciones asignadas en formato rds
saveRDS(txi_condicion_pe, "/export/storage/users/andreavg/transcriptomica/deseq/results/salmon/paired_end/txi/txi_PE_condicion.rds")
saveRDS(txi_condicion_se, "/export/storage/users/andreavg/transcriptomica/deseq/results/salmon/single_end/txi/txi_SE_condicion.rds")