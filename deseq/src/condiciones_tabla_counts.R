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

# Hacemos una función que va a asignarle a cada muestra su condición correspondiente
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

# Vamos a generar la tabla de counts para hisat2 y star, para paired-end y single-end, usando la función creada anteriormente
tabla_condicion_hisat_pe <- asignar_condicion("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_counts.tsv", metadata)
tabla_condicion_hisat_se <- asignar_condicion("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_counts.tsv", metadata)
tabla_condicion_star_pe <- asignar_condicion("/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_counts.tsv", metadata)
tabla_condicion_star_se <- asignar_condicion("/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_counts.tsv", metadata)

# Guardamos las tablas de conteos con las condiciones asignadas en formato tsv 
write.table(tabla_condicion_hisat_pe, "/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_counts_condicion.tsv", sep="\t", quote=F, row.names=F)
write.table(tabla_condicion_hisat_se, "/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_counts_condicion.tsv", sep="\t", quote=F, row.names=F)
write.table(tabla_condicion_star_pe, "/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_counts_condicion.tsv", sep="\t", quote=F, row.names=F)
write.table(tabla_condicion_star_se, "/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_counts_condicion.tsv", sep="\t", quote=F, row.names=F)
