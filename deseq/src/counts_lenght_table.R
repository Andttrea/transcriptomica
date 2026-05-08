library(dplyr)
library(purrr)

# Se creó una función para crear la tabla de conteos y longitudes
# NOTA: Recibira como argumentos la ruta de la carpeta y el patrón de los archivos a leer
counts_lenght_table <- function(folder_path, file_pattern) {
  # Se obtuvo la lista de archivos que coinciden con el patrón
  # NOTA: full.names = TRUE para obtener la ruta completa de los archivos y poder leerlos posteriormente 
  archivos <- list.files(folder_path, pattern = file_pattern, full.names = TRUE)
  # Se leyó cada archivo, saltandonos la primera fila que es de comentario 
  # NOTA: se usó map de purrr para aplicar la función read.table a cada archivo de la lista
  # El .x representa cada elemento de la lista de archivos, es decir, cada ruta de archivo y skip = 1 para saltar la primera fila de cada archivo 
  # check.names = FALSE para evitar que R modifique los nombres de las columnas
  lista_featurecounts <- map(archivos, ~read.table(.x, header = TRUE, skip = 1, check.names = FALSE))

  # Se crea la tabla de longitudes. En este caso, solo tomaremos la primer tabla, debido a que como corresponde al genoma del raton todas las longitudes son iguales
  tabla_length <- lista_featurecounts[[1]] %>%
    # Quitaremos el decimal de los Geneid, por ejemplo: ENSMUSG00000102693.2 -> ENSMUSG00000102693. Esto para que sea más fácil de manejar posteriormente.
    mutate(gene_id = sub("\\..*", "", Geneid)) %>%
    # Se seleccionó solo las columnas Geneid y Length
    select(gene_id, Length)
   
  # Se creó la tabla de conteos.
  tabla_counts <- lista_featurecounts %>%
    # Se seleccionan solo las columnas Geneid y la última columna que corresponde a los conteos
    map(~select(.x, Geneid, last_col())) %>%
    # Se unen todas las tablas de conteos en una sola tabla, usando left_join para unirlas por la columna Geneid
    reduce(left_join, by = "Geneid") %>%
    # Quitamos el decimal 
    mutate(gene_id = sub("\\..*", "", Geneid)) %>% 
    # Seleccionamos la columna gene_id y todas las demás, excluyendo Geneid que tenia el decimal 
    select(gene_id, everything(), -Geneid)
    
  # Limpiamos el nombre de las columnas
  #  Usamos basename para quitar toda la ruta del servidor
  colnames(tabla_counts) <- basename(colnames(tabla_counts))

  # Quitamos los sufijos para dejar el nombre de la muestra únicamente 
  # NOTA: Utilzamos gsub para eliminar las partes no deseadas de los nombres de las columnas, como "_Aligned.out", "_fc_counts.txt" y sufijos para quedarnos solo con el nombre de la muestra 
  colnames(tabla_counts) <- gsub("_Aligned.out.*|_fc_counts\\.txt|_SE.*|_PE.*|\\.bam", "", colnames(tabla_counts)) 

  # Retonamos una lista con la tabla de conteos y la tabla de longitudes 
  return(list(counts = tabla_counts, lengths = tabla_length))

}

# Vamos a generar la tabla de counts para hisat2 y star, para paired-end y single-end, usando la función creada anteriormente
# NOTA: Se le pasa la ruta de la carpeta donde se encuentran los archivos de conteos y el patrón de los archivos a leer, que en este caso son los archivos que terminan con "_fc_counts.txt"
hisat_pe_tables <- counts_lenght_table("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts", "SRR.*_fc_counts\\.txt$")
hisat_se_tables <- counts_lenght_table("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts", "SRR.*_fc_counts\\.txt$")
star_pe_tables <- counts_lenght_table("/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts", "SRR.*_fc_counts\\.txt$")
star_se_tables <- counts_lenght_table("/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts", "SRR.*_fc_counts\\.txt$")

# Se guardan las tablas de conteos y longitudes en formato tsv para su posterior uso en el análisis de DESeq2
# Para la tabla de counts
write.table(hisat_pe_tables$counts, "/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_counts.tsv", sep="\t", quote=F, row.names=F)
write.table(hisat_se_tables$counts, "/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_counts.tsv", sep="\t", quote=F, row.names=F)
write.table(star_pe_tables$counts, "/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_counts.tsv", sep="\t", quote=F, row.names=F)
write.table(star_se_tables$counts, "/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_counts.tsv", sep="\t", quote=F, row.names=F)

# Para la tabla de longitudes
write.table(hisat_pe_tables$lengths, "/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_lengths.tsv", sep="\t", quote=F, row.names=F)
write.table(hisat_se_tables$lengths, "/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_lengths.tsv", sep="\t", quote=F, row.names=F)
write.table(star_pe_tables$lengths, "/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_lengths.tsv", sep="\t", quote=F, row.names=F)
write.table(star_se_tables$lengths, "/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_lengths.tsv", sep="\t", quote=F, row.names=F)