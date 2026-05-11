# Reporte de expresión diferencial

## Planteamiento

El análisis de expresión diferencial permite identificar genes cuya abundancia cambia entre condiciones biológicas distintas. En este trabajo se comparó tejido muscular de *Mus musculus* de 9 meses contra 24 meses, utilizando como punto de partida los archivos de cuantificación generados para tres estrategias de análisis: alineamiento con HISAT2, alineamiento con STAR y pseudoalineamiento con Salmon, los cuales fueron generados anteriormente. 

Nuestro enfoque es preparar correctamente las matrices de conteo y aplicar dos enfoques estadísticos para detectar genes diferencialmente expresados: DESeq2 y edgeR. El objetivo fue evaluar cómo cambian los resultados según el alineador, el modo de secuenciación (single-end o paired-end) y el algoritmo estadístico utilizado.

Los datos, scripts y resultados se encuentran en el directorio `/export/storage/users/andreavg/transcriptomica/deseq`.

## Objetivos

### Objetivo general

Comparar los resultados de expresión diferencial obtenidos a partir de HISAT2, STAR y Salmon mediante los paquetes DESeq2 y edgeR, usando datos de RNA-seq de músculo de ratón en dos edades distintas (9 meses y 24 meses)

### Objetivos particulares

1. Integrar los archivos de `featureCounts` en una matriz global de conteos y una tabla de longitudes por gen para HISAT2 y STAR.
2. Estandarizar los nombres de las muestras sustituyendo los identificadores SRR por sus condiciones biológicas para facilitar la interpretación de los resultados.
4. Ejecutar análisis de expresión diferencial con DESeq2 y edgeR para los datos de HISAT2, STAR y Salmon en modos single-end y paired-end.
5. Comparar los resultados obtenidos mediante PCA, volcano plots y heatmaps de genes significativos.


## Metodología y justificación

### Preparación de tablas de conteo

El primer paso consistió en transformar la salida de `featureCounts`, la cual originalmente estaba separada por muestra, en tablas globales con todas las muestras reunidas. Esto se realizó en el script `counts_lenght_table.R`, donde se unieron los conteos por gen y también se extrajó las longitudes de los genes. Este paso fue necesario para HISAT2 y STAR, ya que para Salmon se trabajó con el objeto `txi` generado previamente por `tximport`. 

Aunque Salmon cuantifica a nivel de transcritos, en este análisis se trabajó a nivel de gen gracias a `tximport`. El archivo `.rds` generado por este flujo ya resuelve el problema de múltiples transcritos por gen, porque la tabla `tx2gene` le indica a `tximport` qué transcritos pertenecen a cada gen. Con ello, los conteos de todas las isoformas se suman en un único registro por GENEID. Como resultado, las matrices de `counts` y `abundance` dentro del objeto `txi` contienen filas identificadas por genes y no por transcritos individuales. Además, `tximport` estima una longitud efectiva por gen considerando qué isoformas se expresaron realmente, lo cual es importante porque la longitud aparente de un gen puede cambiar según las isoformas activas. Por esto, el objeto `txi` es una entrada válida y consistente para el análisis diferencial a nivel de gen.

```bash
# Nos movemos a la carpeta de src
cd /export/storage/users/andreavg/transcriptomica/deseq/src

# Hacemos el script
nano counts_lenght_table.R
```

```r
# ----------------------------------------------- EMPIEZA EL SCRIPT ------------------------------------------------- #

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

# ----------------------------------------------- TERMINA EL SCRIPT ------------------------------------------------- #

```

Posteriormente, el script `condiciones_tabla_counts.R` reemplazó los nombres SRR por condiciones más interpretables, usando la tabla de metadatos del experimento. Para cada muestra se construyó una etiqueta con la edad y el sexo, por ejemplo `mm_9` y `mm_24`, y se añadieron sufijos automáticos para distinguir réplicas.

```bash
nano condiciones_tabla_counts.R
```

```r
# ----------------------------------------------- EMPIEZA EL SCRIPT ------------------------------------------------- #

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

# ----------------------------------------------- TERMINA EL SCRIPT ------------------------------------------------- #
```
```bash
# Se activa el entorno que tiene las librerias
conda activate .main_mamba

# Se corren los scripts
Rscript counts_lenght_table.R
Rscript condiciones_tabla_counts.R
```

### Análisis con DESeq2

El script `Deseq2_hisat_star.R` realizó el análisis para HISAT2 y STAR, mientras que `Deseq2_salmon.R` hizo lo propio para Salmon. En ambos casos el flujo fue similar:

- creación del objeto de DESeq2 desde una matriz de conteos o desde `tximport`,
- filtrado de genes de baja expresión con `filterByExpr`,
- transformación `vst` para construir un PCA,
- ajuste del modelo con `DESeq()`,
- contraste entre `mm_24` y `mm_9`,
- cálculo de una tabla de resultados con log2 fold change y FDR
- generación de volcano plots y heatmaps para los genes significativos.

En estos scripts se usaron umbrales de significancia de `FDR < 0.05` y `|log2FC| > 0.5`. Los genes up y down regulados se guardaron en archivos separados, junto con mapas de TPM en escala log2 para apoyar la visualización.

```bash
# En la carpeta src creamos el script
nano Deseq2_hisat_star.R
```
```r
# ----------------------------------------------- EMPIEZA EL SCRIPT ------------------------------------------------- #

# Cargamos las librerias que utilizaremos 
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(edgeR)

# Hacemos una funsión que realice todo el análisis de DESeq2

Deseq2_analysis <- function(gene_counts, annotation, gene_name_map, output_dir) {

  # Creamos la ruta dinámica de salida para guardar los resultados de cada análisis
  fig_dir <- paste0(output_dir, "figuras/")

  # Para generar la tabla de metadatos primero generamos factores
  # Generamos factores para la edad y el sexo 
  age <- factor(c("mm_24", "mm_9", "mm_24", "mm_9", "mm_24", rep("mm_9", 2)), levels = c("mm_9", "mm_24"))
  sex <- factor(c(rep("m", 7)))
  # Asignamos a la variable sample_names los nombres de las columnas de la tabla de conteos
  sample_names <- colnames(gene_counts)

  # Generamos la tabla de metadatos
  meta_data <- data.frame(sample_names, age, sex)
  # ELiminamos los nombres de las filas y asignamos la columna sample_names como nombres de fila 
  meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var = "sample_names")

  # Checaremos que los nombres de la tabla de metadatos estan en el mismo orden que las columnas de los gene_counts
  k <- all(colnames(gene_counts) == rownames(meta_data))
  # Vemos que si arroja TRUE, entonces esta todo correcto
  print(paste0("Los nombres de las columnas de gene_counts y los nombres de fila de meta_data coinciden: ", k))

  # Cargamos los datos en un objeto DESeq2
  # NOTA: Usamos round() para redondear los valores de conteo a enteros, ya que DESeq2 requiere conteos enteros. 
  # Design lo dejamos como ~ 0 + age para que el modelo no tenga intercepto y podamos comparar directamente las condiciones de edad.
  dds <- DESeqDataSetFromMatrix(countData = round(gene_counts), colData = meta_data, design = ~ 0 + age)
  design <- model.matrix(~ 0 + age)

  # FIltramos los genes que contienen una baja expresión, esto que lo haremos con filterByExpr() de edgeR, que nos permite filtrar los genes que no tienen una expresión suficiente para ser considerados en el análisis diferencial.
  # filterByExpr compara la expresión de cada gen con un umbral de expresión mínima para eliminar los genes que no tienen una expresión suficiente 
  keep <- filterByExpr(dds,design)
  # Vemos cuantos genes se mantienen después del filtrado
  suma_keep <- sum(keep)
  # Guardamos el número de genes que se mantienen después del filtrado en un archivo de texto
  write.table(suma_keep, paste0(output_dir, "genes_filtrados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  # Filtramos el objeto DESeq2 para quedarnos solo con los genes que cumplen con el criterio de expresión suficiente
  dds <- dds[keep,]
  # Borramos keep 
  rm(keep)

  # Generamos el PCA para visualizar la variabilidad de los datos 
  # Primero convertimos los datos a vst, el cual minimiza la variabilidad de los datos y hace que sean más adecuados para el PCA.
  vsd <- vst(dds)
  # Hacemos el PCA, definiendo intgroup = "age" para que el PCA se coloree por la edad
  PCA_plot <- plotPCA(vsd, intgroup = "age") + theme_classic(base_size=25, base_line_size = 1)

  # Guardamos la imagen del PCA
  ggsave(paste0(fig_dir, "PCA_plot.png"), plot = PCA_plot, width = 8, height = 6)

  # Calculamos los factores de normalización, varianza y ajustar el modelo 
  # NOTA: Utilizamos DESeq() para calcular los factores de normalización, la varianza y ajustar el modelo a los datos.
  dds <- DESeq(dds)

  # Ahora calcularemos la TPM (Transcripts Per Million) que calcula cuantos transcritos tendrias si tuvieras un millon de transcritos en total
  # Añadimos la longitud del gen
  # NOTA: mcols es una función que permite añadir metadatos a un objeto, usamos rownames(dds) para asegurarnos que la logitud coincida con el gen correcto
  mcols(dds)$basepairs = annotation[rownames(dds),]

  # Calculamos FPKM y pasamos a log2. FKPM (Fragments Per Kilobase of transcript per Million mapped reads) es una medida de expresión que normaliza por la longitud del gen y por el número total de lecturas mapeadas
  # NOTA: Le añadimos un pseudoconteo de 0.1 para evitar problemas con los logaritmos de cero
  log2_fpkm <- log2(fpkm(dds) + 0.1) 

  # Escribimos la formula para convertir de FKPM a TPM dentro de un epsacio lograrítmico 
  fpkm2tpm_log2 <- function(fpkm) { fpkm - log2(sum(2^fpkm)) + log2(1e6) } 
  # Aplicamos la formula a cada columna de la tabla
  log2_tpm <- apply(log2_fpkm, 2, fpkm2tpm_log2) 

  # Guardamos la tabla 
  gene_names <- gene_name_map[rownames(log2_tpm),]
  write.table(cbind(gene_names, log2_tpm), paste0(output_dir, "TPM_log2-table.txt"), sep="\t", quote=FALSE)

  # Hacemos el  contraste de expresión diferencial 
  # Checamos los nombres 
  resultsNames(dds)
  # Output: 
  # [1] "agemm_9"  "agemm_24"
  # Hacemos el contraste entre nuestras condiciones, en este caso entre mm_24 y mm_9
  contrast <- makeContrasts(m24_vs_m9 = agemm_24 - agemm_9, levels = design) 

  # Hacemos el análisis de expresión diferencial 
  # Extraemos los resultados y los nombres 
  res <- results(dds, contrast=contrast[, "m24_vs_m9"])
  # Añadimos el nombre del gen a la tabla de resultados
  res$Gene_name <- gene_name_map[rownames(res),]

  # Establecemos los tresholds para LFC y FDR
  FDR <- 0.05 # Escogemos un FDR de 0.05 para ser un poco más permisivos 
  LFC <- 0.5

  # Filtramos los resultados para sacar el número de genes up y down regulados, usando los tresholds establecidos
  up <- (res$log2FoldChange > LFC) & (res$padj < FDR)
  # Reemplazamos los valores NA por FALSE, ya que los genes que no cumplen con los criterios de significancia 
  up[which(is.na(up))] = FALSE
  # Guardamos el número de genes up regulados en un archivo de texto
  suma_up <- sum(up)
  write.table(suma_up, paste0(output_dir, "genes_up_regulados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Hacemos lo mismo para los genes down regulados
  down <- (res$log2FoldChange < -LFC) & (res$padj < FDR)
  down[which(is.na(down))] = FALSE
  # Guardamos el número de genes down regulados en un archivo de texto
  suma_down <- sum(down)
  write.table(suma_down, paste0(output_dir, "genes_down_regulados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Guardamos las tablas 
  write.table(res[up,], paste0(output_dir, "deseq-DEG_up_0.05.txt"), sep="\t", quote=FALSE, row.names=TRUE)
  write.table(res[down,], paste0(output_dir, "deseq-DEG_down_0.05.txt"), sep="\t", quote=FALSE, row.names=TRUE)

  # Ahora procederemos a graficar con un volcanoplot
  # Asignamos los colores para las categorías
  vpcolors = c("gray", "#de0e99", "#c17803") 
  names(vpcolors) = c("NO", "DOWN", "UP") 

  # Creamos la columna DE en tu objeto de resultados 'res'
  # Primero marcamos todos como "NO"
  res$DE = "NO" 
  # Luego usamos tus vectores 'up' y 'down' para etiquetar los significativos
  res[up, "DE"] = "UP"
  res[down, "DE"] = "DOWN"

  # Creamos la gráfica con ggplot2
  volcano_plot <- ggplot(data = as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), col = DE)) +
          geom_point(alpha = 0.4, size = 1.5) + 
          labs(title = "Volcano plot: Músculo 24m vs 9m (Escala Auto)", 
               x = "log2 Expression fold change", 
               y = "-log10 FDR") + 
          scale_color_manual(values = vpcolors) +
          geom_vline(xintercept = c(-LFC, LFC), col = "black", linetype = "longdash") +
          geom_hline(yintercept = -log10(FDR), col = "black", linetype = "longdash") +
          theme_classic(base_size = 15, base_line_size = 1)

  
  # Guardamos la imagen del volcano plot
  ggsave(paste0(fig_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)

  # Haremos el heatmap de los genes diferencialmente expresados
  # Primero con los genes más significativos 
  # Filtramos solo los genes que pasan los filtros (Up y Down)
  significant <- res[up | down, ]

  # Ordenamos los genes por su significancia 
  significant_order <- significant[order(significant$padj), ]

  # Tomamos los 2,000 genes más significativos 
  # Usamos rownames para obtener los IDs de Ensembl de esos genes top, dejamos el valor n = 2000 para tener un umbral alto y no estar ajustando entre disitintos alineadores
  top_genes <- head(rownames(significant_order), n = 2000)

  # Calculamos el Z-score a partir de tus datos TPM
  # El Z-score centra la expresión de cada gen: 0 es el promedio, valores positivos son arriba del promedio, negativos abajo.
  zscore_t <- t(scale(t(log2_tpm[top_genes, ])))

  # Ordenamos las columnas para que en el heatmap se muestren en el orden de edad, primero los mm_9 y luego los mm_24
  orden_columnas <- order(meta_data$age)
  zscore_significant <- zscore_t[, orden_columnas]

  # Creamos el Heatmap con ComplexHeatmap
  # km = 2 le pide a R que intente separar los genes en 2 grupos principales (clusters)
  heatmap_plot <- Heatmap(zscore_significant, 
              cluster_rows = T,         
              cluster_columns = F,      # F para que respete el orden de edad que pusimos
              show_row_names = F,       
              name = "Z-score", 
              km = 2,                   
              column_title = "Heatmap de los top genes",
              column_names_gp = gpar(
                  col = "black", 
                  fontsize = 10, 
                  fontface = "bold"
              ),
              top_annotation = HeatmapAnnotation(
                  Edad = meta_data$age[orden_columnas],
                  col = list(Edad = c("mm_9" = "lightpink", "mm_24" = "deeppink"))
              )) 

  
  # Guardamos la imagen del heatmap
  png(paste0(fig_dir, "heatmap_top_genes.png"), width = 11.25, height = 7.5, res = 300, units = "in")
  # Es importante usar draw() para asegurar que todos los elementos se guarden
  draw(heatmap_plot)
  # Cerramos el archivo
  dev.off()


  # Ahora haremos el heatmap para el top 20 genes más significativos
  # Seleccionamos los 20 mejores (los de menor padj)
  top_20_ids <- head(rownames(significant_order), n = 20)

  # Calculamos el Z-score solo para estos 20 genes
  zscore_top <- t(scale(t(log2_tpm[top_20_ids, ])))

  # Aplicamos el orden de las columnas por edad (9m primero, luego 24m)
  zscore_top_ordenado <- zscore_top[, orden_columnas]

  # Creamos el Heatmap con los nombres de los genes visibles
  heatmap_top20 <- Heatmap(zscore_top_ordenado, 
              cluster_rows = T, 
              cluster_columns = F, 
              row_labels = gene_name_map[rownames(zscore_top_ordenado), ], 
              name = "Z-score", 
              km = 2, 
              column_title = "Top 20 Genes Significativos",
              column_names_gp = gpar(
                  col = "black", 
                  fontsize = 10, 
                  fontface = "bold"
              ),
              row_names_gp = gpar(
                  fontsize = 10, 
                  fontface = "italic" # Los nombres de genes siempre van en cursivas
              ),
              top_annotation = HeatmapAnnotation(
                  Edad = meta_data$age[orden_columnas],
                  col = list(Edad = c("mm_9" = "lightpink", "mm_24" = "deeppink"))
              ))

  
  # Guardamos la imagen del heatmap de los top 20 genes
  png(paste0(fig_dir, "heatmap_top_20_genes.png"),  width = 11.25, height = 7.5, res = 300, units = "in")
  # Es importante usar draw() para asegurar que todos los elementos se guarden
  draw(heatmap_top20)
  # Cerramos el archivo
  dev.off()
} 

# Cargamos los archivos que utilizaremos 
# NOTA: Usamos row.names = 1 para que la primera columna, que es gene_id, se use como nombres de fila en lugar de una columna 
# SOlo usaremos el annotation de hisat2 paired-end, ya que es el mismo archivo para todas las alineadores 
annotation <- read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_lengths.tsv", row.names = 1)
#NOTA: Usamos header = FALSE porque el archivo no tiene encabezado, y row.names = 1 para que la primera columna se use como fila 
gene_name_map <- read.delim("/export/storage/users/andreavg/transcriptomica/deseq/data/gencode/mm39-gencode-M36-gene_id-gene_name.txt", header = FALSE, row.names = 1)
# Definimos el directorio base de salida pa guardar los resultados
base <- "/export/storage/users/andreavg/transcriptomica/deseq/results/"

# Aplicamos la función a cada una de las tablas de conteos con los distintos alineadores y tipos de datos (paired-end y single-end)
hisat_paired <- Deseq2_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "hisat2/paired_end/deseq_results/"))
hisat_single <- Deseq2_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "hisat2/single_end/deseq_results/"))
star_paired <- Deseq2_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "star/paired_end/deseq_results/"))
star_single <- Deseq2_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "star/single_end/deseq_results/"))

# ----------------------------------------------- TERMINA EL SCRIPT ------------------------------------------------- #
```

```bash
# Creamos el script de DESeq2 para salmon
nano Deseq2_salmon.R
```
```r
# ----------------------------------------------- EMPIEZA EL SCRIPT ------------------------------------------------- #

# Cargamos las librerias que utilizaremos 
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(edgeR)
library(tximport) # Libreria necesaria para procesar los datos de Salmon

# Hacemos una función que realice todo el análisis de DESeq2 para los datos de Salmon
Deseq2_Salmon_analysis <- function(txi_obj, gene_name_map, output_dir) {

  # Creamos la ruta dinámica de salida para guardar los resultados de cada análisis
  fig_dir <- paste0(output_dir, "figuras/")

  # Para generar la tabla de metadatos primero extraemos los nombres de las muestras del objeto txi
  sample_names <- colnames(txi_obj$counts)

  # Generamos factores para la edad y el sexo 
  # NOTA: Usamos grepl para asignar mm_24 o mm_9 basándonos en el nombre de la condición que ya limpiamos
  age <- age <- factor(c("mm_24", "mm_9", "mm_24", "mm_9", "mm_24", rep("mm_9", 2)), levels = c("mm_9", "mm_24"))
  sex <- factor(c(rep("m", 7)))
  
  # Generamos la tabla de metadatos
  meta_data <- data.frame(sample_names, age, sex)
  # Eliminamos los nombres de las filas y asignamos la columna sample_names como nombres de fila 
  meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var = "sample_names")

  # Checaremos que los nombres de la tabla de metadatos estan en el mismo orden que las columnas del objeto txi
  k <- all(colnames(txi_obj$counts) == rownames(meta_data))
  # Vemos que si arroja TRUE, entonces esta todo correcto
  print(paste0("Los nombres de las columnas de txi_obj$counts y los nombres de fila de meta_data coinciden: ", k))

  # Cargamos los datos en un objeto DESeq2 usando la función específica para tximport
  # NOTA: Usamos DESeqDataSetFromTximport porque utiliza las abundancias y longitudes efectivas calculadas por Salmon
  dds <- DESeqDataSetFromTximport(txi = txi_obj, colData = meta_data, design = ~ 0 + age)
  design <- model.matrix(~ 0 + age)

  # Filtramos los genes que contienen una baja expresión usando filterByExpr() de edgeR
  keep <- filterByExpr(dds, design)
  # Vemos cuantos genes se mantienen después del filtrado
  suma_keep <- sum(keep)
  # Guardamos el número de genes que se mantienen después del filtrado en un archivo de texto
  write.table(suma_keep, paste0(output_dir, "genes_filtrados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  # Filtramos el objeto DESeq2
  dds <- dds[keep,]
  # Borramos keep 
  rm(keep)

  # Generamos el PCA para visualizar la variabilidad de los datos 
  vsd <- vst(dds)
  # Hacemos el PCA coloreado por la edad
  PCA_plot <- plotPCA(vsd, intgroup = "age") + theme_classic(base_size = 25, base_line_size = 1)
  # Guardamos la imagen del PCA
  ggsave(paste0(fig_dir, "PCA_plot.png"), plot = PCA_plot, width = 8, height = 6)

  # Calculamos los factores de normalización, varianza y ajustamos el modelo 
  dds <- DESeq(dds)

  # Ahora calcularemos la TPM (Transcripts Per Million) en escala log2
  # NOTA: Para Salmon extraemos la abundancia directamente del objeto txi y añadimos un pseudoconteo de 0.1
  log2_tpm <- log2(txi_obj$abundance[rownames(dds), ] + 0.1) 

  # Guardamos la tabla de TPMs
  gene_names_tpm <- gene_name_map[rownames(log2_tpm), ]
  write.table(cbind(gene_names_tpm, log2_tpm), paste0(output_dir, "TPM_log2-table.txt"), sep="\t", quote=FALSE)

  # Hacemos el contraste de expresión diferencial entre mm_24 y mm_9
  contrast <- makeContrasts(m24_vs_m9 = agemm_24 - agemm_9, levels = design) 

  # Ejecutamos el análisis y extraemos los resultados
  res <- results(dds, contrast = contrast[, "m24_vs_m9"])
  # Añadimos el nombre del gen a la tabla de resultados
  res$Gene_name <- gene_name_map[rownames(res), ]

  # Establecemos los tresholds para LFC y FDR
  FDR <- 0.05 
  LFC <- 0.5

  # Filtramos los genes Up y Down regulados
  up <- (res$log2FoldChange > LFC) & (res$padj < FDR)
  up[which(is.na(up))] = FALSE
  down <- (res$log2FoldChange < -LFC) & (res$padj < FDR)
  down[which(is.na(down))] = FALSE

  # Guardamos el conteo de genes significativos
  write.table(sum(up), paste0(output_dir, "genes_up_regulados.txt"), row.names = F, col.names = F, quote = F)
  write.table(sum(down), paste0(output_dir, "genes_down_regulados.txt"), row.names = F, col.names = F, quote = F)

  # Guardamos las tablas de resultados finales
  write.table(res[up,], paste0(output_dir, "deseq-Salmon-Up.txt"), sep="\t", quote=FALSE, row.names=TRUE)
  write.table(res[down,], paste0(output_dir, "deseq-Salmon-Down.txt"), sep="\t", quote=FALSE, row.names=TRUE)

  # Hacemos el gráfico de volcan para visualizar los resultados de expresión diferencial
  vpcolors = c("gray", "#057b5a", "#091d8c") 
  names(vpcolors) = c("NO", "DOWN", "UP") 
  res$DE = "NO" 
  res[up, "DE"] = "UP"
  res[down, "DE"] = "DOWN"

  volcano_plot <- ggplot(data = as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), col = DE)) +
          geom_point(alpha = 0.4, size = 1.5) + 
          labs(title = "Volcano plot Salmon") + 
          scale_color_manual(values = vpcolors) +
          geom_vline(xintercept = c(-LFC, LFC), col = "black", linetype = "longdash") +
          geom_hline(yintercept = -log10(FDR), col = "black", linetype = "longdash") +
          theme_classic(base_size = 15, base_line_size = 1)
  ggsave(paste0(fig_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)

  # Cremos los heatmaps para visualizar los patrones de expresión de los genes significativos
  significant <- res[up | down, ]
  significant_order <- significant[order(significant$padj), ]
  
  # Heatmap Top Globales 
  top_genes <- head(rownames(significant_order), n = 2000)
  zscore_t <- t(scale(t(log2_tpm[top_genes, ])))
  orden_columnas <- order(meta_data$age)
  
  heatmap_plot <- Heatmap(zscore_t[, orden_columnas], cluster_rows = T, cluster_columns = F, show_row_names = F, 
                          name = "Z-score", column_title = "Top 2000 Genes (Salmon)",
                          top_annotation = HeatmapAnnotation(Edad = meta_data$age[orden_columnas], # Ordenamos el heatmap según la edad 
                                                            col = list(Edad = c("mm_9" = "#370857", "mm_24" = "#7316b0")))) 
  
  # Guardamos la imagen del heatmap de los top genes 
  png(paste0(fig_dir, "heatmap_top_genes.png"), width = 11.25, height = 7.5, res = 300, units = "in")
  draw(heatmap_plot)
  dev.off() 

  # Heatmap Top 20
  top_20 <- head(rownames(significant_order), n = 20)
  zscore_top <- t(scale(t(log2_tpm[top_20, ])))
  
  # Para los ratones de 9 meses usaremos un deep purple y para los de 24 meses un purple más claro
  heatmap_top20 <- Heatmap(zscore_top[, orden_columnas], cluster_rows = T, cluster_columns = F, 
                           row_labels = gene_name_map[rownames(zscore_top), ], 
                           name = "Z-score", column_title = "Top 20 Genes (Salmon)",
                           top_annotation = HeatmapAnnotation(Edad = meta_data$age[orden_columnas], # Usamos la misma orden de columnas que en el heatmap global para mantener la consistencia visual
                                                             col = list(Edad = c("mm_9" = "#370857", "mm_24" = "#7316b0"))))
  
  # Guardamos la imagen del heatmap de los top 20 genes
  png(paste0(fig_dir, "heatmap_top_20_genes.png"), width = 11.25, height = 7.5, res = 300, units = "in")
  draw(heatmap_top20)
  dev.off()
}

# Cargamos el mapa de nombres de genes 
gene_name_map <- read.delim("/export/storage/users/andreavg/transcriptomica/deseq/data/gencode/mm39-gencode-M36-gene_id-gene_name.txt", header = FALSE, row.names = 1)
base <- "/export/storage/users/andreavg/transcriptomica/deseq/results/salmon/"

# Aplicamos la función a los objetos RDS que procesamos anteriormente
salmon_pe <- Deseq2_Salmon_analysis(readRDS(paste0(base, "paired_end/txi/txi_PE_condicion.rds")), gene_name_map, paste0(base, "paired_end/deseq_results/"))
salmon_se <- Deseq2_Salmon_analysis(readRDS(paste0(base, "single_end/txi/txi_SE_condicion.rds")), gene_name_map, paste0(base, "single_end/deseq_results/")) 

# ----------------------------------------------- TERMINA EL SCRIPT ------------------------------------------------- #
```

```bash
# Activamos entorno 
conda activate .main_mamba 
# Corremos scripts 
Rscript Deseq2_hisat_star.R
Rscript Deseq2_salmon.R
```

### Análisis con edgeR

El script `EdgeR_hisat_star.R` aplicó el mismo planteamiento general, pero usando la estructura de edgeR (`DGEList`) y su flujo estadístico. Para Salmon se empleó `EdgeR_salmon.R`, también sobre el objeto `txi` ya consolidado.

La secuencia de trabajo fue:

- construcción de `DGEList`,
- filtrado con `filterByExpr` y `min.count = 20`,
- normalización TMM con `calcNormFactors`,
- cálculo de `logCPM` y PCA,
- estimación robusta de dispersión,
- ajuste con `glmFit`,
- contraste con `glmLRT`,
- exportación de genes significativos, volcano plots y heatmaps.

```bash
# Creamos en src el script para hisat2 y star
nano EdgeR_hisat_star.R
```

```r
# ----------------------------------------------- EMPIEZA EL SCRIPT ------------------------------------------------- #

# Cargamos las librerias que utilizaremos 
library(edgeR) 
library(limma) 
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)

# Hacemos una función que realice todo el análisis de edgeR
edgeR_analysis <- function(gene_counts, annotation, gene_name_map, output_dir) {

  # Creamos la ruta dinámica de salida para guardar los resultados de cada análisis
  fig_dir <- paste0(output_dir, "figuras/")

  # Para generar la tabla de metadatos primero generamos factores
  # Generamos factores ppara la edad y el sexo 
  age <- factor(c("mm_24", "mm_9", "mm_24", "mm_9", "mm_24", rep("mm_9", 2)), levels = c("mm_9", "mm_24"))
  sex <- factor(c(rep("m", 7)))
  # Asignamos un color para el factor de la edad, en donde mm_9 es lightblue y mm_24 es blue
  sample_colors <- c("blue", "lightblue", "blue", "lightblue", "blue", rep("lightblue", 2)) 
  # Asignamos a la variable sample_names los nombres de las columnas de la tabla de conteos
  sample_names <- colnames(gene_counts)

  # Generamos la tabla de metadatos
  meta_data <- data.frame(sample_names, age, sex)
  # ELiminamos los nombres de las filas y asignamos la columna sample_names como nombres de fila 
  meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var = "sample_names")

  # Checaremos que los nombres de la tabla de metadatos estan en el mismo orden que las columnas de los gene_counts
  k <- all(colnames(gene_counts) == rownames(meta_data))
  # Si arroja TRUE, entonces esta todo correcto
  print(paste0("Los nombres de las columnas de gene_counts y los nombres de fila de meta_data coinciden: ", k))

  # Cargamos los datos en un objeto DGEList de edgeR (Digital Gene Expression List)
  # NOTA: Usamos round() para redondear los valores de conteo a enteros
  # El diseño lo establecemos como ~ 0 + age para no tener intercepto
  dge <- DGEList(counts = round(gene_counts), group = meta_data$age)
  design <- model.matrix(~ 0 + age, data = meta_data)

  # FIltramos los genes que contienen una baja expresión, usando filterByExpr() de edgeR
  # filterByExpr compara la expresión de cada gen con un umbral para eliminar los genes sin expresión suficiente 
  # Hacemos un min.count de 20 para eliminar los genes muy poco expresados y que no se castigue el FDR
  keep <- filterByExpr(dge, design, min.count = 20)
  # Vemos cuantos genes se mantienen después del filtrado
  suma_keep <- sum(keep)
  # Guardamos el número de genes que se mantienen después del filtrado en un archivo de texto
  write.table(suma_keep, paste0(output_dir, "genes_filtrados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Filtramos el objeto DGEList para quedarnos solo con los genes que cumplen el criterio 
  # NOTA: En edgeR es importante indicar keep.lib.sizes=FALSE para que recalcule los tamaños de librería reales tras el filtrado
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  # Borramos keep 
  rm(keep)

  # Calculamos los factores de normalización
  # NOTA: Utilizamos calcNormFactors() que por defecto aplica el método TMM (Trimmed Mean of M-values) y es matemáticamente equivalente a Median of Ratios de DeSeq2
  # Generamos el PCA para visualizar la variabilidad de los datos 
  # En edgeR calculamos los conteos por millón en log2 (logCPM) para estabilizar la varianza
  logCPM <- cpm(dge, log = TRUE, prior.count = 2)
  # Realizamos el PCA con prcomp transponiendo la matriz (genes en columnas, muestras en filas)
  pca_res <- prcomp(t(logCPM), scale. = TRUE)
  
  # Preparamos el data.frame para graficar con ggplot
  pca_data <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], age = meta_data$age)

  # Hacemos el PCA, mapeando el color al factor edad
  PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = age)) +
          geom_point(size = 4) +
          theme_classic(base_size = 25, base_line_size = 1) +
          labs(title = "PCA Plot", x = "PC1", y = "PC2")

  # Guardamos la imagen del PCA
  ggsave(paste0(fig_dir, "PCA_plot.png"), plot = PCA_plot, width = 8, height = 6)

  # Ahora calcularemos la TPM (Transcripts Per Million)
  # Añadimos la longitud del gen al objeto de edgeR para calcular RPKM
  # NOTA: Asignamos las longitudes a dge$genes$Length para que la función rpkm() las reconozca automáticamente
  dge$genes$Length = annotation[rownames(dge), ]

  # Calculamos RPKM (análogo a FPKM) y pasamos a log2. 
  # NOTA: Le añadimos un pseudoconteo de 0.1 para evitar problemas con los logaritmos de cero
  log2_fpkm <- log2(rpkm(dge) + 0.1) 

  # Escribimos la formula para convertir de FKPM a TPM dentro de un espacio lograrítmico 
  fpkm2tpm_log2 <- function(fpkm) { fpkm - log2(sum(2^fpkm)) + log2(1e6) } 
  # Aplicamos la fórmula a cada columna de la tabla
  log2_tpm <- apply(log2_fpkm, 2, fpkm2tpm_log2) 

  # Guardamos la tabla 
  # Agregamos [, 1] para forzar que sea un vector y no un dataframe
  gene_names <- gene_name_map[rownames(log2_tpm), 1]
  write.table(cbind(gene_names, log2_tpm), paste0(output_dir, "TPM_log2-table.txt"), sep="\t", quote=FALSE)

  # Primero estimamos la dispersión de los datos de forma robusta
  # NOTA: robust = TRUE protege las estimaciones de varianza contra muestras atípicas (outliers) 
  dge <- estimateDisp(dge, design, robust = TRUE)
  
  # Ajustamos el modelo usando Likelihood Ratio Test (LRT)
  # NOTA: LRT es el método clásico de edgeR y suele ser menos castigador con el FDR que QLFit cuando hay mucha dispersión
  fit <- glmFit(dge, design, robust = TRUE)

  # Hacemos el contraste entre nuestras condiciones, en este caso entre mm_24 y mm_9
  # NOTA: Los nombres en design son agemm_9 y agemm_24 generados por model.matrix
  contrast <- makeContrasts(m24_vs_m9 = agemm_24 - agemm_9, levels = design) 

  # Ejecutamos la prueba estadística con glmLRT
  lrt <- glmLRT(fit, contrast = contrast[, "m24_vs_m9"])

  # Extraemos los resultados con topTags extrayendo todos los genes (n = Inf)
  res_edgeR <- topTags(lrt, n = Inf)$table

  # Renombramos las columnas de edgeR para que coincidan con las de DESeq2
  res <- as.data.frame(res_edgeR)
  colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
  colnames(res)[colnames(res) == "FDR"] <- "padj"

  # Añadimos el nombre del gen a la tabla de resultados
  res$Gene_name <- gene_name_map[rownames(res), 1]

  # Establecemos los tresholds para LFC y FDR
  FDR <- 0.05 # Escogemos un FDR de 0.05 para ser un poco más permisivos 
  LFC <- 0.5

  # Filtramos los resultados para sacar el número de genes up y down regulados, usando los tresholds establecidos
  up <- (res$log2FoldChange > LFC) & (res$padj < FDR)
  # Reemplazamos los valores NA por FALSE, ya que los genes que no cumplen con los criterios de significancia 
  up[which(is.na(up))] = FALSE
  # Guardamos el número de genes up regulados en un archivo de texto
  suma_up <- sum(up)
  write.table(suma_up, paste0(output_dir, "genes_up_regulados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Hacemos lo mismo para los genes down regulados
  # Aqui buscamos los genes con log2FoldChange menor a -LFC y padj menor a FDR
  down <- (res$log2FoldChange < -LFC) & (res$padj < FDR)
  # Reemplazamos los valores NA por FALSE, ya que los genes que no cumplen con los criterios de significancia
  down[which(is.na(down))] = FALSE
  # Guardamos el número de genes down regulados en un archivo de texto
  suma_down <- sum(down)
  write.table(suma_down, paste0(output_dir, "genes_down_regulados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

  # Guardamos las tablas 
  write.table(res[up,], paste0(output_dir, "edger-DEG_up_0.05.txt"), sep="\t", quote=FALSE, row.names=TRUE)
  write.table(res[down,], paste0(output_dir, "edger-DEG_down_0.05.txt"), sep="\t", quote=FALSE, row.names=TRUE)

  # Ahora procederemos a graficar con un volcanoplot
  # Asignamos los colores para las categorías
  vpcolors = c("gray", "#4411be", "#146f0a") 
  names(vpcolors) = c("NO", "DOWN", "UP") 

  # Creamos la columna DE de resultados 'res'
  # Primero marcamos todos como "NO"
  res$DE = "NO" 
  # Luego usamos los vectores 'up' y 'down' para etiquetar los significativos
  res[up, "DE"] = "UP"
  res[down, "DE"] = "DOWN"

  # Creamos la gráfica con ggplot2
  volcano_plot <- ggplot(data = as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), col = DE)) +
          geom_point(alpha = 0.4, size = 1.5) + 
          labs(title = "Volcano plot: Músculo 24m vs 9m (Escala Auto)", 
               x = "log2 Expression fold change", 
               y = "-log10 FDR") + 
          scale_color_manual(values = vpcolors) +
          geom_vline(xintercept = c(-LFC, LFC), col = "black", linetype = "longdash") +
          geom_hline(yintercept = -log10(FDR), col = "black", linetype = "longdash") +
          theme_classic(base_size = 15, base_line_size = 1)

  
  # Guardamos la imagen del volcano plot
  ggsave(paste0(fig_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)

  # Haremos el heatmap de los genes diferencialmente expresados
  # Primero con los genes más significativos 
  # Filtramos solo los genes que pasan los filtros (Up y Down)
  significant <- res[up | down, ]

  # Ordenamos los genes por su significancia 
  significant_order <- significant[order(significant$padj), ]

    # Si no hay genes significativos, omitimos los heatmaps para evitar errores
    if (nrow(significant_order) == 0) {
        message("No hay genes significativos para generar heatmaps.")
        return(invisible(res))
    }

  # Tomamos los 2,000 genes más significativos 
  # Usamos rownames para obtener los IDs de Ensembl de esos genes top, dejamos el valor n = 2000 para tener un umbral alto y no estar ajustando entre disitintos alineadores
  top_genes <- head(rownames(significant_order), n = 2000)

  # Calculamos el Z-score a partir de tus datos TPM
  # El Z-score centra la expresión de cada gen: 0 es el promedio, valores positivos son arriba del promedio, negativos abajo.
    zscore_t <- t(scale(t(log2_tpm[top_genes, , drop = FALSE])))

  # Ordenamos las columnas para que en el heatmap se muestren en el orden de edad, primero los mm_9 y luego los mm_24
  orden_columnas <- order(meta_data$age)
    zscore_significant <- zscore_t[, orden_columnas, drop = FALSE]

  # Creamos el Heatmap con ComplexHeatmap
  # km = 2 le pide a R que intente separar los genes en 2 grupos principales (clusters)
  heatmap_plot <- Heatmap(zscore_significant, 
              cluster_rows = T,         
              cluster_columns = F,      # F para que respete el orden de edad que pusimos
              show_row_names = F,       
              name = "Z-score", 
              km = 2,                   
              column_title = "Heatmap de los top genes",
              column_names_gp = gpar(
                  col = "black", 
                  fontsize = 10, 
                  fontface = "bold"
              ),
              top_annotation = HeatmapAnnotation(
                  Edad = meta_data$age[orden_columnas],
                  col = list(Edad = c("mm_9" = "lightpink", "mm_24" = "deeppink"))
              )) 

  
  # Guardamos la imagen del heatmap
  png(paste0(fig_dir, "heatmap_top_genes.png"), width = 11.25, height = 7.5, res = 300, units = "in")
  # Es importante usar draw() para asegurar que todos los elementos se guarden
  draw(heatmap_plot)
  # Cerramos el archivo
  dev.off()


  # Ahora haremos el heatmap para el top 20 genes más significativos
  # Seleccionamos los 20 mejores (los de menor padj)
  top_20_ids <- head(rownames(significant_order), n = 20)

  # Calculamos el Z-score solo para estos 20 genes
    zscore_top <- t(scale(t(log2_tpm[top_20_ids, , drop = FALSE])))

  # Aplicamos el orden de las columnas por edad (9m primero, luego 24m)
    zscore_top_ordenado <- zscore_top[, orden_columnas, drop = FALSE]

    if (nrow(zscore_top_ordenado) == 0) {
        message("No hay genes en el top 20 para generar el heatmap.")
        return(invisible(res))
    }

  # Usamos [, 1] para forzar vector. Convertimos a as.character y usamos unname() para limpiar metadatos
  nombres_top20 <- unname(as.character(gene_name_map[rownames(zscore_top_ordenado), 1]))
  
  # Sustituimos NA o espacios vacíos por el ID de Ensembl
  filtro_nulos <- is.na(nombres_top20) | nombres_top20 == "" | nombres_top20 == "NA"
  nombres_top20[filtro_nulos] <- rownames(zscore_top_ordenado)[filtro_nulos]

  # Creamos el Heatmap con los nombres de los genes visibles
  heatmap_top20 <- Heatmap(zscore_top_ordenado, 
              cluster_rows = T, 
              cluster_columns = F, 
              row_labels = nombres_top20, # Ahora es un vector de texto puro garantizado
              name = "Z-score", 
              km = 2, 
              column_title = "Top 20 Genes Significativos",
              column_names_gp = gpar(
                  col = "black", 
                  fontsize = 10, 
                  fontface = "bold"
              ),
              row_names_gp = gpar(
                  fontsize = 10, 
                  fontface = "italic" # Los nombres de genes siempre van en cursivas
              ),
              top_annotation = HeatmapAnnotation(
                  Edad = meta_data$age[orden_columnas],
                  col = list(Edad = c("mm_9" = "lightpink", "mm_24" = "deeppink"))
              ))

  
  # Guardamos la imagen del heatmap de los top 20 genes
  png(paste0(fig_dir, "heatmap_top_20_genes.png"),  width = 11.25, height = 7.5, res = 300, units = "in")
  # Es importante usar draw() para asegurar que todos los elementos se guarden
  draw(heatmap_top20)
  # Cerramos el archivo
  dev.off()
} 

# Cargamos los archivos que utilizaremos 
# NOTA: Usamos row.names = 1 para que la primera columna, que es gene_id, se use como nombres de fila en lugar de una columna 
# SOlo usaremos el annotation de hisat2 paired-end, ya que es el mismo archivo para todas las alineadores 
annotation <- read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_lengths.tsv", row.names = 1)
#NOTA: Usamos header = FALSE porque el archivo no tiene encabezado, y row.names = 1 para que la primera columna se use como fila 
gene_name_map <- read.delim("/export/storage/users/andreavg/transcriptomica/deseq/data/gencode/mm39-gencode-M36-gene_id-gene_name.txt", header = FALSE, row.names = 1)
# Definimos el directorio base de salida pa guardar los resultados
base <- "/export/storage/users/andreavg/transcriptomica/deseq/results/"

# Aplicamos la función edgeR_analysis a cada una de las tablas de conteos 
hisat_paired <- edgeR_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "hisat2/paired_end/edgeR_results/"))
hisat_single <- edgeR_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "hisat2/single_end/edgeR_results/"))
star_paired <- edgeR_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "star/paired_end/edgeR_results/"))
star_single <- edgeR_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "star/single_end/edgeR_results/"))

# ----------------------------------------------- TERMINA EL SCRIPT ------------------------------------------------- #
```

```bash
# Creamos el script para salmon
nano EdgeR_salmon.R
```

```r
# ----------------------------------------------- EMPIEZA EL SCRIPT ------------------------------------------------- #

# Cargamos las librerias que utilizaremos 
library(edgeR) 
library(limma) 
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(tximport) 

# Hacemos una función que realice todo el análisis de edgeR para los datos de Salmon
edgeR_Salmon_analysis <- function(txi_obj, gene_name_map, output_dir) {

  # Creamos la ruta dinámica de salida para guardar los resultados de cada análisis
  fig_dir <- paste0(output_dir, "figuras/")

  # Para generar la tabla de metadatos primero extraemos los nombres de las muestras del objeto txi
  sample_names <- colnames(txi_obj$counts)

  # Generamos factores para la edad y el sexo 
  age <- factor(c("mm_24", "mm_9", "mm_24", "mm_9", "mm_24", rep("mm_9", 2)), levels = c("mm_9", "mm_24"))
  sex <- factor(rep("m", 7))
  
  # Generamos la tabla de metadatos
  meta_data <- data.frame(sample_names, age, sex)
  # Eliminamos los nombres de las filas y asignamos la columna sample_names como nombres de fila 
  meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var = "sample_names")

  # Checaremos que los nombres de la tabla de metadatos estan en el mismo orden que las columnas de los conteos en el objeto txi
  k <- all(colnames(txi_obj$counts) == rownames(meta_data))
  # Si arroja TRUE, entonces esta todo correcto
  print(paste0("Los nombres de las columnas de txi_obj$counts y los nombres de fila de meta_data coinciden: ", k))

  # Cargamos los datos en un objeto DGEList de edgeR
  # NOTA: Usamos round() porque edgeR requiere que los conteos sean enteros, y los extraemos de txi_obj$counts
  dge <- DGEList(counts = round(txi_obj$counts), group = meta_data$age) 
  design <- model.matrix(~ 0 + age, data = meta_data)

  # Filtramos los genes que contienen una baja expresión usando filterByExpr()
  # Hacemos un min.count de 20 para eliminar genes con ruido y mejorar la potencia estadística
  keep <- filterByExpr(dge, design, min.count = 20)
  suma_keep <- sum(keep)
  # Guardamos el número de genes filtrados
  write.table(suma_keep, paste0(output_dir, "genes_filtrados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Filtramos el objeto DGEList recalculando los tamaños de librería
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  rm(keep)

  # Calculamos los factores de normalización usando el método TMM
  dge <- calcNormFactors(dge)

  # Generamos el PCA para visualizar la variabilidad de los datos 
  # Calculamos logCPM para estabilizar la varianza antes del PCA
  logCPM <- cpm(dge, log = TRUE, prior.count = 2)
  pca_res <- prcomp(t(logCPM), scale. = TRUE)
  
  pca_data <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], age = meta_data$age)

  # Hacemos el PCA coloreado por el factor edad
  PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = age)) +
          geom_point(size = 4) +
          theme_classic(base_size = 25, base_line_size = 1) +
          labs(title = "PCA Plot Salmon (edgeR)", x = "PC1", y = "PC2")

  # Guardamos la imagen del PCA
  ggsave(paste0(fig_dir, "PCA_plot.png"), plot = PCA_plot, width = 8, height = 6)

  # Ahora manejaremos la TPM (Transcripts Per Million)
  # NOTA: En Salmon, las abundancias ya vienen calculadas en txi_obj$abundance. Las pasamos a escala log2.
  # Filtramos para que coincidan con los genes que pasaron el filtro de expresión
  log2_tpm <- log2(txi_obj$abundance[rownames(dge), ] + 0.1) 

  # Guardamos la tabla de TPMs
  # Forzamos que los nombres de genes sean un vector
  gene_names_tpm <- gene_name_map[rownames(log2_tpm), 1]
  write.table(cbind(gene_names_tpm, log2_tpm), paste0(output_dir, "TPM_log2-table.txt"), sep="\t", quote=FALSE)

  # Estimamos la dispersión de los datos de forma robusta para protegernos de outliers
  dge <- estimateDisp(dge, design, robust = TRUE)
  
  # Ajustamos el modelo usando Likelihood Ratio Test (LRT)
  fit <- glmFit(dge, design, robust = TRUE)

  # Hacemos el contraste entre mm_24 y mm_9
  contrast <- makeContrasts(m24_vs_m9 = agemm_24 - agemm_9, levels = design) 

  # Ejecutamos la prueba estadística con glmLRT
  lrt <- glmLRT(fit, contrast = contrast[, "m24_vs_m9"])

  # Extraemos los resultados completos
  res_edgeR <- topTags(lrt, n = Inf)$table

  # Renombramos las columnas para que sean consistentes con otros análisis
  res <- as.data.frame(res_edgeR)
  colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
  colnames(res)[colnames(res) == "FDR"] <- "padj"

  # Añadimos el nombre del gen a la tabla de resultados
  res$Gene_name <- gene_name_map[rownames(res), 1]

  # Establecemos los tresholds de significancia
  FDR <- 0.05 
  LFC <- 0.5

  # Obtenemos los genes Up y Down regulados
  up <- (res$log2FoldChange > LFC) & (res$padj < FDR)
  up[which(is.na(up))] = FALSE
  down <- (res$log2FoldChange < -LFC) & (res$padj < FDR)
  down[which(is.na(down))] = FALSE

  # Guardamos los conteos de genes significativos
  write.table(sum(up), paste0(output_dir, "genes_up_regulados.txt"), row.names = F, col.names = F, quote = F)
  write.table(sum(down), paste0(output_dir, "genes_down_regulados.txt"), row.names = F, col.names = F, quote = F)

  # Guardamos las tablas de DEGs
  write.table(res[up,], paste0(output_dir, "edger-DEG_up_0.05.txt"), sep="\t", quote=FALSE, row.names=TRUE)
  write.table(res[down,], paste0(output_dir, "edger-DEG_down_0.05.txt"), sep="\t", quote=FALSE, row.names=TRUE)

  # Procederemos a graficar con un volcanoplot
  vpcolors = c("gray", "#e6160b", "#07a394") 
  names(vpcolors) = c("NO", "DOWN", "UP") 
  res$DE = "NO" 
  res[up, "DE"] = "UP"
  res[down, "DE"] = "DOWN"

  volcano_plot <- ggplot(data = as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), col = DE)) +
          geom_point(alpha = 0.4, size = 1.5) + 
          labs(title = "Volcano plot Salmon (edgeR)", x = "log2 Fold Change", y = "-log10 FDR") + 
          scale_color_manual(values = vpcolors) +
          geom_vline(xintercept = c(-LFC, LFC), col = "black", linetype = "longdash") +
          geom_hline(yintercept = -log10(FDR), col = "black", linetype = "longdash") +
          theme_classic(base_size = 15, base_line_size = 1)
  
  # Guardamos la imagen del volcanoplot
  ggsave(paste0(fig_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)

  # Haremos los heatmaps de los genes diferencialmente expresados
  significant <- res[up | down, ]
  significant_order <- significant[order(significant$padj), ]
  
  # Solo haremos el heatmap si hay genes significativos, para evitar errores en caso de que no se encuentren DEGs
  if (nrow(significant_order) > 0) {
    # Heatmap de los genes más significativos
    # Para los ratones de 9 meses utilizaremos un lightgreen y para los de 24 meses un deep green
    top_genes <- head(rownames(significant_order), n = 2000)
    zscore_t <- t(scale(t(log2_tpm[top_genes, , drop = FALSE])))
    orden_columnas <- order(meta_data$age)
    
    heatmap_plot <- Heatmap(zscore_t[, orden_columnas, drop = FALSE], 
                cluster_rows = T, cluster_columns = F, show_row_names = F, 
                name = "Z-score", km = 2, column_title = "Heatmap Salmon edgeR (Top 2000)",
                top_annotation = HeatmapAnnotation(Edad = meta_data$age[orden_columnas],
                                                  col = list(Edad = c("mm_9" = "#3da528", "mm_24" = "#113706")))) 
    
    # Guardamos la imagen del heatmap
    png(paste0(fig_dir, "heatmap_top_genes.png"), width = 11.25, height = 7.5, res = 300, units = "in")
    draw(heatmap_plot)
    dev.off()

    # Heatmap del Top 20
    top_20_ids <- head(rownames(significant_order), n = 20)
    zscore_top <- t(scale(t(log2_tpm[top_20_ids, , drop = FALSE])))
    nombres_top20 <- unname(as.character(gene_name_map[rownames(zscore_top), 1]))
    
    heatmap_top20 <- Heatmap(zscore_top[, orden_columnas, drop = FALSE], 
                cluster_rows = T, cluster_columns = F, row_labels = nombres_top20,
                name = "Z-score", km = 2, column_title = "Top 20 Genes Salmon (edgeR)",
                row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                top_annotation = HeatmapAnnotation(Edad = meta_data$age[orden_columnas],
                                                  col = list(Edad = c("mm_9" = "#3da528", "mm_24" = "#113706"))))
    
    # Guardamos la imagen del heatmap del Top 20
    png(paste0(fig_dir, "heatmap_top_20_genes.png"), width = 11.25, height = 7.5, res = 300, units = "in")
    draw(heatmap_top20)
    dev.off()
  }
  
} 


# Cargamos el mapa de nombres de genes
gene_name_map <- read.delim("/export/storage/users/andreavg/transcriptomica/deseq/data/gencode/mm39-gencode-M36-gene_id-gene_name.txt", header = FALSE, row.names = 1)
base <- "/export/storage/users/andreavg/transcriptomica/deseq/results/salmon/"

# Aplicamos la función a los objetos RDS de Salmon (Paired-End y Single-End)
# NOTA: Usamos los archivos con las condiciones ya asignadas
salmon_pe_edgeR <- edgeR_Salmon_analysis(readRDS(paste0(base, "paired_end/txi/txi_PE_condicion.rds")), gene_name_map, paste0(base, "paired_end/edgeR_results/"))
salmon_se_edgeR <- edgeR_Salmon_analysis(readRDS(paste0(base, "single_end/txi/txi_SE_condicion.rds")), gene_name_map, paste0(base, "single_end/edgeR_results/"))

# ----------------------------------------------- TERMINA EL SCRIPT ------------------------------------------------- #
```

```bash
# Activamos entorno 
conda activate .main_mamba 
# Corremos scripts 
Rscript EdgeR_hisat_star.R
Rscript EdgeR_salmon.R
```

## Resultados

### Genes retenidos y significativos

Usando los mismos umbrales de corte, los resultados globales fueron los siguientes:

| Alineador | Modo | Método | Genes filtrados | Up | Down | Total DEGs |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| HISAT2 | Paired-end | DESeq2 | 16884 | 62 | 227 | 289 |
| HISAT2 | Paired-end | edgeR | 13657 | 24 | 297 | 321 |
| HISAT2 | Single-end | DESeq2 | 13722 | 65 | 208 | 273 |
| HISAT2 | Single-end | edgeR | 11645 | 28 | 343 | 371 |
| STAR | Paired-end | DESeq2 | 16847 | 59 | 230 | 289 |
| STAR | Paired-end | edgeR | 13582 | 25 | 347 | 372 |
| STAR | Single-end | DESeq2 | 13574 | 59 | 200 | 259 |
| STAR | Single-end | edgeR | 11518 | 27 | 349 | 376 |
| Salmon | Paired-end | DESeq2 | 15387 | 7 | 67 | 74 |
| Salmon | Paired-end | edgeR | 11969 | 4 | 44 | 48 |
| Salmon | Single-end | DESeq2 | 14912 | 30 | 78 | 108 |
| Salmon | Single-end | edgeR | 11844 | 13 | 25 | 38 |


### Comparación visual

#### PCA 

Los análisis de PCA muestran consistentemente una separación por edad (`mm_9` vs `mm_24`) en todos los flujos, pero con diferencias en la compactación de réplicas y en la dispersión general según la estrategia de cuantificación y el modo de secuenciación.

##### Por alineador

![PCA HISAT2 paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/hisat2/paired_end/deseq_results/figuras/PCA_plot.png)
**Figura 1. PCA para HISAT2 (paired-end, DESeq2)**

En la figura 1 se observá réplicas muy próximas dentro de cada edad y una separación clara entre `mm_9` y `mm_24`, lo que indicó baja variabilidad técnica y una señal biológica fuerte.

![PCA STAR paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/star/paired_end/deseq_results/figuras/PCA_plot.png)
**Figura 2. PCA para STAR (paired-end, DESeq2)**

De igual manera, en la figura 2 se mostró un patrón prácticamente superponible al de HISAT2. Las muestras de cada condición forman nubes compactas y bien separadas.

![PCA Salmon paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/salmon/paired_end/deseq_results/figuras/PCA_plot.png)
**Figura 3. PCA para Salmon (paired-end, DESeq2)**

Por su parte, la figura 3 muestrá que la separación por edad persistió, pero las réplicas muestran mayor dispersión que en los alineadores genómicos, lo que sugierió mayor variabilidad atribuible al pseudoalineamiento o a la cuantificación a nivel de transcritos.

##### Por método estadístico: DESeq2 vs edgeR

![PCA STAR single-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/star/single_end/deseq_results/figuras/PCA_plot.png)
**Figura 4. PCA para STAR (single-end, DESeq2)**

![PCA STAR single-end edgeR](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/star/single_end/edgeR_results/figuras/PCA_plot.png)
**Figura 5. PCA para STAR (single-end, edgeR)**

L a figura 4 y 5 muestrá que ambos métodos discriminarón claramente las dos edades. Con DESeq2 las réplicas tendierón a quedar ligeramente más compacta y con edgeR se mostró una dispersión algo mayor pero conservó la estructura de grupos, por lo que la interpretación biológica principal no cambia.

##### Por tipo de librería: single-end vs paired-end

![PCA Salmon single-end edgeR](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/salmon/single_end/edgeR_results/figuras/PCA_plot.png)
**Figura 6. PCA para Salmon (single-end, edgeR)**

![PCA Salmon paired-end edgeR](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/salmon/paired_end/edgeR_results/figuras/PCA_plot.png)
**Figura 7. PCA para Salmon (paired-end, edgeR)**

En las  figuras 6 y 7 se observa que en Salmon se apreció un poco más de sensibilidad al tipo de librería, algunos puntos se alejarón del núcleo del grupo en single-end mientras que en paired-end la agrupación puede ser algo más definida. En HISAT2 y STAR las diferencias entre single- y paired-end son menores y no alteran la separación biológica principal.

En conjunto, HISAT2 y STAR produjó PCAs muy similares, con réplicas compactas y una separación nítida entre edades, tanto en single-end como en paired-end. DESeq2 y edgeR mantuvó la misma estructura general del espacio de muestras; DESeq2 tendió a mostrar clusters ligeramente más compactos, mientras que edgeR presentó dispersión algo mayor sin afectar la interpretación biológica. Salmon conservó la separación por edad pero introdujó mayor dispersión y una mayor dependencia del modo de secuenciación. Por lo tanto, la estrategia de cuantificación/alineamiento (pseudoalineamiento vs alineamiento genómico) tiene un efecto más marcado sobre la forma del PCA que la elección entre DESeq2 y edgeR o el tipo de librería en los alineadores genómicos.

#### Volcano Plot 

##### Por alineador

![Volcano STAR paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/star/paired_end/deseq_results/figuras/volcano_plot.png)
**Figura 8. Volcano plot para STAR (paired-end, DESeq2)**

Como se observa en la Figura 8, STAR mostró una nube amplia de puntos significativos, con predominio de genes diferencialmente expresados hacia valores negativos de log2 fold change y una fracción menor hacia valores positivos. Se identificó una concentración considerable de genes con significancia estadística moderada a alta, lo que sugiere una señal transcriptómica intensa entre las dos edades. La dispersión horizontal de los puntos significativos fue amplia, aunque la mayor densidad se concentró cerca del umbral de cambio de expresión, con algunos genes que alcanzaron valores más extremos.

![Volcano HISAT2 paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/hisat2/paired_end/deseq_results/figuras/volcano_plot.png)
**Figura 9. Volcano plot para HISAT2 (paired-end, DESeq2)**

Como se aprecia en la Figura 9, HISAT2 presentó un patrón visual muy similar al de STAR, con una distribución asimétrica de genes significativos y una acumulación mayor en el lado de downregulation. Se observó una dispersión también amplia, aunque ligeramente más contenida que en STAR, lo que sugiere una señal diferencial algo más compacta. La cantidad de puntos significativos fue elevada, pero la nube de genes extremos pareció menos extendida que en STAR, con una mayor concentración alrededor de cambios de expresión moderados.

![Volcano Salmon paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/salmon/paired_end/deseq_results/figuras/volcano_plot.png)
**Figura 10. Volcano plot para Salmon (paired-end, DESeq2)**

Como se muestra en la Figura 10, Salmon presentó una cantidad menor de genes significativos que en los alineadores genómicos. Se identificó una nube más compacta y una dispersión reducida, lo que sugiere un comportamiento más conservador en la detección de genes diferencialmente expresados. Los cambios de expresión extrema aparecieron en menor número y la mayoría de los puntos significativos se concentró cerca de umbrales moderados de log2 fold change. En consecuencia, Salmon evidenció una señal diferencial más contenida, aunque todavía claramente asociada con el contraste entre edades.

En comparación general, STAR e HISAT2 produjeron volcano plots muy semejantes, con una mayor densidad de genes significativos y una señal más marcada que Salmon. Visualmente, STAR pareció mostrar la nube más amplia de puntos diferenciales, mientras que HISAT2 mantuvo una distribución ligeramente más compacta. Salmon, en cambio, presentó menos genes significativos y una configuración más conservadora. En conjunto, se concluye que STAR y HISAT2 detectaron una señal de expresión diferencial más amplia que Salmon, y que STAR mostró la tendencia visual más expansiva de las tres. 

##### Por método estadístico: DESeq2 vs edgeR

![Volcano STAR paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/star/paired_end/deseq_results/figuras/volcano_plot.png)
**Figura 11. Volcano plot para STAR (paired-end, DESeq2)**

Como se observa en la Figura 11, DESeq2 mostró una nube amplia de genes significativos, con una clara acumulación de puntos por debajo del eje central y una fracción relevante de genes hacia cambios positivos. Se identificó una distribución densa de puntos que cruzaron los umbrales de significancia, lo que sugiere una detección relativamente permisiva de expresión diferencial. La dispersión fue considerable y coexistieron genes con cambios moderados y genes con magnitudes más extremas, especialmente en el lado negativo del eje horizontal.

![Volcano STAR paired-end edgeR](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/star/paired_end/edgeR_results/figuras/volcano_plot.png)
**Figura 12. Volcano plot para STAR (paired-end, edgeR)**

Como se aprecia en la Figura 12, edgeR mantuvo el patrón biológico general, pero con una distribución visualmente más contenida para los genes upregulated  y con un conjunto de genes significativos algo más selectivo. Se observó una concentración marcada de puntos significativos alrededor de cambios de expresión negativos, aunque la nube global pareció menos densa que en DESeq2. A su vez, algunos genes alcanzaron magnitudes de cambio notables, lo que indica que edgeR destacó un subconjunto más acotado pero con valores extremos bien definidos. En este sentido, edgeR mostró una señal más conservadora, pero con contraste pronunciado en genes puntuales.

En la comparación entre ambos métodos, DESeq2 pareció recuperar una mayor cantidad de genes diferencialmente expresados y una nube más extensa de puntos significativos, mientras que edgeR produjo una distribución algo más restringida para los upregulated pero con varios genes de gran magnitud de cambio. Visualmente, DESeq2 ofreció una lectura más densa de la señal diferencial, en tanto que edgeR resaltó menos puntos, aunque algunos aparecieron más extremos. En conclusión, DESeq2 pareció aumentar ligeramente la detección de expresión diferencial, mientras que edgeR se mostró más selectivo.

##### Por tipo de librería: single-end vs paired-end

![Volcano HISAT2 paired-end edgeR](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/hisat2/paired_end/edgeR_results/figuras/volcano_plot.png)
**Figura 13. Volcano plot para HISAT2 (paired-end, edgeR)**

Como se observa en la Figura 13, el análisis paired-end presentó una distribución clara de genes diferencialmente expresados, con predominio de puntos significativos en la zona de downregulation y una fracción menor en la zona de upregulation. Se identificó una nube relativamente compacta, pero con una separación nítida respecto del centro del gráfico, lo que sugiere una señal definida y coherente. La dispersión de los puntos fue moderada y los genes significativos alcanzaron cambios de expresión apreciables, aunque sin una expansión excesiva del eje horizontal.

![Volcano HISAT2 single-end edgeR](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/hisat2/single_end/edgeR_results/figuras/volcano_plot.png)
**Figura 14. Volcano plot para HISAT2 (single-end, edgeR)**

Como se aprecia en la Figura 14, el análisis single-end conservó el patrón general observado en paired-end, pero mostró una dispersión ligeramente mayor y una distribución algo menos compacta de los puntos significativos. Se observó una señal diferencial todavía clara, aunque con una concentración algo menor de genes cruzando los umbrales de significancia. Los cambios de expresión parecieron más conservadores en algunos sectores y, en conjunto, la nube de puntos mostró una definición visual algo menos nítida que en paired-end. Por ello, el modo paired-end pareció ofrecer una estructura más consistente y fácilmente interpretable.

En la comparación entre ambas estrategias de librería, se observó que el patrón biológico principal se mantuvo, pero paired-end produjo una distribución ligeramente más definida y con una señal algo más compacta. Single-end presentó una dispersión mayor y una delimitación menos uniforme de los genes significativos, sin alterar por completo la interpretación global. En consecuencia, el tipo de librería pareció influir muy poco en la detección de genes diferencialmente expresados. 

#### Heatmaps

![Heatmap top 20 genes STAR paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/star/paired_end/deseq_results/figuras/heatmap_top_20_genes.png)
**Figura 15. Heatmap de los 20 genes diferencialmente expresados más significativos para STAR (paired-end, DESeq2)**

La visualización mostró que hay una separación clara entre `mm_9` y `mm_24`, con bloques de expresión recíprocos que reflejarón el contraste biológico entre las edades. 

- Genes observados en azul (relativamente subexpresados en `mm_24`): Serpinh1, Nrep, Tfrc, Smox, Col4a1, Asb1, Ythdf3, Thbd, Phka1, Ibtk, Egfl6, Hif1an, Ulk2, Arl5a, Ppp1r3c, Rimoc1, Clic5.
- Genes observados en rojo (relativamente sobreexpresados en `mm_24`): Rasd2, ENSMUSG00000135686, Dbp.

Observando la columna deeppink (señalando a los ratones de 24 meses) se conluyó que un gen aparezca en rojo implica que su expresión relativa es más alta en las réplicas de `mm_24` que en `mm_9`, mientras que el azul indica una expresión relativamente más baja en `mm_24`. En este heatmap dominó un patrón de represión para `mm_24`: la mayoría de los 20 genes seleccionados aparecen en tonos fríos, lo que sugirió una tendencia general hacia la disminución de la expresión génica en la condición `mm_24`. Esto implicó que la respuesta transcriptómica característica de `mm_24` está marcada por una reducción de la abundancia de numerosos transcritos, mientras que solo un pequeño subconjunto (por ejemplo `Rasd2` y `Dbp`) mostró una sobreexpresión marcada.

![Heatmap top 20 genes STAR paired-end edgeR](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/star/paired_end/edgeR_results/figuras/heatmap_top_20_genes.png)
**Figura 16. Heatmap de los 20 genes diferencialmente expresados más significativos para STAR (paired-end, edgeR)**

El heatmap para edgeR reflejó la misma estructura de separación por condición, con réplicas agrupadas por edad y bloques cromáticos invertidos entre `mm_9` y `mm_24`.

- Genes observados en rojo (relativamente sobreexpresados en `mm_24`): Dbp, Rasd2.
- Genes observados en azul (relativamente subexpresados en `mm_24`): Col4a1, Hr, Clic5, Rimoc1, Ppp1r3c, Ythdf3, Asb1, Ibtk, Arl5a, Pdgfra, Egfr, Thbd, Smox, Serpinf1, Nrep, Tfrc, Ogn, Gsn.

Visualmente, edgeR presentó una gran coincidencia con la imagen generada por DESeq2: la mayoría de los genes top aparecieron subexpresados en `mm_24` (azul) y solo un par de genes muestrarón sobreexpresión (rojo). Entre los dos métodos aplicados sobre STAR, coincidierón de forma explícita genes como `Dbp` y `Rasd2` (sobreexpresados en `mm_24`) y un amplio conjunto de genes subexpresados (por ejemplo `Col4a1`, `Clic5`, `Rimoc1`, `Ppp1r3c`, `Ythdf3`, `Asb1`, `Ibtk`, `Arl5a`, `Thbd`, `Smox`, `Nrep`, `Tfrc`). EdgeR pareció detectar una represión más extensa o consistente en el subconjunto top-20 (es decir, incorpora varios genes adicionales en azul) y también destacó genes propios que no fueron incluidos en la lista top de DESeq2 para STAR (por ejemplo `Hr`, `Pdgfra`, `Egfr`, `Serpinf1`, `Gsn`, `Ogn`), lo que sugierió que edgeR puede ser más sensible o más selectivo hacia ciertos patrones de disminución de expresión en este conjunto.

![Heatmap top 20 genes Salmon paired-end DESeq2](https://raw.githubusercontent.com/Andttrea/transcriptomica/refs/heads/main/deseq/results/salmon/paired_end/deseq_results/figuras/heatmap_top_20_genes.png)
**Figura 17. Heatmap de los 20 genes diferencialmente expresados más significativos para Salmon (paired-end, DESeq2)**

El heatmap obtenido a partir de Salmon mantiene la separación entre edades pero presentó una heterogeneidad de intensidades algo mayor que los alineadores genómicos.

- Genes observados en rojo (relativamente sobreexpresados en `mm_24`): Rasd2, Ube2srt.
- Genes observados en azul (relativamente subexpresados en `mm_24`): Ppp1r3c, Thbd, Ulk2, Arl5a, Ythdf3, Asb1, Hif1an, Serpinh1, Loxl2, Ubs, Ifi205, Ms4a4a, Ldha, Cidec, Myh4, Ogn, Gfpt2, Gvin1.

Visualmente, Salmon compartió con STAR varias anotaciones en azul que indican represión en `mm_24` (por ejemplo `Ppp1r3c`, `Thbd`, `Ulk2`, `Arl5a`, `Ythdf3`, `Asb1`, `Hif1an`, `Serpinh1`, `Ogn`), y también coincidió en `Rasd2` como uno de los genes sobreexpresados en `mm_24`. No obstante, Salmon introdujó genes diferenciales que no aparecen en las listas top-20 de STAR (por ejemplo `Ube2srt`, `Loxl2`, `Ubs`, `Ifi205`, `Ms4a4a`, `Ldha`, `Cidec`, `Myh4`, `Gfpt2`, `Gvin1`), lo que sugierió que el pseudoalineamiento detectó con mayor sensibilidad (o simplemente con distinto ranking) ciertos genes que los alineadores genómicos no colocaron entre los top-20.


Se observó un patrón dominante en los tres heatmaps, en donde predominó la represión relativa en `mm_24` (mayor proporción de genes en azul), por lo que la tendencia general apuntó a una disminución de la expresión de los genes top asociados a los ratones de 24 meses.

Por otra parte, hubo genes conservados entre métodos/alineadores, `Rasd2` y `Dbp` aparecierón repetidamente como sobreexpresados en `mm_24` en STAR (tanto DESeq2 como edgeR) y `Rasd2` también se observaron sobreexpresados en Salmon (DESeq2); genes como `Ppp1r3c`, `Thbd`, `Ulk2`, `Arl5a`, `Ythdf3`, `Asb1`, `Hif1an`, `Serpinh1`, `Clic5`, `Col4a1`, `Rimoc1`, `Ibtk`, `Nrep` y `Tfrc` se conservarón repetidamente en la categoría de subexpresión en `mm_24` entre STAR y Salmon o entre DESeq2 y edgeR aplicados sobre STAR.
Tambien hubo genes que cambiarón según método/alineador, edgeR detectó varios genes en azul que no figuraron en el top-20 de DESeq2 para STAR (por ejemplo `Hr`, `Pdgfra`, `Egfr`, `Serpinf1`, `Gsn`), mientras que Salmon puso en rojo `Ube2srt` y resaltó genes en azul exclusivos (por ejemplo `Loxl2`, `Ubs`, `Ifi205`, `Ms4a4a`, `Ldha`, `Cidec`, `Myh4`, `Gfpt2`, `Gvin1`), indicando selecciones dependientes del flujo.

Entre los tres heatmaps analizados se observó que STAR (especialmente cuando se comparan DESeq2 y edgeR) mostró los patrones más consistentes y definidos, con muchos genes coincidentes y una visualización compacta. edgeR tendió a ampliar la lista de genes subexpresados en `mm_24`, sugiriendo detecciones más extensas de represión en el subconjunto top. Salmon mantuvo la tendencia global hacia la represión en `mm_24`, pero mostró un perfil levemente diferente y detectó varios genes adicionales que no figuraban en las listas top de los alineadores genómicos.


## Discusión

**Efectos del alineador y la estrategia de cuantificación**

El análisis comparativo de los tres alineadores revela patrones sistemáticos en la detección de expresión diferencial. Los alineadores genómicos, HISAT2 y STAR, produjeron resultados notablemente similares, tanto en el número total de genes diferencialmente expresados como en los patrones visuales de PCA y volcano plots. Esto es consistente con la literatura que indica que ambas herramientas, aunque implementan algoritmos de alineamiento distintos, logran recuperar la mayoría de la señal biológica subyacente cuando se aplican a datos de calidad similar. En el presente estudio, HISAT2 detectó entre 259 y 371 genes diferencialmente expresados dependiendo del método estadístico y el tipo de librería, mientras que STAR identificó entre 259 y 376, lo que representa una variación menor al 5%. Esta concordancia sugiere que la elección entre HISAT2 y STAR no constituye un factor determinante en la calidad o cantidad de genes detectados cuando ambos se aplican correctamente al mismo conjunto de datos.

Por el contrario, el pseudoalineador Salmon exhibió un comportamiento marcadamente diferente, detectando entre 38 y 108 genes diferencialmente expresados, lo que representa una reducción del 65-90% en comparación con los alineadores genómicos. Aunque todos los métodos mantuvieron una separación clara entre edades en el espacio de PCA, Salmon mostró mayor dispersión de las réplicas dentro de cada grupo de edad, particularmente en modo paired-end. Este comportamiento puede atribuirse a varias razones técnicas: primero, Salmon cuantifica a nivel de transcritos y luego se agrega a nivel de gen mediante `tximport`, lo que introduce una capa adicional de incertidumbre al consolidar múltiples isoformas. Segundo, los pseudoalineadores emplean algoritmos de k-mers y búsqueda rápida que pueden ser más sensibles a la complejidad de la biblioteca de transcriptomas, especialmente cuando hay ambigüedad en la asignación de reads entre genes o transcritos homólogos. Tercero, la ausencia de información de alineamiento a escala nucleotídica puede resultar en una menor precisión en casos de sobreposición génica o regiones repetidas. Sin embargo, es importante destacar que los genes detectados por Salmon incluyen varios de los denominados "top" detectados por HISAT2 y STAR, como `Rasd2` y varios genes subexpresados en `mm_24` (por ejemplo, `Ppp1r3c`, `Thbd`, `Ulk2`, `Arl5a`, `Ythdf3`, `Asb1`, `Hif1an`, `Serpinh1`), lo que sugiere que aunque Salmon es menos sensible en términos de cantidad de genes, la señal biológica principal se conserva.

**Influencia de los métodos estadísticos**

DESeq2 y edgeR representan dos enfoques complementarios para el análisis de expresión diferencial basados en modelos de distribución binomial negativa. En este estudio, ambos métodos mantuvieron la estructura biológica fundamental del conjunto de datos, como se evidencia en los análisis de PCA donde la separación por edad es igualmente clara con ambos métodos. Sin embargo, se observaron diferencias cuantitativas consistentes: DESeq2 detectó un promedio de 19 genes adicionales por comparación en relación a edgeR, principalmente en la categoría de genes upregulated. En los casos de HISAT2 paired-end, DESeq2 identificó 62 genes upregulated frente a 24 en edgeR; para STAR paired-end, la cifra fue 59 versus 25; y para Salmon paired-end, 7 versus 4.

Estas diferencias se derivan de los distintos marcos de normalización y estimación de dispersión. DESeq2 utiliza factores de tamaño geométricos y un método de estimación de dispersión que combina información empírica con información paramétrica, lo que puede resultar en umbrales de significancia ligeramente más permisivos. Por su parte, edgeR aplica normalización TMM (Trimmed Mean of M-values), que es más robusta a genes altamente expresados, y en este estudio se configuró con un filtrado más estricto (`min.count = 20`), lo que redujo el número de genes candidatos antes del análisis y posiblemente aumentó la selectividad. Los volcano plots revelan que DESeq2 tiende a producir nubes más densas de puntos significativos, mientras que edgeR genera distribuciones más compactas pero con algunos puntos más extremos. Esto es particularmente evidente en comparaciones de STAR, donde el volcano plot de DESeq2 muestra una dispersión más amplia que el de edgeR, aunque ambos mantienen la asimetría característica hacia la downregulation en `mm_24`.

La selección entre ambos métodos implica un balance entre sensibilidad y especificidad. DESeq2 aparenta ser ligeramente más sensible a cambios de expresión moderados, mientras que edgeR tiende a ser más selectivo, priorizando cambios de mayor magnitud. Ambos enfoques resultan en listas de genes con considerable sobreposición (como lo demuestran los heatmaps, donde genes como `Dbp`, `Rasd2` y la mayoría de genes downregulated aparecen en ambos métodos), pero con diferencias en el ranking y la sensibilidad. Por tanto, la recomendación metodológica sería aplicar ambos métodos cuando sea posible y priorizar los genes consenso como más robustos.

**Efecto del tipo de librería: paired-end vs single-end**

El modo de secuenciación (paired-end frente a single-end) mostró un impacto diferenciado según el alineador empleado. Para HISAT2 y STAR, las diferencias fueron comparativamente menores: el número de genes detectados varió entre 15-20 genes en DESeq2 y 27-50 genes en edgeR, representando cambios relativos menores al 10% en la mayoría de casos. Los análisis de PCA evidenciaron que ambos modos de secuenciación mantenían patrones de agrupamiento muy similares, aunque los volcano plots de single-end con edgeR mostraron ligeramente mayor dispersión y menos definición de la nube de puntos significativos.

En contraste, Salmon demostró una sensibilidad mayor al tipo de librería. En DESeq2, los datos paired-end resultaron en 74 genes diferencialmente expresados frente a 108 en single-end, una tendencia inversa a la observada en los alineadores genómicos. En edgeR, la diferencia fue aún más marcada: 48 genes en paired-end contra 38 en single-end. Adicionalmente, los análisis de PCA para Salmon mostraron diferencias visuales más pronunciadas entre modos de librería, con mayor dispersión en single-end. Esta diferencia puede explicarse por la mayor complejidad del pseudoalineamiento con paired-end, donde la información adicional de distancia entre reads puede aumentar la ambigüedad en la asignación cuando los fragmentos son cortos o cuando existen genes similares. En single-end, la ausencia de información de pareja podría paradójicamente aumentar la eficiencia del pseudoalineamiento en este conjunto de datos específico, aunque esto debe interpretarse con cautela ya que single-end típicamente reduce la sensibilidad.

En general, el tipo de librería tuvo un efecto menor sobre los alineadores genómicos, sugiriendo que HISAT2 y STAR son robustos a esta variación. Sin embargo, el pseudoalineador Salmon mostró una dependencia mayor del modo de librería, lo que refuerza la observación anterior de que las características técnicas particulares de Salmon interactúan de manera más compleja con los parámetros experimentales.

**Concordancia entre métodos como indicador de robustez**

Un hallazgo importante es el alto grado de concordancia entre métodos para los genes "top". Los genes identificados consistentemente (como `Rasd2`, `Dbp`, `Serpinh1`, `Ppp1r3c`, `Thbd`) a través de diferentes alineadores y métodos estadísticos pueden considerarse como señales biológicas altamente robustas. Esta concordancia es especialmente evidente en los heatmaps de STAR con DESeq2 y edgeR, donde la mayoría de los 20 genes principales son los mismos o están presentes en ambos análisis. El hecho de que Salmon, a pesar de detectar menos genes totales, también recupere muchos de estos mismos genes "top" refuerza la confianza en que estos representan verdaderos cambios de expresión asociados con el envejecimiento.

Inversamente, genes que difieren entre métodos (por ejemplo, `Hr`, `Pdgfra`, `Egfr`, `Serpinf1`, `Gsn` detectados principalmente por edgeR en STAR; o `Loxl2`, `Ubs`, `Ifi205`, `Ms4a4a`, `Ldha`, `Cidec`, `Myh4`, `Gfpt2`, `Gvin1` únicos a Salmon) pueden representar señales más débiles o dependientes del método, y podrían requerir validación experimental adicional.


## Conclusión

El análisis de expresión diferencial realizado mediante múltiples pipelines de bioinformática ha demostrado que, a pesar de las variaciones técnicas introducidas por diferentes alineadores, métodos estadísticos y tipos de librería, existe un patrón biológico coherente y robusto asociado con el envejecimiento del músculo esquelético en ratones. Los alineadores genómicos HISAT2 y STAR produjeron resultados notablemente concordantes, detectando entre 250 y 380 genes diferencialmente expresados, mientras que el pseudoalineador Salmon identificó un subconjunto más pequeño pero enriquecido en genes altamente significativos. Los métodos estadísticos DESeq2 y edgeR mantuvieron la estructura biológica fundamental aunque con diferencias en sensibilidad, siendo DESeq2 ligeramente más permisivo y edgeR más selectivo. El tipo de librería (paired-end vs single-end) tuvo un impacto menor en los alineadores genómicos, pero mostró una influencia mayor en Salmon, indicando que la robustez de cada estrategia varía según el contexto metodológico.

El consenso entre métodos en genes como `Rasd2`, `Dbp`, `Ppp1r3c`, `Thbd`, `Serpinh1` y otros proporciona mayor confianza en que estos representan verdaderos cambios transcriptómicos asociados con el envejecimiento. La tendencia dominante observada es la represión génica en animales de 24 meses, sugiriendo una disminución de la capacidad anabólica y procesos reparadores en el músculo envejecido, consistente con los mecanismos conocidos de sarcopenia relacionada con la edad. Esta pauta se mantiene consistente a través de PCA, volcano plots y heatmaps, proporcionando múltiples líneas de evidencia que convergen hacia la misma conclusión biológica.

La comparación de múltiples pipelines ha permitido evaluar la robustez y reproducibilidad del análisis, demostrando que aunque existen diferencias cuantitativas en el número de genes detectados, las tendencias biológicas principales se conservan. Esto subraya la importancia de las herramientas bioinformáticas para el análisis de RNA-seq, ya que la elección de alineador, normalización y método estadístico impacta significativamente los resultados cuantitativos. Sin embargo, cuando múltiples enfoques independientes convergen en conclusiones similares, la confianza en la validez biológica de los hallazgos aumenta considerablemente.


## Referencias

1. Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. *Genome Biology*.
2. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*.
3. Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*.
4. Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*.
