# Cargamos las librerias que utilizaremos 
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(edgeR)
library(limma)

# Hacemos una funsión que realice todo el análisis con edgeR

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
  all(colnames(gene_counts) == rownames(meta_data))
  # Si arroja TRUE, entonces esta todo correcto

    # Cargamos los datos en un objeto edgeR
    # NOTA: Usamos round() para redondear los valores de conteo a enteros.
    # DGEList es el contenedor base en edgeR para conteos, factores y metadatos.
    # Design lo dejamos como ~ 0 + age para que el modelo no tenga intercepto y podamos comparar directamente las condiciones de edad.
    dge <- DGEList(counts = round(gene_counts))
    design <- model.matrix(~ 0 + age)

  # FIltramos los genes que contienen una baja expresión, esto que lo haremos con filterByExpr() de edgeR, que nos permite filtrar los genes que no tienen una expresión suficiente para ser considerados en el análisis diferencial.
  # filterByExpr compara la expresión de cada gen con un umbral de expresión mínima para eliminar los genes que no tienen una expresión suficiente 
  # NOTA: Le pasamos el diseño para que el filtro sea consistente con las condiciones que vamos a comparar.
  keep <- filterByExpr(dge, design)
  # Vemos cuantos genes se mantienen después del filtrado
  suma_keep <- sum(keep)
  # Guardamos el número de genes que se mantienen después del filtrado en un archivo de texto
  write.table(suma_keep, paste0(output_dir, "genes_filtrados.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  # Filtramos el objeto DESeq2 para quedarnos solo con los genes que cumplen con el criterio de expresión suficiente
    dge <- dge[keep, , keep.lib.sizes = FALSE]
  # Borramos keep 
  rm(keep)

  # Generamos el PCA para visualizar la variabilidad de los datos 
  # Calculamos logCPM (conteos por millon en log2) para estabilizar la varianza y hacer comparables las muestras.
  # prior.count evita log2(0) y reduce el efecto de genes con muy bajo conteo.
  log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
  # prcomp espera muestras en filas, por eso transponemos la matriz
  pca_res <- prcomp(t(log_cpm))
  pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], age = meta_data$age)
  # Hacemos el PCA, definiendo age para colorear por edad
  PCA_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = age)) +
    geom_point(size = 3) +
    theme_classic(base_size = 25, base_line_size = 1)

  # Guardamos la imagen del PCA
  ggsave(paste0(fig_dir, "PCA_plot.png"), plot = PCA_plot, width = 8, height = 6)

    # Calculamos los factores de normalización, dispersión y ajustamos el modelo 
    # calcNormFactors aplica la normalización TMM para corregir diferencias en profundidad de secuenciación.
    dge <- calcNormFactors(dge)
    # estimateDisp calcula la dispersión por gen (biológica) necesaria para el modelo NB de edgeR.
    dge <- estimateDisp(dge, design)
    # glmQLFit ajusta el modelo con quasi-likelihood, más robusto para tamaños de muestra pequeños.
    fit <- glmQLFit(dge, design)

  # Ahora calcularemos la TPM (Transcripts Per Million) que calcula cuantos transcritos tendrias si tuvieras un millon de transcritos en total
    # Añadimos la longitud del gen
    # NOTA: Usamos rownames(dge) para asegurarnos que la logitud coincida con el gen correcto
    gene_length <- annotation[rownames(dge), 1]

    # Calculamos FPKM y pasamos a log2. FKPM (Fragments Per Kilobase of transcript per Million mapped reads) es una medida de expresión que normaliza por la longitud del gen y por el número total de lecturas mapeadas
    # NOTA: En edgeR usamos rpkm() y pasamos la longitud del gen.
    # Le añadimos un pseudoconteo de 0.1 para evitar problemas con los logaritmos de cero.
    log2_fpkm <- log2(rpkm(dge, gene.length = gene_length) + 0.1) 

  # Escribimos la formula para convertir de FKPM a TPM dentro de un epsacio lograrítmico 
  fpkm2tpm_log2 <- function(fpkm) { fpkm - log2(sum(2^fpkm)) + log2(1e6) } 
  # Aplicamos la formula a cada columna de la tabla
  log2_tpm <- apply(log2_fpkm, 2, fpkm2tpm_log2) 

  # Guardamos la tabla 
  gene_names <- gene_name_map[rownames(log2_tpm),]
  write.table(cbind(gene_names, log2_tpm), paste0(output_dir, "TPM_log2-table.txt"), sep="\t", quote=FALSE)

  # Hacemos el  contraste de expresión diferencial 
  # Checamos los nombres 
  colnames(design)
  # Hacemos el contraste entre nuestras condiciones, en este caso entre mm_24 y mm_9
  contrast <- makeContrasts(m24_vs_m9 = agemm_24 - agemm_9, levels = design) 

  # Hacemos el análisis de expresión diferencial 
  # glmQLFTest aplica el contraste sobre el modelo ajustado.
  qlf <- glmQLFTest(fit, contrast = contrast[, "m24_vs_m9"])
  # topTags devuelve una tabla con logFC, FDR, etc. Usamos n = Inf para no truncar resultados.
  res <- topTags(qlf, n = Inf)$table
  # Renombramos para mantener el resto del pipeline igual que en DESeq2
  res$log2FoldChange <- res$logFC
  res$padj <- res$FDR
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
  vpcolors = c("gray", "#6204a9", "#509a05") 
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
          labs(title = "Volcano plot", 
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
                  col = list(Edad = c("mm_9" = "lightblue", "mm_24" = "blue"))
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
                  col = list(Edad = c("mm_9" = "lightblue", "mm_24" = "blue"))
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
hisat_paired <- edgeR_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/paired_end/featurecounts/global/hisat_pe_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "hisat2/paired_end/edgeR_results/"))
hisat_single <- edgeR_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/hisat2/single_end/featurecounts/global/hisat_se_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "hisat2/single_end/edgeR_results/"))
star_paired <- edgeR_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/featurecounts/global/star_pe_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "star/paired_end/edgeR_results/"))
star_single <- edgeR_analysis(read.delim("/export/storage/users/andreavg/transcriptomica/deseq/results/star/single_end/featurecounts/global/star_se_counts_condicion.tsv", row.names = 1), annotation, gene_name_map, paste0(base, "star/single_end/edgeR_results/"))
