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