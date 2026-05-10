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