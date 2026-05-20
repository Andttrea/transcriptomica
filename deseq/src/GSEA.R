# Cargar librerías necesarias
library(DESeq2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot) # Necesaria para gseaplot2 (gráficos más estéticos)
library(ggplot2)    # Para guardar las gráficas

# Definimos rutas de entrada y salida
input_file <- "/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/deseq_results/DESeq2_results.rds"
output_dir <- "/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/deseq_results/GSEA"

# Cargamos los resultados de DESeq2
res <- readRDS(file = input_file)

# Convertir a data.frame por si viene como objeto DESeqResults
res_df <- as.data.frame(res)

# Preparamos  la lista de genes (gene_list) y filtramos los valores NA en la columna 'stat'. 
# NOTA: DESeq2 suele generar NAs por conteos bajos o valores atípicos, y si pasamos un NA a clusterProfiler, la función fallará
res_df <- res_df[!is.na(res_df$stat), ]

# Ordenamos de mayor a menor basado en el estadístico (stat)
# Esto es esencial para GSEA, ya que evalúa si los genes de una vía se acumulan en los extremos (up o down regulados)
res_df <- res_df[order(-res_df$stat), ]

# Extraer el vector numérico y asignarle los nombres de los genes (Ensembl IDs)
gene_list <- res_df$stat
names(gene_list) <- rownames(res_df)

# 4. Ejecutar el Gene Set Enrichment Analysis (GSEA)
gse <- gseGO(
  geneList     = gene_list,
  ont          = "BP",            # Biological Process
  keyType      = "ENSEMBL",
  OrgDb        = org.Mm.eg.db,    # Base de datos para Mus musculus
  eps          = 1e-300,
  pvalueCutoff = 0.05             # Umbral de significancia
)

# Generamos el gseaplot2 múltiple
# NOTA: Seleccionamos los geneSetID c(1, 43, 392) porque, en conjunto, capturan la fisiopatología 
# de la sarcopenia en la transición de 9 a 24 meses. Se eligieron estas 3 vías ortogonales para evitar redundancia:
# - ID 1 (extracellular matrix organization): Indica el colapso del andamiaje e integridad estructural del tejido.
# - ID 43 (muscle cell proliferation): Refleja la pérdida de la capacidad regenerativa de los mioblastos.
# - ID 392 (response to hypoxia): Muestra la incapacidad del músculo viejo para adaptarse al estrés celular y metabólico.
# Juntas, estas vías (con enriquecimiento negativo) cuentan la historia del declive funcional del músculo.

# Extraemos el nombre del pathway para el título
top_pathway_title <- gse$Description[1] 
p1 <- gseaplot2(gse, geneSetID = c(1, 43,392), title = top_pathway_title, pvalue_table = TRUE)

# Guardar el plot de GSEA en tu carpeta
ggsave(
  filename = file.path(output_dir, "GSEA_top1_pathway.png"),
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white" 
)



