library(dplyr)
library(ggplot2)
library(stringr)

# Se creó una función para generar y guardar el dot plot de enriquecimiento funcional
# NOTA: Recibirá la ruta del archivo de DAVID, la ruta de salida de la imagen, el título del gráfico 
generar_go_dotplot <- function(file_path, output_path, plot_title) {
  
  # Leemos el archivo reportado por DAVID
  datos_go <- read.table(file_path, header = TRUE, sep = ",")
  
  # Procesamiento y limpieza de los datos
  datos_procesados <- datos_go %>%
    # Calculamos el -log10 del P.Value para que las vías más significativas tengan un valor alto y visualmente sea más intuitivo
    mutate(log10_pval = -log10(`P.Value`)) %>%
    # Ordenamos por `P.Value` (de menor a mayor, es decir, más significativos primero)
    arrange(`P.Value`) %>%
    # Seleccionamos solo el top 20 de vías para evitar sobrecargar la gráfica y que sea legible
    head(20)
  
  # Creación del gráfico con ggplot2
  plot <- ggplot(datos_procesados, aes(x = Fold.Enrichment, y = reorder(Term, Fold.Enrichment))) +
    # Dibujamos los puntos. El tamaño dependerá del conteo de genes y el color del -log10(P.Value)
    geom_point(aes(size = Count, color = log10_pval)) +
    # Configuración de colores: un gradiente donde los colores más encendidos (rojos/azules) indiquen mayor significancia
    scale_color_gradient(low = "blue", high = "red") +
    # Títulos y etiquetas de los ejes
    labs(
      title = plot_title,
      x = "Fold Enrichment",
      y = "Términos de GO (Biological Process)",
      size = "Conteo de Genes",
      color = "-log10(P.Value)"
    ) +
    # Aplicamos un tema limpio e ideal para publicaciones académicas
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 10, face = "plain"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold")
    )
  
  # Guardamos la gráfica en la ruta especificada
  # NOTA: Se define un tamaño estándar
  ggsave(output_path, plot = plot,  width = 11.25, height = 7.5, dpi = 600)
  

}

base <- "/export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/deseq_results/GO"

# Ejecutamos la función para los archivos up y down
# NOTA: Usamos file.path() para concatenar la base con las subcarpetas y nombres de archivos de forma segura

# Graficamos para los genes Upregulated 
generar_go_dotplot(
 file_path   = file.path(base, "DAVIDChartReport_deseq-DEG_up_clean_2026-05-19.csv"), 
  output_path = file.path(base, "figuras", "dotplot_GO_upregulated.png"), 
  plot_title  = "UPREGULATED_GO_DOTPLOT"
)

# Graficamos para los genes Downregulated 
generar_go_dotplot(
 file_path   = file.path(base, "DAVIDChartReport_deseq-DEG_down_clean_2026-05-19.csv"), 
 output_path = file.path(base, "figuras", "dotplot_GO_downregulated.png"), 
  plot_title  = "DOWNREGULATED_GO_DOTPLOT"
)

generar_go_dotplot(
  file_path   = file.path(base, "DAVIDChartReport_deseq-DEG_down_clean_2026-05-20_CC.csv"), 
  output_path = file.path(base, "figuras", "dotplot_GO_downregulated_CC.png"), 
  plot_title  = "DOWNREGULATED_GO_DOTPLOT_CC"
)

generar_go_dotplot(
  file_path   = file.path(base, "DAVIDChartReport_deseq-DEG_up_clean_2026-05-20_CC.csv"), 
  output_path = file.path(base, "figuras", "dotplot_GO_upregulated_CC.png"), 
  plot_title  = "UPREGULATED_GO_DOTPLOT_CC"
)