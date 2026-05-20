# Reporte de anotación funcional

## Planteamiento

El análisis de anotación funcional permite interpretar listas de genes diferencialmente expresados a partir de categorías biológicas más amplias, como procesos biológicos, componentes celulares y redes de interacción entre proteínas. En este trabajo se partió del contraste entre músculo de *Mus musculus* de 9 meses y 24 meses, con el propósito de identificar qué funciones biológicas y qué estructuras celulares se asocian con el cambio transcripcional observado durante el envejecimiento del tejido de la extremidad.

Para este reporte se tomó como base el resultado obtenido con STAR en modo paired-end y el análisis estadístico realizado con DESeq2. Esa combinación se seleccionó porque en el reporte de expresión diferencial mostró una separación clara entre las edades, una distribución compacta de las réplicas y una lista de genes significativa suficientemente amplia para alimentar análisis de enriquecimiento y redes. En este contexto, la lista generada por DESeq2 resultó útil para explorar procesos biológicos con mayor poder interpretativo sin perder coherencia con la biología observada previamente.

El flujo de trabajo se concentró en tres herramientas principales. Primero, se utilizó DAVID para realizar el enriquecimiento por Gene Ontology. Después, se empleó STRING para construir redes de interacción proteína-proteína. Finalmente, se ejecutó un análisis GSEA para evaluar si conjuntos de genes completos se concentraban en los extremos de la lista ordenada por estadística diferencial.

Los datos, scripts y figuras se encuentran en el directorio `/export/space3/users/andreavg/transcriptomica/deseq/results/star/paired_end/deseq_results`.

## Objetivos

### Objetivo general

Realizar un análisis de anotación funcional sobre los genes diferencialmente expresados obtenidos con STAR paired-end y DESeq2, con el fin de identificar procesos biológicos, componentes celulares y redes funcionales asociadas al envejecimiento del músculo de la extremidad en *Mus musculus*.

### Objetivos particulares

1. Preparar listas limpias de identificadores Ensembl para su uso en DAVID y STRING.
2. Analizar el enriquecimiento funcional de los genes upregulated y downregulated mediante Gene Ontology en DAVID.
3. Explorar la localización celular y la coherencia biológica de los genes con términos CC_direct y BP_direct.
4. Visualizar los resultados de DAVID con gráficos tipo dot plot generados en R.
5. Construir redes de interacción proteica con STRING para los genes diferencialmente expresados.
6. Ejecutar un análisis GSEA sobre la lista ordenada de resultados de DESeq2 para identificar rutas enriquecidas a nivel global.

## Metodología y justificación

### Selección del dataset de trabajo

La primera decisión metodológica consistió en fijar STAR paired-end con DESeq2 como punto de partida para el análisis funcional. Esa selección se apoyó en el comportamiento observado en el reporte de expresión diferencial, donde ese flujo conservó una separación biológica clara entre mm_9 y mm_24, además de una señal robusta en PCA, volcano plots y heatmaps. A nivel práctico, esa combinación produjo una lista de genes significativa que resultó suficiente para enriquecer términos funcionales sin volver el análisis excesivamente restrictivo.

### Preparación de las listas de genes para DAVID

Las herramientas web de anotación funcional, como DAVID y STRING, necesitan listas puras de identificadores para poder mapear correctamente los genes a sus bases de datos internas. Por esa razón no se subieron tablas completas con columnas de `log2FoldChange`, `pvalue`, `padj` o promedios de expresión, sino archivos reducidos a una sola columna con identificadores Ensembl.

El trabajo de limpieza se realizó con comandos `awk` para conservar únicamente la primera columna de los archivos generados por DESeq2. De esta manera se construyó un archivo de background y dos listas separadas para genes upregulated y downregulated.

```bash
# Nos movemos a la carpeta correspondiente
cd /export/storage/users/andreavg/transcriptomica/deseq/results/star/paired_end/deseq_results

# Hacemos un awk para extraer únicamente la primera columna (los Ensembl IDs)
# NR>1 indica que se salte la primera línea
awk -F'\t' 'NR>1 {print $1}' TPM_log2-table.txt > background.txt

# Realizamos lo mismo para los archivos up y down regulated
awk -F'\t' 'NR>1 {print $1}' deseq-DEG_up_0.05.txt > deseq-DEG_up_clean.txt
awk -F'\t' 'NR>1 {print $1}' deseq-DEG_down_0.05.txt > deseq-DEG_down_clean.txt
```

Este paso fue necesario porque DAVID trabaja mejor cuando recibe una lista simple de IDs, ya que de ese modo puede mapear cada gen de forma directa a sus anotaciones funcionales y evitar errores de formato.

### Análisis de Gene Ontology en DAVID

DAVID por Database for Annotation, Visualization and Integrated Discovery, es una plataforma web diseñada para interpretar listas de genes mediante anotación funcional, enriquecimiento de categorías biológicas y visualización de términos. En este trabajo se utilizó para obtener términos de Gene Ontology sobre los conjuntos upregulated y downregulated.

Una vez cargados el `background` y las listas de genes diferenciados en la interfaz de DAVID, se descargaron los reportes de Gene Ontology correspondientes a `BP_direct` y `CC_direct`. El primero permitió recuperar procesos biológicos asociados al cambio transcripcional, mientras que el segundo sirvió para verificar si los genes eran coherentes con el tejido analizado y con estructuras celulares relevantes para músculo.

Para los genes upregulated, se relajó el `EASE Threshold` a `0.25` en el análisis de `CC_direct` porque la lista contenía solo 59 genes y con el threshold que ya se tenía no permitió encontrar nada. Esa decisión permitió explorar mejor las tendencias dominantes de localización celular sin imponer un filtro demasiado estricto que ocultara señales biológicamente informativas.

Los reportes descargados de DAVID se guardaron en la carpeta `figuras` dentro de `results/star/paired_end/deseq_results/GO`.

### Visualización de GO en R

Posteriormente se utilizó el script `src/plot_Go.R` para transformar los reportes de DAVID en gráficos tipo dot plot. Ese script lee los archivos CSV exportados por DAVID, calcula `-log10(P.Value)` para representar la significancia de forma más intuitiva, ordena los términos por valor de enriquecimiento y selecciona los 20 términos más relevantes para evitar sobrecargar la figura.

El gráfico se construye con `ggplot2`, donde el eje X representa el `Fold.Enrichment`, el eje Y los términos de GO, el tamaño del punto representa el número de genes asociados y el color refleja la significancia estadística. Finalmente, el script guarda las imágenes en la carpeta `figuras` correspondiente al análisis de GO.

En términos prácticos, este script funcionó como una capa de visualización sobre los resultados de DAVID y permitió comparar de forma compacta los términos más informativos para los genes upregulated, downregulated y para las categorías celulares exploradas con `CC_direct`.

### Análisis de STRING

STRING es una plataforma web que integra interacciones proteína-proteína conocidas y predichas, y permite visualizar redes funcionales, agrupamientos y conexiones entre genes de interés. En este trabajo se utilizó para evaluar si los genes diferencialmente expresados formaban módulos o subredes compatibles con funciones biológicas compartidas.

Se cargaron tres entradas distintas: la lista de genes upregulated, la lista de genes downregulated y el conjunto combinado de ambos. Además, se generaron imágenes en dos resoluciones para algunos conjuntos, en particular para `STRING_downregulated` y `STRING_all`, porque en la vista de alta densidad los clusters aparecían muy agrupados. La versión de resolución media y la de alta resolución permitieron conservar conexiones significativas y hacer más legible la estructura global de la red.

Las figuras de STRING se guardaron en la carpeta `STRING` dentro de `results/star/paired_end/deseq_results`.

### Análisis GSEA

El análisis GSEA, o Gene Set Enrichment Analysis, evalúa si conjuntos de genes predefinidos se concentran en los extremos de una lista ordenada por una métrica continua, en lugar de depender únicamente de un umbral de significancia por gen. Esa estrategia resulta útil porque aprovecha la información de todo el ranking generado por DESeq2 y no solo de los genes que superan un corte arbitrario.

El script `src/GSEA.R` cargó el archivo `DESeq2_results.rds`, lo convirtió a `data.frame`, eliminó los genes con valores `NA` en la columna `stat` y ordenó los genes de mayor a menor según dicho estadístico. Después construyó un vector nombrado con identificadores Ensembl y ejecutó `gseGO()` del paquete `clusterProfiler` usando `org.Mm.eg.db` y la ontología `BP`.

El resultado permitió identificar rutas enriquecidas a partir del ranking completo de DESeq2. El script también generó una figura resumen con `gseaplot2`, enfocada en las vías más representativas, y guardó la imagen en la carpeta `GSEA` dentro de `results/star/paired_end/deseq_results`.

### Organización de archivos y reproducibilidad

Para mantener la reproducibilidad del análisis, los scripts se conservaron en la carpeta `src` y las figuras generadas se almacenaron dentro de subdirectorios específicos del flujo STAR paired-end. Esa organización permitió separar claramente los insumos, los objetos intermedios y las salidas finales de GO, STRING y GSEA.

La estructura de trabajo quedó apoyada en cuatro elementos principales: la tabla de genes diferenciales generada por DESeq2, los archivos limpios con identificadores Ensembl, el script de visualización `plot_Go.R` y el script de enriquecimiento `GSEA.R`. Con ello, el análisis funcional quedó trazado de forma clara, desde la preparación de los datos hasta la generación de las figuras finales.

