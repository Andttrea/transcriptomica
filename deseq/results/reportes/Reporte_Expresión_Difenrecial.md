# Reporte de expresión diferencial

## Planteamiento

El análisis de expresión diferencial permite identificar genes cuya abundancia cambia entre condiciones biológicas distintas. En este trabajo se comparó tejido muscular de *Mus musculus* de 9 meses contra 24 meses, utilizando como punto de partida los archivos de cuantificación generados para tres estrategias de análisis: alineamiento con HISAT2, alineamiento con STAR y pseudoalineamiento con Salmon.

A diferencia del reporte de pseudoalineamientos, aquí el énfasis no está en obtener la cuantificación inicial, sino en preparar correctamente las matrices de conteo y aplicar dos enfoques estadísticos para detectar genes diferencialmente expresados: DESeq2 y edgeR. El objetivo fue evaluar cómo cambian los resultados según el alineador, el modo de secuenciación (single-end o paired-end) y el algoritmo estadístico utilizado.

Los datos, scripts y resultados se encuentran en el directorio `/export/storage/users/andreavg/transcriptomica/deseq`.

## Objetivos

### Objetivo general

Comparar los resultados de expresión diferencial obtenidos a partir de HISAT2, STAR y Salmon mediante los paquetes DESeq2 y edgeR, usando datos de RNA-seq de músculo de ratón en dos edades biológicas.

### Objetivos particulares

1. Integrar los archivos de `featureCounts` en una matriz global de conteos y una tabla de longitudes por gen para HISAT2 y STAR.
2. Estandarizar los nombres de las muestras sustituyendo los identificadores SRR por sus condiciones biológicas para facilitar la interpretación de los resultados.
3. Justificar y reutilizar el objeto `txi` de Salmon como entrada para el análisis a nivel de gen.
4. Ejecutar análisis de expresión diferencial con DESeq2 y edgeR para los datos de HISAT2, STAR y Salmon en modos single-end y paired-end.
5. Comparar los resultados obtenidos mediante PCA, volcano plots y heatmaps de genes significativos.

## Metodología y justificación

### Preparación de tablas de conteo

El primer paso consistió en transformar la salida de `featureCounts`, que originalmente está separada por muestra, en tablas globales con todas las muestras reunidas. Esto se realizó en el script [counts_lenght_table.R](../../src/counts_lenght_table.R), donde se unieron los conteos por gen y también se extrajeron las longitudes de los genes. Este paso fue necesario para HISAT2 y STAR, ya que para Salmon se trabajó con el objeto `txi` generado previamente por `tximport`.

Posteriormente, el script [condiciones_tabla_counts.R](../../src/condiciones_tabla_counts.R) reemplazó los nombres SRR por condiciones más interpretables, usando la tabla de metadatos del experimento. Para cada muestra se construyó una etiqueta con la edad y el sexo, por ejemplo `mm_9` y `mm_24`, y se añadieron sufijos automáticos para distinguir réplicas.

### Justificación del uso de `txi` en Salmon

Aunque Salmon cuantifica a nivel de transcritos, en este análisis se trabajó a nivel de gen gracias a `tximport`. El archivo `.rds` generado por este flujo ya resuelve el problema de múltiples transcritos por gen, porque la tabla `tx2gene` le indica a `tximport` qué transcritos pertenecen a cada gen. Con ello, los conteos de todas las isoformas se suman en un único registro por GENEID.

Como resultado, las matrices de `counts` y `abundance` dentro del objeto `txi` contienen filas identificadas por genes y no por transcritos individuales. Además, `tximport` estima una longitud efectiva por gen considerando qué isoformas se expresaron realmente, lo cual es importante porque la longitud aparente de un gen puede cambiar según las isoformas activas. Por eso, el objeto `txi` es una entrada válida y consistente para el análisis diferencial a nivel de gen.

### Análisis con DESeq2

El script [Deseq2_hisat_star.R](../../src/Deseq2_hisat_star.R) realizó el análisis para HISAT2 y STAR, mientras que [Deseq2_salmon.R](../../src/Deseq2_salmon.R) hizo lo propio para Salmon. En ambos casos el flujo fue similar:

- creación del objeto de DESeq2 desde una matriz de conteos o desde `tximport`;
- filtrado de genes de baja expresión con `filterByExpr`;
- transformación `vst` para construir un PCA;
- ajuste del modelo con `DESeq()`;
- contraste entre `mm_24` y `mm_9`;
- cálculo de una tabla de resultados con log2 fold change y FDR;
- generación de volcano plots y heatmaps para los genes significativos.

En estos scripts se usaron umbrales de significancia de `FDR < 0.05` y `|log2FC| > 0.5`. Los genes up y down regulados se guardaron en archivos separados, junto con mapas de TPM en escala log2 para apoyar la visualización.

### Análisis con edgeR

El script [EdgeR_hisat_star.R](../../src/EdgeR_hisat_star.R) aplicó el mismo planteamiento general, pero usando la estructura de edgeR (`DGEList`) y su flujo estadístico. Para Salmon se empleó [EdgeR_salmon.R](../../src/EdgeR_salmon.R), también sobre el objeto `txi` ya consolidado.

La secuencia de trabajo fue:

- construcción de `DGEList`;
- filtrado con `filterByExpr` y `min.count = 20`;
- normalización TMM con `calcNormFactors`;
- cálculo de `logCPM` y PCA;
- estimación robusta de dispersión;
- ajuste con `glmFit`;
- contraste con `glmLRT`;
- exportación de genes significativos, volcano plots y heatmaps.

### Diferencias entre DESeq2 y edgeR

Ambos paquetes modelan conteos con distribuciones binomiales negativas y sirven para detectar expresión diferencial, pero no implementan exactamente el mismo flujo:

- DESeq2 estima factores de tamaño y usa su propio marco de normalización y ajuste del modelo.
- edgeR usa normalización TMM y, en este trabajo, se evaluó con LRT tras estimar dispersión robusta.
- edgeR fue configurado aquí con un filtrado más estricto (`min.count = 20`), por lo que retuvo menos genes antes del ajuste.

En la práctica, esto hace que edgeR y DESeq2 no devuelvan exactamente la misma lista de genes, aunque sí suelen recuperar una fracción grande de la misma señal biológica.

## Resultados

### Resumen de genes retenidos y significativos

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

En los PCA, los tres enfoques separan claramente las muestras por edad, lo que indica que el efecto biológico de 9 meses vs 24 meses domina la variación principal. En HISAT2 y STAR, la separación es muy consistente entre métodos y modos, mientras que en Salmon la dispersión entre réplicas es más evidente y la separación depende algo más del algoritmo usado.

En los volcano plots, HISAT2 y STAR muestran una nube amplia de genes significativamente down regulados en 24 meses, con un grupo menor de genes up regulados. Salmon conserva el patrón general de diferenciación por edad, pero con menos puntos significativos en el caso de DESeq2 y una lista más compacta en edgeR.

Los heatmaps confirman esta tendencia. En HISAT2 y STAR se repiten genes como `Dbp`, `Rasd2`, `Col4a1`, `Egfr`, `Thbd`, `Smox`, `Ogn`, `Nrep`, `Ppp1r3c`, `Asb1`, `Ythdf3`, `Hif1an` y `Arl5a`, con un patrón claro de expresión opuesta entre las dos edades. En Salmon aparecen varios de estos genes, pero también cambia parte del top 20, con más genes poco caracterizados o específicos del cuantificador, lo cual sugiere una sensibilidad distinta del pseudoalineamiento.

### Similitud entre listas de genes

Al comparar las listas de genes significativos, se observó que:

- DESeq2 y edgeR comparten una fracción moderada de genes dentro del mismo alineador y modo, con Jaccard entre aproximadamente 0.52 y 0.63.
- HISAT2 y STAR producen listas muy parecidas entre sí, con Jaccard superior a 0.90 en DESeq2 y cercano a 0.82-0.95 en edgeR.
- Salmon se separa mucho más de HISAT2 y STAR, con Jaccard cercano a 0.16-0.17 frente a los análisis con alineamiento genómico.
- La comparación paired-end vs single-end en HISAT2 y STAR muestra una similitud alta, alrededor de 0.75-0.78.
- En Salmon, el cambio entre paired-end y single-end afecta más, con similitud más baja, alrededor de 0.30-0.37.

## Discusión

Usando los mismos cutoffs, el número de genes obtenidos cambia de forma clara entre métodos. En HISAT2 y STAR se obtienen alrededor de 259-376 DEGs según el algoritmo y el modo, mientras que en Salmon el número baja notablemente en DESeq2 y edgeR, especialmente en paired-end. Esto muestra que el origen de la cuantificación influye mucho más que pequeñas variaciones en el flujo estadístico.

La lista de genes es similar, pero no idéntica. Entre DESeq2 y edgeR existe una coincidencia apreciable, aunque edgeR tiende a recuperar más genes down regulados y DESeq2 algo más de genes up regulados. Esto es coherente con las diferencias de modelo, normalización y estrategia de filtrado entre ambos paquetes. Además, edgeR se aplicó con `min.count = 20`, lo que hizo el filtrado más estricto que en DESeq2.

Si se compara el efecto de cada factor, el que más altera la lista final de genes es la estrategia de cuantificación/alineamiento, especialmente cuando se contrasta Salmon con los alineadores genómicos HISAT2 y STAR. HISAT2 y STAR son muy consistentes entre sí, por lo que el alineador genómico específico no parece ser el factor dominante. En segundo lugar aparece el algoritmo de expresión diferencial, porque DESeq2 y edgeR sí cambian una porción importante de la lista. El tipo de biblioteca, paired-end o single-end, tiene un efecto menor en HISAT2 y STAR y más visible en Salmon, pero en general no supera al cambio entre pseudoalineamiento y alineamiento genómico.

En términos biológicos, el resultado más robusto es la separación entre 9 y 24 meses, visible en todos los PCA y volcano plots. Eso indica que el envejecimiento del músculo deja una huella transcriptómica fuerte y reproducible. Los genes que más se repiten en los métodos de alineamiento genómico están asociados con remodelación tisular, matriz extracelular, señalización y cambios del estado celular, lo que es consistente con un tejido muscular envejecido.

## Conclusión

El análisis de expresión diferencial permitió confirmar una señal clara de envejecimiento en músculo de ratón y mostró que HISAT2 y STAR generan resultados muy cercanos entre sí. DESeq2 y edgeR recuperan una base común importante, pero no exactamente la misma lista de genes. Salmon conserva la tendencia biológica principal, aunque produce listas más distintas y menos compartidas con los alineadores genómicos.

En conjunto, el factor que más impactó los genes identificados fue la estrategia de cuantificación/alineamiento, seguido por el algoritmo estadístico, mientras que el efecto paired-end vs single-end fue menor en HISAT2 y STAR y más variable en Salmon.

## Referencias

1. Anders, S., & Huber, W. (2010). Differential expression analysis for sequence count data. *Genome Biology*.
2. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*.
3. Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*.
4. Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*.
