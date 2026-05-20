# Reporte de anotación funcional

## Selección del dataset para GO y STRING

Para el análisis de enriquecimiento funcional en GO y la construcción de redes de interacción en STRING, se eligió como dataset principal el resultado obtenido con STAR en modo paired-end, priorizando la lista de genes diferencialmente expresados generada con DESeq2.

La elección se fundamentó en que este flujo mostró una separación biológica clara entre las condiciones de 9 y 24 meses, mantuvo réplicas compactas en el PCA y conservó una señal transcriptómica robusta en los volcano plots y heatmaps. Además, STAR y HISAT2 presentaron un comportamiento muy similar en el reporte de expresión diferencial, pero STAR ofreció una representación visual ligeramente más amplia y consistente de los genes significativos, lo que resulta útil para un análisis funcional posterior.

## Razones para elegir STAR paired-end con DESeq2

Se eligió STAR paired-end porque combinó tres ventajas relevantes para GO y STRING. En primer lugar, el modo paired-end mostró una estructura más consistente que el single-end, con una agrupación más definida de las muestras y una menor dispersión entre réplicas. En segundo lugar, STAR conservó una separación por edad tan clara como la observada con HISAT2, pero con una lectura visual más expansiva de la señal diferencial, lo que sugiere una recuperación más informativa de genes candidatos.

En tercer lugar, DESeq2 fue preferido como base principal porque recuperó una señal diferencial amplia y estable, con una cantidad de genes suficiente para realizar enriquecimiento funcional con mayor potencia estadística. Para GO y STRING resulta ventajoso partir de una lista de genes robusta, pero no excesivamente restrictiva, ya que esto permite detectar procesos biológicos y nodos de interacción con mayor capacidad interpretativa. En este contexto, DESeq2 ofreció una lista adecuada para ese propósito sin perder la coherencia biológica observada en los análisis previos.

## Razones para no elegir los otros datasets

### HISAT2 paired-end

HISAT2 paired-end también mostró resultados muy cercanos a los de STAR y constituyó una alternativa válida. Sin embargo, al compararlo con STAR, la señal visual de los genes diferencialmente expresados fue ligeramente menos expansiva y la interpretación global no aportó una ventaja clara sobre STAR. Por esa razón, aunque su desempeño fue sólido, no se seleccionó como opción principal.

### Los datasets single-end

Los flujos single-end fueron descartados porque presentaron una dispersión mayor y una estructura menos compacta en el PCA. Aunque conservaron la separación biológica entre edades, su menor definición visual y su menor consistencia entre réplicas los hicieron menos convenientes para un análisis funcional que requiere una lista de genes lo más estable posible.

### Salmon

Salmon no fue elegido como dataset principal para GO y STRING porque produjo una cantidad considerablemente menor de genes diferencialmente expresados y mostró mayor dispersión en las réplicas, especialmente en la comparación paired-end. Aunque recuperó varios genes biológicamente relevantes, su comportamiento fue más conservador y menos uniforme que el de los alineadores genómicos. Para un análisis funcional basado en enriquecimiento y redes, una lista demasiado reducida puede limitar la detección de términos y conexiones significativas.

### edgeR como lista principal

edgeR tampoco se seleccionó como base principal, no porque sus resultados fueran inconsistentes, sino porque tendió a ser más selectivo. Aunque en algunos casos recuperó un número total de DEGs más alto, lo hizo a costa de identificar muy pocos genes upregulated, lo que refleja un comportamiento más restrictivo. Esa selectividad es útil como validación complementaria, pero para GO y STRING resulta más conveniente iniciar con una lista que retenga una señal más amplia y balanceada de genes diferencialmente expresados. En ese sentido, DESeq2 ofreció un punto de partida más equilibrado para el análisis funcional.

## Conclusión

El dataset más apropiado para realizar GO y STRING fue STAR paired-end analizado con DESeq2. Esta elección se apoyó en su equilibrio entre robustez biológica, separación clara de las condiciones, buena compactación de réplicas y una señal diferencial suficientemente amplia para alimentar un análisis de enriquecimiento funcional. Los demás datasets conservaron valor como comparación y validación, pero STAR paired-end con DESeq2 ofreció la mejor combinación de estabilidad, interpretabilidad y potencia para el objetivo funcional posterior.

