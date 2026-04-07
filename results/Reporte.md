# Análisis de Alineamientos de RNA-seq

## Introducción

Este reporte documenta el análisis bioinformático de secuencias de RNA-seq, incluyendo la descarga de datos, control de calidad, limpieza y alineamiento del genoma. Todos los datos y resultados se encuentran ubicados en el servidor chaac de la licenciatura, en el directorio `/export/storage/users/andreavg/transcriptomica`.

## Estructura del Repositorio

```
transcriptomica/
├── data/
│   ├── clean/                   # Datos limpios después de procesamiento con fastp
│   ├── fastqc/                  # Reportes de control de calidad inicial
│   │   └── multiqc_report_data/ # Datos consolidados de MultiQC
│	│		└── images			 # Plots de multiqc
│   ├── metadata/                # Metadatos y scripts auxiliares
│   └── SRRs/                    # Secuencias raw descargadas
├── results/
│   ├── hisat2/                  # Alineamientos con HISAT2
│   │   ├── paired_end/          # Alineamientos paired-end
│   │   └── single_end/          # Alineamientos single-end
│   ├── star/                    # Alineamientos con STAR
│   │   ├── paired_end/          # Alineamientos paired-end
│   │   └── single_end/          # Alineamientos single-end
│   └── Reporte.md               # Este reporte
└── src/                         # Scripts de análisis
    ├── fastp_limpieza.sh        # Script de limpieza de datos
    ├── hisat2_alineamientos.sh  # Script de alineamiento con HISAT2
    └── star_alineamientos.sh    # Script de alineamiento con STAR
```

## 1. Descarga de Datos

Los datos de secuencias se encuentran en la carpeta `data/SRRs/`, en formato FASTQ descargados de repositorios públicos, con el uso de `fasterq-dump` y `prefetch`


## 2. Control de Calidad Inicial

Usando `fastqc` se procesaron los archivos FASTQ descargados, dando como resultado reportes sobre la calidad de los mismos.

### 2.1 Ejecución de FastQC

Se ejecutó FastQC sobre todas las secuencias raw para evaluar su calidad:

```bash 
fastqc -o fastqc SRRs/*.fastq 
```

Para consolidar los resultados, se ejecutó MultiQC en el directorio fastqc:

```bash 
multiqc . -n multiqc_report.html  
```

### 2.2 Análisis de Resultados

Los reportes HTML de FastQC y MultiQC revelan que los datos presentan buena calidad general. A continuación se detallan los hallazgos principales:

**Calidad de Secuencias:**
- No contienen bases `N` ambiguas
- Todas las secuencias tienen una longitud estándar de 100 pb (Illumina)
- Buena calidad promedio de bases por posición

![multiqc per_sequence_quality_scores_plot](/data/fastqc/multiqc_report_data/images/fastqc_per_base_sequence_quality_plot.png)
Fig 1. Per Sequence Quality Scores

**Contenido de Bases:**
Se observó variación en el contenido de bases en las primeras 15 posiciones de las secuencias, después de las cuales se estabiliza (Fig. 2). Este patrón es típico en librerías de Illumina y puede resolverse eliminando las primeras 15-18 pb durante la limpieza.

![Sequence_content.png](/data/fastqc/multiqc_report_data/images/fastqc_per_base_sequence_content_plot.png)
Fig. 2. Per base sequence content

**Duplicación de Secuencias:**
Se detectan niveles moderados a altos de duplicación de secuencias (Fig. 3). Sin embargo, como estos son datos de RNA-seq, las secuencias duplicadas se conservan, ya que pueden contener información biológica importante para el análisis diferencial.

![fastqc_sequence_duplication_levels_plot.png](/data/fastqc/multiqc_report_data/images/fastqc_sequence_duplication_levels_plot.png)
Fig. 3. Sequence Duplication Levels

**Secuencias Sobrerepresentadas:**
En el análisis de secuencias sobrerepresentadas, la mayoría mostró "No hit", indicando que son señales biológicas y no artefactos de secuenciación.

**Contenido de Adaptadores:**
La muestra `SRR9126963` presenta contaminación por adaptadores Nextera Transposase. Este adaptador debe ser removido durante la limpieza de los datos.

![fastqc_adapter_content_plot.png](/data/fastqc/multiqc_report_data/images/fastqc_adapter_content_plot.png)
Fig. 4. Adapter Content

## 3. Limpieza de Datos con fastp

A pesar de que solo una muestra presenta contaminación por adaptadores, todas las muestras se procesaron con idénticos parámetros de limpieza para garantizar homogeneidad en el análisis posterior. Esto asegura que todos los datos transcriptómicos estén en igualdad de condiciones.

### 3.1 Script de Limpieza 

El script de limpieza usa el programa `fastp`, el código del mismo puede ser encontrado en el Apéndice, Sección-I.

**Tiempo de ejecución:** 10 minutos 30 segundos para todas las muestras.

**Parámetros de fastp utilizados:**
- `-f 15 -F 15`: Recorte de 15 pb desde el inicio de ambos pares de lectura para evitar problemas de calidad
- `--detect_adapter_for_pe`: Detección automática y eliminación de adaptadores paired-end
- `-l 50`: Longitud mínima de 50 pb después del recorte (elimina secuencias muy cortas)
- `--trim_poly_g`: Recorte de colas de poliguanina (útil para NovaSeq)
- Archivos de salida comprimidos en formato .fastq.gz

**Comparación de los estatus**

Antes de la limpieza:
![Estatus fastqc crudos](/data/fastqc/multiqc_report_data/images/fastqc-status-check-heatmap_raw.png)

Después de la limpieza:
![Estuatus fastqc procesados](/data/fastqc/multiqc_report_data/images/fastqc-status-check-heatmap_procesado.png)

* Se puede observar el estatus a color verde en adaptadores y contenido de base por secuencia, la longitud paso a estatus color amarillo pero era esperable debido al recorte.

## 4. Alineamiento de Secuencias

Los alineamientos se realizaron utilizando dos herramientas diferentes: HISAT2 y STAR. Ambas se ejecutaron en modo paired-end y single-end para comparar resultados. Los resultados se guardan en la carpeta `results/` con estructura separada por alineador y tipo de alineamiento.

### 4.1 Alineamiento con HISAT2 

Script disponible en el Apéndice, Sección-II.

**Parámetros de HISAT2 utilizados:**
- `-p 8`: Uso de 8 hilos para procesamiento paralelo eficiente
- `--no-unal`: Excluye lecturas no alineadas del archivo de salida para reducir su tamaño
- `--no-mixed`: Solo alinea pares completos (paired-end)
- `--no-discordant`: Solo alinea pares con orientación correcta y distancia razonable
- Archivos de salida en formato SAM

El tiempo de ejecución de cada alineamiento:
| | SRR9126694 | SRR9126930 | SRR9126963 | SRR9127454 | SRR9127455 | SRR9127457 | SRR9127584 | 
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| single_end | 4m8.636s | 3m25.969s | 2m51.976s | 3m30.032 | 2m38.712s | 3m22.183s | 1m47.139s
| paired_end | 8m49.477s | 7m51.402s | 6m53.150s | 6m11.756s | 5m31.339s | 7m7.690s | 3m54.522s

Total single end: 21 minutos con 44.647 segundos
Total paired end: 46 minutos con 19.336 segundos
Total global: 1 hora, 8 minutos y 3.983 segundos

### 4.2 Alineamiento con STAR

Script disponible en el Apéndice, Sección-III.

**Parámetros de STAR utilizados:**
- `--runThreadN 8`: Uso de 8 hilos para procesamiento paralelo eficiente
- `--readFilesCommand zcat`: Especifica descompresión de archivos gzip en tiempo real
- `--outSAMunmapped None`: Excluye lecturas no alineadas del archivo de salida
- `--outSAMtype SAM`: Formato de salida SAM (compatible con HISAT2 para comparación)
- Archivos de salida incluyen: alineamientos SAM, logs detallados y splice junctions

El tiempo de ejecución de cada alineamiento:
| | SRR9126694 | SRR9126930 | SRR9126963 | SRR9127454 | SRR9127455 | SRR9127457 | SRR9127584 | 
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| single_end | 3m9.891s | 2m14.195s | 2m3.863s | 1m53.746s | 1m35.023s | 1m44.683s | 1m11.076s
| paired_end | 5m15.661s | 4m14.119s | 3m50.596s | 3m12.379s | 3m0.588s | 3m16.733s | 2m33.589s

Total single end: 13 minutos con 52.477 segundos
Total paired end: 25 minutos con 23.665 segundos
Total global: 39 minutos con 16.142 segundos

## 5. Resultados

Los archivos de alineamiento se encuentran organizados en la carpeta `results/` con la siguiente estructura:

- **HISAT2**: Carpetas `paired_end/` y `single_end/` contienen archivos SAM por muestra
- **STAR**: Carpetas `paired_end/` y `single_end/` contienen archivos SAM, logs de ejecución e información de splice junctions

Ambos alineadores proporcionan archivos de resumen que incluyen estadísticas de mapeo, permitiendo comparar la calidad de los alineamientos entre métodos.

### 5.1 Análisis Detallado por Alineador

#### STAR

Single-ended

| SRRs | Number of input reads | % uniquely mapped reads | % of reads mapped to multiple loci | % of reads mapped to too many loci | % of reads unmapped: too many mismatches | % of reads unmapped: too short | % of reads unmapped: other | % of chimeric reads |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SRR9126694 | 26796394 | 73.60% | 14.77% | 0.83% | 0.00% | 10.47% | 0.34% | 0.00% |
| SRR9126930 | 25037181 | 75.17% | 15.58% | 0.31% | 0.00% | 8.88% | 0.06% | 0.00% |
| SRR9126963 | 19724778 | 73.67% | 17.20% | 0.58% | 0.00% | 8.36% | 0.20% | 0.00% |
| SRR9127454 | 18834843 | 80.67% | 13.34% | 0.81% | 0.00% | 4.89% | 0.29% | 0.00% |
| SRR9127455 | 17768667 | 74.63% | 16.04% | 0.40% | 0.00% | 8.83% | 0.10% | 0.00% |
| SRR9127457 | 23080779 | 76.63% | 15.46% | 0.34% | 0.00% | 7.52% | 0.05% | 0.00% |
| SRR9127584 | 10827984 | 78.27% | 14.32% | 1.77% | 0.00% | 5.17% | 0.47% | 0.00% |

Paired-ended

| SRRs | Number of input reads | % uniquely mapped reads | % of reads mapped to multiple loci | % of reads mapped to too many loci | % of reads unmapped: too many mismatches | % of reads unmapped: too short | % of reads unmapped: other | % of chimeric reads |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| SRR9126694 | 26796394 | 67.40% | 5.89% | 0.61% | 0.00% | 25.75% | 0.35% | 0.00% |
| SRR9126930 | 25037181 | 69.33% | 6.42% | 0.23% | 0.00% | 24.00% | 0.02% | 0.00% |
| SRR9126963 | 19724778 | 65.95% | 7.93% | 0.40% | 0.00% | 25.53% | 0.19% | 0.00% |
| SRR9127454 | 18834843 | 77.82% | 5.72% | 0.61% | 0.00% | 15.60% | 0.25% | 0.00% |
| SRR9127455 | 17768667 | 70.45% | 6.55% | 0.27% | 0.00% | 22.66% | 0.07% | 0.00% |
| SRR9127457 | 23080779 | 71.88% | 6.52% | 0.25% | 0.00% | 21.33% | 0.02% | 0.00% |
| SRR9127584 | 10827984 | 75.03% | 7.78% | 1.38% | 0.00% | 15.23% | 0.58% | 0.00% |

#### HISAT2

Single-ended

| SRRs | Number of reads | Aligned 0 times | Aligned exactly 1 time | Aligned >1 times | Overall alignment rate |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SRR9126694 | 26796394 | 20.90% | 66.34% | 12.76% | 79.10% |
| SRR9126930 | 25037181 | 18.09% | 68.26% | 13.64% | 81.91% |
| SRR9126963 | 19724778 | 16.73% | 68.18% | 15.09% | 83.27% |
| SRR9127454 | 18834843 | 11.87% | 75.87% | 12.26% | 88.13% |
| SRR9127455 | 17768667 | 18.16% | 68.15% | 13.69% | 81.84% |
| SRR9127457 | 23080779 | 15.91% | 70.46% | 13.63% | 84.09% |
| SRR9127584 | 10827984 | 12.20% | 74.19% | 13.61% | 87.80% |

Paired-ended

| SRRs | Number of reads | Aligned 0 times | Aligned exactly 1 time | Aligned >1 times | Overall alignment rate |
| :--- | :--- | :--- | :--- | :--- | :--- |
| SRR9126694 | 26796394 | 36.08% | 58.61% | 5.31% | 63.92% |
| SRR9126930 | 25037181 | 33.48% | 60.57% | 5.95% | 66.52% |
| SRR9126963 | 19724778 | 35.23% | 57.80% | 6.97% | 64.77% |
| SRR9127454 | 18834843 | 24.75% | 69.72% | 5.53% | 75.25% |
| SRR9127455 | 17768667 | 32.62% | 61.71% | 5.67% | 67.38% |
| SRR9127457 | 23080779 | 30.45% | 63.41% | 6.13% | 69.55% |
| SRR9127584 | 10827984 | 25.08% | 67.94% | 6.98% | 74.92% |

### Gráfica

Para una mejor observación general de los datos se generó una gráfica.
El código fuente esta disponible en el Apéndice, Sección-IV.

![Gráfica estadísticas](/results/tasa_mapeo_barras_puntos.png)

## 6. Conclusión

El análisis de los alineamientos muestra un comportamiento consistente de las herramientas HISAT2 y STAR al procesar los datos en modos single-end y paired-end.

**Single-end vs. Paired-end:** El modo *single-end* presentó sistemáticamente tasas de mapeo más altas (arriba del 75%) en comparación con el *paired-end* (~65-75%). Esta disminución en el abordaje pareado es esperada debido a las mayores restricciones impuestas necesarias para alinear el fragmento completo (orientación, concordancia y longitud), lo cual descarta lecturas que individualmente podrían alinear pero no satisfacen los criterios genómicos completos.

**Desempeño de Herramientas y Muestras:** Ambos alineadores demostraron alta capacidad de procesamiento; STAR resultó especialmente sensible al descartar fragmentos por ser demasiado cortos en el modo paired-end, mientras que HISAT2 presentó un aumento notable de lecturas no alineadas por incompatibilidad del par. Se destaca de manera pronunciada la muestra `SRR9127454`, la cual logró consistentemente las mejores métricas de mapeo único y menor proporción de datos no alineados en la totalidad de las estrategias bajo prueba.

Los alineamientos obtenidos, particularmente aquellos con alto porcentaje de mapeo único y de forma específica usando el protocolo *paired-end* para minimizar falsos positivos, garantizan datos confiables aptos para utilizarlos con seguridad en subsecuentes estudios de expresión génica diferencial.

**Comparación de tiempos de ejecución:** Relizando la comparación, el modo *paired-end* fue más lento que *single-end* (HISAT2: 46m19.336s vs 21m44.647s; STAR: 25m23.665s vs 13m52.477s), reflejando el costo computacional adicional de imponer concordancia entre pares. Al comparar alineadores, STAR fue consistentemente más rápido que HISAT2 tanto en *single-end* (13m52.477s vs 21m44.647s, 7m52.170s menos) como en *paired-end* (25m23.665s vs 46m19.336s, 20m55.671s menos), con una reducción total de 28m47.841s (39m16.142s vs 1h8m3.983s) al procesar el conjunto completo.

# Apéndice

## Sección-I
```bash
#!/bin/bash

# Extraemos el nombre de la base de la muestra

for R1 in /export/storage/users/andreavg/transcriptomica/data/SRRs/*_1.fastq; do

 base=$(basename "$R1" _1.fastq)

# Definimos el archivo Read2 correspondiente

 R2="/export/storage/users/andreavg/transcriptomica/data/SRRs/${base}_2.fastq"

# Definimos los archivos de salida

out_R1="/export/storage/users/andreavg/transcriptomica/clean/${base}_1_clean.fastq.gz"

out_R2="/export/storage/users/andreavg/transcriptomica/clean/${base}_2_clean.fastq.gz"

echo "Procesando $base"

# Ejecutamos fastp

echo "Ejecutando fastp para $base"
# Flags:
# -i : archivo de entrada para las lecturas del primer par, en este caso el archivo R1 que corresponde a la muestra actual
# -I : archivo de entrada para las lecturas del segundo par, en este caso el archivo R2 que corresponde a la muestra actual
# -o : archivo de salida para las lecturas del primer par, en este caso el archivo out_R1 que corresponde a la muestra actual
# -O : archivo de salida para las lecturas del segundo par, en este caso el archivo out_R2 que corresponde a la muestra actual
# -f : número de bases a recortar desde el inicio de las lecturas, en este caso 15 para ser laxos y evitar problemas de calidad al inicio de las lecturas
# -F : número de bases a recortar desde el inicio de las lecturas del segundo par, en este caso 15 para ser laxos y evitar problemas de calidad al inicio de las lecturas
# --detect_adapter_for_pe : permite a fastp detectar y eliminar adaptadores. En este caso para eliminar a Nextera Transposase Sequence,
# el cúal observamos en la calidad de las lecturas.
# -l : longitud mínima de las lecturas después del recorte, en este caso 50 para eliminar secuencias muy cortas que podrían introducir ruido en el alineamiento
# --trim_poly_g : permite a fastp recortar las colas de poliguanina, lo cual es útil para evitar errores relacionados con la tecnología NovaSeq,
# -h : archivo de salida para el reporte HTML generado por fastp.
 fastp -i "$R1" -I "$R2" \
    -o "$out_R1" -O "$out_R2" \
	-f 15 -F 15 \
	--detect_adapter_for_pe \
	-l 50 \
	--trim_poly_g \
	-h "/export/storage/users/andreavg/transcriptomica/clean/${base}_fastp_report.html"


done

echo "Limpieza completada para todas las muestras"
```

## Sección-II
```bash
#!/bin/bash

# Asignamos las variables para nuestro indice de referencia

index="/export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.hisat/mm39.gencode.M36.hisat"

time_file_SE="/export/storage/users/andreavg/transcriptomica/results/hisat2/single_end/Hisat2_time_SE.txt"

time_file_PE="/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end/Hisat2_time_PE.txt"

# Iniciamos nuestro ciclo

for R1 in /export/storage/users/andreavg/transcriptomica/data/clean/*_1_clean.fastq.gz; do

# Asignamos el nombre de la muestra a partir del nombre del archivo
	base=$(basename $R1 _1_clean.fastq.gz)

# Definimos donde se encuentran los archivos R2
	R2="/export/storage/users/andreavg/transcriptomica/data/clean/${base}_2_clean.fastq.gz"

	echo "Procesando muestra $base"
# Ejecutamos el comando de alineamiento con hisat2 con single-end
	echo "Single-end para $base"
# Flags:
# -p : número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
# -x : especifica el índice de referencia, guardado en la variable "index"
# -U : es el archivo de entrada para single-end, en este caso el archivo R1 que corresponde a la muestra actual
# --no-unal : evita que se escriban las lecturas no alineadas en el archivo de salida, para reducir el tamaño del archivo resultante
# --summary-file : especifica el archivo donde se guardará el resumen del alineamiento
# -S : especifica el archivo de salida para el alineamiento en formato SAM
# time se utiliza para medir el tiempo de ejecución del comando hisat2, y se redirige la salida de error estándar (que es donde se imprime el tiempo) al archivo correspondiente para single-end.
# usando las {} se asegura que el tiempo se mida solo para el comando hisat2 y no para otros comandos que puedan estar en el script.
	{ time hisat2 -p 8 \
		-x "$index" \
		-U "$R1" \
		--no-unal \
		--summary-file "/export/storage/users/andreavg/transcriptomica/results/hisat2/single_end/summary/${base}.summary.txt" \
		-S "/export/storage/users/andreavg/transcriptomica/results/hisat2/single_end/${base}_SE.sam" ; } 2>> "$time_file_SE"


# Ejecutamos hisat2 con paired-end
	echo "Paired-end para $base"
#Flags:
# -p : número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
# -x : especifica el índice de referencia, guardado en la variable "index"
# -1 : es el archivo de entrada para las lecturas del primer par, en este caso el archivo R1 que corresponde a la muestra actual
# -2 : es el archivo de entrada para las lecturas del segundo par, en este caso el archivo R2 que corresponde a la muestra actual
# --no-unal : evita que se escriban las lecturas no alineadas en el archivo de salida, para reducir el tamaño del archivo resultante
# --no-mixed : asegura que solo se alineen los pares completos, esto se hace con el fin de obtener resultados más confiables, evitando alineamientos parciales que podrían ser menos precisos
# --no-discordant : asegura que solo se alineen los pares que están correctamente orientados y a una distancia razonable, esto ayuda a mejorar la calidad del alineamiento
# --summary-file : especifica el archivo donde se guardará el resumen del alineamiento
# -S : especifica el archivo de salida para el alineamiento en formato SAM
# time se utiliza para medir el tiempo de ejecución del comando hisat2, y se redirige la salida de error estándar (que es donde se imprime el tiempo) al archivo correspondiente para single-end.
# usando las {} se asegura que el tiempo se mida solo para el comando hisat2 y no para otros comandos que puedan estar en el script.
	{ time hisat2 -p 8 \
		-x "$index" \
		-1 "$R1" \
		-2 "$R2" \
		--no-unal \
		--no-mixed \
		--no-discordant \
		--summary-file "/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end/summary/${base}.summary.txt" \
		-S "/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end/${base}_PE.sam" ; } 2>> "$time_file_PE"

  

done

echo "Alineamiento completado para todas las muestras"
```

## Sección-III
```bash
#!/bin/bash
# Asignamos las variables para nuestro indice de referencia

index="/export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.star/"
Star_time_file_SE="/export/storage/users/andreavg/transcriptomica/results/star/single_end/STAR_time_SE.txt"
Star_time_file_PE="/export/storage/users/andreavg/transcriptomica/results/star/paired_end/STAR_time_PE.txt"

# Iniciamos nuestro ciclo
for R1 in /export/storage/users/andreavg/transcriptomica/data/clean/*_1_clean.fastq.gz; do

# Asignamos el nombre de la muestra a partir del nombre del archivo

	base=$(basename $R1 _1_clean.fastq.gz)

# Definimos donde se encuentran los archivos R2

	R2="/export/storage/users/andreavg/transcriptomica/data/clean/${base}_2_clean.fastq.gz"

	echo "Procesando muestra: $base"
# Ejecutamos el comando de alineamiento con hisat2 con single-end
	echo "Single-end para $base"
#Flags:
# --runThreadN : número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
# --genomeDir : especifica el directorio del índice de referencia, guardado en la variable "index"
# --readFilesIn : es el archivo de entrada para single-end, en este caso el archivo R1 que corresponde a la muestra actual
# --readFilesCommand : especifica el comando para descomprimir los archivos de entrada, en este caso usamos "zcat" para manejar archivos comprimidos
# --outSAMunmapped None : evita que se escriban las lecturas no alineadas en el archivo de salida, para reducir el tamaño del archivo resultante
# --outSAMtype SAM : especifica el formato de salida para el alineamiento, en este caso SAM para mantener la consistencia con los resultados de hisat2
# --outFileNamePrefix : especifica el prefijo para los archivos de salida
# time se utiliza para medir el tiempo de ejecución del comando STAR y se redirige la salida de error estándar, que es donde se imprime el tiempo, al archivo correspondiente para single-end.
# usando las {} se asegura que el tiempo se mida solo para el comando STAR y no para otros comandos que puedan estar en el script.
	{ time STAR --runThreadN 8 \
		--genomeDir "$index" \
		--readFilesIn "$R1" \
		--readFilesCommand zcat \
		--outSAMunmapped None \
		--outSAMtype SAM \
		--outFileNamePrefix "/export/storage/users/andreavg/transcriptomica/results/star/single_end/${base}_SE_" ; } 2>> "$Star_time_file_SE"

# Ejecutamos hisat2 con paired-end
	echo "Paired-end para $base"
#Flags:
# --runThreadN : número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
# --genomeDir : especifica el directorio del índice de referencia, guardado en la variable "index"
# --readFilesIn : es el archivo de entrada para las lecturas del primer par y el segundo par, en este caso los archivos R1 y R2 que corresponden a la muestra actual
# --readFilesCommand : especifica el comando para descomprimir los archivos de entrada, en este caso usamos "zcat" para manejar archivos comprimidos
# --outSAMunmapped None : evita que se escriban las lecturas no alineadas en el archivo de salida, para reducir el tamaño del archivo resultante
# --outSAMtype SAM : especifica el formato de salida para el alineamiento, en este caso SAM para mantener la consistencia con los resultados de hisat2
# --outFileName Prefix : especifica el prefijo para los archivos de salida
# time se utiliza para medir el tiempo de ejecución del comando STAR y se redirige la salida de error estándar, que es donde se imprime el tiempo, al archivo correspondiente para single-end.
# usando las {} se asegura que el tiempo se mida solo para el comando STAR y no para otros comandos que puedan estar en el script.
	{ time STAR --runThreadN 8 \
		--genomeDir "$index" \
		--readFilesIn "$R1" "$R2" \
		--readFilesCommand zcat \
		--outSAMunmapped None \
		--outSAMtype SAM \
			--outFileNamePrefix "/export/storage/users/andreavg/transcriptomica/results/star/paired_end/${base}_PE_" ; } 2>> "$Star_time_file_PE"

done

echo "Alineamiento completado para todas las muestras"
```

## Sección-IV
```R
#!/usr/bin/env Rscript

# Suprimir el stdin sobre el cargado de datos
suppressPackageStartupMessages({
	library(ggplot2)
	library(dplyr)
})

# Nombre de los SRR 
srr <- c(
	"SRR9126694", "SRR9126930", "SRR9126963", "SRR9127454",
	"SRR9127455", "SRR9127457", "SRR9127584"
)

# STAR: tasa total mapeada = uniquely + multiple + too many loci
# Tomando todos los datos de STAR
# Formando un data.frame
star_single <- data.frame(
	SRR = srr,
	uniquely = c(73.60, 75.17, 73.67, 80.67, 74.63, 76.63, 78.27),
	multiple = c(14.77, 15.58, 17.20, 13.34, 16.04, 15.46, 14.32),
	too_many = c(0.83, 0.31, 0.58, 0.81, 0.40, 0.34, 1.77)
) |>
	mutate(
		Alineador = "STAR",
		Modo = "Single-end",
		Mapeo = uniquely + multiple + too_many
	) |>
	select(SRR, Alineador, Modo, Mapeo)

star_paired <- data.frame(
	SRR = srr,
	uniquely = c(67.40, 69.33, 65.95, 77.82, 70.45, 71.88, 75.03),
	multiple = c(5.89, 6.42, 7.93, 5.72, 6.55, 6.52, 7.78),
	too_many = c(0.61, 0.23, 0.40, 0.61, 0.27, 0.25, 1.38)
) |>
	mutate(
		Alineador = "STAR",
		Modo = "Paired-end",
		Mapeo = uniquely + multiple + too_many
	) |>
	select(SRR, Alineador, Modo, Mapeo)

# HISAT2: usar la tasa global de alineamiento reportada
hisat_single <- data.frame(
	SRR = srr,
	Alineador = "HISAT2",
	Modo = "Single-end",
	Mapeo = c(79.10, 81.91, 83.27, 88.13, 81.84, 84.09, 87.80)
)

hisat_paired <- data.frame(
	SRR = srr,
	Alineador = "HISAT2",
	Modo = "Paired-end",
	Mapeo = c(63.92, 66.52, 64.77, 75.25, 67.38, 69.55, 74.92)
)

# Organizar todos los datos
datos <- bind_rows(star_single, star_paired, hisat_single, hisat_paired) |>
	mutate(
		Grupo = factor(
			paste(Alineador, Modo, sep = " - "),
			levels = c(
				"STAR - Single-end", "STAR - Paired-end",
				"HISAT2 - Single-end", "HISAT2 - Paired-end"
			)
		)
	)

# Obtener el promedio de los datos
promedios <- datos |>
	group_by(Grupo) |>
	summarise(Media = mean(Mapeo), .groups = "drop")

# Formar el gráfico
grafica <- ggplot() +
	geom_col(
		data = promedios,
		aes(x = Grupo, y = Media, fill = Grupo),
		width = 0.70,
		alpha = 0.80,
		show.legend = FALSE
	) +
	geom_point(
		data = datos,
		aes(x = Grupo, y = Mapeo),
		position = position_jitter(width = 0.12, height = 0),
		size = 2.6,
		alpha = 0.95,
		color = "black"
	) +
	stat_summary(
		data = datos,
		aes(x = Grupo, y = Mapeo),
		fun.data = mean_se,
		geom = "errorbar",
		width = 0.15,
		linewidth = 0.5,
		color = "black"
	) +
	scale_y_continuous(
		limits = c(0, 100),
		expand = expansion(mult = c(0, 0.03))
	) +
	labs(
		title = "Tasa de mapeo por alineador y modo",
		subtitle = "Barras: promedio por grupo | Puntos: muestras SRR individuales",
		x = "Grupo",
		y = "Lecturas mapeadas (%)"
	) +
	theme_minimal(base_size = 12) +
	theme(
		axis.text.x = element_text(angle = 20, hjust = 1),
		panel.grid.major.x = element_blank(),
		plot.title = element_text(face = "bold")
	)

print(grafica)

ggsave(
	filename = "results/tasa_mapeo_barras_puntos.png",
	plot = grafica,
	width = 10,
	height = 6,
	dpi = 300,
    bg = "white"
)
```


