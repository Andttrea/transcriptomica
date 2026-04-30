# Reporte de pseudoalineamientos 

## Planteamiento 

El análisis de la expresión génica diferencial depende de una cuantificación precisa de las lecturas secuenciadas. Tras obtener los alineamientos iniciales con las herramientas de STAR e HISAT2, es necesario transformar esos datos en matrices de conteo utilizables. 

Además de los métodos de alineamiento tradicionales, se explora el uso de pseudoalineadores como Salmon, los cuales ofrecen una alternativa rápida y eficiente al asignar lecturas directamente al transcriptoma. En este trabajo, se analizan muestras de tejido muscular de Mus musculus provenientes del dataset Tabula Muris Senis, comparando ratones jóvenes (9 meses) contra ratones viejos (24 meses).

Todos los datos, scripts y resultados se encuentran ubicados en el servidor chaac de la licenciatura, en el directorio `/export/storage/users/andreavg/transcriptomica`.

## Objetivos 

### Objetivo general

Generar y comparar matrices de cuantificación de expresión génica mediante alineamiento tradicional (STAR/HISAT2) y pseudoalineamiento (Salmon) utilizando datos de RNA-seq de músculo de ratón a distintas edades.

### Objetivos particulares

1. Realizar el procesamiento post-alineamiento (conversión, filtrado y ordenamiento) de archivos SAM a BAM.
2. Generar archivos de visualización BigWig (pileups) para el análisis de cobertura.
3. Cuantificar la abundancia de genes mediante featureCounts para los alineadores STAR e HISAT2 en modo Single-End y Paired-End.
4. Ejecutar el pseudoalineamiento con Salmon y convertir las abundancias de transcritos a nivel de gen utilizando la herramienta tximport.


## Estructura del repostorio 

```text
transcriptomica/
├── README.md                      
├── data/                          # Directorio de datos de entrada crudos y pre-procesados
│   ├── bed_file/                  # Archivos BED con regiones de exclusión (mm39)
│   ├── clean/                     
│   ├── fastqc/                   
│   ├── gencode/                   # Archivos de referencia 
│   ├── metadata/                  
│   └── SRRs/                      
├── results/                       # Resultados generados en los análisis
│   ├── figuras/                   
│   ├── hisat2/                    # Archivos de resultados de HISAT2 (alineamientos BAM, pileups, conteos)
│   ├── reportes/                  # Ubicación de los reportes markdown
│   │   ├── Reporte_Hisat2_Star.md 
│   │   └── Reporte_Pseudoalineamientos.md # Reporte actual (cuantificación y Salmon)
│   ├── salmon/                    # Archivos de cuantificación generados por Salmon
│   └── star/                      # Archivos de resultados de STAR (alineamientos BAM, conteos)
└── src/                           # Códigos fuente y scripts (Bash y R) 
    ├── abundancia_tx2gene.R       # Script en R para convertir abundancias de transcritos a genes (tximport)
    ├── fastp_limpieza.sh          
    ├── feature_counts_hisat2.sh   # Conteo de lecturas sobre genes mapeadas con HISAT2
    ├── feature_counts_star.sh     # Conteo de lecturas sobre genes mapeadas con STAR
    ├── generar_pileups.sh         # Script para generar archivos BigWig de cobertura temporal
    ├── graficas.R                 
    ├── hisat2_alineamientos.sh    
    ├── rds_to_tsv.R               # Script para convertir salidas de tximport (.rds) a matrices TSV
    ├── salmon_alineamiento.sh     # Indexación de transcriptoma y cuantificación (pseudoalineamiento) con Salmon
    ├── sam_to_bam.sh              # Script reportado arriba: conversión de SAM a BAM filtrado y ordenado
    └── star_alineamientos.sh      
```

## Metodología y Justificación

### Conversión de SAM a BAM 

Como se discutió en clase, los archivos SAM son pesados y difíciles de procesar. Se procedió a su conversión a formato binario (BAM), filtrando por calidad de mapeo ($MAPQ < 10$) para asegurar que solo los alineamientos confiables se utilicen en la cuantificación. Los archivos fueron ordenados por coordenadas para permitir la generación de índices.

```bash
 # Nos movemos a la carpeta donde crearemos el script
 cd /export/storage/users/andreavg/transcriptomica/src
 # Creamos el script para hacer la conversión de SAM a BAM
 nano sam_to_bam.sh 
 ## ============================================== SCRIPT ============================================== ##
 #!/bin/bash
 
 # Asignamos la variable para el directorio de Hisat2
 SAM_Dir=("/export/storage/users/andreavg/transcriptomica/results/hisat2/single_end"
 "/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end"
 "/export/storage/users/andreavg/transcriptomica/results/star/single_end"
 "/export/storage/users/andreavg/transcriptomica/results/star/paired_end")
 
 # Recorremos cada directorio de la variable SAM_Dir y convertiremos los archivos SAM a BAM
 for DIR in "${SAM_Dir[@]}"; do
	 
     echo "Procesando directorio: $DIR"
	 
     # Recorremos cada archivo SAM en el directorio actual
	 for SAM in "$DIR"/*.sam ; do
    # Declaramos los archivos de salida BAM y el nombre base para cada archivo SAM
        base="${SAM%.sam}"
		BAM="${base}.bam"
	 
		echo "Convirtiendo $(basename $SAM)"
	 
	    #Conversión de SAM a BAM y utilizamos samtools sort para ordenar el BAM resultante. La ordenación es 
        #importante para la indexación posterior y para mejorar el rendimiento en análisis posteriores.
        # Flags:
        # -bSq 10: Convertir a BAM y filtrar por calidad de mapeo (MAPQ >= 10)
        # --threads 5: Usar 5 hilos para la conversión
        # -T : Usar un archivo temporal para la ordenación
        # -o : Especificar el archivo de salida BAM
        samtools view -bSq 10 "$SAM" | samtools sort -@ 5 -T "${base}_temp" -o "$BAM"
	 
	 # Indexamos 
		samtools index "$BAM"
	 
		echo "Se convirtió con éxito $(basename $BAM)"
	 
	 done 
 done 
 
 echo "Conversión completada para todos los directorios"
 ## ============================================== SCRIPT ============================================== ##

 # Le damos permisos al script
 chomd +x sam_to_bam.sh
 # Ejecutamos el script
 ./sam_to_bam.sh

 # Después de obtener los archivos BAM debemos verificar que los archivos no estén vacíos y obtener el número exacto de lecturas que lograron alinearse #   exitosamente al genoma de referencia.
for BAM in /export/storage/users/andreavg/transcriptomica/results/{hisat2,star}/{single_end,paired_end}/*.bam; do
    # Asignamos las variables para asignar nuestro archivo de salida
    out_dir="${BAM%/*}"
    bam_name="${BAM##*/}"
    out_file="${out_dir}/n.align_reads_${bam_name}.txt"
    
    # Iniciamos el ciclo
    if [[ "$BAM" =~ single_end ]]; then
        # Si es single no utiliza la flag 64
        samtools view -@ 25 -c "$BAM" > "$out_file"
    elif [[ "$BAM" =~ paired_end ]]; then
        samtools view -@ 25 -f 64 -c "$BAM" > "$out_file"
    fi
done 

# Después de comprobar que los archivos BAM si tuvieran reads,borramos los archivos BAM.
rm *.sam 
```


### Generación de pileups(BigWig)

Para la visualización de los datos, se generaron archivos BigWig (.bw). Los archivos BigWig (bw) representan el "pileup" o la acumulación de lecturas a lo largo del genoma.

Se utilizó la herramienta bamCoverage de deeptools, aplicando la flag 64 para contar únicamente la primera lectura del par en datos Paired-End, evitando así la duplicación artificial de la señal. Además se aplicó normalización BPM (Bins Per Million) que ajusta la señal según el tamaño de la librería, permitiendo comparar visualmente las muestras. 

bamCoverage necesita un archivo BED  con el fin de evitar coordenadas genómicas conflictivas en análisis genómicos. En este caso utlizamos `mm39.excluderanges`, el cual se encontró en [excluderanges](https://dozmorovlab.github.io/excluderanges/#excluderanges-genomic-ranges-of-problematic-genomic-regions). 
El archivo se encuentra en el siguiente path: `/export/storage/users/andreavg/transcriptomica/data/bed_file`

Solo se corrió bamCoverage para paired-end para optimizar el tiempo de procesamiento en el servidor.

```bash
# Primero creamos la carpeta donde guardaremos los archivos generados por bamCoverage
# Para hisat2
cd /export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end
mkdir pileups
# Para star
cd ../../star/paired_end
mkdir pileups

# Nos movemos a la carpeta src y creamos el script
cd /export/storage/users/andreavg/transcriptomica/src
nano generar_pileups.sh

## ============================================== SCRIPT ============================================== ##

#!/bin/bash

# Asignamos las variales para las rutas de los archivos paired-end, tanto para HIsat2 y para STAR
DIR_HISAT2="/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end"
DIR_STAR="/export/storage/users/andreavg/transcriptomica/results/star/paired_end"

# Asignamos la ruta para el archivo de blacklist, el cual contiene regiones del genoma que deben ser excluidas del análisis 
# debido a que pueden generar señales de fondo no deseadas. 
FILE_BLACKLIST="/export/storage/users/andreavg/transcriptomica/data/bed_file/mm39.excluderanges.bed"


# Iniciamos nuestro ciclo para los archivos BAM de HISAT2
echo "Generando pileups para HISAT2..."

for BAM in "$DIR_HISAT2"/*bam; do 

    # Asignamos el nombre de la muestra a partir del nombre del archivo
    base=$(basename "$BAM" .bam)
    
    # Ejecutamos el comando bamCoverage para generar archivos bigwig a partir de los archivos BAM de HISAT2
    echo "Procesando muestra: $base"
    # Flags:
    # -p 8 : número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
    # -b : especifica el archivo BAM de entrada, que corresponde al archivo BAM de la muestra actual
    # -o : especifica la ruta de salida para el archivo bigwig generado
    # -bs 20 : establece el tamaño de las bins a 20 bases, lo que significa que el conteo de lecturas se realizará en intervalos de 20 bases a lo largo del genoma
    # --blackListFileName : especifica el archivo de blacklist que contiene regiones del genoma a excluir del análisis
    # --normalizeUsing BPM : normaliza los conteos de lecturas utilizando BPM (Reads Per Million), lo que permite comparar la cobertura entre diferentes muestras de manera más justa
    # --skipNAs : omite las posiciones del genoma que no tienen lecturas alineadas, ayudando a reducir el tamaño del archivo de salida y enfocarse en las regiones con datos relevantes
    # --ignoreDuplicates : ignora las lecturas duplicadas, lo que es importante para evitar sesgos en la cobertura debido a PCR
    # --samFlagInclude 64 : incluye solo las lecturas que están alineadas al primer extremo (read 1) de los pares
    bamCoverage -p 8 \
        -b "$BAM" \
        -o "/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end/pileups/${base}.bw" \
        -bs 20 --blackListFileName "$FILE_BLACKLIST" \
        --normalizeUsing BPM --skipNAs --ignoreDuplicates \
        --samFlagInclude 64
done 

# Procederemos a iniciar nuestro ciclo para los archivos BAM de STAR
echo "Generando pileups para STAR..."
for BAM in "$DIR_STAR"/*bam; do

    # Asignamos el nombre de la muestra a partir del nombre del archivo
    base=$(basename "$BAM" .bam)
    
    # Ejecutamos el comando bamCoverage para generar archivos bigwig a partir de los archivos BAM de STAR
    echo "Procesando muestra: $base"
    # Flags:
    # -p 8 : número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
    # -b : especifica el archivo BAM de entrada, que corresponde al archivo BAM de la muestra actual
    # -o : especifica la ruta de salida para el archivo bigwig generado
    # -bs 20 : establece el tamaño de las bins a 20 bases, lo que significa que el conteo de lecturas se realizará en intervalos de 20 bases a lo largo del genoma
    # --blackListFileName : especifica el archivo de blacklist que contiene regiones del genoma a excluir del análisis
    # --normalizeUsing BPM : normaliza los conteos de lecturas utilizando BPM (Reads Per Million), lo que permite comparar la cobertura entre diferentes muestras de manera más justa
    # --skipNAs : omite las posiciones del genoma que no tienen lecturas alineadas, ayudando a reducir el tamaño del archivo de salida y enfocarse en las regiones con datos relevantes
    # --ignoreDuplicates : ignora las lecturas duplicadas, lo que es importante para evitar sesgos en la cobertura debido a PCR
    # --samFlagInclude 64 : incluye solo las lecturas que están alineadas al primer extremo (read 1) de los pares
    bamCoverage -p 8 \
        -b "$BAM" \
        -o "/export/storage/users/andreavg/transcriptomica/results/star/paired_end/pileups/${base}.bw" \
        -bs 20 --blackListFileName "$FILE_BLACKLIST" \
        --normalizeUsing BPM --skipNAs --ignoreDuplicates \
        --samFlagInclude 64
done 


echo "Proceso de generación de pileups completado"

## ============================================== SCRIPT ============================================== ##

# Habilitamos los permisos, activamos el entorno y corremos el script
chmod +x generar_pileups.sh
conda activate deeptools
./generar_pileups.sh

# Desactivamos el entorno
conda deactivate
```

### Cuantificación por feautureCounts


Se deben convertir las coordenadas de alineamiento en una matriz de conteos digitales. Para esto se utilizó featureCounts, una herramienta altamente eficiente para asignar lecturas a características genómicas (exones/genes) basadas en una anotación externa. 
Se generaron matrices de conteo de lecturas para STAR e HISAT2. Para los datos paired-end, se utilizaron los parámetros -p (indicando fragmentos emparejados) y -B (exigiendo que ambos extremos estén mapeados). Se utilizó el archivo de anotación de Gencode M36, el cual se encontró en Gencode y se encuentra en el siguiente path: `/export/storage/users/andreavg/transcriptomica/data/gencode`

```bash
# Nos movemos a la carpeta de hisat2 single_end y paired_end
cd /export/storage/users/andreavg/transcriptomica/results/hisat2/single_end
# Creamos la carpeta para guardar los archivos generados por featureCounts
mkdir featurecounts
cd ../paired_end
mkdir featurecounts

# Hacemos lo mismo para la carpeta de star
cd ../../star/single_end  
mkdir featurecounts
cd ../paired_end
mkdir featurecounts

# Nos movemos a la carpeta src para crear el script
cd /export/storage/users/andreavg/transcriptomica/src
nano feature_counts_hisat2.sh

## ============================================== SCRIPT ============================================== ##
#!/bin/bash 

# Asignamos las variales para las rutas de los archivos single_end y paired_end, para Hisat2
DIR_HISAT2_SE="/export/storage/users/andreavg/transcriptomica/results/hisat2/single_end"
DIR_HISAT2_PE="/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end"

# Asignamos la variable para la ruta del archivo de anotación GTF
GFF_FILE="/export/storage/users/andreavg/transcriptomica/data/gencode/gencode.vM36.chr_patch_hapl_scaff.annotation.gff3.gz"


# Iniciamos nuestro ciclo para los archivos BAM de HISAT2 single-end
for BAM in "$DIR_HISAT2_SE"/*bam; do
    # Asignamos el nombre de la muestra a partir del nombre del archivo
    base=$(basename "$BAM" .bam)
    
    # Ejecutamos el comando featureCounts para generar archivos de conteo a partir de los archivos BAM de HISAT2 
    echo "Procesando muestra: $base"
    #Flags: 
    # -T: número de hilos a usar, en este caso 8 para optimizar el proceso
    # -a: especifica el archivo de anotación GTF o GFF, que es necesario para asignar las lecturas a las características genómicas correctas
    # --largestOverlap: asigna una lectura a la característica con la que tiene la mayor superposición
    # "$BAM": inidiica el archivo BAM de entrada que se va a procesar
    featureCounts -o "/export/storage/users/andreavg/transcriptomica/results/hisat2/single_end/featurecounts/${base}_fc_counts.txt" \
        -T 8 \
        -a "$GFF_FILE" \
        --largestOverlap \
        "$BAM"
done

# Iniciamos nuestro ciclo para los archivos BAM de HISAT2 paired-end
for BAM in "$DIR_HISAT2_PE"/*bam; do
    # Asignamos el nombre de la muestra a partir del nombre del archivo
    base=$(basename "$BAM" .bam)
    # Ejecutamos el comando featureCounts para generar archivos de conteo a partir de los archivos BAM de HISAT2
    echo "Procesando muestra: $base"
     #Flags: 
    # -T: número de hilos a usar, en este caso 8 para optimizar el proceso
    # -a: especifica el archivo de anotación GTF o GFF, que es necesario para asignar las lecturas a las características genómicas correctas
    # --largestOverlap: asigna una lectura a la característica con la que tiene la mayor superposición
    # -p: indica que los datos son de secuenciación paired-end, lo que permite a featureCounts manejar correctamente las lecturas emparejadas
    # -B: requiere que ambas lecturas de un par estén correctamente alineadas para ser contadas. mejora la precisión 
    # "$BAM": inidiica el archivo BAM de entrada que se va a procesar
    featureCounts -o "/export/storage/users/andreavg/transcriptomica/results/hisat2/paired_end/featurecounts/${base}_fc_counts.txt" \
        -T 8 \
        -a "$GFF_FILE" \
        --largestOverlap \
        -p \
        -B \
        "$BAM"
done
## ============================================== SCRIPT ============================================== ##

# Hacemos un script para star
nano feature_counts_hisat2.sh

## ============================================== SCRIPT ============================================== ##
#!/bin/bash 

# Asignamos las variales para las rutas de los archivos single_end y paired_end, para Hisat2
DIR_STAR_SE="/export/storage/users/andreavg/transcriptomica/results/star/single_end"
DIR_STAR_PE="/export/storage/users/andreavg/transcriptomica/results/star/paired_end"

# Asignamos la variable para la ruta del archivo GFF
GFF_FILE="/export/storage/users/andreavg/transcriptomica/data/gencode/gencode.vM36.chr_patch_hapl_scaff.annotation.gff3.gz"


# Iniciamos nuestro ciclo para los archivos BAM de STAR single-end
for BAM in "$DIR_STAR_SE"/*bam; do
    # Asignamos el nombre de la muestra a partir del nombre del archivo
    base=$(basename "$BAM" .bam)
    
    # Ejecutamos el comando featureCounts para generar archivos de conteo a partir de los archivos BAM de STAR 
    echo "Procesando muestra: $base"
    #Flags: 
    # -T: número de hilos a usar, en este caso 8 para optimizar el proceso
    # -a: especifica el archivo de anotación GTF o GFF, que es necesario para asignar las lecturas a las características genómicas correctas
    # --largestOverlap: asigna una lectura a la característica con la que tiene la mayor superposición
    # "$BAM": inidiica el archivo BAM de entrada que se va a procesar
    featureCounts -o "/export/storage/users/andreavg/transcriptomica/results/star/single_end/featurecounts/${base}_fc_counts.txt" \
        -T 8 \
        -a "$GFF_FILE" \
        --largestOverlap \
        "$BAM"
done

# Iniciamos nuestro ciclo para los archivos BAM de STAR paired-end
for BAM in "$DIR_STAR_PE"/*bam; do
    # Asignamos el nombre de la muestra a partir del nombre del archivo
    base=$(basename "$BAM" .bam)
    # Ejecutamos el comando featureCounts para generar archivos de conteo a partir de los archivos BAM de STAR
    echo "Procesando muestra: $base"
     #Flags: 
    # -T: número de hilos a usar, en este caso 8 para optimizar el proceso
    # -a: especifica el archivo de anotación GTF o GFF, que es necesario para asignar las lecturas a las características genómicas correctas
    # --largestOverlap: asigna una lectura a la característica con la que tiene la mayor superposición
    # -p: indica que los datos son de secuenciación paired-end, lo que permite a featureCounts manejar correctamente las lecturas emparejadas
    # -B: requiere que ambas lecturas de un par estén correctamente alineadas para ser contadas. mejora la precisión 
   # "$BAM": inidiica el archivo BAM de entrada que se va a procesar
    featureCounts -o "/export/storage/users/andreavg/transcriptomica/results/star/paired_end/featurecounts/${base}_fc_counts.txt" \
        -T 8 \
        -a "$GFF_FILE" \
        --largestOverlap \
        -p \
        -B \
        "$BAM"
done
## ============================================== SCRIPT ============================================== ##

# Les damos permiso de ejeción
chmod +x feature_counts_hisat2.sh
chmod +x feature_counts_star.sh
# Activamos el entorno
conda activate subread
# Los ejecutamos
./feature_counts_hisat2.sh
./feature_counts_star.sh

# Desactivamos el entornp
conda deactivate
```

### Pseudoalineamiento con Salmon y tximport

Se ejecutó Salmon sobre los archivos FASTQ limpios. Se utilizó el parámetro -l A para la detección automática del tipo de librería.

```bash
# Nos movemos a la carpeta result donde crearemos la carpeta para guardar los archivos generados por salmon y nos movemos a la carpeta 
cd /export/storage/users/andreavg/transcriptomica/results
mkdir salmon
cd salmon

# Dentro de salmon crearemos las carpetas para paired-end y single-end
mkdir single-end
mkdir paired-end

# Nos movemos a la carpeta de src
cd /export/storage/users/andreavg/transcriptomica/src
# Crearemos el script para alinear con salmon tanto para paired-end y single-end.
nano salmon_alineamientos.sh

## ============================================== SCRIPT ============================================== ##
#!/bin/bash

# Asignamos las variables para nuestro indice de referencia
index="/export/space3/users/silvanac/transcriptomica_2026/indexes/mm39.gencode.M36.salmon/"
# Asignamos nuestras variables para nuestros archivos de tiempo finlaes
Salmon_time_file_SE="/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/Salmon_time_SE.txt"
Salmon_time_file_PE="/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/Salmon_time_PE.txt"

# Iniciamos nuestro ciclo
for R1 in /export/storage/users/andreavg/transcriptomica/data/clean/*_1_clean.fastq.gz; do

# Asignamos el nombre de la muestra a partir del nombre del archivo

	base=$(basename $R1 _1_clean.fastq.gz)

# Definimos donde se encuentran los archivos R2

	R2="/export/storage/users/andreavg/transcriptomica/data/clean/${base}_2_clean.fastq.gz"

	echo "Procesando muestra: $base"
# Ejecutamos el comando de alineamiento con salmon para single-end
	echo "Single-end para $base"
#Flags:
# -i: Nos indica en que lugar se encuentra el index
# -l: Se utiliza este parámetro para definir las propiedades de la biblioteca de secuenciación, en este caso utilizamos `A` para que Salmon analice una parte de las lecturas para decidir por sí mismo cuál es la configuración correcta
# -p: Número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
# -r: Nos indica que las lecturas son unpaired y especifica 
# -o: Indicamos donde guardaremos los resultados
# time se utiliza para medir el tiempo de ejecución del comando salmon y se redirige la salida de error estándar, que es donde se imprime el tiempo, al archivo correspondiente para single-end.
# usando las {} se asegura que el tiempo se mida solo para el comando salmon y no para otros comandos que puedan estar en el script.
	{ time salmon quant -i "$index" \
		-l A -p 5 \
		-r "$R1" \
		-o "/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/${base}_SE_" ; } 2>> "$Salmon_time_file_SE" \

# Ejecutamos hisat2 con paired-end
	#echo "Paired-end para $base"
#Flags:
# -i: Nos indica en que lugar se encuentra el index
# -l: Se utiliza este parámetro para definir las propiedades de la biblioteca de secuenciación, en este caso utilizamos `A` para que Salmon analice una parte de las lecturas para decidir por sí mismo cuál es la configuración correcta
# -p: Número de hilos a usar, en este caso usaremos 8 para manejar la carga de trabajo de manera eficiente
# -1: indica el archivo que contiene las lecturas hacia adelante 
# -2: Indica el archivo de las lecturas reversas
# -o: Indicamos donde guardaremos los resultados
# time se utiliza para medir el tiempo de ejecución del comando salmon y se redirige la salida de error estándar, que es donde se imprime el tiempo, al archivo correspondiente para single-end.
# usando las {} se asegura que el tiempo se mida solo para el comando salmon y no para otros comandos que puedan estar en el script.
	{ time salmon quant -i "$index" \
		-l A -p 5 \
		-1 "$R1" \
		-2 "$R2" \
		-o "/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/${base}_PE_" ; } 2>> "$Salmon_time_file_PE"
			

done

echo "Alineamiento completado para todas las muestras"
## ============================================== SCRIPT ============================================== ##

# Le damos permisos de ejecución, activamos entorno y ejecutamos
chomd +x salmon_alineamientos.sh
conda activate salmon
./salmon_alineamientos.sh

# Desactivamos el entorno
conda deactivate
```

Dado que Salmon cuantifica a nivel de transcrito, se empleó tximport en R junto con una tabla de correspondencia tx2gene (generada desde el GFF) para consolidar los conteos a nivel de gen.
Para esto se creo un script en R.

```bash
# Creamos las carpetas en donde guardaremos los archivos de salida de tximport
cd /export/storage/users/andreavg/transcriptomica/results/salmon/single_end
mkdir txi
cd ../paired_end
mkdir txi


# Creamos el script desde la carpeta src, el script tiene el nombre de `abudancia_tx2gene.R`

## ============================================== SCRIPT ============================================== ##
# Cargamos las librerias necesarias
library(GenomicFeatures)
library(tximport)
library(txdbmaker)

# Para convertir las abundancias de transcritos a genes, primero se generará una tabla de mapeo (tx2gene) 
# que relaciona cada transcrito con su gen correspondiente usando la anotación de ratón Gencode vM36.

# Primero creamos la base de datos de TxDb a partir del archivo GFF
txdb <- makeTxDbFromGFF("/export/storage/users/andreavg/transcriptomica/data/gencode/gencode.vM36.chr_patch_hapl_scaff.annotation.gff3.gz")

# Extraemos los IDs de transcritos que existen en txdb
# NOTA: keytype = "TXNAME" nos indica que solo queremos los nombres de los transcritos
tx_ids <- keys(txdb, keytype = "TXNAME")

# Creamos un data frame con los IDs de transcritos y sus correspondientes genes
# NOTA: keys = tx_ids nos indica que solo queremos los IDs de transcritos que existen en txdb
# columns = "GENEID" nos indica que solo queremos los IDs de genes correspondientes a cada transcrito
tx2gene <- select(txdb, keys = tx_ids, keytype = "TXNAME", columns = "GENEID")

# Quitamos la versión de los IDs de transcritos en tx2gene
# Esto es necesario porque los quant.sf de Salmon no incluyen la versión
# NOTA: sub("\\..*", "", tx2gene$TXNAME) nos indica que se reemplace cualquier cosa después del punto (incluyendo el punto) por una cadena vacía
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)

# Vamos a crear una varible con la ruta de los archivos de abundancia de transcritos generados por Salmon
# Creamos una lista de las muestras que queremos procesar
muestras <- c("SRR9126694", "SRR9126930", "SRR9126963", "SRR9127454", "SRR9127455", "SRR9127457", "SRR9127584")

# Construimos las rutas de los archivos single-end
files_SE <- file.path("/export/storage/users/andreavg/transcriptomica/results/salmon/single_end", paste0(muestras, "_SE_"), "quant.sf")
# Asignamos los nombres de las muestras a los archivos
names(files_SE) <- muestras


# Construimos las rutas de los archivos paired-end
files_PE <- file.path("/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end", paste0(muestras, "_PE_"), "quant.sf")
# Asignamos los nombres de las muestras a los archivos
names(files_PE) <- muestras

# Ejecutamos tximport para convertir las abundancias de transcritos a genes utilizando el mapeo tx2gene

# Para los archivos single-end
# NOTA: ignoreTxVersion = TRUE nos indica que se ignore la versión de los transcritos al hacer el mapeo, 
# lo cual es útil si los nombres de los transcritos en los archivos de Salmon no incluyen la versión.
txi_SE <- tximport(files_SE, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Para los archivos paired-end
# NOTA: ignoreTxVersion = TRUE nos indica que se ignore la versión de los transcritos al hacer el mapeo, 
# lo cual es útil si los nombres de los transcritos en los archivos de Salmon no incluyen la versión.
txi_PE <- tximport(files_PE, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Guardamos los resultados en archivos
saveRDS(txi_SE, file = "/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/txi/txi_SE.rds")
saveRDS(txi_PE, file = "/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/txi/txi_PE.rds")
## ============================================== SCRIPT ============================================== ##

# Debido que el script necesita ciertas librerias que se necesitan. Por ello necesitamos un entorno en donde podamos descargar las librerias, para ello tenemos uno propio. 
mamba activate /home/andreavg/.conda/envs/.main_mamba

#Corremos el script
Rscript abundancia_tx2gene.R

# Desactivamos el entorno
mamba deactivate 
```

Guardamos el archivo como un RDS para poder, sin embargo, este formato es esclusivo de R. Para permite que la información sea accesible fuera del entorno de programación de R y para que que los datos puedan ser procesados o integrados fácilmente en otros flujos de trabajo, decidimos exportarlos a formato TSV.
El archivo RDS obtenido a través de tximport funciona como un contenedor que organiza toda la información de la cuantificación de los genes. De toda la estructura de datos disponible, decidimos enfocarnos únicamente en las columnas de Counts y Abundance. Esto se debe a que los Counts (conteos) son el insumo esencial para realizar los cálculos estadísticos de expresión diferencial, mientras que la Abundance (abundancia en TPM) nos da una visión ya normalizada, lo que facilita comparar los niveles de expresión entre distintas muestras de manera más directa y comprensible. 

```bash
# El script para convertir de formato RDS a TSV se encuentra en la carpeta src, guardado con el nombre de `rds_to_tsv.R`

## ============================================== SCRIPT ============================================== ##
# Exportamos las tablas guardadas anteriormente
txi_SE <- readRDS("/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/txi/txi_SE.rds")
txi_PE <- readRDS("/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/txi/txi_PE.rds")

# Hacemos una función para convertir las tablas a formato TSV solamente de counts y abundance
rds_to_tsv <- function(txi, base_rute) {
    # Extraemos la base del nombre del archivo sin la extensión .rds
    base <- sub("\\.rds$", "", base_rute)
    # Escribimos las tablas de counts y abundance en formato TSV
    write.table(txi$counts, file = paste0(base, "_counts.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA )
    write.table(txi$abundance, file = paste0(base, "_abundance.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}

# Convertimos las tablas a formato TSV
rds_to_tsv(txi_SE,"/export/storage/users/andreavg/transcriptomica/results/salmon/single_end/txi/txi_SE.rds")
rds_to_tsv(txi_PE,"/export/storage/users/andreavg/transcriptomica/results/salmon/paired_end/txi/txi_PE.rds")
## ============================================== SCRIPT ============================================== ##

# Ejecutamos el script
Rscript rds_to_tsv.R
```
# Discusión

El análisis comparativo de la cuantificación y el pseudoalineamiento en las muestras de músculo de Mus musculus permite extraer conclusiones técnicas valiosas sobre la eficiencia y precisión de las herramientas empleadas.

Al analizar los resultados de featureCounts (Tabla 1), se observa una diferencia marcada entre los protocolos de secuenciación. Las muestras procesadas en modo Paired-End (PE) obtuvieron sistemáticamente porcentajes de asignación más altos (hasta un 84.48% con HISAT2) en comparación con el modo Single-End (SE).

Esta superioridad de los datos PE se debe a que el par de lecturas proporciona un contexto genómico doble, lo que reduce drásticamente la incertidumbre al asignar fragmentos a genes específicos. Es notable que, aunque STAR es altamente preciso en el alineamiento, HISAT2 mostró una tasa de asignación de lecturas a genes ligeramente mayor en este dataset, lo que sugiere una mayor sensibilidad en la detección global de la señal. 

| Muestra    | Alineador | Modo        | Assigned   | No Features | Ambiguity | Total Reads | % Assigned |
|------------|-----------|-------------|------------|-------------|-----------|-------------|------------|
| SRR9126694 | HISAT2    | Single-End  | 11,405,414 | 5,876,556   | 494,162   | 17,776,132  | 64.16%     |
| SRR9126930 | HISAT2    | Single-End  | 13,282,665 | 3,182,195   | 626,679   | 17,091,539  | 77.71%     |
| SRR9126963 | HISAT2    | Single-End  | 9,990,139  | 3,044,754   | 413,221   | 13,448,114  | 74.29%     |
| SRR9127454 | HISAT2    | Single-End  | 10,446,870 | 3,365,046   | 478,517   | 14,290,433  | 73.10%     |
| SRR9127455 | HISAT2    | Single-End  | 9,326,793  | 2,360,821   | 421,883   | 12,109,497  | 77.02%     |
| SRR9127457 | HISAT2    | Single-End  | 12,817,991 | 2,815,358   | 629,982   | 16,263,331  | 78.82%     |
| SRR9127584 | HISAT2    | Single-End  | 5,393,779  | 2,394,778   | 244,285   | 8,032,842   | 67.15%     |
| SRR9126694 | HISAT2    | Paired-End  | 21,753,713 | 8,690,721   | 968,832   | 31,413,266  | 69.25%     |
| SRR9126930 | HISAT2    | Paired-End  | 25,415,581 | 3,695,601   | 1,219,330 | 30,330,512  | 83.80%     |
| SRR9126963 | HISAT2    | Paired-End  | 17,903,453 | 4,153,309   | 746,158   | 22,802,920  | 78.51%     |
| SRR9127454 | HISAT2    | Paired-End  | 20,295,788 | 5,024,370   | 943,410   | 26,263,568  | 77.28%     |
| SRR9127455 | HISAT2    | Paired-End  | 18,097,975 | 3,003,622   | 828,769   | 21,930,366  | 82.52%     |
| SRR9127457 | HISAT2    | Paired-End  | 24,730,140 | 3,312,224   | 1,230,360 | 29,272,724  | 84.48%     |
| SRR9127584 | HISAT2    | Paired-End  | 10,397,323 | 3,837,409   | 478,360   | 14,713,092  | 70.67%     |
| SRR9126694 | STAR      | Single-End  | 11,664,604 | 7,564,213   | 492,479   | 19,721,296  | 59.15%     |
| SRR9126930 | STAR      | Single-End  | 13,560,046 | 4,635,578   | 625,745   | 18,821,369  | 72.05%     |
| SRR9126963 | STAR      | Single-End  | 10,178,531 | 3,942,084   | 410,788   | 14,531,403  | 70.05%     |
| SRR9127454 | STAR      | Single-End  | 10,628,628 | 4,088,525   | 477,322   | 15,194,475  | 69.95%     |
| SRR9127455 | STAR      | Single-End  | 9,459,415  | 3,380,610   | 420,059   | 13,260,084  | 71.34%     |
| SRR9127457 | STAR      | Single-End  | 13,045,590 | 4,012,544   | 628,011   | 17,686,145  | 73.76%     |
| SRR9127584 | STAR      | Single-End  | 5,473,944  | 2,758,043   | 243,259   | 8,475,246   | 64.59%     |
| SRR9126694 | STAR      | Paired-End  | 23,814,842 | 11,318,262  | 985,864   | 36,118,968  | 65.93%     |
| SRR9126930 | STAR      | Paired-End  | 27,400,885 | 6,066,992   | 1,247,603 | 34,715,480  | 78.93%     |
| SRR9126963 | STAR      | Paired-End  | 19,374,278 | 5,885,362   | 758,926   | 26,018,566  | 74.46%     |
| SRR9127454 | STAR      | Paired-End  | 21,839,029 | 6,519,628   | 955,619   | 29,314,276  | 74.50%     |
| SRR9127455 | STAR      | Paired-End  | 19,523,682 | 4,671,164   | 839,758   | 25,034,604  | 77.99%     |
| SRR9127457 | STAR      | Paired-End  | 26,502,749 | 5,426,839   | 1,251,596 | 33,181,184  | 79.87%     |
| SRR9127584 | STAR      | Paired-End  | 11,171,478 | 4,590,958   | 486,206   | 16,248,642  | 68.75%     |
Tabla 1. Tabla comparativa del summary de featureCounts

En cuanto a la velocidad del procesamiento entre alineadores (STAR y HISAT2) contra los pseudoalineamientos (SALMON). Mientras que el flujo anterior requirió tiempos considerables (aprox. 39 minutos para STAR y 68 minutos para HISAT2), Salmon completó la cuantificación en una fracción del tiempo, con promedios de entre 1.5 y 3 minutos por muestra.

Esta eficiencia del pseudoalineamiento se explica porque Salmon no intenta reconstruir la posición base por base en el genoma, sino que asigna directamente los fragmentos al transcriptoma. Sin embargo, como se discutió en clase, la rapidez de Salmon requiere el paso adicional con tximport en R para colapsar los resultados a nivel de gen y permitir una comparación directa con los métodos basados en alineamiento genómico.

# Conclusiones

Tras completar las seis cuantificaciones requeridas por muestra, se concluye que la alta consistencia en las tasas de asignación confirma el éxito de los procesos de limpieza y filtrado previos. 
El empleo del modo paired-end resultó fundamental, ya que además de incrementar la precisión en el conteo de lecturas, permitió generar visualizaciones de cobertura exentas de duplicaciones artificiales gracias al uso de filtros de paridad. Si bien el pseudoalineamiento con Salmon destaca como una alternativa superior en términos de velocidad y optimización de recursos , presenta la limitante técnica de depender de un transcriptoma de referencia ya generado, lo que restringe el descubrimiento de nuevos transcritos que el mapeo genómico tradicional sí permite identificar. Por esta razón, el flujo de trabajo convencional sigue siendo una herramienta de diagnóstico esencial cuando se requiere una inspección visual detallada o el análisis de la arquitectura del genoma. En definitiva, la elección entre alineamiento y pseudoalineamiento no debe considerarse excluyente; la coordinación de herramientas como STAR e HISAT2 con Salmon , integrada mediante tximport, fortalece la validez técnica de los resultados y proporciona una base robusta para una interpretación biológica profunda.


# Referencias

1. GENCODE - Mouse Release M36. (s. f.). https://www.gencodegenes.org/mouse/release_M36.html 
2. Bash Reference Manual. (s. f.). https://www.gnu.org/software/bash/manual/bash.html
3. Genomic coordinates of problematic genomic regions. (s. f.). https://dozmorovlab.github.io/excluderanges/index.html. 
4. tximport. (s. f.). Bioconductor. https://bioconductor.org/packages/release/bioc/html/tximport.html
5. txdbmaker. (s. f.). Bioconductor. https://bioconductor.org/packages/release/bioc/html/txdbmaker.html 