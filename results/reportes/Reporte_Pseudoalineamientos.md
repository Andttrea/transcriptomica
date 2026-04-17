# Reporte de pseudoalineamientos 

## Introducción 

Este reporte documenta el análisis bioinformático, incluyendo la 

Todos los datos y resultados se encuentran ubicados en el servidor chaac de la licenciatura, en el directorio `/export/storage/users/andreavg/transcriptomica`.

## Estructura del repostorio 



## Conversión de SAM a BAM 

Debido a que los archivos BAM generados por los alineadores Hisat2 y STAR pesan demasiado, lo que debemos de hacer es convertirlos a SAM. Para realizar esta tarea, utilizaremos el siguiente script. Que fue generado en la carpeta src

```bash
 # Creamos el script para hacer la conversión de SAM a BAM
 nano sam_to_bam.sh 
```

El script es el siguiente: 

```bash
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
```

```bash
# Le damos permisos al script
chomd +x sam_to_bam.sh
# Ejecutamos el script
./sam_to_bam.sh
```

Nos debemos de asegurar que los archivos bam contengan algo y ver cuantas reads alineadas se obtuvieron por  lo que hacemos esto:

```bash
# Recorremos las carpetas hisat2 y star, tanto para single_end y paired_end
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

# Después de comprobar que los archivos BAM si tuvieran reads, entonces borramos los archivos BAM.
rm *.sam 
```

### Generación de pileups(BigWig)

Para la visualización de los datos en navegadores como IGV, se generaron archivos BigWig (.bw). Se utilizó la herramienta bamCoverage de deeptools, aplicando la flag 64 para contar únicamente la primera lectura del par en datos Paired-End, evitando así la duplicación artificial de la señal. Se aplicó normalización BPM (Bins Per Million) para permitir la comparación visual entre muestras.

bamCoverage necesita un archivo BED  con el fin de evitar coordenadas genómicas conflictivas en análisis genómicos. En este caso utlizamos `mm39.excluderanges`, el cual se encontró en [excluderanges](https://dozmorovlab.github.io/excluderanges/#excluderanges-genomic-ranges-of-problematic-genomic-regions). 
El archivo se encuentra en el siguiente path: `/export/storage/users/andreavg/transcriptomica/data/bed_file`

Solo se corrió bamCOverage para paired-end. 

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

# Habilitamos los permisos y corremos el script
chmod +x generar_pileups.sh
./generar_pileups.sh
```

### Cuantificación por feautureCounts

Se generaron matrices de conteo de lecturas para STAR e HISAT2. Para los datos Paired-End, se utilizaron los parámetros -p (indicando fragmentos emparejados) y -B (exigiendo que ambos extremos estén mapeados). Se utilizó el archivo de anotación de Gencode M36 (.gff3), el cual se encontró en Gencode y se encuentra en el siguiente path: `/export/storage/users/andreavg/transcriptomica/data/gencode`

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

# Primero creamos la base de datos de TxDb a partir del archivo GFF, el cual tambien utilizamos para los featurecounts
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

# Debido que el script necesita ciertas librerias que se necesitan. Se debe de activar un entorno para poder correrlo
conda activate .main_mamba

#Corremos el script
Rscript abundancia_tx2gene.R



```