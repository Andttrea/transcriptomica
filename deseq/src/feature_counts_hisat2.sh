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