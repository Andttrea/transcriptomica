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
