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
