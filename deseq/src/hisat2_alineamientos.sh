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
