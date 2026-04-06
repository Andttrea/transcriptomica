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